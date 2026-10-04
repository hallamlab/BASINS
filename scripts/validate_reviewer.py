#!/usr/bin/env python3
"""Check the bundled SI reviewer inputs and key scientific outputs (stdlib only)."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def table(path):
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f'Missing or empty table: {path}')
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t' if path.suffix == '.tsv' else ','))


def validate(output):
    data = ROOT/'examples/reviewer/data'
    provenance = json.loads((data/'provenance.json').read_text())
    for name, entry in provenance['files'].items():
        path = data/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError(f'Reviewer input checksum mismatch: {name}')
        if len(table(path)) != entry['rows']:
            raise ValueError(f'Reviewer input row count mismatch: {name}')
    merged = table(output/'modules/biochemical_processing/tables/01_merged_nearest_depth.tsv')
    if len(merged) != 754:
        raise ValueError(f'Expected all 754 chemistry rows after merging, found {len(merged)}')
    def observation_key(row):
        return tuple(float(row[k]) for k in ('Cruise', 'Year', 'Month', 'Day', 'Depth'))
    if Counter(map(observation_key, merged)) != Counter(map(observation_key, table(data/'chemistry.csv'))):
        raise ValueError('Chemistry observation identities changed during merge')
    matches = sum(r.get('CTD_Depth_Used', '').lower() not in ('', 'na', 'nan') for r in merged)
    if matches < 0.9 * len(merged):
        raise ValueError(f'Too few CTD depth matches: {matches}/{len(merged)}')
    print(f'Chemistry: {len(merged)} observations retained; {matches} matched to CTD')
    required = {
        'PCA': 'environmental_pca/tables/eigenvectors_scores.csv',
        'GMM': 'gmm_compartments/tables/compartments_assignments_smoothed.csv',
        'Oxygen': 'oxygen_compartments/tables/o2_compartments_assignments_smoothed.csv',
        'Hybrid': 'hybrid_compartments/tables/compartments_assignments_hybrid.csv',
        'Cruise EOF': 'eof_analysis/tables/eof_eigenvectors_scores_by_cruise.csv',
        'Cruise groups': 'eof_states/tables/gmm_selected_assignments.tsv',
    }
    pca_ids = None
    for label, relative in required.items():
        rows = table(output/'modules'/relative)
        if len(rows) < (10 if label.startswith('Cruise') else 100):
            raise ValueError(f'{label} unexpectedly has only {len(rows)} rows')
        if label == 'PCA':
            pca_ids = {r['cruise_year_month_depth'] for r in rows}
            for row in rows:
                for column in ('PC1', 'PC2'):
                    if not math.isfinite(float(row[column])):
                        raise ValueError(f'Non-finite PCA coordinate: {column}')
        if label in ('GMM', 'Oxygen', 'Hybrid') and {r['cruise_year_month_depth'] for r in rows} != pca_ids:
            raise ValueError(f'{label} observation IDs differ from the PCA cohort')
        print(f'{label}: {len(rows)} rows')
    sparse_path = output/'modules/eof_analysis/tables/eof_sparse_feature_pc_spearman.csv'
    with sparse_path.open(newline='') as stream:
        columns = csv.DictReader(stream).fieldnames or []
    if not {'feature', 'PC', 'spearman_r', 'n_cruises_used', 'coverage'} <= set(columns):
        raise ValueError('EOF sparse correlation table lost its schema (empty tables still need headers)')
    trace = table(output/'logs/nextflow_trace.tsv')
    if not trace or any(r['status'] not in ('COMPLETED', 'CACHED') or r['exit'] != '0' for r in trace):
        raise ValueError('Nextflow trace contains an incomplete or unsuccessful task')
    stages = {r['name'].split(' (')[0] for r in trace}
    for stage in ('BIOCHEM_MERGE', 'BIOCHEM_EIGENVECTORS', 'BIOCHEM_SELECTK',
                  'BIOCHEM_GMM', 'BIOCHEM_O2_SOFT', 'BIOCHEM_HYBRID',
                  'BIOCHEM_STATE_TRANSITIONS', 'BIOCHEM_SUCCESSION_GRAPH',
                  'BIOCHEM_EOF_PIPELINE', 'BIOCHEM_EOF_STATE_CLUSTER',
                  'BIOCHEM_WITHIN_GMM_HDBSCAN', 'BIOCHEM_CONTINUOUS_SECTIONS', 'MASTER_SUMMARY'):
        if stage not in stages:
            raise ValueError(f'Expected reviewer stage missing from trace: {stage}')
    report = output/'summary/report/BASIN_run_report.html'
    if not report.is_file() or report.stat().st_size < 1000:
        raise ValueError('Missing or empty HTML report')
    print('PASS: input integrity, chemistry row preservation, key analyses and report')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    try:
        validate(args.output)
    except (ValueError, KeyError, FileNotFoundError) as exc:
        parser.exit(1, f'Reviewer validation failed: {exc}\n')
