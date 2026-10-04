#!/usr/bin/env python3
"""Reproduce the SI 2009–2012 reviewer subset from the cleaned source CSVs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(chemistry, ctd, output):
    output.mkdir(parents=True, exist_ok=True)
    tables = []
    keys = ('Cruise', 'Year', 'Month', 'Day')
    for path in (chemistry, ctd):
        with path.open(newline='') as stream:
            reader = csv.DictReader(stream)
            tables.append((reader.fieldnames, list(reader)))
    candidates = [r for r in tables[0][1] if 2009 <= int(r['Year']) <= 2012]
    ctd_profiles = {tuple(r[k] for k in keys) for r in tables[1][1]}
    excluded = sorted({tuple(r[k] for k in keys) for r in candidates} - ctd_profiles)
    selected = [r for r in candidates if tuple(r[k] for k in keys) in ctd_profiles]
    profiles = {tuple(r[k] for k in keys) for r in selected}
    matched = [r for r in tables[1][1] if tuple(r[k] for k in keys) in profiles]
    if not selected or {tuple(r[k] for k in keys) for r in matched} != profiles:
        raise ValueError('Every selected chemistry cruise/date must have a CTD profile')
    products = {}
    for name, source, fields, rows in [
        ('chemistry.csv', chemistry, tables[0][0], selected),
        ('ctd.csv', ctd, tables[1][0], matched),
    ]:
        target = output/name
        with target.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        products[name] = dict(rows=len(rows), bytes=target.stat().st_size,
                              sha256=sha256(target), source_filename=source.name,
                              source_sha256=sha256(source))
    manifest = dict(selection='All chemistry rows from 2009–2012 cruises with a matching CTD profile; all CTD rows with matching Cruise, Year, Month, Day. Original field values retained; CSV line endings normalized.',
                    excluded_profiles_without_ctd=[dict(zip(keys, key)) for key in excluded],
                    cruises=len(profiles), chemistry_depths=sorted({float(r['Depth']) for r in selected}),
                    article_doi='10.1038/sdata.2017.159', dataset_doi='10.5061/dryad.nh035',
                    source_note='Subset of the cleaned SI tables supplied by the BASINS authors; not a byte-for-byte extract of the archival Dryad files.',
                    files=products)
    (output/'provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chemistry', type=Path, required=True)
    parser.add_argument('--ctd', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.chemistry, args.ctd, args.output)
