#!/usr/bin/env python3
"""Publish completed BASINS staging outputs in an ASPIRE-style layout."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import shutil
import tempfile
from pathlib import Path


PLOT_SUFFIXES = {".svg", ".pdf", ".png", ".jpg", ".jpeg", ".html"}
PUBLIC_DIRS = ("modules", "intermediates", "references", "summary", "logs")
MODULES = {
    "biochem_processing": "biochemical_processing",
    "env_missingness_sensitivity": "missingness_sensitivity",
    "env_pca": "environmental_pca",
    "env_continuous_sections": "continuous_time_depth_sections",
    "env_compartments_selectk": "compartment_selection",
    "env_compartments_gmm": "gmm_compartments",
    "env_o2_soft_compartments": "oxygen_compartments",
    "env_hybrid_soft_compartments": "hybrid_compartments",
    "env_compare_compartments": "compartment_comparison",
    "env_o2_split_by_gmm": "oxygen_gmm_subcompartments",
    "env_stratification_index": "stratification",
    "env_state_transitions": "state_transitions",
    "env_succession_graphs": "succession_graphs",
    "env_compartment_feature_assoc": "feature_associations",
    "eof_pca": "eof_analysis",
    "eof_states": "eof_states",
    "eof_plots": "eof_modes",
    "gapseq_media": "gapseq_media",
    "genome_modeling": "genome_modeling",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def move_tree(source: Path, module: Path) -> None:
    if not source.is_dir():
        return
    for path in sorted((item for item in source.rglob("*") if item.is_file()), key=lambda item: len(item.parts), reverse=True):
        relative = path.relative_to(source)
        structural = {"tables", "plots", "figures", "results", "data"}
        parts = [part for part in relative.parts[:-1] if part.lower() not in structural]
        bucket = "plots" if path.suffix.lower() in PLOT_SUFFIXES else "tables"
        destination = module / bucket / Path(*parts, path.name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(destination))
    shutil.rmtree(source, ignore_errors=True)
    (module / "tables").mkdir(parents=True, exist_ok=True)
    (module / "plots").mkdir(parents=True, exist_ok=True)


def write_manifests(root: Path) -> list[tuple[str, int, int]]:
    modules = root / "modules"
    summary = root / "summary"
    table_dir = summary / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    counts = []
    for module in sorted(path for path in modules.iterdir() if path.is_dir()):
        tables = list((module / "tables").rglob("*")) if (module / "tables").is_dir() else []
        plots = list((module / "plots").rglob("*")) if (module / "plots").is_dir() else []
        table_count = sum(path.is_file() for path in tables)
        plot_count = sum(path.is_file() for path in plots)
        counts.append((module.name, table_count, plot_count))
        for path in sorted(item for item in module.rglob("*") if item.is_file()):
            relative = path.relative_to(root)
            output_type = "plots" if "plots" in relative.parts else "tables"
            rows.append((module.name, output_type, str(relative), path.stat().st_size, sha256(path)))
    with (table_dir / "module_output_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("module", "output_type", "relative_path", "size_bytes", "sha256"))
        writer.writerows(rows)
    with (table_dir / "module_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("module", "table_files", "plot_files"))
        writer.writerows(counts)
    return counts


def write_report(root: Path, counts: list[tuple[str, int, int]]) -> None:
    summary = root / "summary"
    (summary / "plots").mkdir(parents=True, exist_ok=True)
    (summary / "report").mkdir(parents=True, exist_ok=True)
    maximum = max((tables + plots for _, tables, plots in counts), default=1)
    bars = []
    for index, (name, tables, plots) in enumerate(counts):
        y = 55 + index * 28
        tw = 580 * tables / maximum
        pw = 580 * plots / maximum
        bars.append(f'<text x="190" y="{y + 12}" text-anchor="end">{html.escape(name)}</text>')
        bars.append(f'<rect x="205" y="{y}" width="{tw:.1f}" height="10" fill="#0072B2"/>')
        bars.append(f'<rect x="{205 + tw:.1f}" y="{y}" width="{pw:.1f}" height="10" fill="#D55E00"/>')
        bars.append(f'<text x="810" y="{y + 10}">{tables} tables, {plots} plots</text>')
    height = max(150, 90 + len(counts) * 28)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="920" height="{height}">
<style>text{{font:12px "Times New Roman",serif;fill:#202020}}</style><rect width="100%" height="100%" fill="white"/>
<text x="20" y="28" style="font-size:20px;font-weight:bold">BASINS module outputs</text>{''.join(bars)}</svg>'''
    (summary / "plots" / "module_output_summary.svg").write_text(svg)
    cards = "".join(
        f'<article><h2>{html.escape(name.replace("_", " ").title())}</h2><p>{tables} tables; {plots} plots.</p>'
        f'<a href="../../modules/{html.escape(name)}/tables/">Tables</a> '
        f'<a href="../../modules/{html.escape(name)}/plots/">Plots</a></article>'
        for name, tables, plots in counts
    )
    log_links = "".join(
        f'<li><a href="../../logs/{name}">{label}</a></li>'
        for name, label in (
            ("nextflow_report.html", "Nextflow execution report"),
            ("nextflow_timeline.html", "Execution timeline"),
            ("nextflow_trace.tsv", "Task trace"),
            ("nextflow_dag.html", "Workflow DAG"),
            ("controller.log", "Controller log"),
            ("launch_command.txt", "Launch command"),
            ("nextflow_version.txt", "Nextflow version"),
        ) if (root / "logs" / name).is_file()
    )
    report = f'''<!doctype html><html><head><meta charset="utf-8"><title>BASINS run report</title>
<style>body{{font-family:"Times New Roman",serif;max-width:1100px;margin:2rem auto;color:#202020}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:1rem}}article{{border:1px solid #ddd;border-radius:8px;padding:1rem}}a{{color:#006a8e;margin-right:1rem}}img{{max-width:100%}}</style></head>
<body><h1>BASINS run report</h1><p>Non-interpretive inventory of environmental processing, compartment, stratification, EOF, and diagnostic outputs.</p>
<h2>Output inventory</h2><img src="../plots/module_output_summary.svg" alt="BASINS module output counts"><div class="grid">{cards}</div>
<h2>Master summary</h2><p><a href="../tables/basin_run_overview.tsv">Run overview</a> <a href="../tables/basin_key_outputs.tsv">Key outputs</a> <a href="../tables/basin_module_inventory.tsv">Module inventory</a></p>
<h2>Nextflow run details</h2><ul>{log_links}</ul><p>Exact paths and checksums are in <a href="../tables/module_output_manifest.tsv">module_output_manifest.tsv</a>.</p></body></html>'''
    (summary / "report" / "BASIN_run_report.html").write_text(report)


def organize(root: Path, config: Path | None) -> None:
    scientific = root
    modules = root / "modules"
    summary = root / "summary"
    for directory in (modules, root / "intermediates", root / "references", summary / "tables", summary / "plots", summary / "report", root / "logs"):
        directory.mkdir(parents=True, exist_ok=True)
    for source_name, module_name in MODULES.items():
        move_tree(scientific / source_name, modules / module_name)
    master = scientific / "master_summary"
    if master.is_dir():
        for path in master.iterdir():
            if path.is_file():
                shutil.move(str(path), str(summary / "tables" / path.name))
        shutil.rmtree(master, ignore_errors=True)
    if config and config.is_file():
        shutil.copy2(config, summary / "tables" / "run_config.yml")
    counts = write_manifests(root)
    write_report(root, counts)


def publish(staging: Path, output: Path, config: Path | None) -> None:
    temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.organizing-", dir=output.parent))
    try:
        shutil.copytree(staging, temporary, dirs_exist_ok=True)
        organize(temporary, config)
        backup = Path(tempfile.mkdtemp(prefix=f".{output.name}.backup-", dir=output.parent))
        installed, moved = [], []
        try:
            for name in PUBLIC_DIRS:
                source, destination = temporary / name, output / name
                if destination.exists() or destination.is_symlink():
                    shutil.move(str(destination), str(backup / name)); moved.append(name)
                if source.exists() or source.is_symlink():
                    shutil.move(str(source), str(destination)); installed.append(name)
        except Exception:
            for name in installed:
                target = output / name
                if target.is_dir() and not target.is_symlink(): shutil.rmtree(target)
                elif target.exists() or target.is_symlink(): target.unlink()
            for name in moved:
                shutil.move(str(backup / name), str(output / name))
            raise
        finally:
            shutil.rmtree(backup, ignore_errors=True)
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()
    args.output_dir.resolve().mkdir(parents=True, exist_ok=True)
    publish(args.staging_dir.resolve(), args.output_dir.resolve(), args.config.resolve() if args.config else None)


if __name__ == "__main__":
    main()
