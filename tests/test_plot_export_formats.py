from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_export_module():
    path = ROOT / "processes" / "shared_plot_export.py"
    spec = importlib.util.spec_from_file_location("shared_plot_export", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_shared_export_writes_png_pdf_and_svg(tmp_path):
    module = load_export_module()
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])

    outputs = module.save_figure_all_formats(
        figure, tmp_path / "diagnostic.png", dpi=72
    )
    plt.close(figure)

    assert {path.suffix for path in outputs} == {".png", ".pdf", ".svg"}
    assert all(path.is_file() and path.stat().st_size > 0 for path in outputs)


def test_eof_mode_heatmap_reports_all_written_formats(tmp_path):
    processes_dir = ROOT / "processes"
    sys.path.insert(0, str(processes_dir))
    try:
        module = load_module(
            processes_dir / "eof_mode_plots" / "eof_mode_plots.py",
            "eof_mode_plots",
        )
    finally:
        sys.path.pop(0)

    outputs = module.plot_heatmap(
        A=np.array([[0.5, -0.25], [-0.1, 0.3]]),
        var_names=["oxygen", "nitrate"],
        depth_grid=np.array([0.0, 10.0]),
        title="EOF test",
        out_path=str(tmp_path / "eof_mode_EOF1"),
        figsize=(4.0, 3.0),
        dpi=72,
    )

    assert {path.suffix for path in outputs} == {".png", ".pdf", ".svg"}
    assert all(path.is_file() and path.stat().st_size > 0 for path in outputs)


def test_single_path_plotters_use_shared_export_contract():
    relative_paths = (
        "eigenvectors/env_eigenvectors.py",
        "eof_mode_plots/eof_mode_plots.py",
        "eof_pipeline/env_eof_pipeline.py",
        "eof_state_cluster/eof_state_clustering.py",
        "eof_state_cluster/cruise_group_season_benchmark.py",
        "feature_assoc/env_compartment_feature_assoc.py",
        "gmm/env_compartments_gmm.py",
        "hybrid/env_hybrid_compartment_builder.py",
        "o2_soft/env_compartments_o2_soft.py",
        "selectk/env_compartments_selectk.py",
        "state_transitions/env_state_transition_analysis.py",
        "succession_graph/env_succession_graph.py",
        "within_gmm_hdbscan/env_within_gmm_hdbscan.py",
    )
    for relative in relative_paths:
        source = (ROOT / "processes" / relative).read_text()
        assert "from shared_plot_export import save_figure_all_formats" in source
        assert "save_figure_all_formats(" in source


def test_configurable_plotters_default_to_all_static_formats():
    split_source = (
        ROOT / "processes" / "split_o2_by_gmm" / "env_split_o2_by_gmm.py"
    ).read_text()
    interpretation_source = (
        ROOT
        / "processes"
        / "eof_state_cluster"
        / "eof_state_interpretation.py"
    ).read_text()

    assert 'default="png,pdf,svg"' in split_source
    assert 'default="pdf,png,svg"' in interpretation_source


def test_compartment_biplot_writes_all_static_formats(tmp_path):
    path = (
        ROOT
        / "processes"
        / "compare_compartments"
        / "env_compare_compartments.py"
    )
    name = "env_compare_compartments_biplot_test"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(name, None)

    frame = pd.DataFrame(
        {
            "PC1": [-1.0, -0.5, 0.5, 1.0],
            "PC2": [0.5, -0.5, -0.25, 0.75],
            "component": [0, 0, 1, 1],
        }
    )
    loadings = pd.DataFrame(
        {
            "feature": ["Oxygen", "Nitrate"],
            "PC1": [-0.8, 0.6],
            "PC2": [0.2, -0.5],
        }
    )
    sparse = pd.DataFrame(
        {
            "feature": ["Methane", "Methane"],
            "PC": ["PC1", "PC2"],
            "spearman_r": [0.3, 0.7],
        }
    )
    cfg = type(
        "PlotConfig",
        (),
        {"save_pdf": True, "save_svg": True, "save_png": True, "png_dpi": 72},
    )()
    output = tmp_path / "compartment_biplot"
    module.plot_compartment_biplot(
        df=frame,
        label_col="component",
        palette={"0": "#1f77b4", "1": "#ff7f0e"},
        loadings_df=loadings,
        sparse_corr_df=sparse,
        title="Compartment biplot",
        out_base=str(output),
        cfg=cfg,
    )

    outputs = [output.with_suffix(suffix) for suffix in (".png", ".pdf", ".svg")]
    assert all(path.is_file() and path.stat().st_size > 0 for path in outputs)
    assert module.hybrid_o2_gmm_label("hyb_C0_G2") == "oxic-GMM2"
    assert module.hybrid_o2_gmm_label("hyb_C3_G4") == "anoxic-GMM4"
    hybrid_palette = module.build_hybrid_parent_palette(
        pd.Series(["hyb_C0_G0", "hyb_C0_G1", "hyb_C1_G0"])
    )
    assert set(hybrid_palette) == {"hyb_C0_G0", "hyb_C0_G1", "hyb_C1_G0"}
    assert hybrid_palette["hyb_C0_G0"] != hybrid_palette["hyb_C0_G1"]
