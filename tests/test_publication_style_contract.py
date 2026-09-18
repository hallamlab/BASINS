import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT = Path(__file__).parents[1]
PROCESSES = PROJECT / "processes"
SPEC = importlib.util.spec_from_file_location(
    "shared_plot_style", PROCESSES / "shared_plot_style.py"
)
STYLE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(STYLE)


PLOT_MARKERS = ("matplotlib", "pyplot", "seaborn", "savefig", "plotly")


def test_every_matplotlib_renderer_installs_contract():
    missing = []
    excluded = {"shared_plot_style.py", "shared_plot_export.py"}
    for path in sorted(PROCESSES.rglob("*.py")):
        if path.name in excluded:
            continue
        text = path.read_text(errors="ignore")
        if any(marker in text for marker in PLOT_MARKERS):
            if "install_publication_style" not in text:
                missing.append(str(path.relative_to(PROJECT)))
    assert not missing, f"Plot renderers missing publication contract: {missing}"


def test_no_renderer_declares_a_conflicting_font():
    conflicts = []
    for path in sorted(PROCESSES.rglob("*.py")):
        if path.name == "shared_plot_style.py":
            continue
        text = path.read_text(errors="ignore").lower()
        if (
            "dejavu sans" in text
            or "source sans" in text
            or '"font.family": "sans-serif"' in text
            or "font-family:sans-serif" in text
        ):
            conflicts.append(str(path.relative_to(PROJECT)))
    assert not conflicts, f"Renderers declare non-contract fonts: {conflicts}"


def test_required_font_resolves_to_times_new_roman_file():
    resolved = STYLE.require_font()
    assert resolved.is_file()
    assert "times" in resolved.name.lower()


def test_save_hook_reapplies_times_new_roman(tmp_path):
    STYLE.install_publication_style()
    fig, ax = plt.subplots(figsize=(4, 3))
    title = ax.set_title("BASINS", fontfamily="DejaVu Sans")
    ax.set_xlabel("Environmental score", fontfamily="DejaVu Sans")
    ax.plot([0, 1], [0, 1])
    fig.savefig(tmp_path / "contract.pdf")

    assert title.get_fontfamily()[0] == STYLE.FONT_FAMILY
    for text in fig.findobj(match=matplotlib.text.Text):
        assert text.get_fontfamily()[0] == STYLE.FONT_FAMILY
    plt.close(fig)


def test_basin_color_maps_are_not_rewritten(tmp_path):
    STYLE.install_publication_style()
    fig, ax = plt.subplots(figsize=(4, 3))
    image = ax.imshow([[0, 1], [1, 0]], cmap="viridis")
    fig.savefig(tmp_path / "colors.pdf")
    assert image.get_cmap().name == "viridis"
    plt.close(fig)


def test_biplot_labels_are_opaque_and_collision_free():
    STYLE.install_publication_style()
    fig, ax = plt.subplots(figsize=(8, 6))
    anchors = [
        (-1.2, 0.20),
        (-0.8, 0.25),
        (-0.5, 0.30),
        (0.5, 0.20),
        (0.8, 0.25),
        (1.2, 0.30),
    ]
    texts = [
        ax.text(
            x,
            y,
            f"Feature {index}",
            color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.2},
        )
        for index, (x, y) in enumerate(anchors)
    ]
    ax.scatter([-2, 2], [-1, 1], alpha=0)
    STYLE.layout_biplot_labels(ax, texts, anchors)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [text.get_window_extent(renderer=renderer) for text in texts]

    for index, left in enumerate(boxes):
        assert texts[index].get_bbox_patch().get_alpha() == 1.0
        for right in boxes[index + 1 :]:
            assert not left.overlaps(right)

    origin_x = ax.transData.transform((0, 0))[0]
    left_tip_x = min(ax.transData.transform(anchor)[0] for anchor in anchors)
    right_tip_x = max(ax.transData.transform(anchor)[0] for anchor in anchors)
    for text, anchor, box in zip(texts, anchors, boxes):
        if anchor[0] < 0:
            assert box.x1 < left_tip_x < origin_x
        else:
            assert box.x0 > right_tip_x > origin_x
    plt.close(fig)


def test_shared_biplot_vectors_match_pca_aesthetics():
    STYLE.install_publication_style()
    cloud_scale, vector_scale = STYLE.calculate_biplot_vector_scale(
        [-2, -1, 0, 1, 2],
        [-1, -0.5, 0, 0.5, 1],
        [0.5, 1.0],
    )
    assert cloud_scale > 0
    assert vector_scale == cloud_scale

    fig, ax = plt.subplots(figsize=(4, 3))
    STYLE.draw_biplot_vector(
        ax, 1.0, 0.5, "red", "Core loading", cloud_scale
    )
    assert len(ax.patches) == 2
    assert [patch.get_linewidth() for patch in ax.patches] == [3.2, 2.2]
    assert ax.patches[0].get_facecolor()[:3] == (1.0, 1.0, 1.0)

    STYLE.draw_biplot_vector(
        ax, -1.0, 0.5, "blue", "Sparse correlation", cloud_scale
    )
    assert len(ax.lines) == 2
    assert [line.get_linewidth() for line in ax.lines] == [3.2, 2.2]
    assert all(line.get_linestyle() == "--" for line in ax.lines)
    assert len(ax.collections) == 1
    plt.close(fig)
