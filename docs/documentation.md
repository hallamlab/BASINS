# Build and host this guide

The README is the short entry point; `docs/` contains the full user guide. Scientific software and reference databases are not needed to build the documentation.

## Local preview

```bash
mamba create -n basins-docs -c conda-forge python=3.11 pip
mamba activate basins-docs
python -m pip install -r docs/requirements.txt
python -m sphinx -W --keep-going -b html docs _build/html
python -m http.server 8765 --bind 127.0.0.1 --directory _build/html
```

Open `http://localhost:8765`. `-W` makes documentation warnings fail the build.

## Read the Docs administration

1. Sign in to Read the Docs and connect GitHub.
2. Import `hallamlab/BASINS`. If the preferred project name is taken, choose an available project slug; the displayed title can remain BASINS.
3. Use `.readthedocs.yaml` at the repository root as the configuration path.
4. Select a branch that contains this configuration and the guide. If it is missing from the dropdown, resynchronize repository versions and activate that branch.
5. Build the selected version and inspect its build log. The default repository branch will not build this guide until these files are merged there.
6. After review and merge, point the production/latest version at the intended main development branch and activate stable tags as needed.
7. Add the confirmed public documentation URL to the repository description and README after the project builds successfully.

No credentials belong in this repository. The YAML installs only documentation dependencies. The explicit page allow-list in `docs/conf.py` excludes manuscript drafts, manuscript tables, and local analysis notes. Add new public pages to both that list and the appropriate toctree.

## Workflow figures

The public user-guide diagrams are in `docs/assets/`; Mermaid sources are in `docs/diagrams/`. They describe software behavior and are independent of unpublished manuscript figures. Main workflow SVG/PDF assets and white-background Mermaid previews stay readable in light and dark viewers. The secondary diagrams share a common canvas scale and are centered; click a preview to zoom.

To rebuild diagrams, install `docs/diagram-requirements.txt`, install Chromium with `python -m playwright install chromium`, then run `python scripts/render_workflow_diagrams.py`. The renderer downloads the pinned Mermaid library; ordinary documentation builds use the committed previews without a browser or network renderer.

## Workflow figures

The main SVG and PDF are generated from `docs/diagrams/main-workflow.json` by
`scripts/render_main_workflow.py`. The renderer uses the MP/ASPIRE conventions:
Times typography, numbered modules, computational diamonds, data circles,
blue inputs, green outputs, black arrows and an opaque white background.
The detailed workflow chapter maps the displayed stages to the public code.

After changing the JSON or Mermaid sources, install
`docs/diagram-requirements.txt`, install Playwright Chromium, and run
`python scripts/render_workflow_diagrams.py` from the repository root.
Review both the generated SVG and PDF. Main figures use their own tight canvas
at full page width; supporting diagrams share a centered 2240-unit canvas.
Do not copy unpublished manuscript figures or numerical results into these
public software diagrams.

The landing page and README use the brief workflow (`docs/diagrams/brief-workflow.json`); the workflow guide retains the detailed version (`main-workflow.json`). The renderer rebuilds both SVG and PDF versions.
