# Software and method citations

Record the BASINS repository revision or release with your methods, and cite the tools used by the enabled stages. The [workflow page](workflow.md) shows where they enter the analysis. This guide does not assign an unverified manuscript DOI to BASINS.

| Software | Use | Citation / project |
|---|---|---|
| HDBSCAN | Density-based clustering | McInnes, Healy & Astels (2017), *hdbscan: Hierarchical density based clustering*, JOSS 2(11):205. [DOI](https://doi.org/10.21105/joss.00205); [project](https://github.com/scikit-learn-contrib/hdbscan). |
| UMAP | Low-dimensional embeddings | [Project and citation guidance](https://github.com/lmcinnes/umap). |
| scikit-learn | Statistical learning and preprocessing | [Project and citation guidance](https://scikit-learn.org/stable/about.html). |
| NumPy, SciPy, pandas | Numerical arrays, scientific routines, and tables | [NumPy](https://numpy.org/citing-numpy/), [SciPy](https://scipy.org/citing-scipy/), [pandas](https://pandas.pydata.org/about/citing.html). |

| Nextflow | Workflow orchestration | Di Tommaso et al. (2017), *Nextflow enables reproducible computational workflows*, Nature Biotechnology 35:316–319. [DOI](https://doi.org/10.1038/nbt.3820); [project](https://github.com/nextflow-io/nextflow). |
| GSW / TEOS-10 | Seawater density and physical calculations | [GSW-Python](https://github.com/TEOS-10/GSW-Python); [TEOS-10](https://www.teos-10.org/). |
| Matplotlib / seaborn / NetworkX | Figures and graph analysis | [Matplotlib](https://matplotlib.org/), [seaborn](https://seaborn.pydata.org/), [NetworkX](https://networkx.org/). |

Mamba manages software environments; see the [Mamba project](https://github.com/mamba-org/mamba). Record environment specifications and actual versions alongside scientific citations.
