# Optional media and genome modeling

The environmental workflow can run without genomes. Enable this branch after validating the environmental matrix and group assignments.

## Configure inputs

Copy `examples/genomes_manifest.tsv.example`, supply your own genome FASTAs, and set `genome_modeling.genomes_manifest`. Enable `gapseq_media.enabled` and `genome_modeling.enabled` in the complete run configuration. The local genome manifest referenced by the Saanich Inlet example is not distributed.

Use the field-by-field [configuration reference](CONFIGURATION.md) for basal media, genome metadata, DNA/RNA abundance tables, normalization, solver, flux thresholds, FVA, and counterfactual settings. Keep paths accessible to every execution node. Reconstructed models require gapseq's reference resources as well as the software environment; complete their setup before running on offline nodes.

## What is compared

BASINS layers group-specific measured chemistry onto one common basal medium and records provenance, exclusions, and recipe audits. It reconstructs and gap-fills a model for each genome with gapseq 2.1.0, then applies supported media to clean copies of that fixed model. This avoids confusing changes in model reconstruction with changes in the tested medium. DIAMOND is the default aligner. MMseqs2 requires the explicit external executable setting described in the configuration reference.

COBRApy-based comparisons produce growth predictions, exchange flux variability, pairwise differences, and nutrient counterfactuals. Combined tables summarize results across genomes. These are conditional model predictions, not measurements of realized growth or proof that a genome occupies a niche. Consult the medium recipe and reconstruction QC when a model does not grow.

A missing nutrient-importance value means the nutrient was not supplied/tested, not that its importance is zero. [Output interpretation](outputs.md) defines the reported scores. Cite [gapseq and COBRApy](citations.md) when using this branch.
