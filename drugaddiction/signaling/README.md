# Subtype proportions and CellChat

Pool-level subtype proportion models and expression-based neuronal communication analysis.

## Inputs and software

Use `stage3/01_inputs/neuron_metadata.csv` for proportion models and the excitatory and inhibitory `IEGhigh_payload.rds` files for CellChat. The model packages are `lme4`, `jsonlite` and `digest`. Orthology preparation requires `httr2`, `xml2`, `jsonlite` and network access to the Ensembl 116 archive. CellChat inference uses CellChat 1.6.1 and its human ligand–receptor database.

Run examples from `drugaddiction/` and replace the paths. The `stage3` directory is the run root for these commands. Use `--help` for command arguments.

## Subtype proportion models

```bash
Rscript signaling/fit_binomial_models.R /path/to/work/stage3/01_inputs/neuron_metadata.csv /path/to/work/stage3
```

This script fits binomial mixed models using condition as a fixed effect and pooled sample as a random intercept. Results are written under `04_proportions_GLMM/` in the run root. Set `GLMM_WORKERS` to choose the number of bootstrap workers on Unix-like systems.

Use the proportion figure command in [figures/README.md](../figures/README.md) with the finalized excitatory and inhibitory statistics tables.

## CellChat

```bash
Rscript signaling/prepare_orthology.R --root /path/to/work/stage3
Rscript signaling/fit_cellchat.R --root /path/to/work/stage3
Rscript signaling/summarize_joint_cellchat.R /path/to/work/stage3
Rscript figures/render_cellchat.R /path/to/work/stage3 /path/to/figure-output
```

| Script | Purpose |
| --- | --- |
| `prepare_orthology.R` | Obtain verified rat–human one-to-one orthologs from Ensembl 116 |
| `fit_cellchat.R` | Fit abundance-weighted neuronal communication models |
| `summarize_joint_cellchat.R` | Summarize Ret-Ext versus Ext outgoing strength from the joint model |
| `render_cellchat.R` | Draw Fig. 3P and Fig. S10A–E from the saved summary tables |

CellChat files are saved under `06_CellChat/`. The figure renderer reads `06_CellChat/manuscript_tables/`. Receiver scopes select other excitatory subtypes or all other included neuronal subtypes from the same joint model.
