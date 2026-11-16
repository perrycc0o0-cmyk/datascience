# Subtype proportions and CellChat

Pool-level binomial mixed models and expression-based neuronal communication analysis.

## Inputs and packages

Use `neuron_metadata.csv` for proportion models and the excitatory and inhibitory `IEGhigh_payload.rds` files for CellChat. The model packages include lme4, jsonlite and digest. CellChat uses version 1.6.1. Orthology preparation requires httr2, xml2, jsonlite and network access to the Ensembl 116 archive.

## Usage

Open the required script, edit the input paths, `work_dir`, `result_dir` and `script_path` where present, and check its sample table. Run the complete R file using **Source** in RStudio.

| Script | Purpose |
| --- | --- |
| `fit_binomial_models.R` | Fit subtype proportion models and save statistics under `04_proportions_GLMM/` |
| `prepare_orthology.R` | Obtain rat–human one-to-one orthologs |
| `fit_cellchat.R` | Fit abundance-weighted neuronal communication models |
| `summarize_joint_cellchat.R` | Summarize Ret-Ext versus Ext outgoing strength |

For GLMM calculations, set `bootstrap_workers` near the top of the script. The model uses condition as a fixed effect and pooled sample as a random intercept.

For CellChat, run orthology preparation, model fitting and summary generation in that order. Outputs are saved under `06_CellChat/`. Use the resulting `manuscript_tables/` with `figures/render_cellchat.R`.
