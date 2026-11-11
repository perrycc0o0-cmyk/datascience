# Pooled-sample recruitment models and CellChat

Run the commands from the project directory. Replace `/path/to/work` with the private output directory passed to the upstream workflow. The downstream run root is its `stage3` directory.

## Subtype recruitment

```bash
Rscript signaling/fit_binomial_models.R /path/to/work/stage3/01_inputs/neuron_metadata.csv /path/to/work/stage3
```

The model uses the IEG-high and total nuclear counts for each pool and subtype. It fits a binomial logit model with condition as a fixed effect and an unrestricted pooled-sample random intercept. The primary condition test is a null-parametric-bootstrap likelihood-ratio test with 1,999 simulations and 15-point adaptive Gauss–Hermite quadrature. The code reports Benjamini–Hochberg adjusted values separately for the excitatory and inhibitory subtype families and retains unavailable results as missing. The R dependencies are `lme4`, `jsonlite` and `digest`. Set `GLMM_WORKERS` to control the number of bootstrap workers on Unix-like systems.

For the excitatory composition contrast, each subtype's numerator is the number of IEG-high nuclei of that subtype and its denominator is the full IEG-high excitatory count in the corresponding pool. For inhibitory recruitment the denominator is the full within-subtype count in that pool. The model script is separate from `figures/render_proportions.R`, which renders the finalized statistics supplied to it as input. The plotting script does not fit or alter the model or treat nuclei as independent biological samples.

## Inferred neuronal signaling

```bash
Rscript signaling/prepare_orthology.R --root /path/to/work/stage3
Rscript signaling/fit_cellchat.R --root /path/to/work/stage3
Rscript signaling/summarize_joint_cellchat.R /path/to/work/stage3
Rscript figures/render_cellchat.R /path/to/work/stage3 /path/to/figure-output
```

- `prepare_orthology.R` queries the Ensembl 116 archive for verified rat–human one-to-one orthologs. It requires network access to Ensembl and the `httr2`, `xml2` and `jsonlite` R packages.
- `fit_cellchat.R` loads the upstream normalized IEG-high neuronal payloads in `01_inputs/Excitatory` and `01_inputs/Inhibitory`. It uses CellChat 1.6.1 and its human ligand–receptor database after orthology mapping, triMean expression, population-size weighting, the native `min.cells = 10` filter, 100 permutations and seed 42. The fit writes joint and separate lineage models; the reported comparisons use the **joint** models only.
- `summarize_joint_cellchat.R` compares Ret-Ext and Ext from the joint model. The sender is excluded from its receiver set. An excitatory receiver scope is a subset of the joint model, not a separate excitatory-only fit. The all-other-subtypes scope includes eligible inhibitory receivers. Groups not eligible for native CellChat inference are missing rather than zero. A tested and filtered interaction can contribute zero; an untested pathway remains missing.
- `render_cellchat.R` renders the selected total-strength comparisons (Fig. 3P) and the supplementary panels (Fig. S10A–E). The positive and negative pathway panels rank signed differences from the 13 candidate pathways; they are descriptive. The heatmap displays the 13 selected pathways for the same five comparisons as Fig. 3P, not an additive decomposition of total outgoing strength.

CellChat results are expression-based predictions. Native ligand–receptor P values filter edges within a fitted network; they are not a test of the condition contrast across independent pools.
