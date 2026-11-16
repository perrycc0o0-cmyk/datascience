# Single-nucleus RNA-seq processing

This module runs the study's upstream rat infralimbic cortex workflow from nine Cell Ranger libraries. It uses the fixed quality-control settings, per-library ambient RNA correction and doublet detection, sample-only Harmony integration, and neuronal annotations used for the current analysis. The scripts are specific to this cohort; cluster numbers and count checks should not be carried over to unrelated datasets.

## Inputs and runtime

Supply the original filtered and raw HDF5 matrices in separate sample directories beneath one input directory. Use the source sample names in [`config/sample_naming.json`](config/sample_naming.json): `Con1`–`Con3`, `Ext1`–`Ext3`, and `Ret_ext1`–`Ret_ext3`. Every directory must contain `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5`. The workflow maps the source names to `WD1`–`WD3`, `Ext1`–`Ext3`, and `Ret-Ext1`–`Ret-Ext3` without modifying the matrices. Each library represents one biological pool of five rats; nuclei within a library are not independent biological replicates.

Use Python 3 and R with the package versions in [`config/runtime_versions.json`](config/runtime_versions.json). The recorded analysis environment was R 4.5.2 on Linux. Install packages in a local R library and expose it through `R_LIBS_USER` if necessary. The plot scripts require Cairo and Fontconfig (`fc-match`), with Arial available locally. An optional `--fontconfig-file` may be supplied to the `preprocess` and `neurons` steps. Matrices, analysis objects, per-nucleus tables, and font files are not distributed in this repository.

## Run the stages

Each R script runs directly, without a Python launcher or a prepared R session. Keep the source tree intact; `script_io.R` provides path and configuration utilities and never invokes an analysis. Use a private work directory outside the code folder. `--root` identifies the stage's working directory, not the source directory. Each command creates only its own outputs; it does not automatically execute the next stage. All executables accept `--help`.

Run the following commands individually from `drugaddiction/`. Absolute script paths also work from another current directory.

```bash
# Library processing and broad-cell annotation
Rscript upstream/stage1/00_runtime_preflight.R --root /path/to/work/stage1 --input-dir /path/to/matrices
Rscript upstream/stage1/01_process_libraries.R --root /path/to/work/stage1 --input-dir /path/to/matrices
Rscript upstream/stage1/02_global_celltype.R --root /path/to/work/stage1
Rscript upstream/stage1/03_plot_global_review.R --root /path/to/work/stage1
Rscript upstream/stage1/06_match_six_major_classes.R --root /path/to/work/stage1
Rscript upstream/stage1/07_plot_six_class_review.R --root /path/to/work/stage1

# After reviewing the broad-cell marker results
Rscript upstream/stage2/00_runtime_preflight.R --root /path/to/work/stage2 --stage1 /path/to/work/stage1
Rscript upstream/stage2/01_recluster_neurons.R --root /path/to/work/stage2 --stage1 /path/to/work/stage1
Rscript upstream/stage2/02_annotate_neuron_subtypes.R --root /path/to/work/stage2
Rscript upstream/stage2/03_plot_neuron_review.R --root /path/to/work/stage2

# After reviewing the neuronal subtype markers
Rscript upstream/stage3/00_apply_accepted_subtype_override.R /path/to/work/stage2
Rscript upstream/stage3/01_prepare_analysis_inputs.R --root /path/to/work/stage3 --stage2 /path/to/work/stage2
```

The first command for a stage initializes its configuration from the unchanged templates in `config/` and records local paths in the private stage directory. No scripts are copied into that work directory. Later scripts reuse that saved configuration. Use `--config /path/to/config.json` for an explicit configuration. Library processing requires `--input-dir` only when the saved configuration has no input path. Neuronal steps can infer a sibling `stage1` directory; input preparation can infer a sibling `stage2` directory. Use the explicit `--stage1` and `--stage2` options shown above for other layouts. If Arial needs a local Fontconfig file, pass the same `--fontconfig-file /path/to/fonts.conf` when initializing the processing or neuronal stage. Existing configurations with different paths or values are rejected rather than silently overwritten.

Review the global marker panels before neuronal reclustering and the neuronal subtype panels before input preparation. Resume failed per-library processing by calling `01_process_libraries.R` with the same arguments plus `--resume`; use the same input matrices, configuration and scripts. Completed analysis outputs are not overwritten. The original marker-review renderers are included here; presentation-only layout revisions are separate from the processing and annotation stages. The direct entry points do not change QC thresholds, correction methods, annotation rules, seeds, analysis order, or plotting calculations.

| Standalone script | Required saved input | Output location within its stage directory |
| --- | --- | --- |
| `stage1/00_runtime_preflight.R` | Original raw and filtered matrices | `qa/` environment checks |
| `stage1/01_process_libraries.R` | Original raw and filtered matrices | `01_ambient/`, `02_qc_doublet/` |
| `stage1/02_global_celltype.R` | Per-library singlet objects | `03_global/` |
| `stage1/03_plot_global_review.R` | Saved global marker payload | `04_celltype_review/` |
| `stage1/06_match_six_major_classes.R` | Saved global object and marker payload | `05_six_class_review/` annotations |
| `stage1/07_plot_six_class_review.R` | Six-class annotations and global payload | `05_six_class_review/` figures |
| `stage2/00_runtime_preflight.R` | Accepted six-class object | `qa/` environment checks |
| `stage2/01_recluster_neurons.R` | Accepted six-class object | `01_Exc/`, `02_Inh/` |
| `stage2/02_annotate_neuron_subtypes.R` | Reclustering objects and completion record | `03_annotation/` |
| `stage2/03_plot_neuron_review.R` | Saved subtype marker payloads | `04_marker_review/` |
| `stage3/00_apply_accepted_subtype_override.R` | Inhibitory annotation CSV | Mapping in the supplied stage-2 directory |
| `stage3/01_prepare_analysis_inputs.R` | Neuronal objects, annotations and accepted mapping | `01_inputs/` |

### Stage 1: Library processing and broad classes

`stage1/01_process_libraries.R` reads matched raw and filtered HDF5 matrices, corrects ambient RNA using SoupX, and preserves both corrected integer counts and observed counts for retained nuclei. Initial gene and nucleus eligibility uses Seurat's `min.cells = 3` and `min.features = 300`. The fixed filters use observed counts: more than 700 and fewer than 5,000 genes, more than 500 and fewer than 30,000 UMIs, mitochondrial percentage below 5, and `log10(genes) / log10(UMIs)` above 0.75. Ribosomal percentage and corrected-count metrics are retained for review, not used to rescue or exclude nuclei. scDblFinder runs on the corrected counts of eligible nuclei within each library. The doublet-rate prior is `0.004 * (original called barcodes / 1000)`. The study retained 131,480 singlets after the filters and doublet removal. The script writes per-library checkpoints, the called-barcode QC audit, and an attrition summary.

`stage1/02_global_celltype.R` merges the singlets, normalizes counts and selects 2,000 variable genes. It computes PCA, selects PCs using the study's variance rule, and runs Harmony with `sample_id` as the only integration variable. Clustering uses Louvain resolution 1.0 for the primary result. The global UMAP and marker-based broad-class proposals are written without importing old cluster labels.

`stage1/06_match_six_major_classes.R` assigns initially mixed or unresolved clusters to their nearest of six broad-class reference profiles using library-balanced correlations. Original labels and uncertainty flags are preserved. This is a constrained assignment, not independent biological confirmation of the uncertain clusters. The analysis object for the next stage is `stage1/05_six_class_review/global_celltype_six_class.rds`.

### Stage 2: Neuronal subtype annotation

`stage2/01_recluster_neurons.R` separately reclusters excitatory and inhibitory nuclei from the six-class object, using the existing corrected counts without another QC filter. The neuronal configurations retain 2,000 variable genes, sample-only Harmony, and resolution 0.5. `stage2/02_annotate_neuron_subtypes.R` uses the marker rules in [`config/subtype_markers_and_style.json`](config/subtype_markers_and_style.json). It keeps original cluster and lineage-review fields. Annotation and reclustering do not perform a between-condition expression test. The resulting objects are `stage2/01_Exc/neuron_reclustered.rds` and `stage2/02_Inh/neuron_reclustered.rds`; marker review tables are generated under `stage2/03_annotation`.

### Stage 3: IEG-high definition and downstream inputs

`stage3/00_apply_accepted_subtype_override.R` applies the accepted reassignment of inhibitory subcluster 6 from `GABA_Lamp5` to `GABA_Etv1`. It checks the original label and the cohort's expected counts before generating the mapping and leaves other clusters unchanged. No per-nucleus mapping is distributed.

`stage3/01_prepare_analysis_inputs.R` aligns the accepted annotations to the corrected RNA layer. It sums the log-normalized expression of `Arc`, `Bdnf`, `Btg2`, `Fos`, `Fosl2`, `Homer1`, `Npas4`, and `Nr4a1`, defining IEG-high at or above the 95th percentile within each neuronal class across all conditions. It aggregates corrected integer counts by pool, neuronal class and IEG group, preserving observed counts separately for audit. The downstream input directory is `stage3`; its main files are:

- `01_inputs/pseudobulk_inputs.rds`: corrected count matrices and metadata for all, IEG-high and IEG-low populations.
- `01_inputs/neuron_metadata.csv`: one row per annotated nucleus with pool, condition, subtype, IEG-high status and preserved annotation flags.
- `01_inputs/Excitatory/IEGhigh_payload.rds` and `01_inputs/Inhibitory/IEGhigh_payload.rds`: corrected counts, normalized expression and aligned metadata for each IEG-high population.
- `01_inputs/<class>/IEGhigh_corrected_pseudobulk_counts.csv` and `01_inputs/<class>/IEGhigh_pool_metadata.csv`: sample-level counts and pool metadata for the corresponding class.

Keep generated outputs outside version control. Package checks, completion markers, parameter records and input audits are written to the private run directory.
