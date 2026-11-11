# Single-nucleus RNA-seq processing

This module runs the study's upstream rat infralimbic cortex workflow from nine Cell Ranger libraries. It uses the fixed quality-control settings, per-library ambient RNA correction and doublet detection, sample-only Harmony integration, and neuronal annotations used for the current analysis. The scripts are specific to this cohort; cluster numbers and count checks should not be carried over to unrelated datasets.

## Inputs and runtime

Supply the original filtered and raw HDF5 matrices in separate sample directories beneath one input directory. Use the source sample names in [`config/sample_naming.json`](config/sample_naming.json): `Con1`–`Con3`, `Ext1`–`Ext3`, and `Ret_ext1`–`Ret_ext3`. Every directory must contain `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5`. The workflow maps the source names to `WD1`–`WD3`, `Ext1`–`Ext3`, and `Ret-Ext1`–`Ret-Ext3` without modifying the matrices. Each library represents one biological pool of five rats; nuclei within a library are not independent biological replicates.

Use Python 3 and R with the package versions in [`config/runtime_versions.json`](config/runtime_versions.json). The recorded analysis environment was R 4.5.2 on Linux. Install packages in a local R library and expose it through `R_LIBS_USER` if necessary. The plot scripts require Cairo and Fontconfig (`fc-match`), with Arial available locally. An optional `--fontconfig-file` may be supplied to the `preprocess` and `neurons` steps. Matrices, analysis objects, per-nucleus tables, and font files are not distributed in this repository.

## Run the stages

Use a private output directory outside the code folder. Each command writes the native stage-specific files below that directory.

```bash
python3 upstream/run.py preprocess --input-dir /path/to/matrices --run-dir /path/to/work
python3 upstream/run.py global --run-dir /path/to/work
python3 upstream/run.py neurons --run-dir /path/to/work
python3 upstream/run.py prepare --run-dir /path/to/work
```

Review the global marker panels before proceeding to `neurons`, and the neuronal subtype panels before `prepare`. A failed per-library run can be resumed with the same `preprocess` command plus `--resume`; use the same input matrices, configuration and scripts. Completed analysis outputs are not overwritten. The original marker-review renderers are included here; presentation-only layout revisions are separate from the processing and annotation stages.

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
