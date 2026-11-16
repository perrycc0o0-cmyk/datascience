# Single-nucleus RNA-seq processing

Library QC, ambient RNA correction, doublet removal, cell annotation, neuronal subclustering and IEG-high input preparation.

## Inputs and packages

Place each sample's `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5` in a sample folder under `input_dir`. Use the source folder names `Con1`–`Con3`, `Ext1`–`Ext3` and `Ret_ext1`–`Ret_ext3`. The sample table maps these to the WD, Ext and Ret-Ext labels.

The processing environment used R 4.5.2 on Linux, Seurat 5.3.1, SeuratObject 5.2.0, SoupX 1.6.2, scDblFinder 1.24.10 and Harmony 1.2.4. Plotting requires Arial, Cairo and Fontconfig. Each script lists its packages at the beginning.

## Usage

Edit `input_dir`, `work_dir`, `result_dir` and `script_path` in the required file. Set `script_path` to that file's local path. Check the sample table and inline parameter lists, then run the complete script using **Source** in RStudio. Use the same sample definitions and processing settings across stages.

| Order | Script | Purpose |
| --- | --- | --- |
| 1 | `stage1/00_runtime_preflight.R` | Check packages, input files and fonts |
| 2 | `stage1/01_process_libraries.R` | Correct ambient RNA, apply QC and remove doublets |
| 3 | `stage1/02_global_celltype.R` | Integrate libraries and cluster broad cell classes |
| 4 | `stage1/03_plot_global_review.R` | Plot broad-cell marker evidence |
| 5 | `stage1/06_match_six_major_classes.R` | Apply six-class annotation |
| 6 | `stage1/07_plot_six_class_review.R` | Plot six-class results |
| 7 | `stage2/00_runtime_preflight.R` | Check neuronal-analysis inputs |
| 8 | `stage2/01_recluster_neurons.R` | Recluster excitatory and inhibitory nuclei |
| 9 | `stage2/02_annotate_neuron_subtypes.R` | Annotate neuronal subtypes |
| 10 | `stage2/03_plot_neuron_review.R` | Plot neuronal subtype markers |
| 11 | `stage3/00_apply_accepted_subtype_override.R` | Apply the accepted inhibitory subtype mapping |
| 12 | `stage3/01_prepare_analysis_inputs.R` | Prepare IEG-high counts and metadata |

Review the broad-cell markers before stage 2 and neuronal subtype markers before stage 3. The annotation rules and expected counts are specific to this cohort. To resume interrupted library processing, set `resume <- TRUE` with the same paths, script and settings. Set `selected_samples` to a character vector to process selected libraries, or leave it `NULL` for all libraries.

## Outputs

Stage 1 saves corrected matrices, QC summaries, singlet objects and broad-cell plots. Stage 2 saves neuronal objects, subtype annotations and marker plots. Stage 3 writes the following files under `01_inputs/`:

- `neuron_metadata.csv`
- `Excitatory/IEGhigh_payload.rds`
- `Inhibitory/IEGhigh_payload.rds`
- `pseudobulk_inputs.rds`
- Class-specific pseudobulk count and sample metadata CSVs

These files are the inputs for the expression, signaling and figure scripts. Parameter snapshots are saved with the results as JSON.
