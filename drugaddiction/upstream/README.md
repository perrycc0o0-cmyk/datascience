# Single-nucleus RNA-seq processing

Library processing, broad-cell annotation, neuronal subclustering and IEG-high input preparation for the rat infralimbic cortex dataset.

## Inputs and software

Place each library's `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5` in its own directory under the matrix input directory. Use the source names in [config/sample_naming.json](config/sample_naming.json): `Con1`–`Con3`, `Ext1`–`Ext3` and `Ret_ext1`–`Ret_ext3`. The scripts map these to `WD1`–`WD3`, `Ext1`–`Ext3` and `Ret-Ext1`–`Ret-Ext3`.

Use the R packages listed in [config/runtime_versions.json](config/runtime_versions.json). The recorded processing environment used R 4.5.2 on Linux. Set `R_LIBS_USER` for a local package library if needed. Plotting requires Cairo, Fontconfig (`fc-match`) and Arial.

## Commands

Run the commands from `drugaddiction/`, replacing the example paths. Use a work directory outside the code folder. `--root` specifies the output directory for that stage.

```bash
# Library processing and broad-cell annotation
Rscript upstream/stage1/00_runtime_preflight.R --root /path/to/work/stage1 --input-dir /path/to/matrices
Rscript upstream/stage1/01_process_libraries.R --root /path/to/work/stage1 --input-dir /path/to/matrices
Rscript upstream/stage1/02_global_celltype.R --root /path/to/work/stage1
Rscript upstream/stage1/03_plot_global_review.R --root /path/to/work/stage1
Rscript upstream/stage1/06_match_six_major_classes.R --root /path/to/work/stage1
Rscript upstream/stage1/07_plot_six_class_review.R --root /path/to/work/stage1

# After reviewing the broad-cell markers
Rscript upstream/stage2/00_runtime_preflight.R --root /path/to/work/stage2 --stage1 /path/to/work/stage1
Rscript upstream/stage2/01_recluster_neurons.R --root /path/to/work/stage2 --stage1 /path/to/work/stage1
Rscript upstream/stage2/02_annotate_neuron_subtypes.R --root /path/to/work/stage2
Rscript upstream/stage2/03_plot_neuron_review.R --root /path/to/work/stage2

# After reviewing the neuronal subtype markers
Rscript upstream/stage3/00_apply_accepted_subtype_override.R /path/to/work/stage2
Rscript upstream/stage3/01_prepare_analysis_inputs.R --root /path/to/work/stage3 --stage2 /path/to/work/stage2
```

The initial command for each stage records its configuration under that stage's `config/` directory. Use `--config /path/to/config.json` to supply a configuration file. For a local Arial configuration, add `--fontconfig-file /path/to/fonts.conf` when initializing stage 1 or stage 2. Keep the source directory structure intact.

To resume interrupted library processing, repeat the processing command with the same inputs and settings plus `--resume`. Existing outputs and differing saved configurations are protected from overwrite. Use `--help` for command arguments.

## Script inputs and outputs

| Script | Input | Output within the stage directory |
| --- | --- | --- |
| `stage1/00_runtime_preflight.R` | Raw and filtered matrices | `qa/` environment checks |
| `stage1/01_process_libraries.R` | Raw and filtered matrices | `01_ambient/`, `02_qc_doublet/` |
| `stage1/02_global_celltype.R` | Per-library singlet objects | `03_global/` |
| `stage1/03_plot_global_review.R` | Global marker payload | `04_celltype_review/` |
| `stage1/06_match_six_major_classes.R` | Global object and marker payload | `05_six_class_review/` annotations |
| `stage1/07_plot_six_class_review.R` | Six-class annotations and global payload | `05_six_class_review/` figures |
| `stage2/00_runtime_preflight.R` | Accepted six-class object | `qa/` environment checks |
| `stage2/01_recluster_neurons.R` | Accepted six-class object | `01_Exc/`, `02_Inh/` |
| `stage2/02_annotate_neuron_subtypes.R` | Reclustering objects | `03_annotation/` |
| `stage2/03_plot_neuron_review.R` | Subtype marker payloads | `04_marker_review/` |
| `stage3/00_apply_accepted_subtype_override.R` | Inhibitory annotation CSV | Mapping under the supplied stage-2 directory |
| `stage3/01_prepare_analysis_inputs.R` | Neuronal objects, annotations and mapping | `01_inputs/` |

The annotation rules and expected counts are specific to this cohort. The accepted inhibitory mapping reassigns subcluster 6 from `GABA_Lamp5` to `GABA_Etv1`.

## Downstream input files

Input preparation writes the following files under `stage3/01_inputs/`:

- `neuron_metadata.csv`: nuclear annotations, pool and condition labels, subtypes and IEG-high status.
- `Excitatory/IEGhigh_payload.rds` and `Inhibitory/IEGhigh_payload.rds`: corrected counts, normalized expression, feature annotations and nuclear metadata.
- `pseudobulk_inputs.rds`: aggregated counts and pool metadata.
- `<class>/IEGhigh_corrected_pseudobulk_counts.csv` and `<class>/IEGhigh_pool_metadata.csv`: class-specific count tables and sample information.

Use these files with the commands in the [expression](../expression/README.md), [signaling](../signaling/README.md) and [figure](../figures/README.md) instructions.
