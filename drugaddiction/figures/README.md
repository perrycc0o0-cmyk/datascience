# Manuscript figures

UMAP, marker expression, cell composition, subtype proportion, module heatmap and CellChat figures.

## Requirements

Use Arial, a Cairo-enabled R build, `ggplot2`, `ragg` and `svglite`. Run examples from `drugaddiction/` and replace the paths. Use `--help` to view command arguments. The upstream marker and subtype plots are described in [upstream/README.md](../upstream/README.md).

## Scripts and inputs

| Script | Panels | Inputs |
| --- | --- | --- |
| `render_major_markers.R` | Fig. 3B, C and Fig. S7D | Accepted six-class Seurat object |
| `render_composition.R` | Fig. 3D, K and Fig. S9A | Global coordinate/annotation CSV and neuronal metadata CSV |
| `render_proportions.R` | Fig. 3L and Fig. S9B | Neuronal metadata and finalized excitatory/inhibitory statistics CSVs |
| `render_module_heatmap.R` | Fig. 3I | Excitatory IEG-high payload |
| `render_cellchat.R` | Fig. 3P and Fig. S10A–E | Joint CellChat summary tables |

## Commands

```bash
Rscript figures/render_major_markers.R \
  /path/to/work/stage1/05_six_class_review/global_celltype_six_class.rds /path/to/figure-output

Rscript figures/render_composition.R \
  /path/to/work/stage1/05_six_class_review/source_data/global_cell_annotations_UMAP.csv.gz \
  /path/to/work/stage3/01_inputs/neuron_metadata.csv /path/to/figure-output

Rscript figures/render_proportions.R /path/to/work/stage3/01_inputs/neuron_metadata.csv \
  /path/to/excitatory_statistics.csv /path/to/inhibitory_statistics.csv /path/to/figure-output

Rscript figures/render_module_heatmap.R \
  /path/to/work/stage3/01_inputs/Excitatory/IEGhigh_payload.rds /path/to/figure-output

Rscript figures/render_cellchat.R /path/to/work/stage3 /path/to/figure-output
```

Use a separate output directory for each renderer. Figures are exported as PDF, SVG and 600-dpi PNG files.

Each proportion statistics CSV must have `subtype`, `P`, `BH_adjusted_P` and `Label`, with one row per subtype. Use the metadata subtype IDs with underscores. For an unavailable test, leave the three result fields empty. The renderer checks the supplied values and labels before plotting.

Expression figure commands are listed in [expression/README.md](../expression/README.md). Prepare CellChat summary tables using the [signaling instructions](../signaling/README.md).
