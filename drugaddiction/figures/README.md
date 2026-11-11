# Manuscript figure rendering

Use Arial, a Cairo-enabled R build, `ggplot2`, `ragg` and `svglite`. The scripts write PDF, SVG and 600 dpi PNG to the specified output directory. Upstream marker and subtype review scripts are listed in the [upstream documentation](../upstream/README.md).

| Script | Panels | Inputs |
| --- | --- | --- |
| `render_major_markers.R` | Fig. 3B, C and Fig. S7D | Accepted six-class Seurat object |
| `render_composition.R` | Fig. 3D, K and Fig. S9A | Six-class global coordinate/annotation CSV and neuronal metadata CSV |
| `render_proportions.R` | Fig. 3L and Fig. S9B | Neuronal metadata CSV, finalized excitatory statistics CSV and finalized inhibitory statistics CSV |
| `render_module_heatmap.R` | Fig. 3I | Excitatory IEG-high payload RDS |
| `render_cellchat.R` | Fig. 3P and Fig. S10A–E | Joint CellChat summary tables under `stage3/06_CellChat/manuscript_tables` |

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

Each proportion statistics CSV must contain `subtype`, `P`, `BH_adjusted_P` and `Label`, one row per subtype. Use the metadata subtype IDs with underscores. Leave all three result fields empty for an unavailable test. Brackets use the supplied adjusted values and labels; the script verifies consistency with the P values and does not fit a model. A missing or non-significant test is not given a significance bracket. Pool percentages with a zero denominator are omitted.

The major-class renderer uses the accepted coordinates, labels and log-normalized expression, without fitting another embedding. It uses Neurod6 in the marker dot plot and Aqp4 in the marker feature grid. The composition stacks are descriptive, with all nuclei pooled within a condition for display. The recruitment panels instead show each observed pool separately. The CellChat renderer requires the original CellChat plotting API; see the [signaling documentation](../signaling/README.md) for computation and receiver-scope definitions.
