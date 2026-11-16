# Manuscript figures

UMAP, marker expression, composition, subtype proportion, module heatmap and CellChat plots.

## Usage

Install the packages listed at the top of the required script. Figure export uses Arial, Cairo, ragg and svglite. Edit the input paths, `figure_dir` and sample table, then run the complete file using **Source** in RStudio. Choose a separate output directory for each renderer.

| Script | Figures | Input |
| --- | --- | --- |
| `render_major_markers.R` | Fig. 3B/C and Fig. S7D | Six-class Seurat object |
| `render_composition.R` | Fig. 3D/K and Fig. S9A | Global annotations and neuronal metadata |
| `render_proportions.R` | Fig. 3L and Fig. S9B | Neuronal metadata and finalized statistics |
| `render_module_heatmap.R` | Fig. 3I | Excitatory IEG-high payload |
| `render_cellchat.R` | Fig. 3P and Fig. S10A–E | Joint CellChat summary tables |

Proportion statistics CSVs must contain `subtype`, `P`, `BH_adjusted_P` and `Label`, with one row per subtype. Use subtype IDs matching the nuclear metadata. Leave the result fields empty for unavailable tests.

The heatmap gene list and color palette are defined near the top of `render_module_heatmap.R`. Figures are saved as PDF, SVG and 600-dpi PNG.

Upstream marker plots are described in [upstream/README.md](../upstream/README.md). Volcano, trajectory, GSEA and network figures are described in [expression/README.md](../expression/README.md).
