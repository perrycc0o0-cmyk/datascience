# Rat infralimbic cortex snRNA-seq analysis

Analysis code and figure-generation scripts for comparing withdrawal alone (WD), extinction (Ext), and retrieval-extinction (Ret-Ext) in the rat infralimbic cortex. The study includes three pooled biological samples per condition, with tissue from five rats per pool.

## Contents

| Directory | Contents |
| --- | --- |
| [upstream/](upstream/README.md) | Library QC, SoupX correction, doublet removal, Harmony integration, cell annotation, neuronal subclustering and IEG-high classification |
| [expression/](expression/README.md) | DESeq2 pseudobulk differential expression, Monocle 2 trajectories, GSEA and gene–category networks |
| [signaling/](signaling/README.md) | Subtype proportion models and CellChat analysis |
| [figures/](figures/README.md) | UMAP, marker, composition, proportion, heatmap and communication plots |

## Requirements

Use R and the package versions listed in the module documentation. The main packages are Seurat 5.3.1, SeuratObject 5.2.0, SoupX 1.6.2, scDblFinder 1.24.10, Harmony 1.2.4, DESeq2 1.50.2, Monocle 2.38.0, DDRTree 0.1.5, fgsea 1.36.2 and CellChat 1.6.1. The gene–category network builder uses Python 3. Figure export requires Arial and the graphics packages specified in each module.

## Usage

1. Prepare the raw and filtered 10x HDF5 matrices for each library as described in [upstream/README.md](upstream/README.md).
2. Run library processing, broad-cell annotation and marker review.
3. Review the neuronal subtype annotations and prepare the IEG-high analysis inputs.
4. Use the expression, signaling and figure commands with the resulting input files.

Run commands from `drugaddiction/`. Replace example paths with your local input and output paths. Keep analysis outputs outside the code directory and use a new output directory for a new run. Use `--help` to view a script's arguments.

```bash
Rscript upstream/stage1/01_process_libraries.R --root /path/to/work/stage1 --input-dir /path/to/10x_libraries
Rscript upstream/stage1/02_global_celltype.R --root /path/to/work/stage1
```

Continue with the annotation and input-preparation commands in the upstream instructions. The resulting `stage3/01_inputs/` directory supplies the downstream analyses:

```bash
Rscript expression/pseudobulk_deseq2.R /path/to/work/stage3/01_inputs/Excitatory/IEGhigh_payload.rds /path/to/work/stage3/01_inputs/Inhibitory/IEGhigh_payload.rds /path/to/work/deseq2
Rscript expression/plot_deseq2_volcanoes.R /path/to/work/deseq2 /path/to/work/volcanoes
```

The module instructions list the required reference files, output tables and figure commands. Sequencing matrices, analysis objects and generated results are not included.

## License

The original analysis code and accompanying documentation in `drugaddiction/` are released under the [MIT License](LICENSE). Third-party software and external data or reference resources remain subject to their own terms.
