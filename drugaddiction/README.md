# Rat infralimbic cortex snRNA-seq analysis

Scripts for comparing withdrawal alone (WD), extinction (Ext) and retrieval-extinction (Ret-Ext) in the rat infralimbic cortex. The dataset contains three pooled samples per condition, with tissue from five rats per pool.

## Contents

| Directory | Analysis |
| --- | --- |
| [upstream/](upstream/README.md) | QC, SoupX correction, doublet removal, Harmony integration, cell annotation and IEG-high classification |
| [expression/](expression/README.md) | DESeq2 pseudobulk analysis, Monocle 2 trajectories, GSEA and gene–category networks |
| [signaling/](signaling/README.md) | Subtype proportion models and CellChat |
| [figures/](figures/README.md) | Marker, UMAP, composition, proportion, heatmap and communication plots |

## Usage

1. Install the packages listed at the beginning of the required script.
2. Open the script and edit its paths, including `script_path` where present.
3. Check `sample_cols`, `metadata` and the parameters near the top of the file.
4. Run the complete R script in RStudio using **Source**. Run the Python network script in a Python editor.
5. Read the saved tables and figures from the specified output directories.

Start with library processing and cell annotation, then prepare the neuronal inputs for the downstream analyses. The module instructions list the script order and input files. Keep results outside the code directory and use a new output directory for a new analysis. Sequencing matrices and external reference datasets are supplied separately.

The original code and documentation are available under the [MIT License](LICENSE).
