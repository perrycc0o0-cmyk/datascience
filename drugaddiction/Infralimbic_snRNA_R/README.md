# Rat infralimbic cortex snRNA-seq

The notebooks cover library processing and cell annotation, followed by analyses of IEG-high neuronal populations. Tables and figures are written to the output paths defined in each notebook.

## Usage

### Library processing and annotation

Set `input_dir` to the directory containing one folder per library. Each library folder must contain `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5`. The source folder names are `Con1`–`Con3`, `Ext1`–`Ext3` and `Ret_ext1`–`Ret_ext3`; the metadata map these to the WD, Ext and Ret-Ext pool labels. Use the same `work_dir`, sample definitions and settings across the notebooks. The processing environment used R 4.5.2 on Linux, Seurat 5.3.1, SeuratObject 5.2.0, SoupX 1.6.2, scDblFinder 1.24.10 and Harmony 1.2.4. The plotting notebooks require Arial; the upstream checks also use Cairo and Fontconfig.

| Notebooks | Steps and outputs |
| --- | --- |
| [01](01_QC.Rmd) | Correct ambient RNA, filter nuclei and save per-library singlet objects under `stage1/`. |
| [02](02_Processing.Rmd), [03](03_Global_marker_review.Rmd), [04](04_Major_cell_annotation.Rmd), [05](05_Major_cell_review.Rmd) | Integrate and cluster, review markers, assign the six major cell types and plot the annotations. The six-class object is `stage1/05_six_class_review/global_celltype_six_class.rds`. |
| [07](07_Neuron_subclustering.Rmd), [08](08_Neuron_annotation.Rmd), [09](09_Neuron_marker_review.Rmd), [10](10_Final_neuron_annotation.Rmd) | Use the six-class object to recluster excitatory and inhibitory nuclei, review the subtype markers and apply the accepted inhibitory mapping. |
| [11](11_IEG_high_inputs.Rmd) | Write `stage3/01_inputs/neuron_metadata.csv`, both class-specific `IEGhigh_payload.rds` files, `pseudobulk_inputs.rds` and the count and pool metadata tables. |
| [12](12_Major_cell_figures.Rmd), [13](13_Composition_figures.Rmd) | Render the annotated UMAP and marker panels (Fig. 3B/C and Fig. S7D) and composition panels (Fig. 3D/K and Fig. S9A). |

Review the broad-class and neuronal subtype markers before proceeding to downstream analysis. Use a new output directory for a new run. In notebook 01, leave `selected_samples` as `NULL` to process all libraries, or provide a character vector to select libraries. To resume an interrupted run, use the same notebook, paths and settings as when it started; checkpoints from the earlier `.R` layout retain the original script checksum.

Optional diagnostics: use the [library check](checks/00_Preflight.Rmd) before notebook 01 to check matrices, packages and fonts, or the [neuronal check](checks/06_Neuron_preflight.Rmd) before notebook 07 to check the saved input and Seurat environment. Set their paths as in the analysis notebooks. These diagnostic reports are not required by the analysis.

### Downstream analysis and figures

Open each notebook in RStudio and run its chunks in order after setting its input and output paths. Use fresh R sessions for the different analyses. The expression environment used R 4.5.1, DESeq2 1.50.2, edgeR 4.8.2, Matrix 1.7.4, Monocle 2.38.0, DDRTree 0.1.5, igraph 2.2.1, dplyr 1.1.4 and fgsea 1.36.2. CellChat uses version 1.6.1. The figure notebooks use Arial and their listed export packages (including Cairo, ragg or svglite as applicable).

| Notebooks or script | Inputs and outputs |
| --- | --- |
| [14](14_Subtype_proportions.Rmd), [15](15_Proportion_figures.Rmd) | Fit binomial mixed models from `neuron_metadata.csv`. Specify the finalized excitatory and inhibitory statistics tables to draw the proportion panels (Fig. 3L and Fig. S9B). |
| [16](16_Pseudobulk_DESeq2.Rmd), [17](17_Differential_expression_figures.Rmd) | Use the excitatory and inhibitory payloads to generate DESeq2 tables and volcano plots (Fig. 3H and supplementary contrasts). |
| [18](18_Signature_trajectories.Rmd), [19](19_Trajectory_figures.Rmd) | Use the excitatory payload and the Ret-Ext versus Ext DEG table to fit trajectories, summarize by pool and plot the selected signatures (Fig. 3F/G). |
| [20](20_Module_heatmap.Rmd) | Use the excitatory payload and the notebook's gene modules to render the expression heatmap (Fig. 3I). |
| [21](21_Candidate_GSEA.Rmd), [22](22_Enrichment_figures.Rmd) | Use the DESeq2 results and frozen mapping tables to run candidate GSEA and plot pathways selected at nominal P < 0.05 with FDR and NES (Fig. 3J and supplementary plots). |
| [23](23_Gene_category_network.py), [24](24_Gene_category_figures.Rmd) | Run the Python script on the matching DESeq2 and GSEA outputs. Then use its node and edge tables to plot gene-category membership (Fig. S8A). |
| [25](25_Orthology.Rmd), [26](26_CellChat.Rmd), [27](27_Outgoing_signaling.Rmd), [28](28_CellChat_figures.Rmd) | Obtain one-to-one orthologs from the Ensembl 116 archive, fit the abundance-weighted CellChat models with the two payloads, and summarize and plot joint-model outputs (Fig. 3P and Fig. S10). |

The proportion renderer takes finalized CSVs with `subtype`, `P`, `BH_adjusted_P` and `Label` columns, using subtype IDs from the nuclear metadata. Leave the result fields empty for unavailable tests.

For GSEA, supply separate rat gene maps with `SYMBOL,ENTREZID` columns and a pathway-membership table with `ID,Entrez` columns. The source maps used org.Rn.eg.db 3.22.0 and AnnotationDbi 1.72.0. The tables must cover the GO, KEGG and Reactome identifiers listed in notebook 21. For the network figure, `down_layout` and `up_layout` can be set to archived layout CSVs or left as `NULL` to calculate the layouts.
