# Expression analysis and figures

Pseudobulk differential expression, signature-guided trajectories, candidate GSEA and gene–category networks for IEG-high neurons.

## Inputs and packages

Use the excitatory and inhibitory `IEGhigh_payload.rds` files from stage 3. Each contains corrected counts, nuclear metadata and feature annotations.

The expression environment used R 4.5.1, DESeq2 1.50.2, edgeR 4.8.2, Matrix 1.7.4, Monocle 2.38.0, DDRTree 0.1.5, igraph 2.2.1, dplyr 1.1.4 and fgsea 1.36.2. Figure export requires Arial, Cairo, ragg and svglite. The Python network builder uses Python 3.

## Usage

Open a script, edit the paths and sample table near the top, and run the complete file. Use **Source** in RStudio for R files and a Python editor for the Python file. Choose an empty directory for each analysis or figure output.

| Script | Input | Output |
| --- | --- | --- |
| `pseudobulk_deseq2.R` | Excitatory and inhibitory payloads | DESeq2 all-gene and BH-significant tables |
| `plot_deseq2_volcanoes.R` | DESeq2 result directory | Volcano plots |
| `signature_trajectories.R` | Excitatory payload and Ret-Ext versus Ext DEG table | Trajectories, pool summaries and test results |
| `plot_signature_trajectories.R` | Trajectory result directory | Fig. 3F/G |
| `candidate_gsea.R` | DESeq2 results and reference tables | Candidate GSEA results and leading-edge genes |
| `plot_candidate_gsea.R` | GSEA result directory | GSEA plots |
| `gene_category_network.py` | DESeq2 and GSEA results | Gene–category node and edge tables |
| `plot_gene_category_network.R` | Network result directory | Gene–category network figures |

The trajectory renderer uses the 40-gene signatures. The GSEA renderer displays pathways with nominal P < 0.05 and shows FDR and NES. Network edges represent gene–category membership.

For GSEA, supply matching frozen rat gene-mapping tables with columns `SYMBOL,ENTREZID` for each neuronal class, and a pathway-membership table with columns `ID,Entrez`. The source mappings used org.Rn.eg.db 3.22.0 and AnnotationDbi 1.72.0. Pathway IDs and category definitions are listed in the scripts; the reference tables must cover those GO, KEGG and Reactome entries.

For archived network layouts, set `down_layout` and `up_layout` to the matching CSV paths. Leave both `NULL` to calculate the layout. Figures are saved as PDF, SVG and 600-dpi PNG.
