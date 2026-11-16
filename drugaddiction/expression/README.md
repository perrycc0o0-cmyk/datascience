# Expression analysis and figures

Pseudobulk differential expression, signature-guided trajectories, candidate GSEA and gene–category networks for IEG-high neurons.

## Inputs and software

Use the excitatory and inhibitory `IEGhigh_payload.rds` files from `stage3/01_inputs/`. Each payload contains corrected integer counts, nuclear metadata and feature annotations. Pool IDs are `WD1`–`WD3`, `Ext1`–`Ext3` and `Ret-Ext1`–`Ret-Ext3`.

The recorded expression environment used R 4.5.1, DESeq2 1.50.2, edgeR 4.8.2, Matrix 1.7.4, Monocle 2.38.0, DDRTree 0.1.5, igraph 2.2.1, dplyr 1.1.4 and fgsea 1.36.2. The network builder uses the Python 3 standard library. Figure export requires Arial Regular and Bold, Cairo, `ggplot2`, `ggrepel`, `ragg`, `svglite` and `systemfonts`.

Run the examples from `drugaddiction/` and replace the paths. Keep the helper and configuration files in their supplied locations. Use an empty output directory and `--help` to view command arguments.

## Differential expression and volcano plots

```bash
Rscript expression/pseudobulk_deseq2.R /path/to/work/stage3/01_inputs/Excitatory/IEGhigh_payload.rds /path/to/work/stage3/01_inputs/Inhibitory/IEGhigh_payload.rds /path/to/work/deseq2
Rscript expression/plot_deseq2_volcanoes.R /path/to/work/deseq2 /path/to/work/volcanoes
```

The analysis compares conditions using pool-level counts and DESeq2 with independent filtering and BH adjustment. Output tables are saved under `datasets/<class>__<contrast>/tables/`: `ALL_genes.csv` contains the complete feature table, and `BH005_genes.csv` contains genes with adjusted P < 0.05. The volcano renderer reads these results.

## Signature-guided trajectories

```bash
Rscript expression/signature_trajectories.R /path/to/work/stage3/01_inputs/Excitatory/IEGhigh_payload.rds /path/to/work/deseq2/datasets/Excitatory__RetExt_vs_Ext/tables/ALL_genes.csv /path/to/work/trajectories
Rscript expression/plot_signature_trajectories.R /path/to/work/trajectories /path/to/work/trajectory_figures
```

The first script fits the Monocle 2 trajectories and writes pool summaries and comparison results. The figure renderer uses the 40-gene signatures. Pseudotime describes relative transcriptional position.

## GSEA and gene–category networks

Prepare the following external reference tables:

| File | Required columns | Contents |
| --- | --- | --- |
| Excitatory gene mapping | `SYMBOL,ENTREZID` | Rat gene mappings for the excitatory analysis |
| Inhibitory gene mapping | `SYMBOL,ENTREZID` | Rat gene mappings for the inhibitory analysis |
| Pathway membership | `ID,Entrez` | One row per candidate pathway–gene pair |

Use the matching frozen annotation snapshots for exact reproduction. The source gene mappings used `org.Rn.eg.db` 3.22.0 and `AnnotationDbi` 1.72.0. The pathway table requires the GO, KEGG and Reactome memberships for the 43 candidates in `config/candidate_pathways.csv`, including `R-RNO-399719`. These reference snapshots are not included in the repository.

```bash
Rscript expression/candidate_gsea.R /path/to/work/deseq2 /path/to/references/excitatory_symbol_entrez.csv /path/to/references/inhibitory_symbol_entrez.csv /path/to/references/rat_candidate_pathway_membership.csv /path/to/work/gsea
Rscript expression/plot_candidate_gsea.R /path/to/work/gsea /path/to/work/gsea_figures
python expression/gene_category_network.py /path/to/work/deseq2 /path/to/work/gsea /path/to/work/gene_category_network
Rscript expression/plot_gene_category_network.R /path/to/work/gene_category_network /path/to/work/network_figures
```

The GSEA renderer displays pathways with nominal P < 0.05 and encodes FDR and NES. The network builder uses the selected excitatory pathways and category assignments in `config/network_categories.csv`. Network links indicate gene–category membership.

To use archived network coordinates, append the downregulated and upregulated `*_nodes_and_layout.csv` files to the network renderer command. The files must match the current network tables.

## Figure outputs

| Figure | Renderer |
| --- | --- |
| Fig. 3F/G | `plot_signature_trajectories.R` |
| Fig. 3H and supplementary volcanoes | `plot_deseq2_volcanoes.R` |
| Fig. 3J and inhibitory GSEA | `plot_candidate_gsea.R` |
| Excitatory gene–category network | `plot_gene_category_network.R` |

Renderers export PDF, SVG and 600-dpi PNG files.
