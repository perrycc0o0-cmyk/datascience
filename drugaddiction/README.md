# Drug addiction snRNA-seq analysis

Analysis and figure-generation code for single-nucleus RNA sequencing of the rat infralimbic cortex (IL). The workflow covers per-library quality control, ambient RNA correction, doublet detection, sample-batch integration, cell annotation, neuronal subclustering, IEG-high classification, pseudobulk differential expression, marker-program scoring, pathway enrichment, proportion models, and inferred cell–cell communication.

## Experimental design

The study includes three conditions: WD, Ext, and Ret–Ext. Each condition has three independent pooled biological samples, each comprising IL tissue from five rats. One sequencing library represents each pooled sample, giving nine libraries in total. The biological replicate number is **n = 3 per condition**; nuclei within a library are not independent animal replicates.

Harmony uses `sample` as the batch variable, not experimental condition. Differential expression uses unnormalized integer UMI counts from the RNA assay, aggregated by sample and analyzed with DESeq2. Marker-program scores retain nucleus-level distributions. Proportion models include condition as a fixed effect and sample as a random intercept.

## Directory structure

```text
drugaddiction/
├── config/
│   ├── analysis.json               Analysis parameters and input/output paths
│   ├── samples.csv                 Libraries, conditions, and 10x file paths
│   ├── major_annotations.csv       Cluster-to-major-cell-class mapping
│   ├── subtype_annotations.csv     Neuronal cluster-to-subtype mapping
│   ├── markers.json                Cell-class and neuronal subtype markers
│   ├── rat_human_orthologs.tsv     One-to-one orthologs for CellChat
│   ├── heatmap_genes.csv           Genes, modules, and heatmap order
│   └── enrichment_terms.csv        Term-name patterns for enrichment plots
├── R/                             Processing, statistical, and plotting helpers
├── scripts/                       Analysis scripts 01–12 and annotation utility
│   └── README.md                  Per-script documentation and parameters
└── README.md
```

Input data and analysis results are excluded from the code directory. Results default to `results/`; use `--output` to specify another results root.

## R environment and dependencies

Install the following packages in a compatible R/Bioconductor environment:

| Function | R packages |
|---|---|
| Sparse matrices and single-nucleus objects | Seurat, SeuratObject, Matrix |
| Ambient RNA correction | **SoupX 1.6.2** |
| Doublet detection | **scDblFinder 1.24.10**, SingleCellExperiment, SummarizedExperiment |
| Batch integration | harmony |
| Differential expression | DESeq2 |
| Proportion models | lme4 |
| GO/KEGG enrichment | clusterProfiler, AnnotationDbi, org.Rn.eg.db, KEGGREST |
| Cell–cell communication | CellChat |
| Configuration and plotting | jsonlite, ggplot2, ggrepel, scales, igraph |

The preprocessing script explicitly requires the SoupX and scDblFinder versions listed above. Each analysis stage saves `sessionInfo.txt` to record its software environment. KEGG analysis requires access to KEGG services, and communication analysis requires the human CellChat database. PNG and PDF figures use Arial. PDF export requires Cairo support in R; ensure Arial is available to the runtime.

## Prepare the inputs

### 1. Library table

Edit [`config/samples.csv`](config/samples.csv), with one row per library and the following columns:

```csv
sample,condition,filtered_h5,raw_h5
Con1,WD,data/Con1/filtered_feature_bc_matrix.h5,data/Con1/raw_feature_bc_matrix.h5
```

`sample` must be unique. Use `WD`, `Ext`, and `Ret-Ext` for `condition`, with three rows per condition. The file identifiers `Con1–Con3` correspond to WD1–WD3 in the manuscript, and `Ret_ext1–Ret_ext3` correspond to Ret–Ext1–Ret–Ext3. These are file identifiers and display names for the same libraries, not additional samples.

`filtered_h5` is the cell-called 10x matrix. `raw_h5` is the corresponding raw matrix, including empty droplets, from the same library. Features and barcodes must be aligned. Paths may be relative to the project working directory or absolute. Input feature names should be rat gene symbols.

### 2. Major cell-class mapping

[`config/major_annotations.csv`](config/major_annotations.csv) contains `cluster,cell_type` columns and mappings for 38 clusters. Cluster identifiers are specific to their clustering object; matching numbers alone do not justify transferring a mapping to another object. Review the assignments against `02_annotation/cluster_markers.csv` and marker plots.

Major classes include `Excitatory`, `Inhibitory`, `Oligodendrocytes`, `OPC`, `Astrocytes`, and `Microglia`. Other clearly supported classes may retain their own names. Clusters without sufficient evidence should not be forced into these six classes. Unassigned clusters remain `Unknown`; `annotated.rds`, used for neuronal analysis, is saved only after all clusters are annotated.

### 3. Neuronal subtype mapping

Edit [`config/subtype_annotations.csv`](config/subtype_annotations.csv) with `neuron_type,cluster,subtype` columns. Use `Excitatory` or `Inhibitory` for `neuron_type`. Cluster identifiers must come from the corresponding independently reclustered neuronal object.

Script 03 exports cluster markers, marker DotPlots/UMAPs, and `subtype_annotations_template.csv` for each lineage. Use these to complete the manual mapping; major-cell-class cluster identifiers are not interchangeable with neuronal subcluster identifiers. Apply the reviewed mapping using `apply_annotations.R --level subtype`. Complete subtype annotation is required for downstream IEG, differential-expression, proportion, and communication analyses.

### 4. Rat–human ortholog table

CellChat uses the human database and requires [`config/rat_human_orthologs.tsv`](config/rat_human_orthologs.tsv), with these tab-separated columns:

```text
rat_symbol	human_symbol	orthology_type
```

Supply genuine one-to-one ortholog relationships. Use `one_to_one` or `ortholog_one2one` for `orthology_type`. Retained rat and human symbols must each be unique. Uppercasing rat symbols is not a substitute for ortholog mapping. Communication analysis stops if a complete usable mapping is unavailable.

### 5. Markers and heatmap genes

[`config/markers.json`](config/markers.json) contains cell-class and neuronal subtype marker panels. [`config/heatmap_genes.csv`](config/heatmap_genes.csv) contains `gene_order,module_id,module,gene` columns specifying the fixed six-module heatmap groups and display order. These are configured functional groups, not newly discovered clustering modules.

## Run order

Run commands from `drugaddiction/`. Use the same configuration and results root across stages. Execute the following commands individually:

```bash
Rscript scripts/01_preprocess_libraries.R --config config/analysis.json --output results
Rscript scripts/02_integrate_annotate.R --config config/analysis.json --output results
```

Review the cluster markers and plots from script 02, then complete or confirm `config/major_annotations.csv`. If `02_annotation/annotated.rds` has not been created, apply the complete mapping without reclustering:

```bash
Rscript scripts/apply_annotations.R --level major
```

After confirming the major-cell-class mapping, recluster neurons:

```bash
Rscript scripts/03_neuron_subclustering.R --config config/analysis.json --output results
```

Review subtype markers and plots, complete `config/subtype_annotations.csv`, and apply the mapping:

```bash
Rscript scripts/apply_annotations.R --level subtype
```

This creates `03_neurons/Excitatory_annotated.rds` and `03_neurons/Inhibitory_annotated.rds`. Downstream IEG analysis and cell-type plotting preferentially read these reviewed objects. Once both lineage mappings are complete, run:

```bash
Rscript scripts/04_ieg_definition.R --config config/analysis.json --output results
Rscript scripts/05_pseudobulk_DE.R --config config/analysis.json --output results
Rscript scripts/06_marker_program_scores.R --config config/analysis.json --output results
Rscript scripts/07_GO_KEGG_enrichment.R --config config/analysis.json --output results
Rscript scripts/08_subtype_GLMM.R --config config/analysis.json --output results
Rscript scripts/09_CellChat_analysis.R --config config/analysis.json --output results
Rscript scripts/10_plot_cell_types.R --config config/analysis.json --output results
Rscript scripts/11_plot_expression.R --config config/analysis.json --output results
Rscript scripts/12_plot_communication_proportions.R --config config/analysis.json --output results
```

### Annotation and results directories

Numbered scripts create their own stage directories and **refuse to overwrite nonempty stage directories**. When changing processing, clustering, or statistical parameters, choose a new results root, for example `--output results_reviewed`.

For initial annotation, first obtain markers and plots, then complete the mapping and apply it to the same results root using `apply_annotations.R`. There is no need to delete the original objects or rerun script 01. The annotation utility also refuses to overwrite existing annotated objects. To revise a saved complete annotation, preserve the current results and use the required upstream objects with the revised mapping under a new results root. Do not mix downstream results from different annotation versions.

To repeat only a downstream stage, copy the retained upstream stage directories to the same locations under the new root. For example, script 02 requires `01_preprocessing/merged.rds`, and script 03 requires `02_annotation/annotated.rds`, within that root. Changing `--output` alone does not make a script read upstream objects from another root.

Annotation templates and marker results support manual review; scripts do not fill missing labels from expected names or nucleus counts. The annotation commands above use the default `config/analysis.json` and its output path. Keep configuration and output paths consistent when changing the results root.

## Scripts and manuscript analyses

| Script | Analysis | Methods / figure content |
|---|---|---|
| `01_preprocess_libraries.R` | Fixed QC, SoupX, scDblFinder, library merging | snRNA-seq preprocessing and QC |
| `02_integrate_annotate.R` | PCA, sample Harmony, global clustering, major-cell-class annotation | Major-cell-class UMAPs, markers, supplementary QC |
| `03_neuron_subclustering.R` | Independent Exc/Inh clustering, subtype annotation, marker plots | Neuronal subtypes in Fig. 3 and related Fig. S6 panels |
| `apply_annotations.R` | Apply reviewed cell-class or subtype mappings | Cell identity annotation |
| `04_ieg_definition.R` | Eight-gene expression sum and class-specific 95th percentile | IEG-high definition, thresholds, distributions |
| `05_pseudobulk_DE.R` | Sample-by-neuronal-class/IEG-class counts, DESeq2, BH | Condition-dependent expression in Fig. 3 and Fig. S7 |
| `06_marker_program_scores.R` | IEG-high signatures and nucleus-level scores | Fig. 3F/G score distributions |
| `07_GO_KEGG_enrichment.R` | Direction-specific GO/KEGG enrichment and backgrounds | Fig. 3 pathways and related supplementary results |
| `08_subtype_GLMM.R` | Binomial GLMM with condition fixed effect, sample random intercept, BH | Fig. 3L and Fig. S8B |
| `09_CellChat_analysis.R` | One-to-one ortholog mapping and IEG-high subtype communication | Predicted communication in Fig. 3P and Fig. S9 |
| `10_plot_cell_types.R` | Cell-class/subtype UMAPs, DotPlots, marker UMAPs | Cell identity plots in Fig. 3 and Fig. S6 |
| `11_plot_expression.R` | Volcanoes, enrichment, scores, fixed six-module heatmap | Expression and pathway plots in Fig. 3 and Fig. S7 |
| `12_plot_communication_proportions.R` | Subtype proportions, GLMM annotations, communication plots | Fig. 3L, Fig. S8B, Fig. S9, communication summaries |

See [`scripts/README.md`](scripts/README.md) for per-script inputs, outputs, parameters, and statistical interpretation. Its “Implementation settings” section collects settings not individually expanded in the manuscript, including SoupX preclustering, DESeq2 filtering and fitting, GLMM inference and optimization, ortholog policy, and plot selection rules.

The scripts cover the computational components of Fig. 3 (B–D, F–L, P) and related supplementary figures. Wet experiments in Fig. 3M–O are outside this workflow. `config/analysis.json` is the analysis-parameter entry point. CellChat results are expression-based communication inferences, not direct measurements of neuronal projections or connection strength.
