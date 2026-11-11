# Rat infralimbic cortex snRNA-seq workflow

Analysis code for the comparison of withdrawal alone (WD), extinction (Ext), and retrieval-extinction (Ret-Ext) in the rat infralimbic cortex. This directory contains the computational workflow used for the current snRNA-seq analyses and their figures. It does not contain sequencing matrices, serialized study objects, or generated results.

## Experimental design

There are three independent pooled biological samples per condition and nine libraries in total. Each pooled sample comprises infralimbic tissue from five rats. The pooled sample, not a nucleus or a rat within a pool, is the biological replicate for sample-level comparisons. Use `WD`, `Ext`, and `Ret-Ext` as condition names. Original library identifiers beginning with `Con` map to WD.

## Workflow

| Directory | Analysis | Scope |
|---|---|---|
| `upstream/` | Per-library quality control, ambient RNA correction, doublet removal, sample-batch integration, cell annotation, neuronal subclustering, and IEG-high classification | Fig. 3B–E and the snRNA-seq quality-control and annotation panels |
| `expression/` | Pairwise pseudobulk differential expression, signature-guided pseudotime, preranked GSEA, and gene–category association | Fig. 3F–J and related supplementary expression panels |
| `signaling/` | Neuronal proportions and CellChat inference and summaries | Fig. 3K, L and P, and related supplementary panels |
| `figures/` | Figure assembly and plotting from saved analysis tables | Computational figure panels |

Read each module's instructions before running it. Start with the per-library workflow, complete the reviewed cell annotations, and pass the resulting analysis objects and pseudobulk inputs to the downstream modules. The upstream and downstream modules keep intermediate objects separate from source code. Use a new output directory for a new analysis run.

From the `drugaddiction/` directory, the upstream entry point is:

```bash
python upstream/run.py preprocess --input-dir /path/to/10x_libraries --run-dir /path/to/work/study
python upstream/run.py global --run-dir /path/to/work/study
python upstream/run.py neurons --run-dir /path/to/work/study
python upstream/run.py prepare --run-dir /path/to/work/study
```

Keep the work directory outside `drugaddiction/`. The `prepare` step writes the handoff files under the work directory's `stage3/01_inputs/` subdirectory. Each downstream module documents which files it reads and which reference inputs must be supplied separately.

## Input data and software

The upstream workflow requires paired raw and cell-called 10x feature-barcode matrices from each of the nine libraries. The raw matrices provide the background droplets for ambient RNA estimation. The documented reference inputs include gene annotation, cell-cluster mappings, gene sets, and rat-to-human orthologs where applicable. Input paths and output paths are supplied locally; do not put credentials, private storage paths, or study matrices into the repository.

Use the package versions listed in the module instructions. The source analyses used SoupX 1.6.2, scDblFinder 1.24.10, Seurat 5.3.1, SeuratObject 5.2.0, Harmony 1.2.4, DESeq2 1.50.2, Monocle 2.38.0, DDRTree 0.1.5, fgsea 1.36.2, and CellChat 1.6.1. Use an R/Bioconductor environment compatible with the versions required by a given module. Figure export uses Arial when available and requires the indicated graphics devices.

## Analysis definitions

### Preprocessing and neuronal states

Ambient RNA correction and doublet detection are performed independently within each library. The workflow uses fixed quality-control criteria and retains corrected integer RNA counts. Harmony integrates sample batches; experimental condition is not an integration covariate. The reported analysis retained 131,480 singlet nuclei. This is an observed result, not a target to which quality-control parameters are tuned.

IEG-high nuclei are identified using the summed log-normalized expression of eight immediate-early genes. The 95th-percentile threshold is calculated separately within excitatory and inhibitory neuronal classes across the study conditions.

### Differential expression and pseudotime

Pairwise differential expression uses corrected integer counts summed within each pooled sample and neuronal class. The primary DESeq2 contrasts use the six pools from the two conditions being compared, with default independent filtering and Benjamini–Hochberg adjustment.

The expression trajectories use a separate nine-pool DESeq2 fit. For the Ret-Ext versus Ext contrast, genes are ranked by nominal P within each fold-change direction. The displayed signatures use 40 genes per direction; this is a ranking rule, not a requirement that 40 genes pass a significance threshold. The Monocle 2 workflow fits the trajectories to IEG-high excitatory nuclei and summarizes pseudotime within each pooled sample. Pseudotime is a relative coordinate along a fitted transcriptional trajectory. The signatures and comparisons come from the same study samples, so these displays do not constitute independent validation of the signatures.

### Enrichment and functional associations

Preranked GSEA uses the signed DESeq2 Wald statistic. The current plots display candidate gene sets with nominal P < 0.05 and report candidate-set BH-FDR separately. A gene–category network links selected genes to functional categories. Its edges indicate category membership, not protein–protein interactions.

### Subtype proportions and predicted communication

Subtype composition within the IEG-high population and the IEG-high fraction within a subtype have different denominators. The proportion-model code and the plotting of final statistical annotations are separate steps. The plot scripts can read final statistical tables without fitting a model.

The CellChat comparison uses the abundance-weighted joint neuronal model. The receiver scopes in the outgoing-strength summary are subsets of that same network: other excitatory subtypes or all other included neuronal subtypes. These estimates are inferred from expression, not direct measurements of projections or synaptic strength.

## Method references

- Young and Behjati. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience* 9, giaa151 (2020).
- Germain et al. Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research* 10, 979 (2021; version 2, 2022).
- Hochgerner et al. Neuronal types in the mouse amygdala and their transcriptional response to fear conditioning. *Nature Neuroscience* 26, 2237–2249 (2023).
- Qiu et al. Reversed graph embedding resolves complex single-cell trajectories. *Nature Methods* 14, 979–982 (2017).
- Berdugo-Vega et al. Cognitive rejuvenation through partial reprogramming of engram cells. *Neuron* 114, 1102–1116.e7 (2026). doi:10.1016/j.neuron.2025.11.028.

The workflow covers computational snRNA-seq analyses. Histology, physiology, and behavioral measurements are outside its scope.
