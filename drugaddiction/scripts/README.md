# Script documentation

This directory organizes preprocessing, annotation, IEG-high classification, differential expression, marker-program scoring, enrichment, subtype proportion models, CellChat analysis, and plotting for rat IL snRNA-seq. Parameters are centralized in `../config/analysis.json`. Markers, subtype mappings, heatmap genes, and pathway display rules are stored in their respective configuration files.

Each condition contains three independent pooled biological samples, each comprising IL tissue from five rats. Between-condition differential expression uses the pooled sample as the replicate. Nucleus-level marker distributions and sample-level condition effects are handled by distinct statistical procedures.

## Usage and file relationships

Call scripts from the project root. All paths are relative to that root, not to `scripts/`. For example:

```bash
Rscript scripts/01_preprocess_libraries.R --config config/analysis.json --output results
Rscript scripts/02_integrate_annotate.R --config config/analysis.json --output results
```

Numbered scripts share two optional arguments. The annotation utility also uses `--level`:

| Argument | Default | Purpose |
|---|---|---|
| `--config` | `config/analysis.json` | Analysis configuration |
| `--output` | Configured `paths.output`, default `results` | Shared results root for all stages |
| `--level` | `major` | Annotation level for `apply_annotations.R`: `major` or `subtype` |

The global random seed is 42. Outputs are organized by stage. A numbered script stops if its target stage directory is nonempty, preserving existing objects and figures. Use a new results root when changing processing settings or saved annotations, and retain the upstream stage files required for that analysis.

The workflow is:

```text
01 Library preprocessing
  -> 02 Global integration and major-cell-class annotation
     -> 03 Independent neuronal subclustering
        -> 04 IEG-high classification
           -> 05 Pseudobulk DE -> 07 Enrichment
           -> 06 Marker-program scores
           -> 08 Proportion GLMMs
           -> 09 CellChat

10 Plot cell classes, QC, and neuronal subtypes
11 Plot DE, scores, six-module heatmap, and enrichment
12 Plot proportions and communication
```

`config/samples.csv` must specify nine unique samples, with conditions `WD`, `Ext`, and `Ret-Ext` and three samples per condition. Each row identifies the filtered and raw 10x HDF5 files. Major-cell-class and subtype assignments come from `config/major_annotations.csv` and `config/subtype_annotations.csv`. Cluster identifiers must match the corresponding object. Unassigned clusters remain `Unknown` and cannot enter downstream analyses requiring complete annotation.

After clustering and marker review, apply the manually confirmed mappings without reclustering:

```bash
Rscript scripts/apply_annotations.R --level major --output results
Rscript scripts/apply_annotations.R --level subtype --output results
```

Apply the major mapping after step 02 to create `02_annotation/annotated.rds`. Apply subtype mappings after step 03 to create `03_neurons/Excitatory_annotated.rds` and `Inhibitory_annotated.rds`. Downstream readers preferentially use `*_annotated.rds`. If step 02 already saved a complete annotated object, do not repeat the major annotation utility. Existing target objects are never overwritten.

## Scripts and figure correspondence

| File | Main analysis | Corresponding content |
|---|---|---|
| `01_preprocess_libraries.R` | Fixed QC, SoupX, scDblFinder | Methods: preprocessing, ambient RNA, doublets; QC data |
| `02_integrate_annotate.R` | Global PCA, sample Harmony, clustering, major classes | Fig. 3B–D; integration, clustering, annotation |
| `03_neuron_subclustering.R` | Independent Exc/Inh clustering and subtype annotation | Fig. S6E–J; neuronal subtypes |
| `04_ieg_definition.R` | Class-specific 95th-percentile IEG-high definition | Basis for Fig. 3E and downstream IEG-high analyses |
| `05_pseudobulk_DE.R` | Sample-level DESeq2 | Basis for Fig. 3H/I and Fig. S7 volcanoes |
| `06_marker_program_scores.R` | Two marker-program signatures and nucleus-level scores | Fig. 3F/G and related score data |
| `07_GO_KEGG_enrichment.R` | Up/down GO and KEGG enrichment | Fig. 3J and other conditions/classes |
| `08_subtype_GLMM.R` | Binomial model: condition fixed effect, sample random intercept | Fig. 3L, Fig. S8B |
| `09_CellChat_analysis.R` | IEG-high communication inference and summaries | Communication figures and Fig. S9 data |
| `10_plot_cell_types.R` | Cell classes, QC, subtype UMAPs and markers | Fig. 3B–D, Fig. S6A–J |
| `11_plot_expression.R` | Scores, volcanoes, module heatmap, enrichment | Fig. 3F–J, Fig. S7A–C, separate result figures |
| `12_plot_communication_proportions.R` | Composition, GLMM proportions, communication summaries | Fig. 3K/L/P, Fig. S8A/B, Fig. S9A–E |
| `apply_annotations.R` | Apply confirmed cluster mappings to existing objects | Cell-class/subtype annotation without recomputing embeddings or clusters |

These snRNA-seq scripts do not generate the experimental schematic in Fig. 3A or the wet-experiment images and animal-level results in Fig. 3M–O.

## Implementation settings

The table below collects computational choices not individually expanded in Methods. Complete inputs, outputs, and parameters follow in the individual script sections. These are implementation settings, not additional biological groups or predetermined results.

| Stage | Implementation |
|---|---|
| Library QC and SoupX | QC uses uncorrected observed counts after initial Seurat filtering. SoupX estimates contamination in these initially retained called nuclei; the recorded fixed-QC mask is then applied to corrected counts. Preclustering uses LogNormalize 10,000, 2,000 HVGs, up to 30 PCs, and resolution 0.5, distinct from final major-class resolution 1.0. |
| Global integration | PC selection uses two criteria based on percentages of PC standard deviations, as detailed in step 02, not variance percentages based on squared standard deviations. Graph k=20; resolutions 0.1, 0.2, 0.4, 0.6, 0.8, 1.0; Harmony `project.dim=FALSE`. Only positive cluster markers are exported; the Wilcoxon test remains two-sided. |
| Neuronal subtypes | Harmony receives up to 30 PCs; graphs and UMAPs use the first 19 available Harmony dimensions. Subset scaling does not repeat UMI/mitochondrial regression. |
| IEG-high | R `quantile(type=7)`; include ties satisfying `score >= cutoff`. Use available genes from the eight-gene panel and export availability. |
| DESeq2 | Total counts ≥10 and positive counts in ≥3 samples; fit all nine libraries jointly with `~ condition`. Wald test, parametric dispersion, ratio size factors, independent filtering enabled, automatic outlier replacement disabled, default Cook's distance handling retained. |
| Marker programs | Merge Exc/Inh for within-condition IEG-high/low signature discovery. Seurat `min.pct=0.1` and `logfc.threshold=0`. Score the same two signatures separately in All/Exc/Inh scopes with `nbin=24`, `ctrl=100`, seed=42. Distribution tests use the unpaired two-sided Wilcoxon normal approximation with continuity correction. |
| GO/KEGG | Background comprises count-prefiltered genes mappable to each annotation resource; `maxGSSize=500`. Analyze GO BP/CC/MF separately. Retain one-to-many Entrez mappings and deduplicate identifiers. KEGG uses online data at execution time and saves gene–pathway membership. |
| Binomial GLMM | Ext reference; lme4 `glmer`, `nAGQ=1`, `optimizer=bobyqa`, two-sided Wald z inference, ORs and Wald 95% CIs. Singularity tolerance 1e-4. Failed/nonconverged fits have P=NA. Planned BH family sizes are six and eight. Estimate the full sample random-intercept variance without extra weighting. |
| CellChat | Genuine one-to-one rat–human mapping, `triMean` estimator, and Overall/WD/Ext/Ret–Ext analyses within each All/Exc/Inh scope. Exclude self-edges from cross-subtype outgoing summaries, save autocrine quantities separately, and set relative changes to NA when the Ext baseline is zero. |
| Visualization | Arial, base size 14, 600-dpi PNG and Cairo PDF. Up to 20 volcano gene labels, 12 enrichment terms per direction, and 20 communication pathways. See steps 10–12 and plotting helpers for ranking, colors, and canvas sizes. |

### Data and statistical definitions

- `RNA/counts` contains SoupX-corrected, unnormalized integer UMI counts, aggregated by sample for DESeq2. `observed/counts` retains counts before ambient correction. Neither represents Harmony coordinates or scaled expression. “Unnormalized” does not mean “uncorrected by SoupX.”
- Enrichment foregrounds require DE BH-adjusted P≤0.05 and the relevant direction, without the volcano's 1.15-fold magnitude filter. Complete term tables retain all returned results. `*_significant.csv` follows clusterProfiler's P/adjusted-P/q-value filtering, whereas plotting selects from complete tables using actual adjusted P≤0.05 and term-name patterns. Plot selections and significant tables may therefore differ and are not interchangeable definitions.
- Score-plot stars use raw nucleus-level Wilcoxon P values. DESeq2, enrichment, and GLMM results use their respective BH-adjusted P values. Signature discovery and Seurat cluster-marker `p_val_adj` columns use Bonferroni adjustment. Similar column names do not make these adjustments interchangeable.
- The heatmap averages expression within Ext/Ret–Ext and then Z-standardizes the two means for each gene. With two columns, nonconstant genes have Z values of ±1/√2. Color therefore primarily conveys direction, not comparable effect magnitudes across genes. Original condition means are saved separately.
- `Con1–Con3` and `Ret_ext1–Ret_ext3` are file identifiers corresponding to manuscript WD1–WD3 and Ret–Ext1–Ret–Ext3; the design still contains nine libraries.
- Retained nucleus counts, IEG-high cutoffs, DEG counts, P values, and communication percentages are computed from the inputs, not enforced as target values.

### Figure scope

Communication scope `All` follows the main Methods and includes all annotated Exc and Inh subtypes. Scope `Excitatory` includes only Exc subtypes. If `CrossType` in an earlier supplementary legend denotes selected subtypes, its inputs differ from `All`; changing the label alone does not make them equivalent. Proportion figures use step 08 GLMM/BH results. Animal-level t tests for wet experiments are outside that model.

The supplied supplementary Fig. S7 legend explicitly lists volcano panels A–C. Other condition/class volcanoes and enrichment figures are retained as separate outputs, without assigning supplementary panel identifiers in scripts or filenames.

## 01 — `01_preprocess_libraries.R`

### Purpose and inputs

Process nine libraries independently. `config/samples.csv` identifies `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5`. Raw matrices provide the ambient profile; filtered matrices define initially called nuclei. For multimodal HDF5 files, only `Gene Expression` is used. Features and barcodes must be unique.

### Processing and parameters

1. Initialize Seurat objects with genes detected in at least three nuclei and nuclei with at least 300 detected genes.
2. Compute detected genes, UMIs, mitochondrial and ribosomal percentages, and `log10(nFeature_RNA) / log10(nCount_RNA)` from observed counts.
3. Apply strict fixed-QC inequalities: `700 < nFeature_RNA < 5000`, `500 < nCount_RNA < 30000`, `percent.mt < 5`, complexity ratio `> 0.75`. Mitochondrial and ribosomal patterns are `^Mt-` and `^Rp[sl]`. Gene symbols must match these patterns; the script stops if no mitochondrial features match.
4. Run SoupX 1.6.2 per library. Preliminary clustering uses LogNormalize 10,000, 2,000 vst HVGs, ScaleData, up to 30 PCs, a PCA neighbor graph, and resolution 0.5. These clusters are not final biological annotations. Use `autoEstCont(forceAccept = FALSE)` and `adjustCounts(roundToInt = TRUE)`.
5. Run scDblFinder 1.24.10 on fixed-QC-passing nuclei with `clusters = FALSE` and expected doublet fraction `0.004 × original_called / 1000`. Use seed 42 separately for each library and retain predicted singlets only.
6. Store corrected integer counts in `RNA` and uncorrected counts for the same nuclei in `observed`. Prefix barcodes with sample identifiers and retain sample, condition, original QC, and doublet metadata.

### Outputs

Under `01_preprocessing/`:

- `libraries/{sample}.rds`: retained nuclei for each library.
- `merged.rds`: merged object for step 02.
- `QC/{sample}_nuclei.csv`: per-nucleus observed QC, filter flags, ambient fractions, and doublet results.
- `QC/{sample}_thresholds.png/.pdf`: four QC histograms, 60 bins, dashed fixed thresholds, 10×7-inch canvas.
- `ambient_models/{sample}.rds`: SoupX contamination-fit information.
- `library_summary.csv`: initial counts, removals, retained singlets, and contamination summaries.
- `objects.csv` and `sessionInfo.txt`: object paths and software environment.

Retention counts are outputs of the inputs and filtering steps, not filtering parameters.

## 02 — `02_integrate_annotate.R`

### Purpose and inputs

Read `01_preprocessing/merged.rds` for global integration, clustering, UMAP, and major-class marker review. Annotations come from `config/major_annotations.csv`, not a fixed target nucleus count.

### Processing and parameters

- RNA LogNormalize, scale factor 10,000.
- Select 2,000 vst HVGs and scale them while regressing `nCount_RNA` and `percent.mt`.
- Target 50 PCs, limited by nucleus and variable-feature counts.
- PC selection uses each PC standard deviation as a percentage of the sum of standard deviations. Candidate 1 is the first PC with cumulative percentage >90% and individual percentage <5%. Candidate 2 is the last adjacent percentage drop >0.1 percentage points, plus one. Take the smallest finite candidate, the 30-PC cap, and available PC count, with a minimum of two PCs. Export both candidates.
- Harmony uses `sample` only, the selected PCs, and `project.dim = FALSE`; `condition` is not a batch variable.
- Harmony neighbor graph `k.param = 20`; retain clustering columns for resolutions 0.1, 0.2, 0.4, 0.6, 0.8, and 1.0, with final identity at 1.0.
- UMAP uses the same Harmony dimensions, `n.neighbors = 30`, `min.dist = 0.3`, `metric = cosine`, seed 42.
- `FindAllMarkers` exports positive markers only, with `min.pct = 0.25`, `logfc.threshold = 0.25`, up to 1,000 nuclei per cluster, seed 42.

### Outputs and annotation

Under `02_annotation/`, save `all_nuclei.rds`, `cluster_markers.csv`, `PCA_contributions.csv`, `PCA_selection.csv`, `UMAP_coordinates.csv`, `major_annotations.csv`, and the annotation template. Once every cluster is annotated, also save `annotated.rds` for step 03 and final major-class plots.

Generate a PCA elbow plot, cluster UMAP, and marker DotPlots/FeaturePlots by cluster. When annotation is complete, also plot by cell type.

| Class | Markers |
|---|---|
| Excitatory | Slc17a7, Satb2, Slc30a3 |
| Inhibitory | Gad1, Gad2, Dlx1 |
| Oligodendrocytes | Mog, Mbp, St18 |
| OPC | Cacng4, Pdgfra, Vcan |
| Astrocytes | Slco1c1, Gfap, Aqp4 |
| Microglia | P2ry12, Tmem119, Cx3cr1 |

Diagnostic UMAPs use 9×7 inches and point size 0.3; the elbow plot uses 8×6 inches. Marker diagnostic settings are defined in `R/preprocessing.R`. Step 10 generates manuscript-summary plots.

## 03 — `03_neuron_subclustering.R`

### Purpose and inputs

Read the fully annotated `02_annotation/annotated.rds` and independently recluster Excitatory and Inhibitory neurons. Repeat normalization, variable-feature selection, PCA, Harmony, and graph construction for each lineage. Read `neuron_type,cluster,subtype` mappings from `config/subtype_annotations.csv`.

### Parameters

| Step | Excitatory | Inhibitory |
|---|---|---|
| LogNormalize scale factor | 10,000 | 10,000 |
| vst HVGs | 2,000 | 2,000 |
| ScaleData | No additional regression | No additional regression |
| Target PCA dimensions | Up to 100, limited by available data | Up to 100, limited by available data |
| Harmony | First 30 available PCs, sample-only | First 30 available PCs, sample-only |
| Graph/UMAP dimensions | First 19 available Harmony dimensions | First 19 available Harmony dimensions |
| `k.param` | 20 | 30 |
| Final resolution | 0.5 | 0.5 |
| UMAP | neighbors=10, min.dist=0.8, spread=6, repulsion.strength=1 | Seurat defaults: neighbors=30, min.dist=0.3, spread=1, repulsion.strength=1 |
| Seed | 42 | 42 |

`FindAllMarkers` uses the same min.pct=0.25, logFC threshold=0.25, positive-marker-only setting, and maximum of 1,000 nuclei per cluster as step 02.

### Markers and outputs

Exc markers: GLU_L6_CLA—Gng2; GLU_L5_6—Syt6, Foxp2, Tle4; GLU_L2_L3—Glis3; GLU_L5—Parm1; GLU_L4_5—Il1rapl2, Rorb; GLU_ADGM—Rerg, Arhgap6.

Inh markers: GABA_Amb—Slc17a6; GABA_Etv1—Etv1; GABA_Lamp5—Lamp5; GABA_MSN—Gng7; GABA_Sst—Sst; GABA_Hs3st2—Hs3st2; GABA_Pvalb—Slit2; GABA_Vip—Vip.

Under `03_neurons/`, `Excitatory.rds` and `Inhibitory.rds` retain all nuclei in each lineage. Lineage subdirectories contain cluster markers, annotations, per-nucleus metadata, parameters, UMAP coordinates, cluster/subtype UMAPs, and marker plots. Diagnostic UMAPs use 8×7 inches and point size 0.5. `subtype_annotations_template.csv` defines the mapping structure. Complete subtype annotation before subtype-dependent downstream analyses.

## 04 — `04_ieg_definition.R`

### Purpose and inputs

Read Excitatory and Inhibitory objects from `03_neurons/`, preferentially using `*_annotated.rds` from `apply_annotations.R`. Define IEG-high separately within each neuronal class, not from a combined-class percentile.

### Parameters and definition

- Join RNA layers, verify nonnegative integer counts, and LogNormalize with scale factor 10,000.
- Prespecified IEGs: `Arc, Bdnf, Btg2, Fos, Fosl2, Homer1, Npas4, Nr4a1`. Use genes available in the RNA matrix; export availability per lineage and stop if none is available.
- `ieg_sum` is the sum of available defining IEG log-normalized expression values within a nucleus, not a raw UMI sum or Z-score.
- Compute the 95th percentile separately for Exc and Inh, pooling WD, Ext, and Ret–Ext within each lineage, using `quantile(type = 7)`.
- Classify `ieg_sum >= cutoff` as IEG-high. Boundary ties can yield slightly more than 5%. Do not redefine cutoffs by sample, condition, or subtype.

Outputs under `04_ieg/` include both full lineage objects, `*_IEG_availability.csv`, `nucleus_metadata.csv`, `IEG_cutoffs.csv`, and `IEG_counts.csv`. Metadata contain Cell, sample, condition, neuron_type, subtype, ieg_sum, IEGhigh, and IEGclass. Cutoff tables record the genes and number used. Cutoffs are data-derived, not fixed to historical values or target counts.

## 05 — `05_pseudobulk_DE.R`

### Purpose and inputs

Read step 04 objects and retain IEG-high nuclei. Sum SoupX-corrected, unnormalized integer UMI counts from the RNA assay by sample within each neuronal class, forming a gene-by-nine-sample matrix. Each sample must contribute IEG-high nuclei; missing replicates are not replaced with artificial all-zero columns.

### Model and parameters

- Count prefilter: total counts ≥10 across nine samples and positive counts in at least three samples. Configured by `de.min_total_counts` and `de.min_samples`.
- DESeq2 design `~ condition`, with levels WD, Ext, Ret–Ext; fit all nine samples jointly.
- Size factors `ratio`, dispersion fit `parametric`, Wald test, `parallel = FALSE`.
- `minReplicatesForReplace = Inf` disables automatic outlier replacement. Retain default Cook's distance handling in results.
- `results(alpha = 0.05, pAdjustMethod = BH, independentFiltering = TRUE)`. Apply BH within each neuronal class and contrast; export independent-filter thresholds and quantiles.
- Contrasts: `Ext_vs_WD`, `RetExt_vs_WD`, `RetExt_vs_Ext`. Positive log2FoldChange means higher expression in the first-named condition.
- The DEG flag requires padj≤0.05 and `abs(log2FoldChange) >= log2(1.15)`, approximately 0.2016. This magnitude criterion filters outputs; it is not a composite-null Wald threshold.

### Outputs

Under `05_DESeq2/`:

- `counts/{neuron_type}_pseudobulk_counts.csv` and sample metadata.
- `objects/{neuron_type}_DESeq2.rds`: fitted DESeq2 objects.
- `tables/{neuron_type}/{contrast}_all_genes.csv`: all input genes, including tested, total_count, n_samples_expressed, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, contrast, significant, and direction.
- `*_significant.csv`: genes meeting both adjusted-P and fold-change criteria.
- `DE_summary.csv` and `filtering_parameters.csv`: gene counts, up/down counts, and independent-filter settings.

Genes failing count prefiltering remain in all-genes tables with NA test fields. Count-prefiltered genes can also lack padj after independent filtering. Raw and BH-adjusted P values are saved separately and never substituted for one another.

## 06 — `06_marker_program_scores.R`

### Purpose and inputs

Merge Exc/Inh objects from step 04 for signature discovery, retaining the class-specific IEG-high/low labels. Discovery compares nuclei within Ext or within Ret–Ext, not between conditions through DESeq2.

### Signatures and scoring

1. Compare IEG-high versus IEG-low separately within Ext and Ret–Ext using two-sided Seurat Wilcoxon tests: RNA data slot, `min.pct = 0.1`, `logfc.threshold = 0`, `only.pos = FALSE`, base=2, seed=42.
2. Seurat `p_val_adj` is Bonferroni-adjusted. Retain `p_val_adj < 0.05` and `avg_log2FC > 0`, excluding the eight defining IEGs.
3. Rank by adjusted P, then gene symbol, and select up to 500 genes for `Ext_IEGhigh_UP` and `RetExt_IEGhigh_UP`. If fewer qualify, retain the available number without adding nonqualifying genes.
4. Apply the same two signatures separately to All IEG-high, Exc IEG-high, and Inh IEG-high objects. Compute expression-matched controls independently in each scope. AddModuleScore uses RNA data, all RNA genes in the object as the pool, ctrl=100, nbin=24, seed=42, and search=FALSE.
5. Export scores for all three conditions. Plotting and distribution comparisons use Ext versus Ret–Ext only: two-sided, unpaired Wilcoxon, `exact = FALSE`, continuity correction enabled. Nucleus counts are not biological replicate numbers.

### Outputs

Under `06_scores/`, `signature_DE/` contains within-condition marker results, `gene_sets/` contains CSV/TXT signatures, and `signature_summary.csv` records candidate and selected counts. `scores_long.csv` contains per-nucleus scores, scope (All, Excitatory, Inhibitory), program, sample, condition, and subtype. `score_statistics.csv` records raw nucleus-level Wilcoxon P values, nucleus counts, medians, and W statistics. `Prism/` contains Ext/Ret–Ext score columns, padded with NA to equal lengths without adding observations.

Fig. 3F/G use the Excitatory scope for the Ret–Ext/Ext programs. The other scopes are retained as additional data. Scores describe nucleus-level marker-program expression distributions and do not replace sample-level condition DE in step 05.

## 07 — `07_GO_KEGG_enrichment.R`

### Purpose and inputs

Read the three all-genes contrast tables for each neuronal class from step 05. Select padj≤0.05 and split positive/negative log2FoldChange into Up/Down. Do not additionally impose the volcano's 1.15-fold magnitude criterion on enrichment foregrounds. Run GO BP/CC/MF and rat KEGG over-representation analysis separately for each direction.

### Parameters and background

- GO: `org.Rn.eg.db`, SYMBOL keys, separate BP/CC/MF analyses, `pool = FALSE`, without combining their testing families.
- KEGG: SYMBOL→ENTREZID mapping, organism `rno`, keyType `kegg`, current KEGG data, `use_internal_data = FALSE`. Save KEGGREST rno gene–pathway membership.
- Start the background from all count-prefiltered genes in that class/contrast, not the whole genome; intersect it with the corresponding GO ontology or KEGG annotation.
- BH adjustment; `pvalueCutoff = 0.05`, `qvalueCutoff = 0.2`, `minGSSize = 5`, `maxGSSize = 500`. Significant tables additionally require p.adjust≤0.05. Each ontology, direction, neuronal class, and contrast defines a separate GO testing family; KEGG is likewise analyzed separately.
- Retain and export one-to-many SYMBOL→ENTREZID mappings; deduplicate Entrez identifiers for analysis. Unmappable genes are excluded from the corresponding foreground/background.

### Outputs

Under `07_enrichment/tables/{neuron_type}/{contrast}/`:

- `Up_input_genes.csv` and `Down_input_genes.csv`.
- `tested_gene_universe.csv`, SYMBOL/ENTREZID mappings, and separately mapped GO BP/CC/MF and KEGG universes.
- Direction-by-database tables such as `Up_GO_BP.csv` containing all returned test results; `*_significant.csv` contains results meeting configured criteria.
- Tables retain Description, pvalue, p.adjust, qvalue, Count, geneID, direction, and database. Fitted enrichment objects are saved under `objects/`.

Run an enrichment test only when eligible genes and a background exist. Empty results have fixed column definitions; related or nonsignificant terms are not added to fill plots.

## 08 — `08_subtype_GLMM.R`

### Purpose and inputs

Read step 04 nucleus metadata and the sample table to compare Ext versus Ret–Ext. Six pooled samples define the sample random-intercept groups. No extra small weight or partial-variance correction is applied to sample effects.

### Numerators, denominators, and model

| Panel | success | total |
|---|---|---|
| Fig. 3L | IEG-high nuclei of the Exc subtype in that sample | All Exc IEG-high nuclei in the same sample |
| Fig. S8B | IEG-high nuclei of the Inh subtype in that sample | All nuclei of the same Inh subtype in that sample |

Fig. 3L tests subtype composition within Exc IEG-high nuclei. Fig. S8B tests the IEG-high fraction within an Inh subtype. Generate sample-by-subtype counts accordingly. Flag and exclude zero denominators rather than treating their proportions as zero.

The model is:

```r
cbind(success, total - success) ~ condition + (1 | sample)
```

Use the binomial family and logit link, condition fixed effect, sample random intercept, and Ext reference. The coefficient represents Ret–Ext relative to Ext. Defaults are `nAGQ = 1`, optimizer `bobyqa`, and two-sided Wald z inference. Each row contains grouped Bernoulli outcomes from one sample. Disable `check.nobs.vs.nlev` and `check.nobs.vs.nRE` for this grouped-binomial representation without removing or reweighting the random intercept.

Export beta, SE, z, raw P, odds ratio, Wald 95% CI, random-intercept SD, convergence, and singularity flags. Singularity tolerance is 1e-4. Nonconverged fits do not report P values. Require at least two nonzero-denominator samples per condition and four rows overall. Do not fit when all successes or all failures are zero.

Apply BH separately to the two analyses, using six planned Exc tests and eight planned Inh tests. Keep these family sizes even if some subtypes are not estimable. The `analysis` identifiers are `excitatory_IEGhigh_composition` and `inhibitory_IEGhigh_fraction`.

### Outputs

Under `08_GLMM/`: `pool_counts.csv`, `model_statistics.csv`, `model_diagnostics.csv`, and `models.rds`. Plot sample points/ranges from pool_counts and annotate significance from the full GLMM BH results. Code and tables use analysis identifiers; figure numbers appear only in this documentation.

## 09 — `09_CellChat_analysis.R`

### Purpose and inputs

Read step 04 Exc/Inh IEG-high data and `config/rat_human_orthologs.tsv`. Required columns are `rat_symbol,human_symbol,orthology_type`. Accept `one_to_one` or `ortholog_one2one` only, with unique symbols on both sides. Changing symbol capitalization is not a substitute.

### Parameters and network scope

- Use ortholog-mapped, log-normalized RNA data, not Harmony coordinates or scaled residuals.
- Database: `CellChatDB.human`; grouping: `subtype`.
- Steps: subsetData → identifyOverExpressedGenes → identifyOverExpressedInteractions → computeCommunProb → filterCommunication → computeCommunProbPathway → aggregateNet.
- `computeCommunProb(type = triMean, raw.use = TRUE, population.size = TRUE, seed.use = 42)`.
- `filterCommunication(min.cells = 10)`.
- Scopes All, Excitatory, and Inhibitory; analyze Overall, WD, Ext, and Ret–Ext within each scope, giving 12 networks. All uses mapped human genes shared by both lineages.
- Compute `netP` communication centrality when pathways are inferred. Record warnings, errors, and centrality status for each network.
- Exclude autocrine edges from outgoing summaries and save autocrine information separately. Compute Ret–Ext-minus-Ext differences, fold changes, and percentage changes. Set fold/percentage changes to NA when the Ext baseline is zero.

### Outputs

Under `09_CellChat/`: mapping summaries and actual ortholog mappings; per-scope/condition objects, group sizes, nucleus metadata, ligand–receptor communication, pathway communication, outgoing pathways, and count/weight edges. Summary tables include `network_status.csv`, `outgoing_activity_self_excluded.csv`, `RetExt_vs_Ext_outgoing_comparison.csv`, and a GLU_L6_CLA/GLU_L5_6 subset.

These are predicted/inferred communication quantities based on expression and a database. Outgoing weights and condition differences are not measured projection strengths or independent sample-level condition tests.

## 10 — `10_plot_cell_types.R`

### Inputs, figure correspondence, and outputs

Read the complete major-class object from step 02, both neuronal objects from step 03, and `config/markers.json`. Save under `10_celltype_figures/`, organized into `cell_types/`, `QC/`, and `neurons/`. Filenames describe classes, metrics, or markers rather than manuscript panel identifiers. Export both PNG and PDF.

| Figure | Content | Canvas (inches) |
|---|---|---|
| Fig. 3B | Six major cell classes, UMAP | 8×7 |
| Fig. 3C | Major-class marker DotPlot | 12×5 |
| Fig. 3D | Ext/Ret–Ext major-class composition | 6×6 |
| Fig. S6A | Detected-gene violins by library | 10×5 |
| Fig. S6B | percent.ribo violins by library | 10×5 |
| Fig. S6C | UMI violins by library | 10×5 |
| Fig. S6D | Major-class marker FeaturePlots | 6×5.5 per gene |
| Fig. S6E/H | Exc/Inh cluster UMAPs | 7×6 |
| Fig. S6F/I | Exc/Inh subtype UMAPs | 8×6 |
| Fig. S6G/J | Exc/Inh marker DotPlots | 10×5 |

Also export subtype marker UMAPs, Exc/Inh coordinates, and major-class composition data. QC plots preferentially use `nFeature_observed`, `nCount_observed`, and `percent.ribo_observed` where available, otherwise the corresponding RNA QC fields. The first two use pre-ambient-correction counts when available. Composition percentages divide class counts by all class counts within each condition; they are not averages of sample percentages and have no additional condition-effect test.

Summary UMAPs use point size 0.15 and `raster = FALSE`. DotPlot expression is scaled per gene across groups and clipped to [-2,2], with dot.scale=8 and light-blue-to-red colors; size represents the detected fraction. FeaturePlot colors run from `#EEEEEE` to `#753B8E`. Use Arial, base size 14, and 600-dpi PNG.

## 11 — `11_plot_expression.R`

### Inputs and volcanoes

Read step 05 all-genes tables, step 06 score tables, step 07 enrichment tables, the step 04 Exc object, and `config/heatmap_genes.csv`.

| Figure | Comparison / neuronal class |
|---|---|
| Fig. 3H | Exc Ret–Ext vs Ext |
| Fig. S7A | Inh Ret–Ext vs Ext |
| Fig. S7B | Exc Ext vs WD |
| Fig. S7C | Exc Ret–Ext vs WD |
| Separate output | Inh Ext vs WD |
| Separate output | Inh Ret–Ext vs WD |

Volcano x is log2FoldChange; y is the actual `-log10(padj)`. Horizontal and vertical dashed thresholds are padj=0.05 and ±log2(1.15). Plot finite fold changes with valid adjusted P values. Color encodes fold change, clipped to [-2,2]. Point-size range is 0.7–3.2, alpha=0.72. Label up to 20 genes meeting both criteria, ordered by padj, decreasing absolute fold change, then gene symbol. Labels use Arial italic, size=4.5, repel seed=42. Default canvas: 9×7 inches.

### Scores, heatmap, and enrichment

- **Fig. 3F/G:** Ret–Ext/Ext marker-program distributions in the Exc scope, limited to Ext and Ret–Ext. Violin settings: `trim = FALSE`, alpha=0.8; overlaid boxplot width=0.12. The center is the median, the box spans the IQR, and whiskers follow the default 1.5×IQR rule rather than simple IQR error bars. Stars use raw two-sided nucleus-level Wilcoxon P values. Canvas: 5×6 inches. All/Inh scopes are also saved under `scores/`, with scope-specific filenames.
- **Fig. 3I:** Read 44 candidate genes in six configured functional modules. Retain genes present in the object and meeting padj≤0.05 in Exc Ret–Ext vs Ext; export membership. Average log-normalized expression within each condition's Exc IEG-high nuclei, then Z-standardize the two condition means per gene. Set zero-variance rows to zero. Keep `gene_order` and facet by module. Width=9 inches; height=`max(5, number_of_genes × 0.22)`. These are fixed functional groups, not newly discovered modules.
- **Fig. 3J and other enrichment outputs:** Fig. 3J corresponds to Exc Ret–Ext vs Ext. Inh and other contrasts are separate files under `enrichment/`, named by class and contrast. Select actual padj≤0.05 and term names matching `config/enrichment_terms.csv`, up to 12 terms per direction. Patterns cover addiction, synapse, neuron, neurotransmitter, endocannabinoid, glutamate, GABA, dopamine, cAMP, axon, dendrite, plasticity, vesicle, secretory, learning, and memory. Patterns affect display only, not tests or P values. The x-axis is direction-signed `-log10(padj)`; size=Count, color=padj, shape=database. Width=11 inches; height=`max(5, 2 + 0.32 × number_of_terms)`.

Outputs under `11_expression_figures/` are organized into `volcano/`, `enrichment/`, `scores/`, and `heatmap/`. Filenames use class, contrast, or program identifiers, not manuscript panel identifiers. Also save enrichment plot data, heatmap membership, and expression tables. If no genes/terms qualify, export an explanatory empty plot without substituting raw P for adjusted P.

## 12 — `12_plot_communication_proportions.R`

### Inputs and composition/proportion plots

Read step 04 objects, step 08 counts/GLMM statistics, and step 09 CellChat summaries. Save under `12_communication_figures/`, organized into `proportions/` and `communication/`, using descriptive analysis filenames rather than panel identifiers.

- **Fig. 3K / Fig. S8A:** Stacked Ext/Ret–Ext IEG-high subtype composition. Divide each subtype's IEG-high count by the corresponding lineage's total IEG-high count within condition. Canvas: 7×6 inches.
- **Fig. 3L / Fig. S8B:** Each point is an observed pooled-sample proportion; vertical lines span the minimum to maximum of the three samples. Point size=2.7; horizontal jitter width=0.065, seed=42; facet by subtype, three columns. Stars come from step 08 BH-adjusted P values within each analysis. Plotting does not change the model. Canvases: 10×7 and 10×9 inches, respectively, with each analysis's specified denominator.

### CellChat summary plots

| Figure | Content and display rules | Canvas (inches) |
|---|---|---|
| Fig. 3P | GLU_L6_CLA/GLU_L5_6 outgoing-weight percentage changes, Ret–Ext relative to Ext; All/Excitatory facets | 9×5 |
| Fig. S9A | Exc subtype outgoing weights in Ext/Ret–Ext; All/Excitatory scopes | 12×6 |
| Fig. S9B | Directed Exc differential network, excluding self-edges and zero differences | 8×8 |
| Fig. S9C | Increased GLU_L6_CLA outgoing pathways; up to 20 ranked differences | Width 11; height `max(5, 2 + 0.3 × number_of_entries)` |
| Fig. S9D | Decreased GLU_L5/GLU_L5_6 outgoing pathways; up to 20 ranked decreases | Same as above |
| Fig. S9E | Up to 20 pathways with largest absolute changes for these three subtypes; `1000 × difference` | 13×6 |

Differential networks use a circular layout, node size=22, label cex=0.9, and arrow size=0.35. Edge width is `0.4 + 4 × abs(difference) / max(abs(difference))`. Red indicates Ret–Ext increases and blue indicates decreases. Pathway summaries exclude source=target autocrine edges. An edge absent from one condition's network is treated as zero in differences; percentage changes are not constructed from zero baselines.

Also export composition data, network edges, and pathway-difference tables. Communication bars and networks display inferred quantities and differences, without sample-level significance stars for these aggregate quantities.

## `apply_annotations.R` — Apply cell-class or subtype mappings

Read existing clustered objects and configured mappings without repeating QC, PCA, Harmony, graph construction, clustering, or UMAP. With `--level major`, read `02_annotation/all_nuclei.rds` and the major_annotations table, require a cell_type label for every actual cluster, and save `02_annotation/annotated.rds` and `major_annotations_applied.csv`.

With `--level subtype`, read `03_neurons/Excitatory.rds`, `Inhibitory.rds`, and the subtype_annotations table. Both lineages must be completely mapped to the specified six Exc or eight Inh subtype names. Validate both objects before saving `Excitatory_annotated.rds` and `Inhibitory_annotated.rds`. Duplicate/absent cluster identifiers, empty labels, and Unknown assignments are not accepted as complete annotation. Stop if a target RDS already exists.

## Helper files

### `R/common.R`

Shared configuration and I/O utilities. `pipeline_context()` reads config/output/level arguments and JSON settings, then sets seed=42. `new_stage()` prevents overwriting an existing stage. `write_table()` writes CSV without row names and represents missing values as `NA`. `save_object()` writes gzip-compressed RDS and refuses overwriting. `write_session()` saves sessionInfo.

`read_samples()` requires nine unique samples, three per condition, without missing sample, condition, filtered_h5, or raw_h5 entries. `neuron_object_path()` prefers `03_neurons/{type}_annotated.rds`, falling back to `{type}.rds`. `join_rna_layers()` supports Seurat v5 counts layers; `rna_matrix()` supports Assay5 layers and older Assay slots. `integer_counts()` checks nonnegative integers with tolerance 1e-8. `require_neuron_metadata()` requires complete sample/condition/subtype metadata without Unknown. `bh()` applies BH adjustment. `plot_export()` saves white-background PNG and Cairo PDF, default 600 dpi, inch units, `limitsize = FALSE`.

### `R/preprocessing.R`

HDF5 import, version requirements, PC selection, sample Harmony, annotation mapping, and diagnostic marker plots. SoupX/scDblFinder versions must exactly match configuration. HDF5 import uses gene symbols and unique feature names. PC counts are limited to nucleus count−1 and variable-feature count−1.

`run_sample_harmony()` fixes `group.by.vars = sample`, `reduction.use = pca`, `project.dim = FALSE`, and `reduction.save = harmony`, using PCA embeddings from the RNA workflow. `annotation_mapping()` rejects duplicate and absent cluster identifiers. Unfilled labels remain Unknown rather than being assigned to neighboring types.

Marker panels cover the six major classes, six Exc subtypes, and eight Inh subtypes. Diagnostic DotPlots use dot.scale=10, per-gene scaling, lightgrey/blue/red colors, width `max(10, number_of_markers × 0.45)`, and height 7 inches. Marker FeaturePlots use UMAP, point size=0.3, `order = TRUE`, min.cutoff=q10, max.cutoff=q90, and an 8×7-inch canvas. These diagnostic displays are configured separately from summary plots in `R/plotting.R` and do not alter expression or clustering.

### `R/plotting.R`

Shared themes, colors, significance stars, summary UMAPs, DotPlots, volcanoes, and enrichment plots. Font=Arial, base size=14. Condition colors are WD `#9E9E9E`, Ext `#49B8D1`, Ret–Ext `#EF4936`. Named vectors fix class/subtype colors across plots.

| Groups, in order | Colors, in the same order |
|---|---|
| Excitatory, Inhibitory, Oligodendrocytes, OPC, Astrocytes, Microglia | `#3B6FBF`, `#E85D5D`, `#B5895A`, `#9A6FD1`, `#4DB894`, `#F0986D` |
| GLU_L6_CLA, GLU_L5_6, GLU_L2_L3, GLU_L5, GLU_L4_5, GLU_ADGM | `#8DD3C7`, `#FFFFB3`, `#BEBADA`, `#FB8072`, `#80B1D3`, `#FDB462` |
| GABA_Amb, GABA_Etv1, GABA_Lamp5, GABA_MSN, GABA_Sst, GABA_Hs3st2, GABA_Pvalb, GABA_Vip | `#66C2A5`, `#FC8D62`, `#8DA0CB`, `#E78AC3`, `#A6D854`, `#FFD92F`, `#E5C494`, `#FF6347` |

`stars()` uses strict inequalities: P<0.0001 gives ****, P<0.001 ***, P<0.01 **, P<0.05 *, otherwise ns; missing values give NA. Each plot specifies the P-value type supplied to this function. `collect_enrichment()` reads complete term tables, excludes duplicate `*_significant.csv` files, and retains actual padj before display filtering. Term names do not receive database suffixes; database information remains in data columns and shape encoding.

### `R/communication.R`

Strict one-to-one rat–human mapping, IEG-high expression inputs, merging of shared mapped features across lineages, network long tables, and outgoing-condition comparisons. Accept only one_to_one/ortholog_one2one mappings with unique symbols on both sides; capitalization conversion is not ortholog mapping.

`cellchat_input()` checks finite, nonnegative RNA data, unique nucleus names, valid conditions, and retained sample/subtype metadata. `merge_cellchat_inputs()` requires nonoverlapping nucleus identifiers across lineages. `outgoing_network_summary()` saves cross-subtype outgoing and autocrine weights/counts separately. `compare_outgoing_conditions()` calculates Ret–Ext−Ext, Ret–Ext/Ext, and `100 × (Ret–Ext−Ext) / Ext`. For zero baselines, mark zero_baseline and report fold/percentage changes as NA.

## Parameters and interpretation

Configuration contains both analysis thresholds and implementation details. Count prefilters, quantile type, and GLMM optimizer/integration settings define computational behavior; they are not biological findings, historical counts, or guarantees that a particular subtype is significant. Retain configuration, inputs, and outputs when changing parameters so each figure remains linked to its actual statistical table.

The units differ across analyses: DESeq2 uses pooled-sample counts; GLMM sample random intercepts account for within-sample nuclear correlation; marker-program scores describe nucleus-level distributions; CellChat provides expression-supported communication inference. Raw, BH-adjusted, and Bonferroni-adjusted P values are saved and labeled according to their own analysis.
