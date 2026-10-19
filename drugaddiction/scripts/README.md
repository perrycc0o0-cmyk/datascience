# 脚本说明

本目录按分析顺序组织大鼠 IL 单核 RNA 测序的预处理、注释、IEG-high 定义、差异表达、marker-program 评分、通路富集、亚群比例模型、CellChat 和绘图。参数集中在 `../config/analysis.json`；marker、亚群注释、热图基因和通路显示规则分别位于相应配置表中。

每个条件包含 3 个独立 pooled samples，每个 sample 由 5 只大鼠的 IL 组织混合构成。条件差异表达以 sample 为重复单位；单核 marker-program 分布描述与样本层面的条件效应推断采用不同统计步骤。

## 使用方式与文件关系

从项目根目录调用脚本，路径均相对于项目根目录，而非当前 `scripts/` 目录。例如：

```bash
Rscript scripts/01_preprocess_libraries.R --config config/analysis.json --output results
Rscript scripts/02_integrate_annotate.R --config config/analysis.json --output results
```

所有编号脚本接受相同的两个可选参数；另外的`apply_annotations.R`使用`--level`选择注释层级：

| 参数 | 默认值 | 用途 |
|---|---|---|
| `--config` | `config/analysis.json` | 分析配置文件 |
| `--output` | 配置中的 `paths.output`，默认 `results` | 各分析阶段共用的输出根目录 |
| `--level` | `major` | `apply_annotations.R`的注释层级，可选`major`或`subtype` |

全局随机种子为 42。输出按阶段分目录；如果目标阶段目录已有内容，脚本停止，不覆盖已有对象或图表。重新选择配置或注释时，应使用新的输出根目录，并保留该次分析所需的上游阶段文件。

分析顺序为：

```text
01 文库预处理 → 02 全局整合和大类注释 → 03 神经元重聚类
                                             ↓
                                       04 IEG-high
                                      ↙   ↓    ↘
                                 05 DE   06评分  08比例GLMM
                                   ↓             09 CellChat
                                 07富集

10 绘制大类、QC和亚群图
11 绘制DE、评分、六模块热图和富集图
12 绘制比例和通信图
```

`config/samples.csv` 必须提供 9 个唯一 sample，条件标签为 `WD`、`Ext`、`Ret-Ext`，每条件 3 个。每行指定对应的 filtered 与 raw 10x HDF5 文件。大类映射读取 `config/major_annotations.csv`，亚群映射读取 `config/subtype_annotations.csv`。注释表的 cluster 必须对应相同分析对象中的实际 cluster 编号；未分配的群保留为 `Unknown`，不能进入要求完整注释的下游步骤。

完成聚类和marker检查后，可以通过下列命令应用人工确认的映射，无需重新计算聚类：

```bash
Rscript scripts/apply_annotations.R --level major --output results
Rscript scripts/apply_annotations.R --level subtype --output results
```

major层级在步骤02之后应用，生成`02_annotation/annotated.rds`；subtype层级在步骤03之后应用，生成`03_neurons/Excitatory_annotated.rds`和`Inhibitory_annotated.rds`。读取神经元对象的下游步骤优先采用`*_annotated.rds`。若步骤02已直接生成完整注释对象，无需重复应用major注释；已有目标对象不会被覆盖。

## 脚本与图表对应

| 文件 | 主要分析 | 对应内容 |
|---|---|---|
| `01_preprocess_libraries.R` | 固定 QC、SoupX、scDblFinder | Methods：数据预处理、ambient RNA 与 doublet；QC数据 |
| `02_integrate_annotate.R` | 全局 PCA、sample Harmony、聚类、大类注释 | Fig. 3B–D；Methods：整合、聚类与注释 |
| `03_neuron_subclustering.R` | Exc／Inh独立重聚类与亚群注释 | Fig. S6E–J；Methods：神经元亚群 |
| `04_ieg_definition.R` | Exc／Inh各自95%分位数定义IEG-high | Fig. 3E及后续IEG-high分析的数据基础 |
| `05_pseudobulk_DE.R` | sample层面DESeq2差异表达 | Fig. 3H、I及Fig. S7火山图的数据基础 |
| `06_marker_program_scores.R` | 两套marker-program基因集及核级评分 | Fig. 3F、G及对应补充评分 |
| `07_GO_KEGG_enrichment.R` | 上／下调基因GO和KEGG富集 | Fig. 3J、Fig. S7F及其他条件富集 |
| `08_subtype_GLMM.R` | condition固定效应＋sample随机截距二项模型 | Fig. 3L、Fig. S8B |
| `09_CellChat_analysis.R` | IEG-high神经元推断通信与输出汇总 | CellChat相关主图及Fig. S9的数据基础 |
| `10_plot_cell_types.R` | 大类、QC、亚群UMAP与marker图 | Fig. 3B–D、Fig. S6A–J |
| `11_plot_expression.R` | 评分、火山图、模块热图与富集 | Fig. 3F–J、Fig. S7A–F |
| `12_plot_communication_proportions.R` | 组成、GLMM比例与通信汇总图 | Fig. 3K、L、P，Fig. S8A、B，Fig. S9A–E |
| `apply_annotations.R` | 将确认后的cluster映射应用到已有对象 | 大类和亚群注释，不重新计算降维或聚类 |

Fig. 3A的实验设计示意和M–O的湿实验图像／动物级结果不由这些snRNA-seq脚本生成。

## 01 — `01_preprocess_libraries.R`

### 目的和输入

分别处理 9 个测序文库。输入为 `config/samples.csv` 指向的 `filtered_feature_bc_matrix.h5` 与 `raw_feature_bc_matrix.h5`。raw 矩阵用于估计 ambient RNA；filtered 矩阵提供初始 called nuclei。读取多模态 HDF5 时仅使用 `Gene Expression`，feature 和 barcode 必须唯一。

### 处理步骤和参数

1. 创建 Seurat 对象时保留在至少 3 个细胞核中检测到的基因，初始细胞核至少检测到 300 个基因。
2. 在 observed counts 上计算 detected genes、UMIs、线粒体比例、核糖体比例和 `log10(nFeature_RNA) / log10(nCount_RNA)`。
3. 固定 QC 使用严格不等号：`700 < nFeature_RNA < 5000`、`500 < nCount_RNA < 30000`、`percent.mt < 5`、复杂度比值 `> 0.75`。线粒体匹配模式为 `^Mt-`，核糖体匹配模式为 `^Rp[sl]`；参考基因符号必须与模式一致，找不到任何线粒体基因时停止。
4. SoupX 1.6.2 在每个文库内估计和校正污染。用于 SoupX 的预聚类为 LogNormalize 10,000、vst 高变基因 2,000、ScaleData、最多 30 PCs、PCA邻接图、resolution 0.5；此预聚类不作为最终生物学注释。`autoEstCont(forceAccept = FALSE)`，`adjustCounts(roundToInt = TRUE)`。
5. 在通过固定 QC 的核上运行 scDblFinder 1.24.10，`clusters = FALSE`。每文库预期 doublet rate 为 `0.004 × original_called / 1000`；各文库随机种子均为42。仅保留判定为 singlet 的核。
6. 最终 `RNA` assay 保存校正后的整数 counts；`observed` assay保留相同核的未校正 counts。barcode 加 sample 前缀，元数据保留 sample、condition、原始 QC 和 doublet 信息。

### 输出

输出目录 `01_preprocessing/`：

- `libraries/{sample}.rds`：每文库保留核对象。
- `merged.rds`：合并后的对象，供步骤02读取。
- `QC/{sample}_nuclei.csv`：逐核原始QC、各过滤条件、污染比例和doublet结果。
- `QC/{sample}_thresholds.png/.pdf`：4个QC指标的直方图，60 bins，虚线标示固定阈值；画布10×7英寸。
- `ambient_models/{sample}.rds`：SoupX污染拟合信息。
- `library_summary.csv`：初始核数、各步骤去除数量、保留singlets和污染比例汇总。
- `objects.csv`、`sessionInfo.txt`：对象路径与软件环境。

保留数量由实际输入和各步骤结果产生，不作为筛选参数。

## 02 — `02_integrate_annotate.R`

### 目的和输入

读取 `01_preprocessing/merged.rds`，完成全细胞类型整合、聚类、UMAP和大类marker检查。大类注释来自 `config/major_annotations.csv`，不通过固定核数确定。

### 处理步骤和参数

- RNA LogNormalize，scale factor 10,000。
- vst高变基因2,000；在这些基因上ScaleData，并回归 `nCount_RNA` 和 `percent.mt`。
- PCA目标50 PCs，受细胞核数和可用高变基因数限制。
- PC选择：以PC标准差占标准差总和的百分比计算贡献；候选一为累计贡献>90%且当前PC贡献<5%的第一个PC，候选二为相邻贡献下降>0.1个百分点的最后位置再加1；结合最大30 PCs及实际可用PC数，选择候选中的最小值，至少2 PCs。两个候选都保留在输出表中。
- Harmony仅以 `sample` 为分组变量，使用所选PC；不将 `condition` 加入batch变量。`project.dim = FALSE`。
- Harmony邻接图 `k.param = 20`；保存resolution 0.1、0.2、0.4、0.6、0.8、1.0的聚类列，最终身份为resolution 1.0。
- UMAP使用相同Harmony维度，`n.neighbors = 30`、`min.dist = 0.3`、`metric = cosine`、seed 42。
- `FindAllMarkers`：仅正marker、`min.pct = 0.25`、`logfc.threshold = 0.25`、每群最多1,000个核、seed 42。

### 输出和注释

输出目录 `02_annotation/`。主要文件为 `all_nuclei.rds`、`cluster_markers.csv`、`PCA_contributions.csv`、`PCA_selection.csv`、`UMAP_coordinates.csv`、`major_annotations.csv`和注释模板。所有大类完成注释后另外保存 `annotated.rds`；步骤03和正式大类图读取这一对象。

绘制PCA elbow、cluster UMAP及按cluster的marker DotPlot／FeaturePlot；完整注释时再绘制按cell_type的图。大类marker为：

| 类别 | Marker |
|---|---|
| Excitatory | Slc17a7、Satb2、Slc30a3 |
| Inhibitory | Gad1、Gad2、Dlx1 |
| Oligodendrocytes | Mog、Mbp、St18 |
| OPC | Cacng4、Pdgfra、Vcan |
| Astrocytes | Slco1c1、Gfap、Aqp4 |
| Microglia | P2ry12、Tmem119、Cx3cr1 |

本阶段诊断UMAP为9×7英寸、point size 0.3，PCA elbow为8×6英寸。marker图参数由 `R/preprocessing.R` 提供；论文汇总版图由步骤10绘制。

## 03 — `03_neuron_subclustering.R`

### 目的和输入

读取完整大类注释对象 `02_annotation/annotated.rds`，分别取Excitatory和Inhibitory重新聚类。每类单独进行归一化、高变基因选择、PCA、Harmony和邻接图计算。亚群映射读取 `config/subtype_annotations.csv`，字段为 `neuron_type,cluster,subtype`。

### 参数

| 步骤 | Excitatory | Inhibitory |
|---|---|---|
| LogNormalize scale factor | 10,000 | 10,000 |
| vst高变基因 | 2,000 | 2,000 |
| ScaleData | 不额外回归变量 | 不额外回归变量 |
| PCA目标维度 | 100，受可用数据限制 | 100，受可用数据限制 |
| Harmony | 前30个可用PC，sample-only | 前30个可用PC，sample-only |
| 建图和UMAP维度 | 前19个可用Harmony维度 | 前19个可用Harmony维度 |
| `k.param` | 20 | 30 |
| 最终resolution | 0.5 | 0.5 |
| UMAP | neighbors=10，min.dist=0.8，spread=6，repulsion.strength=1 | Seurat `RunUMAP`缺省neighbors=30、min.dist=0.3、spread=1、repulsion.strength=1 |
| 随机种子 | 42 | 42 |

`FindAllMarkers`使用与步骤02相同的0.25表达比例阈值、0.25 logFC阈值、only.pos及每群最多1,000核设置。

### Marker和输出

Exc亚群marker：GLU_L6_CLA—Gng2；GLU_L5_6—Syt6、Foxp2、Tle4；GLU_L2_L3—Glis3；GLU_L5—Parm1；GLU_L4_5—Il1rapl2、Rorb；GLU_ADGM—Rerg、Arhgap6。

Inh亚群marker：GABA_Amb—Slc17a6；GABA_Etv1—Etv1；GABA_Lamp5—Lamp5；GABA_MSN—Gng7；GABA_Sst—Sst；GABA_Hs3st2—Hs3st2；GABA_Pvalb—Slit2；GABA_Vip—Vip。

输出目录 `03_neurons/`。`Excitatory.rds`和`Inhibitory.rds`保存完整谱系核；各谱系子目录保存cluster marker、注释、逐核metadata、参数和UMAP坐标，以及cluster／subtype UMAP和marker图。诊断UMAP画布8×7英寸、point size 0.5。`subtype_annotations_template.csv`提供需填写的映射结构。完整亚群注释后才能继续要求subtype的下游统计。

## 04 — `04_ieg_definition.R`

### 目的和输入

分别读取 `03_neurons/`的Excitatory与Inhibitory对象，优先使用通过`apply_annotations.R`保存的`*_annotated.rds`；在各神经元大类内部定义IEG-high，不能将两大类混在一起计算一个分位数。

### 参数和定义

- 合并RNA layers，检查counts为非负整数，再进行LogNormalize 10,000。
- 预定义八个IEG：`Arc, Bdnf, Btg2, Fos, Fosl2, Homer1, Npas4, Nr4a1`；使用其中在RNA矩阵内可用的基因。逐谱系输出availability；全部缺失时停止。
- `ieg_sum`为同一核可用预定义IEG的log-normalized表达值之和，不是原始UMI之和，也不是z-score。
- 在Exc或Inh各自内部，跨WD、Ext和Ret–Ext共同计算第95百分位，`quantile(type = 7)`。
- `ieg_sum >= cutoff`记为IEG-high。边界存在并列值时，纳入比例可略大于5%。不按每个sample、条件或亚群重新定义阈值。

输出目录 `04_ieg/`：两类完整对象、逐谱系`*_IEG_availability.csv`、`nucleus_metadata.csv`、`IEG_cutoffs.csv`和`IEG_counts.csv`。逐核metadata包含Cell、sample、condition、neuron_type、subtype、ieg_sum、IEGhigh和IEGclass。cutoff表同时记录实际使用的IEG数量和名称。阈值由数据计算，不写死历史阈值或目标核数。

## 05 — `05_pseudobulk_DE.R`

### 目的和输入

读取步骤04的两类对象，仅使用IEG-high核。将RNA counts在每个神经元大类内按sample求和，形成gene×9个pooled samples的整数矩阵。每个sample必须实际包含该类IEG-high核；不会用全零列伪造缺失的重复。

### 模型和参数

- count过滤：全9个sample总counts≥10，且至少3个sample的该基因counts>0。可通过 `de.min_total_counts` 和 `de.min_samples` 调整。
- DESeq2设计为 `~ condition`，condition水平为WD、Ext、Ret–Ext；使用全9个sample估计模型。
- size-factor方法为`ratio`；dispersion fit为`parametric`；Wald检验；`parallel = FALSE`。
- `minReplicatesForReplace = Inf`，不自动替换outlier counts；其余DESeq2结果计算保留函数的默认Cook's distance处理。
- `results(alpha = 0.05, pAdjustMethod = BH, independentFiltering = TRUE)`。BH在每个神经元大类、每个contrast内部进行；独立过滤的阈值及分位数另表记录。
- 对比：`Ext_vs_WD`、`RetExt_vs_WD`、`RetExt_vs_Ext`。log2FoldChange>0表示前者相对后者上调。
- 最终DEG标志为padj≤0.05且`abs(log2FoldChange) >= log2(1.15)`，约0.2016；此幅度阈值用于筛选输出，不是Wald检验的复合零假设阈值。

### 输出

输出目录 `05_DESeq2/`：

- `counts/{neuron_type}_pseudobulk_counts.csv`及sample metadata。
- `objects/{neuron_type}_DESeq2.rds`：拟合后的DESeq2对象。
- `tables/{neuron_type}/{contrast}_all_genes.csv`：所有输入基因；包括tested、total_count、n_samples_expressed、baseMean、log2FoldChange、lfcSE、stat、pvalue、padj、contrast、significant和direction。
- 同目录的`*_significant.csv`：达到padj和幅度双阈值的基因。
- `DE_summary.csv`、`filtering_parameters.csv`：实际可检验数量、上下调数量和独立过滤参数。

未通过预过滤的基因保留在all-genes表中，检验字段为NA。通过counts过滤但被独立过滤的基因也可能没有padj；原始P和BH-adjusted P分别保存，不相互替换。

## 06 — `06_marker_program_scores.R`

### 目的和输入

读取步骤04的Exc和Inh对象，合并用于marker-program基因集发现，但保留已经按各大类定义的IEG-high／low标签。gene-set discovery使用Ext或Ret–Ext各自条件内部的核，不是条件间DESeq2比较。

### 基因集和评分参数

1. 分别在Ext和Ret–Ext中，将IEG-high与IEG-low进行双侧Seurat Wilcoxon比较；RNA data slot、`min.pct = 0.1`、`logfc.threshold = 0`、`only.pos = FALSE`、base=2、seed=42。
2. Seurat的`p_val_adj`为Bonferroni校正；保留`p_val_adj < 0.05`且`avg_log2FC > 0`的基因，移除预定义八个IEG。
3. 按校正P值、基因名依次排序，选择最多500个基因，得到`Ext_IEGhigh_UP`和`RetExt_IEGhigh_UP`。不足500时保留实际候选数量，不补入不满足条件的基因。
4. 相同两套signature分别应用于全部IEG-high、Exc IEG-high、Inh IEG-high三个对象；每个对象独立计算表达匹配背景。AddModuleScore使用RNA data slot、pool为该对象全部RNA基因、ctrl=100、nbin=24、seed=42、search=FALSE。
5. 输出全部三个条件的IEG-high评分；绘图和核级评分分布比较仅使用Ext与Ret–Ext。Wilcoxon为双侧、非配对、`exact = FALSE`、continuity correction开启；不将核数解释为生物学重复数。

### 输出

输出目录 `06_scores/`：`signature_DE/`保存组内marker结果，`gene_sets/`保存基因集CSV和TXT，`signature_summary.csv`记录候选数和所选数。`scores_long.csv`逐核保存score，并含scope（All、Excitatory、Inhibitory）、program、sample、condition和subtype。`score_statistics.csv`记录核级Wilcoxon的P值、核数、中位数和W统计量。`Prism/`为Ext／Ret–Ext双列评分表，列长不同时用NA补齐排版，不添加观测。

Fig. 3F、G使用Excitatory scope的Ret–Ext与Ext programs；另外两个scope保留为对应补充数据。评分用于描述核的marker-program表达分布，不替代步骤05的sample层面条件DE。

## 07 — `07_GO_KEGG_enrichment.R`

### 目的和输入

读取步骤05两类神经元的三个contrast all-genes表，以padj≤0.05筛选，再按log2FoldChange正负分成Up和Down。富集前景不额外要求1.15倍幅度阈值；该幅度阈值用于火山图DEG标志，不能直接替代富集输入定义。对每个方向分别进行GO BP、CC、MF及大鼠KEGG过度代表分析。

### 参数和背景

- GO使用`org.Rn.eg.db`、SYMBOL键；ontology分别为BP、CC、MF，不合并检验族，`pool = FALSE`。
- KEGG通过SYMBOL→ENTREZID映射，organism=`rno`、keyType=`kegg`，使用当前KEGG数据；`use_internal_data = FALSE`。KEGGREST的rno gene–pathway对应表同时保存。
- 背景首先限定为该类、该contrast通过counts预过滤的全部基因，而非全基因组；再分别与GO ontology或KEGG可注释成员取交集。
- 多重校正为BH；`pvalueCutoff = 0.05`、`qvalueCutoff = 0.2`、`minGSSize = 5`、`maxGSSize = 500`；最终显著表另要求p.adjust≤0.05。GO各ontology、方向、神经元类和contrast分别计算其富集检验族，KEGG同样单独计算。
- 一对多SYMBOL→ENTREZID映射保留并导出，分析时ENTREZID去重。不可映射基因不加入对应数据库的输入或背景。

### 输出

输出目录 `07_enrichment/`。`tables/{neuron_type}/{contrast}/`下保存：

- `Up_input_genes.csv`、`Down_input_genes.csv`。
- `tested_gene_universe.csv`、SYMBOL／ENTREZID映射、GO BP／CC／MF及KEGG各自的mapped universe。
- `Up_GO_BP.csv`等方向×数据库表，保存完整返回检验结果；`*_significant.csv`保存达到配置阈值的结果。
- 富集表保留Description、pvalue、p.adjust、qvalue、Count、geneID、direction和database；对象保存在`objects/`。

仅当存在符合输入条件的基因和背景时调用对应富集分析；无结果时输出具有固定字段的空表，不以邻近通路或不显著通路补齐。

## 08 — `08_subtype_GLMM.R`

### 目的和输入

读取步骤04逐核metadata和sample表，比较Ext与Ret–Ext。6个pooled samples为sample随机截距分组，不对sample贡献设定额外小权重，不采用部分方差调整。

### 分子、分母与模型

| Panel | success | total |
|---|---|---|
| Fig. 3L | 该sample、该Exc亚群中的IEG-high核数 | 同一sample全部Exc IEG-high核数 |
| Fig. S8B | 该sample、该Inh亚群中的IEG-high核数 | 同一sample该Inh亚群的全部核数 |

两者分母不同：Fig. 3L检验IEG-high Exc集合内部的亚群组成，Fig. S8B检验同一Inh亚群的IEG-high比例。按上述定义分别生成sample×subtype计数，零分母标记并排除，不将其比例当作0。

模型为：

```r
cbind(success, total - success) ~ condition + (1 | sample)
```

默认参数：binomial family、logit link、condition固定效应、sample随机截距；Ext为参照，系数为Ret–Ext相对Ext；`nAGQ = 1`、optimizer=`bobyqa`、双侧Wald z检验。每一行是同一sample内Bernoulli核的分组二项计数；针对这种分组数据关闭`check.nobs.vs.nlev`和`check.nobs.vs.nRE`检查，不改变模型的随机截距。

输出beta、SE、z、原始P、odds ratio及Wald 95% CI、随机截距SD、收敛状态和singular标志。singular判定容差为1e-4；非收敛模型不报告P。每组至少2个非零分母sample且共至少4行，success或failure全为0时不拟合。

BH按分析类型分别校正，Fig. 3L对应的Exc组成分析计划检验6个亚群，Fig. S8B对应的Inh比例分析计划检验8个亚群；即使部分亚群不可估计，BH的计划检验数仍保持6或8，不因缺失结果缩小。表格的`analysis`字段分别为`excitatory_IEGhigh_composition`和`inhibitory_IEGhigh_fraction`。

### 输出

输出目录 `08_GLMM/`：`pool_counts.csv`、`model_statistics.csv`、`model_diagnostics.csv`和`models.rds`。图中sample点和范围线读取pool_counts；显著性读取完整GLMM的BH结果。代码和表格用分析名称作为标识，图号仅用于本文说明中的对应关系。

## 09 — `09_CellChat_analysis.R`

### 目的和输入

读取步骤04的Exc／Inh IEG-high数据和 `config/rat_human_orthologs.tsv`。该表须包含`rat_symbol,human_symbol,orthology_type`，仅接受`one_to_one`或`ortholog_one2one`映射，rat和human symbol两端均必须唯一。仅改变大小写不能替代此映射表。

### 参数和网络范围

- 使用映射后的RNA log-normalized data，而不是Harmony坐标或scaled residuals。
- 数据库为`CellChatDB.human`；group.by为`subtype`。
- 标准步骤：subsetData → identifyOverExpressedGenes → identifyOverExpressedInteractions → computeCommunProb → filterCommunication → computeCommunProbPathway → aggregateNet。
- `computeCommunProb(type = triMean, raw.use = TRUE, population.size = TRUE, seed.use = 42)`。
- `filterCommunication(min.cells = 10)`。
- 三个scope：All、Excitatory、Inhibitory；每个scope分别分析Overall、WD、Ext、Ret–Ext，共12个网络。All scope使用两类共同映射到的人类基因。
- 存在推断通路时计算`netP`的通信中心性；单个网络的warning、error和centrality状态分别记录。
- outgoing汇总排除自分泌边；autocrine信息另外单独保存。Ret–Ext相对Ext计算差值、fold change及百分比变化；Ext基线为0时fold change／百分比变化记为NA。

### 输出

输出目录 `09_CellChat/`：mapping summary和具体同源映射；每个scope／condition子目录中的对象、group sizes、核metadata、配体受体通信、pathway communication、outgoing pathways以及count／weight edges。总表包括`network_status.csv`、`outgoing_activity_self_excluded.csv`、`RetExt_vs_Ext_outgoing_comparison.csv`和针对GLU_L6_CLA／GLU_L5_6的子表。

这些结果是基于表达和数据库的predicted／inferred communication；outgoing权重及其条件差值不是实测投射强度，也不是sample层面的独立重复检验。

## 10 — `10_plot_cell_types.R`

### 输入、图号和输出

读取步骤02完整大类对象、步骤03两个亚群对象，以及 `config/markers.json`。输出到`10_celltype_figures/`，按`cell_types/`、`QC/`和`neurons/`组织，文件名使用具体细胞类别、指标或marker名称，不含论文图号。所有图同时生成PNG和PDF。

| 图 | 内容 | 画布（英寸） |
|---|---|---|
| Fig. 3B | 六大类UMAP | 8×7 |
| Fig. 3C | 六大类marker DotPlot | 12×5 |
| Fig. 3D | Ext／Ret–Ext大类组成 | 6×6 |
| Fig. S6A | 每文库detected genes小提琴图 | 10×5 |
| Fig. S6B | 每文库percent.ribo小提琴图 | 10×5 |
| Fig. S6C | 每文库UMIs小提琴图 | 10×5 |
| Fig. S6D | 大类marker FeaturePlots | 每基因6×5.5 |
| Fig. S6E、H | Exc、Inh未注释cluster UMAP | 7×6 |
| Fig. S6F、I | Exc、Inh亚群UMAP | 8×6 |
| Fig. S6G、J | Exc、Inh marker DotPlot | 10×5 |

另输出各亚群marker UMAP、Exc／Inh坐标及大类组成数据。大类QC优先使用`nFeature_observed`、`nCount_observed`和可用的`percent.ribo_observed`字段，分别回退至对应RNA QC字段；前两项可用时绘制未污染校正的原始计数。组成百分比是各条件内该大类核数除以该条件所有大类核数，不对sample比例取平均，不附加条件效应检验。

论文UMAP使用point size 0.15、`raster = FALSE`。DotPlot将每基因跨群相对表达裁剪至[-2,2]，dot.scale=8，颜色由浅蓝到红色；点大小表示检测比例。FeaturePlot颜色为`#EEEEEE`到`#753B8E`。统一字体Arial、base size 14、PNG 600 dpi。

## 11 — `11_plot_expression.R`

### 输入与火山图

读取步骤05 all-genes表、步骤06评分表、步骤07富集表、步骤04Exc对象，以及`config/heatmap_genes.csv`。

| 图 | 条件／神经元类别 |
|---|---|
| Fig. 3H | Exc Ret–Ext vs Ext火山图 |
| Fig. S7A | Inh Ret–Ext vs Ext火山图 |
| Fig. S7B | Exc Ext vs WD火山图 |
| Fig. S7C | Exc Ret–Ext vs WD火山图 |
| Fig. S7D | Inh Ext vs WD火山图 |
| Fig. S7E | Inh Ret–Ext vs WD火山图 |

火山图x轴为log2FoldChange，y轴为真实`-log10(padj)`；水平虚线为padj=0.05，垂直虚线为±log2(1.15)。只绘制有限FC和非缺失padj。颜色按FC连续变化，范围[-2,2]，超界截断；点大小范围0.7–3.2，alpha=0.72。最多标注20个满足padj和FC双阈值的基因，依padj、绝对FC、基因名排序；标签Arial italic、size=4.5，repel seed=42。默认画布9×7英寸。

### 评分、热图与富集

- Fig. 3F／G：Exc scope的Ret–Ext／Ext marker-program分布。仅Ext、Ret–Ext，violin为`trim = FALSE`、alpha=0.8，叠加width=0.12的boxplot；双侧核级Wilcoxon原始P转换为星号；画布5×6英寸。All和Inh scope同样写入`scores/`，文件名以scope区分。
- Fig. 3I：读取配置的44个候选基因及六个功能模块；仅保留在Exc Ret–Ext vs Ext中padj≤0.05且对象内存在的基因，输出membership表。分别在两条件的Exc IEG-high核上计算平均log-normalized表达，再对每基因的这两个均值按行进行Z标准化；零方差记为0。按`gene_order`固定排序、按module分面；画布宽9，高为`max(5, 基因数×0.22)`。该步骤是固定功能分组展示，不进行新的模块发现。
- Fig. 3J、Fig. S7F：分别为Exc／Inh Ret–Ext vs Ext富集显示；其他条件另存`Other_enrichment/`。仅选真实padj≤0.05且名称匹配`config/enrichment_terms.csv`的通路；每个方向最多12条。默认关键词包括addiction、synapse、neuron、neurotransmitter、endocannabinoid、glutamate、GABA、dopamine、cAMP、axon、dendrite、plasticity、vesicle、secretory、learning和memory。关键词仅决定显示，不改变富集检验或P值。横轴为带Up／Down方向的`-log10(padj)`，点大小为Count、颜色为padj、形状为数据库。画布宽11，高为`max(5, 2 + 0.32×通路数)`。

输出目录 `11_expression_figures/`，按`volcano/`、`enrichment/`、`scores/`和`heatmap/`组织，文件名使用神经元类别、contrast和program等分析信息，不含论文图号。除PNG／PDF外，另保存富集plot data、热图membership和表达长表。无符合标准的基因或通路时输出说明空图，不使用原始P替代padj。

## 12 — `12_plot_communication_proportions.R`

### 输入和组成／比例图

读取步骤04的两类对象、步骤08计数和GLMM统计、步骤09CellChat汇总表。输出目录为`12_communication_figures/`，按`proportions/`和`communication/`组织，文件名采用分析名称，不含论文图号。

- Fig. 3K、Fig. S8A：Ext与Ret–Ext的IEG-high亚群组成堆叠柱；每条件内各亚群IEG-high核数除以对应谱系IEG-high总核数。画布7×6英寸。
- Fig. 3L、Fig. S8B：每个点为一个pooled sample的实际比例，竖线为三个sample的最小值到最大值，point size=2.7，水平jitter width=0.065、seed=42；按subtype分面，每行3列。星号来自步骤08各分析类型内的BH-adjusted P，不在绘图时更改统计模型。画布分别10×7和10×9英寸；两图采用各自定义的分母。

### CellChat汇总图

| 图 | 内容与显示规则 | 画布（英寸） |
|---|---|---|
| Fig. 3P | GLU_L6_CLA／GLU_L5_6的Ret–Ext相对Ext outgoing权重百分比变化；All与Excitatory scope分面 | 9×5 |
| Fig. S9A | Exc亚群Ext／Ret–Ext outgoing weight柱图；All／Excitatory scope | 12×6 |
| Fig. S9B | Exc scope条件差异有向网络，排除自连接和零差异边 | 8×8 |
| Fig. S9C | GLU_L6_CLA增加的outgoing pathways，按差值取最多20条 | 宽11，高`max(5, 2+0.3×条目数)` |
| Fig. S9D | GLU_L5／GLU_L5_6减少的outgoing pathways，按减少幅度取最多20条 | 同上 |
| Fig. S9E | 上述三个亚群绝对变化最大的最多20条通路矩阵，显示`1000×difference` | 13×6 |

差异网络以圆形布局绘制，节点大小22、标签cex=0.9、箭头size=0.35；边宽为`0.4 + 4×abs(difference)/最大abs(difference)`，Ret–Ext增加为红、减少为蓝。通路汇总均排除source=target的自分泌边；缺失的条件边按该网络中未检测到的0处理，不对零基线构造百分比。

该脚本另外输出组成数据、网络edges及pathway differences表。通信柱图和网络图仅显示推断量及其差值，不给这些聚合量添加sample层面的显著性星号。

## `apply_annotations.R` — 应用大类／亚群注释

该脚本读取已有聚类对象和配置的cluster映射，不重新做QC、PCA、Harmony、邻接图、聚类或UMAP。命令行`--level major`读取`02_annotation/all_nuclei.rds`及major_annotations配置表，要求每个实际cluster都有cell_type标签，保存`02_annotation/annotated.rds`及`major_annotations_applied.csv`。

`--level subtype`分别读取`03_neurons/Excitatory.rds`、`Inhibitory.rds`及subtype_annotations配置表。两谱系都必须完成映射，subtype名称须属于6个Exc或8个Inh的指定名称；两个对象确认完成后，再分别保存`Excitatory_annotated.rds`和`Inhibitory_annotated.rds`。cluster重复、cluster不属于该对象、空标签或Unknown均不允许作为完整注释进入后续流程。已有目标RDS时停止，避免覆盖人工确认结果。

## 辅助文件

### `R/common.R`

编号脚本共同依赖的输入输出工具。`pipeline_context()`读取config、output和level命令行选项、JSON配置并设置seed=42；`new_stage()`防止覆盖已有阶段；`write_table()`写CSV、NA写为`NA`、不导出行名；`save_object()`以gzip保存RDS且拒绝覆盖；`write_session()`保存sessionInfo。

`read_samples()`要求9个唯一sample、每条件3个，并要求sample、condition、filtered_h5、raw_h5无缺失。`neuron_object_path()`优先定位`03_neurons/{type}_annotated.rds`，无该文件时定位`{type}.rds`。`join_rna_layers()`兼容Seurat v5多counts layers；`rna_matrix()`兼容Assay5 layer与旧Assay slot。`integer_counts()`以1e-8容差检查非负整数；`require_neuron_metadata()`要求sample、condition、subtype完整且没有Unknown。`bh()`为BH调整。`plot_export()`同时保存白底PNG与Cairo PDF，默认600 dpi、单位英寸、`limitsize = FALSE`。

### `R/preprocessing.R`

提供HDF5读取、版本要求、PC选择、sample Harmony、注释映射和marker诊断图工具。SoupX及scDblFinder版本严格按配置匹配；HDF5使用gene symbols、唯一feature名称。PC数不超过核数−1及高变基因数−1。

`run_sample_harmony()`固定`group.by.vars = sample`、`reduction.use = pca`、`assay.use = RNA`、`project.dim = FALSE`、`reduction.save = harmony`。`annotation_mapping()`拒绝重复cluster和不属于该对象的cluster；未填写项为Unknown，不自动归入相邻类型。

本文件同时提供六大类和6个Exc／8个Inh的marker panel。诊断DotPlot为dot.scale=10、按基因scale，颜色lightgrey／blue／red；画布宽`max(10, marker数×0.45)`、高7。marker FeaturePlot为UMAP、point size=0.3、`order = TRUE`、min.cutoff=q10、max.cutoff=q90，画布8×7。该诊断显示与`R/plotting.R`的论文汇总显示分别设置，不改变表达矩阵或分群。

### `R/plotting.R`

提供统一绘图主题、配色、星号规则、论文UMAP、DotPlot、火山图和富集图函数。字体Arial，base size=14。条件颜色：WD `#9E9E9E`、Ext `#49B8D1`、Ret–Ext `#EF4936`。大类及亚群配色由同名命名向量固定，保证跨图一致。

| 分组 | 配色（按顺序） |
|---|---|
| Excitatory、Inhibitory、Oligodendrocytes、OPC、Astrocytes、Microglia | `#3B6FBF`、`#E85D5D`、`#B5895A`、`#9A6FD1`、`#4DB894`、`#F0986D` |
| GLU_L6_CLA、GLU_L5_6、GLU_L2_L3、GLU_L5、GLU_L4_5、GLU_ADGM | `#8DD3C7`、`#FFFFB3`、`#BEBADA`、`#FB8072`、`#80B1D3`、`#FDB462` |
| GABA_Amb、GABA_Etv1、GABA_Lamp5、GABA_MSN、GABA_Sst、GABA_Hs3st2、GABA_Pvalb、GABA_Vip | `#66C2A5`、`#FC8D62`、`#8DA0CB`、`#E78AC3`、`#A6D854`、`#FFD92F`、`#E5C494`、`#FF6347` |

`stars()`使用严格小于：P<0.0001为****、P<0.001为***、P<0.01为**、P<0.05为*，其余为ns，缺失为NA；传入的P值类型由各图指定。`collect_enrichment()`读取完整富集表、排除重复的`*_significant.csv`，保留每个条目的真实padj，再由显示函数筛选。通路名称不会添加数据库尾标；数据库保留在数据列和shape映射中。

### `R/communication.R`

负责严格一对一rat–human映射、IEG-high表达输入、两个谱系的共同gene集合合并、网络矩阵长表和outgoing条件比较。映射类型仅接受one_to_one／ortholog_one2one，并要求两端基因符号唯一；不将大小写转换称为同源映射。

`cellchat_input()`检查RNA data非负有限、核名唯一、condition标签合法，保留sample和subtype。`merge_cellchat_inputs()`要求不同谱系的核标识不冲突。`outgoing_network_summary()`分别保存跨亚群outgoing与autocrine权重／交互数。`compare_outgoing_conditions()`计算Ret–Ext−Ext、Ret–Ext／Ext、`100×(Ret–Ext−Ext)/Ext`；基线为0时标记zero_baseline，fold change和percent change为NA。

## 参数与解释边界

配置中包含分析阈值，也包含实现细节。例如DE counts预过滤、分位数type、GLMM优化器和积分设置用于明确代码行为；它们不应与生物学结论、历史核数或特定亚群必然显著相混淆。更改这些参数需要保留对应配置、输入及输出，使图表与其实际统计表一致。

本文各分析使用的单位不同：DESeq2使用pooled sample计数；GLMM以sample随机截距处理同一sample内核的相关性；marker-program评分描述核级表达分布；CellChat为表达支持的通信推断。原始P、BH-adjusted P与Bonferroni-adjusted P按各自分析分别保存和标注。
