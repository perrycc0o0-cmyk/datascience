# Drug addiction snRNA-seq analysis

大鼠下边缘皮层（IL）的单核转录组分析与论文图表代码。分析包括逐文库质控、ambient RNA 校正、doublet 识别、样本批次整合、细胞注释、神经元亚群分析、IEG-high 定义、pseudobulk 差异表达、marker-program 评分、通路富集、比例模型及推断性细胞通信分析。

## 实验设计

包含 WD、Ext 和 Ret–Ext 三个条件，每个条件有 3 个独立 pooled samples。每个 sample 由 5 只大鼠的 IL 组织混合构成，共 9 个测序文库。条件间生物学重复数为 **n = 3/condition**；同一文库中的细胞核不是独立动物重复。

Harmony 的批次变量为 `sample`，不把 experimental condition 作为待消除的批次变量。差异表达以 sample 为单位汇总 RNA assay 中未归一化的整数 UMI counts 后进行 DESeq2 分析；marker-program 评分保留核级分布展示。比例模型以 condition 为固定效应、sample 为随机截距。

## 目录

```text
drugaddiction/
├── config/
│   ├── analysis.json               分析参数和输入输出路径
│   ├── samples.csv                 文库、条件及 10x 文件路径
│   ├── major_annotations.csv       cluster 到主要细胞类别的映射
│   ├── subtype_annotations.csv     神经元 cluster 到亚群的映射
│   ├── markers.json                细胞类别与神经元亚群 marker
│   ├── rat_human_orthologs.tsv     CellChat 使用的一对一同源基因表
│   ├── heatmap_genes.csv           六模块热图的基因与顺序
│   └── enrichment_terms.csv        通路可视化的名称匹配规则
├── R/                             数据处理、统计和作图辅助函数
├── scripts/                       01–12 分析、作图及注释应用脚本
│   └── README.md                  逐脚本说明与详细参数
└── README.md
```

输入数据和分析产生的结果不纳入代码目录。默认结果写入 `results/`，可以使用 `--output` 指定其他结果根目录。

## R 环境与依赖

请在兼容的 R/Bioconductor 环境中安装下列依赖：

| 功能 | R 包 |
|---|---|
| 稀疏矩阵与单核对象 | Seurat、SeuratObject、Matrix |
| Ambient RNA 校正 | **SoupX 1.6.2** |
| Doublet 识别 | **scDblFinder 1.24.10**、SingleCellExperiment、SummarizedExperiment |
| 批次整合 | harmony |
| 差异表达 | DESeq2 |
| 比例模型 | lme4 |
| GO/KEGG 富集 | clusterProfiler、AnnotationDbi、org.Rn.eg.db、KEGGREST |
| 细胞通信 | CellChat |
| 配置及绘图 | jsonlite、ggplot2、ggrepel、scales、igraph |

预处理脚本显式要求上述 SoupX 和 scDblFinder 版本。各分析阶段保存 `sessionInfo.txt`，记录相应软件环境。KEGG 分析需要访问 KEGG 服务；通信分析需要可用的 human CellChat 数据库。PNG/PDF 使用 Arial 字体，PDF 导出依赖 R 的 Cairo 支持；请确保运行环境可找到 Arial。

## 准备输入

### 1. 文库表

编辑 [`config/samples.csv`](config/samples.csv)，每行一个文库，必须包含：

```csv
sample,condition,filtered_h5,raw_h5
Con1,WD,data/Con1/filtered_feature_bc_matrix.h5,data/Con1/raw_feature_bc_matrix.h5
```

`sample` 必须唯一；`condition` 使用 `WD`、`Ext`、`Ret-Ext`。每个条件需要 3 行。配置中的 `Con1–Con3` 对应稿件的 `WD1–WD3`，`Ret_ext1–Ret_ext3` 对应稿件的 `Ret–Ext1–Ret–Ext3`；这些是同一文库的文件标识与展示名称，不是额外样本。`filtered_h5` 为 cell-called 10x matrix，`raw_h5` 为同一文库包含空液滴的 raw matrix，二者的 feature 和 barcode 必须对应。路径相对于项目工作目录，也可以填写绝对路径。输入 feature 名使用大鼠 gene symbols。

### 2. 主要细胞类别映射

[`config/major_annotations.csv`](config/major_annotations.csv) 使用 `cluster,cell_type` 两列，包含 38 个 cluster 的映射。cluster 编号只对相应的聚类对象有意义，不能仅因编号相同就将映射套到其他对象。请结合 `02_annotation/cluster_markers.csv` 和 marker 图审阅映射。

主要类别包括 `Excitatory`、`Inhibitory`、`Oligodendrocytes`、`OPC`、`Astrocytes` 和 `Microglia`。其他明确类别可以保留自己的名称；没有足够证据的 cluster 不应强行并入这六类。映射未完成的 cluster 保留 `Unknown`，完整注释后才生成供神经元分析读取的 `annotated.rds`。

### 3. 神经元亚群映射

编辑 [`config/subtype_annotations.csv`](config/subtype_annotations.csv)，填写 `neuron_type,cluster,subtype` 三列。`neuron_type` 使用 `Excitatory` 或 `Inhibitory`；cluster 以相应神经元独立重聚类后的编号为准。

脚本 03 导出每类神经元的 cluster markers、marker DotPlot/UMAP 和 `subtype_annotations_template.csv`。结合这些证据人工填写映射，不能用主要细胞类别的 cluster 编号替代神经元重聚类编号。使用 `apply_annotations.R --level subtype` 将审阅后的映射应用到对象。完整亚群注释是后续 IEG、差异表达、比例和通信分析的前提。

### 4. Rat–human 同源基因表

CellChat 使用 human 数据库，需要提供 [`config/rat_human_orthologs.tsv`](config/rat_human_orthologs.tsv)，包含以下制表符分隔的列：

```text
rat_symbol	human_symbol	orthology_type
```

填入真实的一对一同源关系；`orthology_type` 使用 `one_to_one` 或 `ortholog_one2one`。保留的 rat symbol 和 human symbol 均须唯一。仅将大鼠符号转为大写不能替代同源映射。未提供完整映射时，通信分析会停止。

### 5. Marker 与热图基因

[`config/markers.json`](config/markers.json) 保存细胞类别和神经元亚群 marker 面板。[`config/heatmap_genes.csv`](config/heatmap_genes.csv) 包含 `gene_order,module_id,module,gene`，用于固定六模块热图的基因分组和显示顺序。热图模块使用该配置，不把显示分组解释为自动发现的聚类模块。

## 运行顺序

在 `drugaddiction/` 目录下运行。所有阶段应使用同一份配置和同一个结果根目录；下列命令可逐条执行。

```bash
Rscript scripts/01_preprocess_libraries.R --config config/analysis.json --output results
Rscript scripts/02_integrate_annotate.R --config config/analysis.json --output results
```

检查脚本 02 的 cluster markers 和 marker 图，填写或确认 `config/major_annotations.csv`。如果 `02_annotation/annotated.rds` 尚未生成，使用以下命令把完整映射应用到已有对象，不需要重新聚类：

```bash
Rscript scripts/apply_annotations.R --level major
```

确认主要细胞类别映射完整，再进行神经元重聚类：

```bash
Rscript scripts/03_neuron_subclustering.R --config config/analysis.json --output results
```

检查脚本 03 的亚群 markers 和 marker 图，填写 `config/subtype_annotations.csv`，再应用映射：

```bash
Rscript scripts/apply_annotations.R --level subtype
```

该步骤生成 `03_neurons/Excitatory_annotated.rds` 和 `03_neurons/Inhibitory_annotated.rds`。后续 IEG 分析及细胞类型绘图优先使用这两个已审阅注释对象。确认兴奋性和抑制性神经元亚群映射完整，再运行下游统计与绘图：

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

### 注释与结果目录

编号脚本创建自己的阶段目录，**拒绝覆盖非空的阶段目录**。修改处理、聚类或统计参数后，需要选择新的输出根目录，例如将命令中的 `--output results` 改为 `--output results_reviewed`。

首次注释可以先生成 cluster markers 和 marker 图，再填写映射并通过 `apply_annotations.R` 应用到同一结果根目录，不需要删除原对象或重跑脚本 01。注释应用脚本也不覆盖已存在的注释对象。如果需要修改已保存的完整注释，先保留当前结果，再在新的输出根目录中使用所需上游对象和修改后的映射，不要混用不同注释版本的下游结果。

若只重做后续处理阶段，可以将保留的上游阶段目录复制到新根目录的相同位置，再从需要重做的阶段继续。例如，重做脚本 02 需要新根目录中已有 `01_preprocessing/merged.rds`；重做脚本 03 需要其中已有 `02_annotation/annotated.rds`。单独改变后续脚本的 `--output` 不会自动读取另一个结果根目录中的上游对象。注释阶段保存的模板和 marker 结果用于人工判断，脚本不会根据名称或预期核数补齐缺失标签。上面的注释应用命令使用默认配置 `config/analysis.json` 和其中的输出路径；更换结果根目录时也应保持该配置中的路径一致。

## 脚本与论文分析的对应关系

| 脚本 | 分析内容 | Methods / 图表内容 |
|---|---|---|
| `01_preprocess_libraries.R` | 固定QC、SoupX、scDblFinder、逐文库对象合并 | snRNA-seq 数据处理和质控 |
| `02_integrate_annotate.R` | PCA、sample Harmony、全局聚类、主要细胞类型注释 | 主要细胞类型 UMAP、marker 及补充QC |
| `03_neuron_subclustering.R` | Exc/Inh 独立聚类、亚群注释和 marker 图 | Fig. 3 神经元亚群及 Fig. S6 相关图 |
| `apply_annotations.R` | 将人工审阅的主要类别或神经元亚群映射写入已聚类对象 | 细胞身份注释 |
| `04_ieg_definition.R` | 八基因表达总和与类别内第95百分位筛选 | IEG-high 定义、阈值和分布 |
| `05_pseudobulk_DE.R` | sample × 神经元类别/IEG class 的整数counts汇总、DESeq2、BH | Fig. 3 与 Fig. S7 条件差异表达 |
| `06_marker_program_scores.R` | IEG-high 相关基因集与核级 marker-program 评分 | Fig. 3F/G 评分分布 |
| `07_GO_KEGG_enrichment.R` | 分方向 GO/KEGG 富集与检验背景导出 | Fig. 3 通路富集及相应补充结果 |
| `08_subtype_GLMM.R` | condition 固定效应、sample 随机截距的二项GLMM，BH校正 | Fig. 3L、Fig. S8B 比例比较 |
| `09_CellChat_analysis.R` | 一对一同源映射、IEG-high 亚群的通信推断 | Fig. 3P 与 Fig. S9 的 predicted communication |
| `10_plot_cell_types.R` | 主要类别/神经元亚群 UMAP、DotPlot 和 marker UMAP | Fig. 3 与 Fig. S6 等细胞身份图 |
| `11_plot_expression.R` | DEG 火山图、富集图、评分及固定六模块热图 | Fig. 3、Fig. S7 表达和通路图 |
| `12_plot_communication_proportions.R` | 亚群比例、GLMM显著性、通信网络及通路图 | Fig. 3L、Fig. S8B、Fig. S9 及通信图 |

逐脚本的输入、输出、参数与统计含义见 [`scripts/README.md`](scripts/README.md)。稿件未逐项展开的计算设置集中列在该文档的“实现参数汇总”中，包括 SoupX 预聚类、DESeq2 预过滤与拟合设置、GLMM 检验及优化参数、同源映射策略和绘图筛选规则。这些脚本对应 Fig. 3 的计算分析部分（B–D、F–L、P）及相关补充图；Fig. 3M–O 的湿实验不属于本代码流程。`config/analysis.json` 是分析参数入口。CellChat 输出为基于表达数据的推断性通信，不等同于直接测得的神经投射或连接强度。
