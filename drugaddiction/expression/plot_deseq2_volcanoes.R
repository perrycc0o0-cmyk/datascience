# plot deseq2 volcanoes

# 1. Packages
library(grid)
library(ggplot2)
library(ggrepel)
library(ragg, include.only = character())
library(svglite, include.only = character())
library(systemfonts, include.only = character())

# 2. Paths
work_dir <- "/path/to/work"
input_dir <- file.path(work_dir, "deseq2")
figure_dir <- file.path(work_dir, "volcanoes")
result_dir <- figure_dir

# 3. Samples and conditions
sample_cols <- c("WD1", "WD2", "WD3", "Ext1", "Ext2", "Ext3", "Ret-Ext1", "Ret-Ext2", "Ret-Ext3")
metadata <- data.frame(
  pool_id = sample_cols,
  source_sample_id = c("Con1", "Con2", "Con3", "Ext1", "Ext2", "Ext3", "Ret_ext1", "Ret_ext2", "Ret_ext3"),
  model = rep("Rat IL", 9L),
  condition = factor(rep(c("WD", "Ext", "Ret-Ext"), each = 3L), levels = c("WD", "Ext", "Ret-Ext")),
  rats_per_pool = rep(5L, 9L), stringsAsFactors = FALSE)

# 5. Helper functions
# Input checks and native DESeq2 calculations for the expression workflows.
require_packages <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install the required R packages: ", paste(absent, collapse = ", "))
}

limit_threads <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
  }
}

new_output_dir <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) {
    stop("Choose an empty output directory: ", path)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

read_table <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-BOM")
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

write_session <- function(out) {
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
}

pool_ids <- sample_cols
pool_conditions <- setNames(as.character(metadata$condition), sample_cols)

read_ieg_payload <- function(path, neuronal_class) {
  require_packages("Matrix")
  q <- readRDS(path)
  stopifnot(all(c("counts", "metadata", "gene_annotation") %in% names(q)))
  x <- q$counts
  md <- q$metadata
  ann <- q$gene_annotation
  stopifnot(inherits(x, "sparseMatrix"),
            all(c("Cell", "pool_id", "condition", "IEGhigh") %in% names(md)),
            all(c("seurat_feature", "gene_id", "original_symbol") %in% names(ann)),
            identical(colnames(x), as.character(md$Cell)),
            identical(rownames(x), as.character(ann$seurat_feature)),
            !anyDuplicated(md$Cell), !anyDuplicated(rownames(x)),
            all(md$IEGhigh), setequal(md$pool_id, pool_ids),
            identical(as.character(md$condition), unname(pool_conditions[as.character(md$pool_id)])),
            all(is.finite(x@x)), all(x@x >= 0), all(x@x == round(x@x)))
  if ("neuron_type" %in% names(md)) stopifnot(all(md$neuron_type == neuronal_class))
  q
}

aggregate_pools <- function(q) {
  group <- Matrix::sparseMatrix(i = seq_len(ncol(q$counts)),
                               j = match(q$metadata$pool_id, pool_ids),
                               x = 1, dims = c(ncol(q$counts), length(pool_ids)))
  counts <- as.matrix(q$counts %*% group)
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts <= .Machine$integer.max))
  dimnames(counts) <- list(rownames(q$counts), pool_ids)
  storage.mode(counts) <- "integer"
  counts
}

make_pool_metadata <- function(pools, levels) {
  data.frame(pool_id = pools, condition = factor(unname(pool_conditions[pools]), levels),
             row.names = pools)
}

expression_filter <- function(counts, metadata) {
  unname(edgeR::filterByExpr(edgeR::DGEList(counts, group = metadata$condition),
                            design = model.matrix(~ condition, metadata)))
}

fit_native_deseq <- function(counts, metadata, eligible, numerator, denominator,
                             allow_convergence_retry = FALSE,
                             contrast = c("condition", numerator, denominator)) {
  stopifnot(identical(colnames(counts), rownames(metadata)), any(eligible),
            length(eligible) == nrow(counts), !anyNA(eligible))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts[eligible, , drop = FALSE],
                                       metadata, design = ~ condition)
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", sfType = "ratio",
                       betaPrior = FALSE, minReplicatesForReplace = 7L,
                       quiet = TRUE, parallel = FALSE)
  initial <- S4Vectors::mcols(dds)$betaConv
  retried <- any(is.na(initial) | !initial)
  if (retried && allow_convergence_retry) {
    dds <- DESeq2::nbinomWaldTest(dds, betaPrior = FALSE, maxit = 1000, quiet = TRUE)
  }
  stopifnot(all(S4Vectors::mcols(dds)$betaConv),
            all(DESeq2::counts(dds) == counts[eligible, , drop = FALSE]),
            !("replaceCounts" %in% SummarizedExperiment::assayNames(dds)))
  res <- DESeq2::results(dds, contrast = contrast,
                         alpha = 0.05, independentFiltering = TRUE,
                         pAdjustMethod = "BH", parallel = FALSE)
  threshold <- as.numeric(S4Vectors::metadata(res)$filterThreshold)
  family <- is.finite(res$pvalue) & res$baseMean >= threshold
  stopifnot(length(threshold) == 1L, is.finite(threshold),
            identical(unname(is.finite(res$padj)), unname(family)))
  if (any(family)) {
    stopifnot(all(abs(res$padj[family] - p.adjust(res$pvalue[family], "BH")) < 1e-12))
  }
  list(model = dds, result = res, filter_threshold = threshold,
       convergence_retry = retried && allow_convergence_retry)
}

native_de_table <- function(fit, counts, annotation, eligible, neuronal_class,
                            numerator, denominator) {
  r <- as.data.frame(fit$result)
  idx <- match(rownames(counts), rownames(r))
  mc <- S4Vectors::mcols(fit$model)
  cooks <- SummarizedExperiment::assays(fit$model)[["cooks"]]
  max_cooks <- apply(cooks, 1, function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_)
  z <- data.frame(
    gene = rownames(counts), gene_id = annotation$gene_id,
    original_symbol = annotation$original_symbol, eligible = eligible,
    tested = eligible & is.finite(r$pvalue[idx]),
    in_BH_family = eligible & is.finite(r$padj[idx]),
    baseMean = r$baseMean[idx], log2FoldChange = r$log2FoldChange[idx],
    lfcSE = r$lfcSE[idx], stat = r$stat[idx], pvalue = r$pvalue[idx],
    padj = r$padj[idx], maxCooks = max_cooks[idx], betaConv = mc$betaConv[idx],
    neuron_type = neuronal_class, numerator = numerator, denominator = denominator,
    method = "DESeq2", n_pools = ncol(counts), n_per_condition = 3L,
    independentFiltering = TRUE, alpha = 0.05, stringsAsFactors = FALSE)
  z$status <- ifelse(!eligible, "filtered_low_expression_before_fit",
                     ifelse(!z$tested, "native_P_withheld",
                            ifelse(!z$in_BH_family, "independent_filter_low_mean", "tested_with_native_BH")))
  z$BH_significant <- z$in_BH_family & z$padj < 0.05
  z$display_DEG <- z$BH_significant & abs(z$log2FoldChange) >= log2(1.15)
  z$DE_direction <- ifelse(z$display_DEG, ifelse(z$log2FoldChange > 0, "Up", "Down"),
                           ifelse(z$in_BH_family, "Not significant", "Not in BH family"))
  z
}

# Arial figure export. Fonts are supplied by the operating system.
load_figure_packages <- function(extra = character()) {
  require_packages(unique(c("grid", "ragg", "svglite", "systemfonts", extra)))
  stopifnot(capabilities("cairo"))
  regular <- systemfonts::match_fonts("Arial")
  bold <- systemfonts::match_fonts("Arial", weight = "bold")
  paths <- tolower(basename(c(regular$path, bold$path)))
  if (length(paths) != 2L || !all(grepl("arial", paths))) {
    stop("Install Arial Regular and Arial Bold before rendering.")
  }
}

export_plot <- function(stem, width, height, draw) {
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width, height, family = "Arial", bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
  ragg::agg_png(paste0(stem, ".png"), width, height, units = "in", res = 600, background = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
  svglite::svglite(paste0(stem, ".svg"), width, height, bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
}

# 6. Analysis and output
# Render the six pairwise DESeq2 contrasts using native BH-adjusted values.
load_figure_packages(c("ggplot2", "ggrepel"))
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(grid) })
out <- new_output_dir(figure_dir)
idx<-data.frame(figure=c('Figure_3',rep('Figure_S8',5)),panel=c('H','A','B','E','C','D'),
 kind=c(rep('Excitatory',3),rep('Inhibitory',3)),contrast=c('RetExt_vs_Ext','Ext_vs_WD','RetExt_vs_WD','RetExt_vs_Ext','Ext_vs_WD','RetExt_vs_WD'),
 width_in=5.1,height_in=c(4.34807692307692,4.24784810126582,4.28035714285714,4.692,4.09243902439024,4.4673417721519),stringsAsFactors=FALSE)
ok <- function(id, label, value) if (!isTRUE(value)) stop(id, ': ', label)
fc<-log2(1.15)
for(i in seq_len(nrow(idx))) {
 r<-idx[i,];id<-paste(r$kind,r$contrast,sep='__');path<-file.path(input_dir,'datasets',id,'tables','ALL_genes.csv')
 z<-read_table(path)
 ok(id,'original gene universe and native method',nrow(z)>0L&&!anyDuplicated(z$gene)&&all(z$method=='DESeq2')&&all(z$independentFiltering)&&all(z$alpha==.05))
 shown<-is.finite(z$padj)&is.finite(z$log2FoldChange);d<-z[shown,,drop=FALSE]
 if (!nrow(d)) stop('No finite native BH-adjusted values for ', id)
 d$display_P<-d$padj;d$minus_log10_display<--log10(pmax(d$padj,.Machine$double.xmin))
 d$display_metric<-'BH-adjusted P';d$display_direction<-ifelse(d$padj<.05&abs(d$log2FoldChange)>=fc,ifelse(d$log2FoldChange>0,'Up','Down'),'Not significant')
 d$label_text<-d$original_symbol;bad<-is.na(d$label_text)|!nzchar(d$label_text);d$label_text[bad]<-d$gene[bad]
 dup<-duplicated(d$label_text)|duplicated(d$label_text,fromLast=TRUE);d$label_text[dup]<-d$gene[dup]
 labels<-do.call(rbind,lapply(c('Up','Down'),function(s){u<-d[d$display_direction==s,,drop=FALSE];head(u[order(u$padj,u$pvalue,-abs(u$log2FoldChange),u$gene),,drop=FALSE],8L)}))
 extent<-max(abs(d$log2FoldChange),fc);ymax<-max(d$minus_log10_display,-log10(.05));exc<-r$kind=='Excitatory'
 title<-paste0(r$kind,' IEGhigh: ',unique(d$numerator),' vs ',unique(d$denominator))
 p<-ggplot(d,aes(log2FoldChange,minus_log10_display,colour=log2FoldChange))+
  geom_point(aes(size=minus_log10_display),alpha=.60,shape=16,stroke=0)+
  geom_vline(xintercept=c(-1,1)*fc,linetype=2,linewidth=.3,colour='black')+
  geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.3,colour='black')+
  scale_colour_gradientn(colours=c('#3288bd','#66c2a5','#ffffbf','#f46d43','#9e0142'),values=c(0,.25,.5,.75,1),limits=c(-extent,extent),name='log2 FC Change')+
  scale_size_continuous(range=c(.45,3.4),limits=c(0,ymax),name='-log10(padj)')+
  scale_x_continuous(limits=c(-extent*1.10,extent*1.10),expand=expansion(mult=.04))+
  scale_y_continuous(limits=c(0,ymax),expand=expansion(mult=c(.015,.18)))+
  labs(title=title,x='log2 Fold Change',y=expression(-log[10]('padj')))+theme_classic(base_size=12,base_family='Arial')+
  theme(text=element_text(face='bold',family='Arial'),axis.line=element_line(linewidth=.65*1.6,colour='black'),
   axis.ticks=element_line(linewidth=.55*1.6,colour='black'),axis.text=element_text(size=8.925*1.15,colour='black',face='bold'),
   axis.title=element_text(size=10.2*1.15,face='bold'),plot.title=element_text(size=(if(exc)9 else 7.5)*1.15,face='bold',hjust=.5),
   legend.title=element_text(size=5.95*1.15,face='bold'),legend.text=element_text(size=5.1*1.15,face='bold'),legend.key.size=unit(2.5,'mm'),
   plot.margin=margin(25,6,7,15),plot.background=element_rect(fill='white',colour=NA))+
  guides(colour=guide_colourbar(title='log2 FC Change',order=1,barheight=unit(if(exc)10 else 12,'mm'),barwidth=unit(if(exc)1.25 else 1.6,'mm')),
   size=guide_legend(title='-log10(padj)',order=2))
 if(nrow(labels))p<-p+geom_label_repel(data=labels,aes(label=label_text),family='Arial',fontface='bold',fill='white',
  size=(if(exc)3.25 else 3.2)*1.3,label.size=.08,label.padding=.06,box.padding=.25,point.padding=.15,min.segment.length=0,
  segment.colour='#999999',segment.size=.12,force=5,force_pull=.3,max.overlaps=Inf,max.time=5,max.iter=100000,seed=20261010,show.legend=FALSE)
 b<-ggplot_build(p)
 ok(id,'all finite native BH values rendered at correct coordinates',nrow(b$data[[1]])==sum(shown)&&identical(as.numeric(b$data[[1]]$x),d$log2FoldChange)&&identical(as.numeric(b$data[[1]]$y),d$minus_log10_display))
 ok(id,'labels only native BH significant FC115 genes',all(labels$padj<.05&abs(labels$log2FoldChange)>=fc))
 ok(id,'actual adjusted P threshold at .05',identical(b$data[[3]]$yintercept,-log10(.05)))
 export_plot(file.path(out, paste0(id, "_DESeq2_BH")), r$width_in, r$height_in, function() print(p))
 sd <- file.path(out, "plot_source_data", id)
 write_table(d, file.path(sd, "points.csv"))
 write_table(labels, file.path(sd, "labels.csv"))
 omitted <- z[!shown, , drop = FALSE]
 omitted$reason <- rep("Native adjusted P or log2 fold change unavailable", nrow(omitted))
 write_table(omitted, file.path(sd, "omitted_unavailable.csv"))
}
