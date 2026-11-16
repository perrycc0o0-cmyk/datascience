# plot gene category network

# 1. Packages
library(grid)
library(igraph)
library(ragg, include.only = character())
library(svglite, include.only = character())
library(systemfonts, include.only = character())

# 2. Paths
work_dir <- "/path/to/work"
input_dir <- file.path(work_dir, "gene_category_network")
figure_dir <- file.path(work_dir, "network_figures")
result_dir <- figure_dir
down_layout <- NULL
up_layout <- NULL

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
layout_files <- c(down_layout, up_layout)
stopifnot(length(layout_files) %in% c(0L, 2L))
# Draw the down/up gene-category network without panel tags or label backgrounds.
load_figure_packages(if (length(layout_files) == 0L) "igraph" else character())
suppressPackageStartupMessages(library(grid))
if (length(layout_files) == 0L) suppressPackageStartupMessages(library(igraph))
out <- new_output_dir(figure_dir)
invisible(Sys.setlocale("LC_COLLATE", "C"))
nodes <- read_table(file.path(input_dir, "display_nodes.csv"))
edges <- read_table(file.path(input_dir, "display_edges.csv"))
stopifnot(all(edges$edge_type == "gene-category membership"),
          !anyDuplicated(paste(nodes$direction, nodes$node_id)),
          !anyDuplicated(paste(edges$direction, edges$category_id, edges$gene)))
types <- c("Downregulated", "Upregulated")
W <- 220; H <- 240; gap <- 6
max_count <- max(nodes$count[nodes$node_type == "category"])
stopifnot(is.finite(max_count), max_count > 0)
fc_limit <- 2
palette <- colorRampPalette(c("#235581", "#F7F7F7", "#A40E20"))(1001)
check <- function(label, value) if (!isTRUE(value)) stop(label)
color <- function(x){
 transformed <- sign(x)*sqrt(pmin(abs(x)/fc_limit,1))
 palette[1+round((transformed+1)*500)]
}
count_breaks <- c(3,10,25,max_count)
module_size <- function(n)8.2*sqrt(n/max_count)
anchors <- rbind(C1=c(-.9,1),C2=c(.9,-1),C3=c(-.4,-.05),C4=c(.9,.1),C5=c(.8,1.1),C6=c(-.6,-1.1))
network_W <- W-35;network_H <- H-28
place_labels <- function(z){
 z$node_x_mm <- 12+(z$x+1)/2*(network_W-24)
 z$node_y_mm <- 12+(z$y+1)/2*(network_H-24)
 z$label_width_mm <- z$label_height_mm <- NA_real_
 for(k in seq_len(nrow(z))){
  gp<-gpar(fontfamily='Arial',fontface=z$fontface[k],fontsize=z$label_size_mm[k]*72/25.4,lineheight=.91)
  gr<-textGrob(z$display_label[k],gp=gp)
  z$label_width_mm[k]<-convertWidth(grobWidth(gr),'mm',valueOnly=TRUE)+.6
  z$label_height_mm[k]<-convertHeight(grobHeight(gr),'mm',valueOnly=TRUE)+.6
 }
 pair_distance<-as.matrix(dist(cbind(z$node_x_mm,z$node_y_mm)))
 density<-rowSums(pair_distance<20)-1
 order_nodes<-c(which(z$node_type=='category')[order(-z$count[z$node_type=='category'])],
  which(z$node_type=='gene')[order(-density[z$node_type=='gene'],z$node_id[z$node_type=='gene'])])
 z$label_x_mm<-z$label_y_mm<-NA_real_
 boxes<-matrix(numeric(),nrow=0,ncol=4)
 radii<-c(seq(2.8,40,by=1.2),seq(42,90,by=2))
 candidates<-expand.grid(radius=radii,angle=seq(0,2*pi,length.out=49)[-49])
 for(k in order_nodes){
  w<-z$label_width_mm[k];h<-z$label_height_mm[k]
  cx<-z$node_x_mm[k]+candidates$radius*cos(candidates$angle)
  cy<-z$node_y_mm[k]+candidates$radius*sin(candidates$angle)
  legal<-cx-w/2>1&cx+w/2<network_W-1&cy-h/2>1&cy+h/2<network_H-1
  cost<-candidates$radius+if(z$node_type[k]=='category') .4*abs(cos(candidates$angle))else 0
  chosen<-NA_integer_
  for(j in order(cost)){
   if(!legal[j])next
   box<-c(cx[j]-w/2-.5,cx[j]+w/2+.5,cy[j]-h/2-.5,cy[j]+h/2+.5)
   if(nrow(boxes)&&any(box[1]<boxes[,2]&box[2]>boxes[,1]&box[3]<boxes[,4]&box[4]>boxes[,3]))next
   nr<-z$point_size_mm/2+.35
   if(any(z$node_x_mm+nr>box[1]&z$node_x_mm-nr<box[2]&z$node_y_mm+nr>box[3]&z$node_y_mm-nr<box[4]))next
   chosen<-j;break
  }
  check(paste(z$direction[k],z$node_id[k],'nonoverlapping fixed text position found'),!is.na(chosen))
  z$label_x_mm[k]<-cx[chosen];z$label_y_mm[k]<-cy[chosen]
  boxes<-rbind(boxes,c(cx[chosen]-w/2-.5,cx[chosen]+w/2+.5,cy[chosen]-h/2-.5,cy[chosen]+h/2+.5))
 }
 check(paste(unique(z$direction),'fixed text boxes do not overlap'),all(vapply(seq_len(nrow(boxes)),function(k){
  if(k==nrow(boxes))return(TRUE)
  b<-boxes[(k+1):nrow(boxes),,drop=FALSE]
  !any(boxes[k,1]<b[,2]&boxes[k,2]>b[,1]&boxes[k,3]<b[,4]&boxes[k,4]>b[,3])
 },logical(1))))
 z
}

make_layout <- function(i, type, z, e) {
 g <- graph_from_data_frame(data.frame(from=e$category_id,to=e$gene),directed=FALSE,
  vertices=data.frame(name=z$node_id))
 stopifnot(ecount(g)==nrow(e),vcount(g)==nrow(z))
 set.seed(20261105L+i)
 initial <- matrix(0,nrow(z),2)
 for(k in seq_len(nrow(z))){
  if(z$node_type[k]=='category')initial[k,]<-anchors[z$node_id[k],]
  else{
   ids<-e$category_id[e$gene==z$node_id[k]]
   initial[k,]<-colMeans(anchors[ids,,drop=FALSE])+rnorm(2,0,.32)
  }
 }
 initial <- initial*sqrt(nrow(z))*.7
 fixed <- which(z$node_type=='category')
 lowerx <- lowery <- rep(-Inf,nrow(z));upperx <- uppery <- rep(Inf,nrow(z))
 lowerx[fixed] <- upperx[fixed] <- initial[fixed,1]
 lowery[fixed] <- uppery[fixed] <- initial[fixed,2]
 xy <- layout_with_fr(g,coords=initial,niter=4000L,start.temp=.5,grid='nogrid',
  minx=lowerx,maxx=upperx,miny=lowery,maxy=uppery)
 xy[,1]<-(xy[,1]-mean(range(xy[,1])))/diff(range(xy[,1]))*2
 xy[,2]<-(xy[,2]-mean(range(xy[,2])))/diff(range(xy[,2]))*2
 z$x<-xy[,1];z$y<-xy[,2]
 z$display_label <- z$label
 gene_ix <- which(z$node_type=='gene')
 z$display_label[gene_ix] <- e$symbol[match(z$node_id[gene_ix],e$gene)]
 stopifnot(!anyNA(z$display_label[gene_ix]),!anyDuplicated(z$display_label[gene_ix]))
 z$display_label[z$node_type=='category'] <- vapply(z$label[z$node_type=='category'],
  function(s)paste(strwrap(s,width=21),collapse='\n'),character(1))
 z$fontface <- ifelse(z$node_type=='category','bold.italic','bold')
 z$label_size_mm <- ifelse(z$node_type=='category',4.1,3.25)
 z$point_size_mm <- ifelse(z$node_type=='category',module_size(z$count),1.9)
 z$node_color <- color(z$log2FC)
 z <- place_labels(z)
 z
}

plot_nodes <- plot_edges <- list()
for (i in seq_along(types)) {
  type <- types[i]
  z <- nodes[nodes$direction == type, , drop = FALSE]
  e <- edges[edges$direction == type, , drop = FALSE]
  stopifnot(nrow(z) > 0L, nrow(e) > 0L, all(c(e$category_id, e$gene) %in% z$node_id))
  if (length(layout_files) == 2L) {
    approved <- read_table(layout_files[i])
    ix <- match(approved$node_id, z$node_id)
    stopifnot(!anyDuplicated(approved$node_id), setequal(z$node_id, approved$node_id),
              all(approved$direction == type), identical(approved$node_type, z$node_type[ix]),
              all(abs(approved$log2FC - z$log2FC[ix]) < 1e-12),
              all(approved$count == z$count[ix]), identical(approved$categories, z$categories[ix]))
    z <- approved
  } else {
    scratch <- tempfile(fileext = ".pdf")
    grDevices::cairo_pdf(scratch, W / 25.4, H / 25.4, family = "Arial", bg = "white")
    z <- tryCatch({
      grid.newpage()
      make_layout(i, type, z, e)
    }, finally = {
      grDevices::dev.off()
      unlink(scratch)
    })
  }
  plot_nodes[[type]] <- z
  plot_edges[[type]] <- e
  write_table(z, file.path(out, paste0(type, "_nodes_and_layout.csv")))
}
draw_network_fixed <- function(type){
 z<-plot_nodes[[type]];e<-plot_edges[[type]]
 pushViewport(viewport(xscale=c(0,network_W),yscale=c(0,network_H)))
 on.exit(popViewport(),add=TRUE)
 U<-function(v)unit(v,'native')
 for(j in seq_len(nrow(e))){
  a<-match(e$category_id[j],z$node_id);b<-match(e$gene[j],z$node_id)
  grid.lines(U(c(z$node_x_mm[a],z$node_x_mm[b])),U(c(z$node_y_mm[a],z$node_y_mm[b])),
   gp=gpar(col='#999999',lwd=1.05))
 }
 for(j in seq_len(nrow(z))){
  dx<-z$label_x_mm[j]-z$node_x_mm[j];dy<-z$label_y_mm[j]-z$node_y_mm[j]
  t<-min(if(dx!=0)z$label_width_mm[j]/2/abs(dx)else Inf,
         if(dy!=0)z$label_height_mm[j]/2/abs(dy)else Inf)
  ex<-z$label_x_mm[j]-dx*t;ey<-z$label_y_mm[j]-dy*t
  if(sqrt(dx^2+dy^2)>4)grid.lines(U(c(z$node_x_mm[j],ex)),U(c(z$node_y_mm[j],ey)),
   gp=gpar(col='#C1C1C1',lwd=.45))
 }
 for(j in seq_len(nrow(z)))grid.circle(U(z$node_x_mm[j]),U(z$node_y_mm[j]),
  r=unit(z$point_size_mm[j]/2,'mm'),gp=gpar(fill=z$node_color[j],col=NA))
 for(j in seq_len(nrow(z))){
  grid.text(z$display_label[j],U(z$label_x_mm[j]),U(z$label_y_mm[j]),
   gp=gpar(fontfamily='Arial',fontface=z$fontface[j],fontsize=z$label_size_mm[j]*72/25.4,lineheight=.91))
 }
}
draw_panel <- function(type){
 pushViewport(viewport(xscale=c(0,W),yscale=c(0,H),gp=gpar(fontfamily='Arial',fontface='bold')))
 on.exit(popViewport(),add=TRUE)
 U<-function(v)unit(v,'native')
 grid.rect(gp=gpar(fill='white',col=NA))
 grid.rect(x=U(W/2),y=U(H-8),width=U(W-4),height=U(12),
  gp=gpar(fill='white',col='#111111',lwd=2.2))
 grid.text(type,x=U(W/2),y=U(H-8),gp=gpar(fontsize=17,fontfamily='Arial'))
 pushViewport(viewport(x=U(1),y=U(4),just=c('left','bottom'),width=U(W-35),height=U(H-28)))
 draw_network_fixed(type)
 popViewport()
 x<-W-21;yt<-H-33
 grid.text('count',x=U(x),y=U(yt),just='left',gp=gpar(fontsize=11.5))
 fill<-if(type=='Downregulated')'#4487B0'else'#B51C2B'
 for(k in seq_along(count_breaks)){
  yy<-yt-8*k
  grid.circle(x=U(x+3),y=U(yy),r=unit(module_size(count_breaks[k])/2,'mm'),gp=gpar(fill=fill,col=NA))
  grid.text(as.character(count_breaks[k]),x=U(x+9),y=U(yy),just='left',gp=gpar(fontsize=10))
 }
 top<-yt-46;bottom<-top-47
 grid.text('log2FC',x=U(x),y=U(top+7),just='left',gp=gpar(fontsize=11.5))
 for(k in 1:200)grid.rect(x=U(x+2.5),y=U(bottom+(k-.5)/200*(top-bottom)),
  width=unit(4,'mm'),height=U((top-bottom)/200+.01),
  gp=gpar(fill=color(-fc_limit+(k-.5)/200*(2*fc_limit)),col=NA))
 ticks<-pretty(c(-fc_limit,fc_limit),4);ticks<-ticks[ticks>= -fc_limit&ticks<=fc_limit]
 for(v in ticks)grid.text(format(v,trim=TRUE),x=U(x+7),y=U(bottom+(v+fc_limit)/(2*fc_limit)*(top-bottom)),
  just='left',gp=gpar(fontsize=10))
}

draw_side <- function(type) {
  grid.newpage()
  draw_panel(type)
}
for (type in types) {
  export_plot(file.path(out, paste0(type, "_gene_category_network")), W / 25.4, H / 25.4,
              function() draw_side(type))
}
draw_combined <- function() {
  grid.newpage()
  for (i in seq_along(types)) {
    pushViewport(viewport(x = unit((i - 1) * (W + gap), "mm"), y = unit(0, "mm"),
                          just = c("left", "bottom"), width = unit(W, "mm"), height = unit(H, "mm")))
    draw_panel(types[i])
    popViewport()
  }
}
export_plot(file.path(out, "Combined_gene_category_network"), (2 * W + gap) / 25.4, H / 25.4, draw_combined)
