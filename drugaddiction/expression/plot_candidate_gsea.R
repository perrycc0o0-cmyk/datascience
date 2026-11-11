#!/usr/bin/env Rscript
# Render nominal-P-selected candidate pathways; ridge fill represents candidate-family FDR.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(here, "workflow_helpers.R"))
source(file.path(here, "figure_helpers.R"))
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: Rscript plot_candidate_gsea.R <gsea_output_dir> <output_dir>")
load_figure_packages()
suppressPackageStartupMessages(library(grid))
out <- new_output_dir(args[2])
types <- c("Excitatory", "Inhibitory")
selected <- leading <- list()
for (type in types) {
  key <- paste0(type, "__RetExt_vs_Ext")
  full <- read_table(file.path(args[1], "tables", paste0(key, "_ALL43_GSEA.csv")))
  z <- read_table(file.path(args[1], "tables", paste0(key, "_P005_displayed_GSEA.csv")))
  finite <- is.finite(full$pvalue)
  expected <- full[finite & full$pvalue < 0.05, , drop = FALSE]
  expected <- expected[order(-expected$NES, expected$candidate_order), , drop = FALSE]
  stopifnot(nrow(full) == 43L, identical(z$ID, expected$ID),
            isTRUE(all.equal(z, expected, check.attributes = FALSE)),
            isTRUE(all.equal(full$BH[finite], p.adjust(full$pvalue[finite], "BH"), tolerance = 1e-12)))
  selected[[type]] <- z
  leading[[type]] <- read_table(file.path(args[1], "source_data", paste0(key, "_ALL43_leading_edge_genes.csv")))
}
all_selected <- do.call(rbind, selected)
values <- unlist(lapply(types, function(type) {
  x <- leading[[type]]
  x$log2FoldChange[x$ID %in% selected[[type]]$ID]
}))
xr <- range(c(0, values[is.finite(values)]))
pad <- max(0.1, 0.08 * diff(xr))
xr <- xr + c(-pad, pad)
logFDR <- -log10(pmax(all_selected$BH, .Machine$double.xmin))
flim <- c(0, max(0.25, ceiling(max(c(0, logFDR)) / 0.25) * 0.25))
ne <- max(0.5, ceiling(max(c(0, abs(all_selected$NES))) / 0.5) * 0.5)
nlim <- c(-ne, ne)
sizes <- if (nrow(all_selected)) unique(as.numeric(quantile(abs(all_selected$NES), c(0, 0.5, 1)))) else numeric()
palette_FDR <- colorRampPalette(c("#F3E7CA", "#9B6B21"))(1001)
palette_NES <- colorRampPalette(c("#588FC0", "#FFFFFF", "#ED7D68"))(1001)
color <- function(value,limits,palette){
 stopifnot(is.finite(value),value>=limits[1]-1e-10,value<=limits[2]+1e-10)
 palette[1L+round(1000*(value-limits[1])/diff(limits))]
}
drawridge <- function(z,lead,title,W,H){
 pushViewport(viewport(xscale=c(0,W),yscale=c(0,H),gp=gpar(fontfamily='Arial',fontface='bold')))
 on.exit(popViewport(),add=TRUE)
 U<-function(x)unit(x,'native')
 line<-function(x,y,...)grid.lines(U(x),U(y),...)
 txt<-function(label,x,y,...){
   dots<-list(...);if(is.null(dots$gp))dots$gp<-gpar()
   dots$gp$fontfamily<-'Arial';dots$gp$font<-2L
   do.call(grid.text,c(list(label=label,x=U(x),y=U(y)),dots))
 }
 grid.rect(gp=gpar(fill='white',col=NA))
 txt(title,W*.5,H-8,gp=gpar(fontsize=11.5,fontface='bold'))
 n<-nrow(z)
 if(!n){
   txt('No candidate pathway has P < 0.05.',W*.5,H*.5,gp=gpar(fontsize=10))
   return(invisible(list(density=data.frame(),modes=data.frame())))
 }
 L<-81;R<-132;X<-function(v)L+(v-xr[1])/diff(xr)*(R-L)
 yy<-if(n==1)(H+10)/2 else seq(H-35,33,length.out=n)
 step<-if(n>1)abs(diff(yy)[1])else min(25,H*.3)
 density_records<-mode_records<-list()
 for(i in seq_len(n)){
   row<-z[i,];d<-lead[lead$ID==row$ID,,drop=FALSE]
   value<-d$log2FoldChange[is.finite(d$log2FoldChange)]
   y<-yy[i];fill<-color(-log10(max(row$BH,.Machine$double.xmin)),flim,palette_FDR)
   label<-paste(strwrap(row$Description,width=38),collapse='\n')
   txt(label,67,y+.22*step,just='right',gp=gpar(fontsize=9.2))
   line(c(L,R),c(y,y),gp=gpar(col='black',lwd=1))
   radius<-2.34*sqrt(abs(row$NES)/ne)
   grid.circle(U(74),U(y+.22*step),r=unit(radius,'mm'),gp=gpar(fill=color(row$NES,nlim,palette_NES),col=NA))
   mode<-'none';bandwidth<-NA_real_
   if(length(unique(value))>=3L){
     bandwidth<-bw.nrd0(value)
     den<-density(value,bw=bandwidth,n=512,from=xr[1],to=xr[2])
     height<-den$y/max(den$y)*step*.67
     grid.polygon(U(c(X(den$x[1]),X(den$x),X(tail(den$x,1)))),U(c(y,y+height,y)),
       gp=gpar(fill=fill,col='black',lwd=.8))
     density_records[[i]]<-data.frame(ID=row$ID,x_log2FC=den$x,density=den$y,bandwidth,
       display_height_mm=height,NES=row$NES,P=row$pvalue,BH=row$BH,display_order=i)
     mode<-'density'
   }else if(length(value)){
     for(v in unique(value)){
       line(rep(X(v),2),c(y,y+step*.25),gp=gpar(lwd=1))
       grid.points(U(X(v)),U(y+step*.27),pch=21,size=unit(1.8,'mm'),gp=gpar(fill=fill,col='black'))
     }
     mode<-'rug_less_than3_unique_values'
   }else{
     txt('No leading-edge values',mean(c(L,R)),y+step*.2,gp=gpar(fontsize=7))
     mode<-'no_leading_edge_values'
   }
   mode_records[[i]]<-data.frame(ID=row$ID,Description=row$Description,display_order=i,
     ES=row$ES,NES=row$NES,P=row$pvalue,BH=row$BH,leading_edge_n=length(value),
     unique_log2FC_n=length(unique(value)),bandwidth,mode,dot_radius_mm=radius)
 }
 line(rep(X(0),2),c(26,H-17),gp=gpar(col='#707070',lty=2,lwd=1))
 line(c(L,R),c(22,22),gp=gpar(lwd=2))
 ticks<-pretty(xr,4);ticks<-ticks[ticks>=xr[1]&ticks<=xr[2]]
 for(t in ticks){
   line(rep(X(t),2),c(22,20),gp=gpar(lwd=2))
   txt(format(signif(t,3)),X(t),15.5,gp=gpar(fontsize=8.3))
 }
 txt(expression('Leading-edge gene log'[2]*'(Fold Change)'),mean(c(L,R)),6.5,gp=gpar(fontsize=9))
 bar<-function(palette,limits,y1,y2,title,ticks){
   txt(title,140,y2+6,just='left',gp=gpar(fontsize=9))
   for(j in 1:100)grid.rect(U(144),U(y1+(j-.5)*(y2-y1)/100),width=unit(2.5,'mm'),
     height=U((y2-y1)/100+.01),gp=gpar(fill=palette[1+round((j-.5)*1000/100)],col=NA))
   for(t in ticks)txt(format(signif(t,2)),148,y1+(t-limits[1])/diff(limits)*(y2-y1),
     just='left',gp=gpar(fontsize=8.3))
 }
 pt<-pretty(flim,3)
 bar(palette_FDR,flim,H*.65,H*.82,'-log10(FDR)',pt[pt>=flim[1]&pt<=flim[2]])
 bar(palette_NES,nlim,H*.38,H*.55,'NES',c(-ne,0,ne))
 txt('abs(NES)',140,H*.29,just='left',gp=gpar(fontsize=9))
 for(i in seq_along(sizes)){
   y<-H*(.24-.06*(i-1))
   grid.circle(U(144),U(y),r=unit(2.34*sqrt(sizes[i]/ne),'mm'),gp=gpar(fill='black',col=NA))
   txt(format(signif(sizes[i],2)),150,y,just='left',gp=gpar(fontsize=8.3))
 }
 invisible(list(density=if(length(density_records))do.call(rbind,density_records)else data.frame(),
   modes=do.call(rbind,mode_records)))
}

W <- 165.1
H <- max(180, max(vapply(selected, nrow, integer(1))) * 9.5 + 38)
gap <- 7
for (type in types) {
  stem <- paste0(type, "_RetExt_vs_Ext_candidate_GSEA")
  draw <- function() {
    grid.newpage()
    drawridge(selected[[type]], leading[[type]], paste0(type, " IEGhigh: Ret-Ext vs Ext"), W, H)
  }
  export_plot(file.path(out, stem), W / 25.4, H / 25.4, draw)
}
draw_combined <- function() {
  grid.newpage()
  for (i in seq_along(types)) {
    type <- types[i]
    pushViewport(viewport(x = unit((i - 1) * (W + gap), "mm"), y = unit(0, "mm"),
                          just = c("left", "bottom"), width = unit(W, "mm"), height = unit(H, "mm")))
    drawridge(selected[[type]], leading[[type]], paste0(type, " IEGhigh: Ret-Ext vs Ext"), W, H)
    grid.text(LETTERS[i], x = unit(6, "mm"), y = unit(H - 8, "mm"), just = "left",
              gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 13))
    popViewport()
  }
}
export_plot(file.path(out, "Combined_candidate_GSEA"), (2 * W + gap) / 25.4, H / 25.4, draw_combined)
