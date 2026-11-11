#!/usr/bin/env Rscript
# Render the six pairwise DESeq2 contrasts using native BH-adjusted values.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(here, "workflow_helpers.R"))
source(file.path(here, "figure_helpers.R"))
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: Rscript plot_deseq2_volcanoes.R <deseq2_output_dir> <output_dir>")
load_figure_packages(c("ggplot2", "ggrepel"))
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(grid) })
out <- new_output_dir(args[2])
idx<-data.frame(figure=c('Figure_3',rep('Figure_S8',5)),panel=c('H','A','B','E','C','D'),
 kind=c(rep('Excitatory',3),rep('Inhibitory',3)),contrast=c('RetExt_vs_Ext','Ext_vs_WD','RetExt_vs_WD','RetExt_vs_Ext','Ext_vs_WD','RetExt_vs_WD'),
 width_in=5.1,height_in=c(4.34807692307692,4.24784810126582,4.28035714285714,4.692,4.09243902439024,4.4673417721519),stringsAsFactors=FALSE)
ok <- function(id, label, value) if (!isTRUE(value)) stop(id, ': ', label)
fc<-log2(1.15)
for(i in seq_len(nrow(idx))) {
 r<-idx[i,];id<-paste(r$kind,r$contrast,sep='__');path<-file.path(args[1],'datasets',id,'tables','ALL_genes.csv')
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
