#!/usr/bin/env Rscript
# Render the selected 40-gene trajectories using the stored pooled means and BH20.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(here, "workflow_helpers.R"))
source(file.path(here, "figure_helpers.R"))
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: Rscript plot_signature_trajectories.R <trajectory_output_dir> <output_dir>")
load_figure_packages("ggplot2")
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
out <- new_output_dir(args[2])
tests <- read_table(file.path(args[1], "tables", "all20_tests.csv"))
pools <- read_table(file.path(args[1], "source_data", "all20_pool_means.csv"))
tt <- subset(tests, N == 40)
pp <- subset(pools, size_N == 40)
stopifnot(nrow(tests) == 20L, nrow(tt) == 2L, nrow(pp) == 18L, all(tt$status == "success"),
          isTRUE(all.equal(p.adjust(tests$Welch_P, "BH", n = 20L), tests$BH20, tolerance = 1e-12)))
write_table(pp, file.path(out, "selected_trajectory_pool_means.csv"))
write_table(tt, file.path(out, "selected_trajectory_statistics.csv"))
cols <- c(Ext='#003BFF',`Ret-Ext`='#FF2222');
for(letter in c('F','G')){
 t <- tt[tt$panel==letter,]; d<-pp[pp$panel==letter & pp$condition %in% names(cols),]
 d$condition <- factor(d$condition,levels=names(cols)); d<-d[order(d$condition,d$pool_id),]
 d$x <- as.numeric(d$condition)+rep(c(-.16,-.05,.16),2)
 s<-data.frame(condition=factor(names(cols),levels=names(cols)),x=1:2,
   mean=c(t$mean_Ext,t$mean_RetExt),SEM=c(t$SEM_Ext,t$SEM_RetExt))
 s$low<-s$mean-s$SEM;s$high<-s$mean+s$SEM
 stopifnot(abs(t.test(sample_mean_pseudotime~condition,data=d)$p.value-t$Welch_P)<1e-12)
 top<-max(d$sample_mean_pseudotime,s$high); upper<-if(letter=='F')2.8 else 3.7
 bracket<-top+.22; star_y<-bracket+.14
 stopifnot(all(d$sample_mean_pseudotime >= 0), star_y < upper)
 stars<-if(t$BH20<.01)'**'else if(t$BH20<.05)'*'else'ns'
 bottom_axis<-grobTree(
   segmentsGrob(x0=unit(10,'pt'),x1=unit(1,'npc'),y0=unit(-7,'pt'),y1=unit(-7,'pt'),gp=gpar(lwd=2.85,lineend='butt')),
   segmentsGrob(x0=unit(c((1-.45)/2.1,(2-.45)/2.1),'npc'),x1=unit(c((1-.45)/2.1,(2-.45)/2.1),'npc'),
     y0=unit(-7,'pt'),y1=unit(-11,'pt'),gp=gpar(lwd=2.85,lineend='butt')))
 p<-ggplot(d,aes(x,sample_mean_pseudotime,colour=condition))+
  geom_point(size=2.9)+
  geom_errorbar(data=s,aes(x=x,ymin=low,ymax=high,colour=condition),inherit.aes=FALSE,width=.16,linewidth=.8)+
  geom_point(data=s,aes(x=x,y=mean,colour=condition),inherit.aes=FALSE,shape=21,fill='white',size=4.9,stroke=1)+
  scale_colour_manual(values=cols)+scale_x_continuous(limits=c(.45,2.55),breaks=1:2,labels=names(cols),expand=expansion(mult=0))+
  scale_y_continuous(limits=c(0,upper),breaks=seq(0,upper,by=if(letter=='F').5 else 1),expand=expansion(mult=0))+
  annotate('segment',x=1,xend=1,y=bracket-.065,yend=bracket,linewidth=.55)+
  annotate('segment',x=1,xend=2,y=bracket,yend=bracket,linewidth=.55)+
  annotate('segment',x=2,xend=2,y=bracket,yend=bracket-.065,linewidth=.55)+
  annotate('text',x=1.5,y=star_y,label=stars,family='Arial',size=9.2,fontface='bold')+
  labs(x=NULL,y=paste0(if(letter=='F')'Ret-Ext-up'else'Ext-up',' signature\npseudotime'))+
  theme_classic(base_family='Arial',base_size=14.5)+theme(text=element_text(face='plain'),
   axis.text=element_text(size=14.5,colour='black',face='bold'),axis.title.y=element_text(size=16,face='bold',margin=margin(r=10)),
   axis.line=element_line(linewidth=1.0),axis.line.x=element_blank(),axis.ticks=element_line(linewidth=1.0),axis.ticks.x=element_blank(),
   axis.ticks.length=unit(4,'pt'),axis.text.x=element_text(margin=margin(t=11)),
   legend.position='none',plot.margin=margin(25,15,22,15),plot.background=element_rect(fill='white',colour=NA))+
  coord_cartesian(clip='off')+annotation_custom(bottom_axis)
 b<-ggplot_build(p);stopifnot(identical(as.numeric(b$data[[1]]$y),d$sample_mean_pseudotime),identical(as.character(b$data[[7]]$label),stars),b$layout$panel_params[[1]]$y.range[1]==0)
 stem <- if (letter == "F") "RetExt_signature_pseudotime" else "Ext_signature_pseudotime"
 export_plot(file.path(out, stem), 3.65, 4.7, function() print(p))
}
