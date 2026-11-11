#!/usr/bin/env Rscript
# Freeze accepted neuronal labels and expression inputs without changing upstream objects.
a <- commandArgs(TRUE); stopifnot(length(a)==1L)
out <- normalizePath(a[1],mustWork=TRUE)
suppressPackageStartupMessages(library(jsonlite))
c0 <- read_json(file.path(out,'config/analysis_contract.json'),simplifyVector=FALSE)
flat <- function(x) as.character(unlist(x,use.names=FALSE))
stage2 <- file.path(c0$base,c0$upstream_stage2)
rt <- read_json(file.path(stage2,'config/stage2_config.json'),simplifyVector=FALSE)$runtime
.libPaths(unique(c(flat(rt$r_libraries),.libPaths())))
Sys.setenv(OMP_NUM_THREADS='4',OMP_THREAD_LIMIT='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(digest)})
pools <- c(paste0('WD',1:3),paste0('Ext',1:3),paste0('Ret-Ext',1:3))
cond <- setNames(rep(c('WD','Ext','Ret-Ext'),each=3),pools); ieg <- flat(c0$IEG_genes)
preview <- file.path(stage2,'07_cluster6_to_Etv1_preview/source_data/preview_cell_mapping.csv.gz')
stopifnot(file.exists(preview))
cfg <- c0; cfg$pools <- as.list(pools);cfg$conditions <- as.list(unname(cond))
cfg$lineage_root <- stage2;cfg$frozen_annotation <- file.path(stage2,'03_annotation')
cfg$accepted_Inh_annotation_resolved <- preview
cfg$primary_counts <- 'RNA: SoupX-corrected integer UMI counts'
cfg$secondary_counts <- 'observed: original UMI counts; audit only'
write_json(cfg,file.path(out,'config/analysis_config.json'),pretty=TRUE,auto_unbox=TRUE)
for(d in c('01_inputs','logs','provenance'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
stamp <- function(...) {cat(format(Sys.time(),tz='UTC',usetz=TRUE),...,'\n');flush.console()}
status <- function(state,stage)write_json(list(status=state,stage=stage,UTC=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'provenance/01_input_status.json'),auto_unbox=TRUE,pretty=TRUE)
atomic_rds <- function(x,path){if(file.exists(path))stop('Refusing completed overwrite: ',path);tmp<-paste0(path,'.incomplete');saveRDS(x,tmp,compress=FALSE);stopifnot(file.rename(tmp,path))}
allmd<-audit<-counts<-observed<-poolmeta<-list();gene_reference<-NULL
tryCatch({
for(type in c('Excitatory','Inhibitory')){
 id<-if(type=='Excitatory')'Exc'else'Inh'; folder<-if(id=='Exc')'01_Exc'else'02_Inh'
 dest<-file.path(out,'01_inputs',type);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 input<-file.path(stage2,folder,'neuron_reclustered.rds');annfile<-file.path(stage2,'03_annotation',paste0(id,'_cell_annotations_UMAP.csv.gz'))
 sources<-c(input,annfile,if(id=='Inh')preview)
 before<-file.info(sources)[,c('size','mtime')]
 status('running',paste('read accepted',type));stamp('READ',type)
 ann<-read.csv(gzfile(annfile),stringsAsFactors=FALSE,check.names=FALSE)
 stopifnot(!anyDuplicated(ann$cell_id))
 ann$original_subtype_proposed<-ann$subtype_proposed
 ann$user_accepted_label_change<-FALSE
 if(id=='Inh'){
  z<-read.csv(gzfile(preview),stringsAsFactors=FALSE,check.names=FALSE);j<-match(ann$cell_id,z$cell_id)
  stopifnot(!anyNA(j),!anyDuplicated(z$cell_id),setequal(z$cell_id,ann$cell_id),all(as.character(z$subcluster[j])==as.character(ann$subcluster)),all(z$sample_id[j]==ann$sample_id),all(z$condition_label[j]==ann$condition_label),all(z$original_label[j]==ann$subtype_proposed),max(abs(z$UMAP_1[j]-ann$UMAP_1))<1e-10,max(abs(z$UMAP_2[j]-ann$UMAP_2))<1e-10)
  ann$subtype_proposed<-z$preview_label[j]
  ann$user_accepted_label_change<-ann$subtype_proposed!=ann$original_subtype_proposed
  stopifnot(sum(ann$user_accepted_label_change)==1377L,all(ann$subcluster[ann$user_accepted_label_change]==6),sum(ann$subtype_proposed=='GABA_Etv1')==8139L,sum(ann$subtype_proposed=='GABA_Lamp5')==1085L)
 }
 o<-readRDS(input);cells<-colnames(o);md<-o[[]];idx<-match(cells,ann$cell_id)
 stopifnot(!anyNA(idx),!anyDuplicated(cells),setequal(cells,ann$cell_id),all(ann$sample_id[idx]==md$sample_id),all(ann$condition_label[idx]==md$condition_label),all(as.character(ann$subcluster[idx])==as.character(md$subcluster)),ncol(o)==if(id=='Exc')62991L else 21693L)
 ann<-ann[idx,,drop=FALSE]
 x<-LayerData(o[['RNA']],layer='data',fast=FALSE);rna<-LayerData(o[['RNA']],layer='counts',fast=FALSE)
 stopifnot(identical(colnames(x),cells),identical(dimnames(x),dimnames(rna)),all(ieg%in%rownames(x)),inherits(x,'sparseMatrix'),inherits(rna,'sparseMatrix'),all(is.finite(rna@x)),all(rna@x>=0),all(rna@x==round(rna@x)))
 genes<-rownames(x); fm<-o[['RNA']][[]]
 fm<-data.frame(seurat_feature=genes,fm[genes,intersect(c('gene_id','original_symbol','duplicate_symbol'),names(fm)),drop=FALSE],row.names=NULL)
 stopifnot(all(c('seurat_feature','gene_id','original_symbol')%in%names(fm)))
 if(is.null(gene_reference)){gene_reference<-fm;write.csv(fm,file.path(out,'01_inputs/gene_annotation.csv'),row.names=FALSE)}else stopifnot(identical(fm,gene_reference))
 score<-as.numeric(Matrix::colSums(x[ieg,,drop=FALSE]));threshold<-as.numeric(quantile(score,c0$IEG_quantile,type=as.integer(c0$IEG_quantile_type)));high<-score>=threshold
 direct<-log1p(sweep(as.matrix(rna[ieg,,drop=FALSE]),2,Matrix::colSums(rna),'/')*10000)
 normalized_max_error<-max(abs(direct-as.matrix(x[ieg,,drop=FALSE])));stopifnot(normalized_max_error<1e-7);rm(direct)
 m<-data.frame(Cell=cells,neuron_type=type,pool_id=as.character(md$sample_id),condition=as.character(md$condition_label),new_cluster=as.character(md$subcluster),subtype=ann$subtype_proposed,biological_annotation=ann$subtype_proposed,passes_lineage_conflict_review=!ann$lineage_identity_warning,annotation_review_required=ann$lineage_identity_warning,ieg_sum=score,IEGhigh=high,UMAP_1=ann$UMAP_1,UMAP_2=ann$UMAP_2,stringsAsFactors=FALSE)
 review_cols<-c('source_cell_id','source_sample_id','initial_major_class','major_class_six','forced_display_only','vascular_review_flag','original_global_cluster','original_global_seurat_clusters','lineage_identity_warning','lineage_warning_reason','subtype_annotation_status','subtype_low_score_margin','Amb_container_reason','original_subtype_proposed','user_accepted_label_change','scDblFinder_class','final_call')
 m<-cbind(m,ann[,intersect(review_cols,names(ann)),drop=FALSE])
 stopifnot(all(m$condition==unname(cond[m$pool_id])),all(is.finite(score)),!anyNA(m[,c('Cell','neuron_type','pool_id','condition','new_cluster','subtype','IEGhigh')]))
 write.csv(m,file.path(dest,'neuron_metadata.csv'),row.names=FALSE)
 write.csv(data.frame(Cell=cells,t(as.matrix(x[ieg,,drop=FALSE])),check.names=FALSE),file.path(dest,'IEG_expression_per_nucleus.csv'),row.names=FALSE)
 write.csv(as.data.frame(table(m$pool_id,m$subtype)),file.path(dest,'accepted_pool_subtype_counts.csv'),row.names=FALSE)
 for(pop in c('All','IEGhigh','IEGlow')){
  keep<-switch(pop,All=rep(TRUE,length(cells)),IEGhigh=high,IEGlow=!high)
  z<-sparseMatrix(i=which(keep),j=match(m$pool_id[keep],pools),x=1,dims=c(length(cells),length(pools)),dimnames=list(cells,pools))
  pb<-as.matrix(rna%*%z);storage.mode(pb)<-'integer';n<-as.integer(table(factor(m$pool_id[keep],levels=pools)))
  pm<-data.frame(pool_id=pools,condition=unname(cond[pools]),n_high=n,total_counts=as.numeric(colSums(pb)),population=pop,neuron_type=type,stringsAsFactors=FALSE)
  stopifnot(all(pb>=0),all(n>0),sum(n)==sum(keep))
  write.csv(data.frame(gene=genes,pb,check.names=FALSE),file.path(dest,paste0(pop,'_corrected_pseudobulk_counts.csv')),row.names=FALSE)
  write.csv(pm,file.path(dest,paste0(pop,'_pool_metadata.csv')),row.names=FALSE)
  counts[[paste(type,pop,sep='__')]]<-pb;poolmeta[[paste(type,pop,sep='__')]]<-pm
 }
 atomic_rds(list(normalized=x[,high,drop=FALSE],counts=rna[,high,drop=FALSE],metadata=m[high,,drop=FALSE],gene_annotation=fm,IEG_threshold=threshold),file.path(dest,'IEGhigh_payload.rds'))
 atomic_rds(x,file.path(dest,'all_neuron_normalized_RNA.rds'))
 obs<-LayerData(o[['observed']],layer='counts',fast=FALSE);stopifnot(identical(dimnames(obs),dimnames(rna)),all(obs@x>=0),all(obs@x==round(obs@x)))
 z<-sparseMatrix(i=which(high),j=match(m$pool_id[high],pools),x=1,dims=c(length(cells),length(pools)))
 pb_obs<-as.matrix(obs%*%z);dimnames(pb_obs)<-list(genes,pools);storage.mode(pb_obs)<-'integer'
 write.csv(data.frame(gene=genes,pb_obs,check.names=FALSE),file.path(dest,'IEGhigh_observed_pseudobulk_counts_audit.csv'),row.names=FALSE);observed[[type]]<-pb_obs
 audit[[type]]<-data.frame(neuron_type=type,n_total=nrow(m),n_high=sum(high),IEG_threshold=threshold,quantile=.95,quantile_type=7,inclusive_comparison='>=',scope='all conditions combined, separately within neuron class',available_IEGs=paste(ieg,collapse=';'),normalized_RNA_max_error=normalized_max_error,input_RDS=input,input_bytes=before$size[1],input_mtime_UTC=format(before$mtime[1],tz='UTC',usetz=TRUE),annotation_csv=annfile,annotation_MD5=unname(tools::md5sum(annfile)),accepted_override_MD5=if(id=='Inh')unname(tools::md5sum(preview))else NA_character_,matrix_RNA_data_SHA256=digest(x,algo='sha256'),matrix_RNA_counts_SHA256=digest(rna,algo='sha256'),stringsAsFactors=FALSE)
 allmd[[type]]<-m;stopifnot(identical(before,file.info(sources)[,c('size','mtime')]))
 stamp('PREPARED',type,nrow(m),sum(high),'threshold',threshold)
 rm(o,x,rna,obs,md,ann,m,z,pb,pb_obs,fm);gc()
}
allcols<-unique(unlist(lapply(allmd,names),use.names=FALSE))
allmd<-lapply(allmd,function(z){for(k in setdiff(allcols,names(z)))z[[k]]<-NA;z[,allcols,drop=FALSE]})
meta<-do.call(rbind,allmd);rownames(meta)<-NULL;stopifnot(nrow(meta)==84684L,!anyDuplicated(meta$Cell))
write.csv(meta,file.path(out,'01_inputs/neuron_metadata.csv'),row.names=FALSE)
write.csv(do.call(rbind,audit),file.path(out,'01_inputs/IEGhigh_definition_audit.csv'),row.names=FALSE)
atomic_rds(list(counts=counts,pool_metadata=poolmeta,observed_IEGhigh_counts=observed),file.path(out,'01_inputs/pseudobulk_inputs.rds'))
write_json(list(status='complete',n_neuron_nuclei=nrow(meta),n_IEGhigh=sum(meta$IEGhigh),no_reclustering=TRUE,no_old_results_overwritten=TRUE,source_stage2=stage2,accepted_Inh_annotation=preview,primary_counts=cfg$primary_counts,session_info=capture.output(sessionInfo())),file.path(out,'provenance/01_input_complete.json'),auto_unbox=TRUE,pretty=TRUE)
status('complete','frozen accepted inputs');stamp('INPUT_PREPARATION_COMPLETE')
},error=function(e){status('failed',conditionMessage(e));stamp('FAILED',conditionMessage(e));quit(status=1)})
