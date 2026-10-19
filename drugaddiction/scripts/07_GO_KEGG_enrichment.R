suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
})
source('R/common.R')

ctx <- pipeline_context()
cfg <- ctx$cfg
out <- new_stage(ctx, '07_enrichment')
orgdb <- org.Rn.eg.db
available_symbols <- AnnotationDbi::keys(orgdb, keytype = 'SYMBOL')
contrasts <- c('Ext_vs_WD', 'RetExt_vs_WD', 'RetExt_vs_Ext')
ontologies <- c('BP', 'CC', 'MF')
kegg_links <- KEGGREST::keggLink('pathway', 'rno')
kegg_membership <- unique(data.frame(
  ENTREZID = sub('^rno:', '', names(kegg_links)),
  pathway = sub('^path:', '', unname(kegg_links)), stringsAsFactors = FALSE
))
write_table(kegg_membership, file.path(out, 'KEGG_gene_pathway_membership.csv'))

empty_terms <- function() data.frame(
  ID = character(), Description = character(), GeneRatio = character(), BgRatio = character(),
  pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(), geneID = character(),
  Count = integer(), direction = character(), database = character(), stringsAsFactors = FALSE
)

term_table <- function(fit, direction, database, all_terms = TRUE) {
  if (is.null(fit)) return(empty_terms())
  x <- if (all_terms) fit@result else as.data.frame(fit)
  if (!nrow(x)) return(empty_terms())
  x$direction <- direction
  x$database <- database
  rownames(x) <- NULL
  x
}

write_terms <- function(fit, table_dir, direction, database) {
  all <- term_table(fit, direction, database, all_terms = TRUE)
  significant <- term_table(fit, direction, database, all_terms = FALSE)
  if (nrow(significant)) {
    significant <- significant[!is.na(significant$p.adjust) &
                                 significant$p.adjust <= as.numeric(cfg$ora$padj), , drop = FALSE]
  }
  prefix <- paste(direction, database, sep = '_')
  write_table(all, file.path(table_dir, paste0(prefix, '.csv')))
  write_table(significant, file.path(table_dir, paste0(prefix, '_significant.csv')))
  invisible(NULL)
}

for (lineage in c('Excitatory', 'Inhibitory')) {
  for (comparison in contrasts) {
    deg <- utils::read.csv(stage_path(ctx, '05_DESeq2', 'tables', lineage,
                                     paste0(comparison, '_all_genes.csv')),
                           stringsAsFactors = FALSE, check.names = FALSE)
    required <- c('gene', 'tested', 'padj', 'log2FoldChange')
    if (!all(required %in% names(deg))) stop('Incomplete differential-expression table')
    table_dir <- file.path(out, 'tables', lineage, comparison)
    background <- unique(as.character(deg$gene[deg$tested %in% TRUE]))
    background_symbols <- intersect(background, available_symbols)
    if (!length(background_symbols)) stop('No tested symbols map to org.Rn.eg.db')
    mappings <- clusterProfiler::bitr(background_symbols, fromType = 'SYMBOL', toType = 'ENTREZID',
                                      OrgDb = orgdb, drop = FALSE)
    mappings <- mappings[!is.na(mappings$ENTREZID) & nzchar(mappings$ENTREZID), , drop = FALSE]
    mappings <- unique(mappings[, c('SYMBOL', 'ENTREZID'), drop = FALSE])
    write_table(data.frame(gene = background), file.path(table_dir, 'tested_gene_universe.csv'))
    write_table(mappings, file.path(table_dir, 'SYMBOL_ENTREZ_mapping.csv'))
    go_mapping <- AnnotationDbi::select(orgdb, keys = background_symbols, keytype = 'SYMBOL',
                                        columns = c('GOALL', 'ONTOLOGYALL'))
    go_universes <- lapply(ontologies, function(ontology) {
      unique(go_mapping$SYMBOL[!is.na(go_mapping$GOALL) & go_mapping$ONTOLOGYALL %in% ontology])
    })
    names(go_universes) <- ontologies
    for (ontology in ontologies) {
      write_table(data.frame(gene = go_universes[[ontology]]),
                  file.path(table_dir, paste0('GO_', ontology, '_mapped_universe.csv')))
    }
    kegg_universe <- intersect(unique(as.character(mappings$ENTREZID)), kegg_membership$ENTREZID)
    write_table(data.frame(ENTREZID = kegg_universe), file.path(table_dir, 'KEGG_mapped_universe.csv'))
    eligible <- !is.na(deg$padj) & deg$padj <= as.numeric(cfg$de$padj) &
      !is.na(deg$log2FoldChange)
    for (direction in c('Up', 'Down')) {
      selected <- eligible & if (direction == 'Up') deg$log2FoldChange > 0 else deg$log2FoldChange < 0
      genes <- unique(as.character(deg$gene[selected]))
      write_table(data.frame(gene = genes), file.path(table_dir, paste0(direction, '_input_genes.csv')))
      for (ontology in ontologies) {
        universe <- go_universes[[ontology]]
        mapped_genes <- intersect(genes, universe)
        fit <- if (length(mapped_genes) && length(universe)) clusterProfiler::enrichGO(
          gene = mapped_genes, OrgDb = orgdb, keyType = 'SYMBOL', ont = ontology,
          universe = universe, pAdjustMethod = 'BH',
          pvalueCutoff = as.numeric(cfg$ora$pvalue_cutoff),
          qvalueCutoff = as.numeric(cfg$ora$qvalue_cutoff),
          minGSSize = as.integer(cfg$ora$min_gs_size), maxGSSize = as.integer(cfg$ora$max_gs_size),
          readable = FALSE, pool = FALSE
        ) else NULL
        write_terms(fit, table_dir, direction, paste0('GO_', ontology))
        if (!is.null(fit)) save_object(fit, file.path(out, 'objects', lineage, comparison,
                                                     paste0(direction, '_GO_', ontology, '.rds')))
      }
      kegg_genes <- unique(as.character(mappings$ENTREZID[mappings$SYMBOL %in% genes]))
      fit <- if (length(kegg_genes) && length(kegg_universe)) clusterProfiler::enrichKEGG(
        gene = kegg_genes, organism = 'rno', keyType = 'kegg', universe = kegg_universe,
        pAdjustMethod = 'BH', pvalueCutoff = as.numeric(cfg$ora$pvalue_cutoff),
        qvalueCutoff = as.numeric(cfg$ora$qvalue_cutoff),
        minGSSize = as.integer(cfg$ora$min_gs_size), maxGSSize = as.integer(cfg$ora$max_gs_size),
        use_internal_data = FALSE
      ) else NULL
      write_terms(fit, table_dir, direction, 'KEGG')
      if (!is.null(fit)) save_object(fit, file.path(out, 'objects', lineage, comparison,
                                                   paste0(direction, '_KEGG.rds')))
    }
  }
}
write_session(file.path(out, 'sessionInfo.txt'))
