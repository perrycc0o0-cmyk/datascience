#!/usr/bin/env Rscript
# Freeze public Ensembl116 rat-human homology independently of model inference.
args <- commandArgs(trailingOnly = TRUE)
i <- match("--root", args)
if (is.na(i) || i == length(args)) stop("Usage: Rscript prepare_orthology.R --root RUN_ROOT")
root <- normalizePath(args[i + 1L], mustWork = TRUE)
out <- file.path(root, "06_CellChat", "orthology")
if (file.exists(file.path(out, "complete.json"))) stop("A frozen orthology export already exists")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(root, "runtime", "biomart_cache")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(BIOMART_CACHE = cache)
options(timeout = 20)
stopifnot(requireNamespace("httr2", quietly = TRUE), requireNamespace("xml2", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
stamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
tryCatch({
  host <- "https://jun2026.archive.ensembl.org"
  endpoint <- paste0(host, "/biomart/martservice")
  fetch <- function(params, filename) {
    req <- do.call(httr2::req_url_query, c(list(httr2::request(endpoint)), params))
    response <- httr2::req_perform(httr2::req_timeout(req, 20))
    bytes <- httr2::resp_body_raw(response)
    writeBin(bytes, file.path(out, filename))
    list(text = rawToChar(bytes), url = httr2::resp_url(response), status = httr2::resp_status(response))
  }
  registry <- fetch(list(type = "registry"), "official_registry.xml")
  registry_xml <- xml2::read_xml(registry$text)
  entry <- xml2::xml_find_all(registry_xml, ".//MartURLLocation[@name='ENSEMBL_MART_ENSEMBL']")
  if (length(entry) != 1L || xml2::xml_attr(entry, "database") != "ensembl_mart_116" || xml2::xml_attr(entry, "displayName") != "Ensembl Genes 116") stop("Official registry is not Ensembl Genes116")
  release <- data.frame(version = 116, url = host, database = xml2::xml_attr(entry, "database"), display_name = xml2::xml_attr(entry, "displayName"))
  write.csv(release, file.path(out, "official_archive_inventory.csv"), row.names = FALSE)
  attribute_response <- fetch(list(type = "attributes", dataset = "rnorvegicus_gene_ensembl", mart = "ENSEMBL_MART_ENSEMBL"), "official_attributes.tsv")
  attributes <- read.delim(text = attribute_response$text, header = FALSE, quote = "", comment.char = "", fill = TRUE, stringsAsFactors = FALSE)
  names(attributes)[1L] <- "name"
  req <- c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_orthology_type")
  if (!all(req %in% attributes$name)) stop("The required official homology attributes are unavailable")
  confidence <- "hsapiens_homolog_orthology_confidence"
  selected <- c(req, confidence[confidence %in% attributes$name])
  query <- paste0('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="rnorvegicus_gene_ensembl" interface="default">',
    paste0('<Attribute name="', selected, '"/>', collapse = ""), '</Dataset></Query>')
  writeLines(query, file.path(out, "official_query.xml"), useBytes = TRUE)
  response <- fetch(list(query = query), "official_raw_response.tsv")
  if (grepl("Query ERROR|<html", response$text, ignore.case = TRUE)) stop("Official BioMart returned an error rather than a gene table")
  raw <- read.delim(text = response$text, header = FALSE, quote = "", comment.char = "", colClasses = "character", check.names = FALSE)
  if (ncol(raw) != length(selected)) stop("Unexpected column count in the official response")
  names(raw) <- selected
  names(raw)[match(req, names(raw))] <- c("rat_gene_id", "rat_symbol", "human_gene_id", "human_symbol", "orthology_type")
  if (confidence %in% names(raw)) names(raw)[names(raw) == confidence] <- "orthology_confidence"
  else raw$orthology_confidence <- NA_real_
  raw <- unique(raw)
  if (!nrow(raw) || !any(grepl("^ENSRNOG", raw$rat_gene_id))) stop("Empty/non-rat official response")
  saveRDS(raw, file.path(out, "rat_human_Ensembl116_raw.rds"), compress = "gzip")
  con <- gzfile(file.path(out, "rat_human_Ensembl116_raw.tsv.gz"), "wt")
  write.table(raw, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(con)
  eligible <- grepl("^ENSRNOG", raw$rat_gene_id) & raw$orthology_type == "ortholog_one2one" & !is.na(raw$human_gene_id) & nzchar(raw$human_gene_id) & !is.na(raw$human_symbol) & nzchar(raw$human_symbol)
  candidate <- unique(raw[which(eligible), , drop = FALSE])
  count_distinct <- function(value, key) vapply(split(value, key), function(v) length(unique(v)), integer(1))
  r_n <- count_distinct(candidate$human_gene_id, candidate$rat_gene_id)
  h_n <- count_distinct(candidate$rat_gene_id, candidate$human_gene_id)
  s_n <- count_distinct(candidate$human_gene_id, candidate$human_symbol)
  r_symbol_n <- count_distinct(candidate$rat_symbol, candidate$rat_gene_id)
  candidate$accepted <- r_n[candidate$rat_gene_id] == 1L & h_n[candidate$human_gene_id] == 1L & s_n[candidate$human_symbol] == 1L & r_symbol_n[candidate$rat_gene_id] == 1L
  write.table(candidate, file.path(out, "one2one_candidate_uniqueness_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  accepted <- candidate[which(candidate$accepted), setdiff(names(candidate), "accepted"), drop = FALSE]
  if (anyDuplicated(accepted$rat_gene_id)) stop("One-to-one evidence has duplicate rat IDs; no arbitrary record selection")
  if (nrow(accepted) < 1000L || anyDuplicated(accepted$rat_gene_id) || anyDuplicated(accepted$human_gene_id) || anyDuplicated(accepted$human_symbol)) stop("Invalid one-to-one export")
  write.table(accepted, file.path(out, "rat_human_Ensembl116_one2one.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  saveRDS(list(raw = raw, accepted = accepted), file.path(out, "orthology_lookup.rds"), compress = "gzip")
  files <- list.files(out, full.names = TRUE)
  info <- list(status = "complete", started_UTC = stamp, completed_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), release = 116L,
    archive = release, actual_mart_host = endpoint, registry_url = registry$url, table_response_url = response$url, dataset = "rnorvegicus_gene_ensembl", selected_attributes = selected,
    raw_rows = nrow(raw), rat_gene_ids = length(unique(raw$rat_gene_id)), one2one_rows = nrow(accepted),
    confidence_preserved_not_posthoc_filtered = TRUE, no_study_expression_uploaded = TRUE,
    source = "Official Ensembl BioMart XML endpoint; whole public rat gene homology table; registry verified Genes116", cache = cache, timeout_seconds_each_request = 20,
    transport = "httr2 public XML query; biomaRt automatic archive discovery was unavailable",
    artifact_md5 = as.list(tools::md5sum(files)), R = R.version.string, httr2 = as.character(utils::packageVersion("httr2")),
    session_info = capture.output(sessionInfo()))
  jsonlite::write_json(info, file.path(out, "complete.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
  cat("Frozen Ensembl116 one-to-one pairs:", nrow(accepted), "\n")
}, error = function(e) {
  failure_path <- if (file.exists(file.path(out, "failure.json"))) file.path(out, paste0("failure_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_", Sys.getpid(), ".json")) else file.path(out, "failure.json")
  jsonlite::write_json(list(status = "failed", started_UTC = stamp, failed_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), error = conditionMessage(e)), failure_path, pretty = TRUE, auto_unbox = TRUE)
  stop(e)
})
