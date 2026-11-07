suppressPackageStartupMessages({ library(SingleCellExperiment); library(httr2) })
source("analysis/config.R")
out <- CFG$annotation_dir
dir.create(out, showWarnings = FALSE, recursive = TRUE)

e <- new.env(); nm <- load(CFG$local_rda, envir = e)
rd <- data.frame(row = rownames(e[[nm[1]]]), SYMBOL = as.character(rowData(e[[nm[1]]])$SYMBOL),
                 ENSEMBL = as.character(rowData(e[[nm[1]]])$ENSEMBL))
rm(e); invisible(gc())

sets <- rbind(
  data.frame(set = "core", symbol = utils::read.delim(CFG$core_file)$gene),
  data.frame(set = "hallmark_ifna", symbol = unique(utils::read.delim(CFG$hallmark_ifna)$gene_symbol)),
  data.frame(set = "hallmark_ifng", symbol = unique(utils::read.delim(CFG$hallmark_ifng)$gene_symbol)))

json <- function(req) tryCatch(
  resp_body_json(req_perform(req_retry(req, max_tries = 4)), simplifyVector = FALSE),
  error = function(e) stop("request failed: ", req$url, " - ", conditionMessage(e), call. = FALSE))

s1 <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x[[1]])
sj <- function(x) if (is.null(x) || !length(x)) "" else paste(unlist(x), collapse = ";")

syms <- unique(sets$symbol)
hgnc <- do.call(rbind, lapply(syms, function(s) {
  Sys.sleep(0.12)
  r <- json(request(paste0("https://rest.genenames.org/fetch/symbol/", utils::URLencode(s, reserved = TRUE))) |>
              req_headers(Accept = "application/json"))
  d <- r$response$docs
  if (length(d) != 1) return(data.frame(symbol = s, hgnc_n = length(d), hgnc_id = NA, approved = NA,
                                        prev_symbol = "", alias_symbol = "", hgnc_ensembl = NA, status = NA))
  d <- d[[1]]
  data.frame(symbol = s, hgnc_n = 1L, hgnc_id = s1(d$hgnc_id), approved = s1(d$symbol),
             prev_symbol = sj(d$prev_symbol), alias_symbol = sj(d$alias_symbol),
             hgnc_ensembl = s1(d$ensembl_gene_id), status = s1(d$status))
}))
hgnc <- merge(sets, hgnc, by = "symbol", sort = FALSE)

cand <- do.call(rbind, lapply(seq_len(nrow(hgnc)), function(i) {
  h <- hgnc[i,]
  names_i <- c(h$symbol, strsplit(h$prev_symbol, ";", fixed = TRUE)[[1]])
  hit <- rd$SYMBOL %in% names_i | (!is.na(h$hgnc_ensembl) & rd$ENSEMBL == h$hgnc_ensembl)
  if (!any(hit)) return(NULL)
  data.frame(set = h$set, symbol = h$symbol, rd[hit,])
}))

ids <- unique(cand$ENSEMBL)
info <- json(request("https://rest.ensembl.org/info/data") |> req_headers(Accept = "application/json"))
asm <- json(request("https://rest.ensembl.org/info/assembly/homo_sapiens") |> req_headers(Accept = "application/json"))
look <- json(request("https://rest.ensembl.org/lookup/id") |>
               req_headers(Accept = "application/json") |>
               req_body_json(list(ids = as.list(ids))))

parse_lookup <- function(look) do.call(rbind, lapply(ids, function(id) {
  x <- look[[id]]
  if (is.null(x)) return(data.frame(ENSEMBL = id, in_current = FALSE, seq_region = NA, biotype = NA,
                                    display_name = NA, ens_hgnc_id = NA, description = NA))
  d <- s1(x$description)
  acc <- regmatches(d, regexpr("Source:HGNC Symbol;Acc:(HGNC:)?[0-9]+", d))
  acc <- if (length(acc)) paste0("HGNC:", sub(".*Acc:(HGNC:)?", "", acc)) else NA
  data.frame(ENSEMBL = id, in_current = TRUE, seq_region = s1(x$seq_region_name), biotype = s1(x$biotype),
             display_name = s1(x$display_name), ens_hgnc_id = acc, description = d)
}))

ens <- parse_lookup(look)
look37 <- json(request("https://grch37.rest.ensembl.org/lookup/id") |>
                 req_headers(Accept = "application/json") |>
                 req_body_json(list(ids = as.list(ids))))
ens37 <- parse_lookup(look37)
names(ens37)[names(ens37) == "in_current"] <- "in_grch37"
ens37$in_current <- ens37$in_grch37
asm37 <- json(request("https://grch37.rest.ensembl.org/info/assembly/homo_sapiens") |> req_headers(Accept = "application/json"))
gone <- ens$ENSEMBL[!ens$in_current]
ens$archive_latest <- NA; ens$archive_is_current <- NA; ens$archive_release <- NA
if (length(gone)) {
  arc <- json(request("https://rest.ensembl.org/archive/id") |>
                req_headers(Accept = "application/json") |>
                req_body_json(list(id = as.list(gone))))
  for (a in arc) {
    k <- ens$ENSEMBL == a$id
    ens$archive_latest[k] <- s1(a$latest); ens$archive_is_current[k] <- s1(a$is_current)
    ens$archive_release[k] <- s1(a$release)
  }
}

utils::write.table(hgnc, file.path(out, "hgnc_sets.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(cand, file.path(out, "candidate_rows.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(ens, file.path(out, "ensembl_rows.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(ens37, file.path(out, "ensembl_rows_grch37.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
chroms <- function(a) unlist(lapply(a$top_level_region, function(r) if (identical(r$coord_system, "chromosome")) r$name))
primary <- chroms(asm)
primary37 <- intersect(chroms(asm37), c(1:22, "X", "Y", "MT"))
writeLines(c(sprintf("retrieved_utc\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
             sprintf("ensembl_release\t%s", paste(unlist(info$releases), collapse = ";")),
             sprintf("ensembl_assembly\t%s", asm$assembly_name),
             sprintf("primary_chromosomes\t%s", paste(primary, collapse = ";")),
             sprintf("grch37_assembly\t%s", asm37$assembly_name),
             sprintf("grch37_primary_chromosomes\t%s", paste(primary37, collapse = ";")),
             "grch37\thttps://grch37.rest.ensembl.org/lookup/id (POST); build reported on GEO GSE96583 (hg19)",
             "hgnc\thttps://rest.genenames.org/fetch/symbol/<symbol>",
             "ensembl\thttps://rest.ensembl.org/lookup/id (POST); archive/id for IDs not in current release"),
           file.path(out, "retrieval.tsv"))
message("hgnc records: ", nrow(hgnc), "; candidate rows: ", nrow(cand), "; ensembl ids: ", length(ids),
        "; not in current release: ", length(gone))
