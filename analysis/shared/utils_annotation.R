map_genes <- function(symbols, row_symbols, row_ids) {
  n <- vapply(symbols, function(g) sum(row_symbols == g, na.rm = TRUE), integer(1))
  row <- vapply(symbols, function(g) {
    i <- which(row_symbols == g)
    if (length(i) == 1) row_ids[i] else NA_character_
  }, character(1))
  data.frame(gene = symbols, n_matches = unname(n), row = unname(row),
             status = ifelse(n == 1, "mapped", ifelse(n == 0, "missing", "ambiguous")))
}

resolve_symbols <- function(cand, hgnc, ens, primary) {
  h <- unique(hgnc[, c("symbol", "hgnc_id", "hgnc_ensembl")])
  x <- merge(unique(cand[, c("symbol", "row", "SYMBOL", "ENSEMBL")]), h, by = "symbol")
  x <- merge(x, ens[, c("ENSEMBL", "in_current", "seq_region", "ens_hgnc_id")], by = "ENSEMBL", all.x = TRUE)
  cur <- x$in_current %in% TRUE
  same <- cur &
    !is.na(x$ens_hgnc_id) &
    x$ens_hgnc_id == x$hgnc_id
  prim <- x$seq_region %in% primary
  hgnc_row <- !is.na(x$hgnc_ensembl) & x$ENSEMBL == x$hgnc_ensembl
  x$class <- ifelse(!cur, "unresolved_not_in_current_ensembl",
                    ifelse(is.na(x$ens_hgnc_id), "unresolved_no_hgnc_xref",
                           ifelse(!same, "distinct_gene",
                                  ifelse(hgnc_row & prim, "selected",
                                         ifelse(!prim, "alt_annotation", "unresolved_same_gene_primary_duplicate")))))
  sym <- do.call(rbind, lapply(split(x, x$symbol), function(s) {
    sel <- s$row[s$class == "selected"]
    peers <- if (length(sel) == 1) s[s$SYMBOL == s$SYMBOL[s$row == sel] & s$row != sel,] else s[0,]
    ok <- length(sel) == 1 && all(peers$class %in% c("alt_annotation", "distinct_gene"))
    data.frame(symbol = s$symbol[1], hgnc_id = s$hgnc_id[1], status = if (ok) "resolved" else "unresolved",
               row = if (ok) sel else NA_character_,
               reason = if (ok) "" else if (length(sel) != 1) "no unique HGNC-ID row on a primary chromosome"
               else paste("unconfirmed same-symbol row(s):",
                          paste(peers$row[!peers$class %in% c("alt_annotation", "distinct_gene")], collapse = ";")))
  }))
  missing <- setdiff(h$symbol, sym$symbol)
  if (length(missing))
    sym <- rbind(sym, data.frame(symbol = missing, hgnc_id = h$hgnc_id[match(missing, h$symbol)],
                                 status = "unresolved", row = NA_character_, reason = "no candidate row in dataset"))
  rownames(x) <- rownames(sym) <- NULL
  list(rows = x[order(x$symbol, x$class != "selected", x$row),], symbols = sym[order(sym$symbol),])
}

linked_rows <- function(res) {
  r <- res$rows
  unique(r$row[r$class %in% c("selected", "alt_annotation") | r$ENSEMBL %in% r$hgnc_ensembl])
}

load_annotation <- function(dir = CFG$annotation_dir) {
  rd <- function(f) utils::read.delim(file.path(dir, f), stringsAsFactors = FALSE)
  meta <- utils::read.delim(file.path(dir, "retrieval.tsv"), header = FALSE, stringsAsFactors = FALSE)
  chr <- function(k) strsplit(meta$V2[meta$V1 == k], ";", fixed = TRUE)[[1]]
  list(hgnc = rd("hgnc_sets.tsv"), cand = rd("candidate_rows.tsv"),
       ens = list(grch38 = rd("ensembl_rows.tsv"), grch37 = rd("ensembl_rows_grch37.tsv")),
       primary = list(grch38 = chr("primary_chromosomes"), grch37 = chr("grch37_primary_chromosomes")),
       meta = meta)
}

resolve_set <- function(ann, set, build = "grch37") {
  resolve_symbols(ann$cand[ann$cand$set == set,], ann$hgnc[ann$hgnc$set == set,],
                  ann$ens[[build]], ann$primary[[build]])
}

hallmark_linked_rows <- function(ann, sets = c("hallmark_ifna", "hallmark_ifng"))
  unique(unlist(lapply(sets, function(s) lapply(c("grch37", "grch38"), function(b)
    linked_rows(resolve_set(ann, s, b))))))
