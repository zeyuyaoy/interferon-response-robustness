for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R",
            "analysis/time_course/scripts/utils_gse226572.R")) source(f)
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "analysis/time_course/scripts"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
D <- "analysis/time_course"; OUT <- "results/time_course"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LABELS_FILE <- "cache/gse226572/cell_labels.tsv.gz"
PRIMARY <- c("0h", "4h", "24h"); MIN_CELLS <- 20; FLAG_CELLS <- 50; FOLD_FLAG <- 10
chk <- list()
add <- function(check, pass, value = NA) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))
out <- function(df, name) write_result(df, name, RUN_ID, dir = OUT)

s <- utils::read.delim(file.path(D, "samples.tsv"), stringsAsFactors = FALSE)
dl <- utils::read.delim(file.path(D, "downloads.tsv"), colClasses = "character")
s <- merge(s, dl[, c("gsm", "file", "sha256")], by = "gsm")
s <- s[order(s$donor, s$time_h, s$library_rep),]
add("24 GEO samples with recorded downloads", nrow(s) == 24, nrow(s))
add("sha256 of every file equals downloads.tsv", all(unname(tools::sha256sum(s$file)) == s$sha256))
add("titles parse to donor, time and library replicate", all(!is.na(s$donor) &
                                                               !is.na(s$time_h) &
                                                               s$library_rep %in% 1:2))
add("3 donors", length(unique(s$donor)) == 3, paste(unique(s$donor), collapse = ";"))
add("untreated samples have no treatment field and stimulated samples do",
    all(is.na(s$treatment_field[s$untreated]) | s$treatment_field[s$untreated] == "") &&
      all(grepl("IFN", s$treatment_field[!s$untreated])))
add("only 0h has a second library", all(s$library_rep[s$time_h != 0] == 1) && all(table(s$donor[s$time_h == 0]) == 2))
add("every donor has 0h, 4h and 24h", all(tapply(s$time, s$donor, function(t) all(PRIMARY %in% t))))

core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
cand <- utils::read.delim("analysis/inputs/annotation/candidate_rows.tsv", stringsAsFactors = FALSE)
cand <- cand[cand$set == "core",]
mk <- unique(unlist(MARKERS))
ifn <- unique(c(core, utils::read.delim("analysis/inputs/hallmark_ifna.tsv")$gene_symbol, utils::read.delim("analysis/inputs/hallmark_ifng.tsv")$gene_symbol))
add("no marker gene is a panel or Hallmark IFN gene", !any(mk %in% ifn), paste(intersect(mk, ifn), collapse = ";"))

G <- NULL; cells <- list(); files <- list(); attrs <- list()
for (i in seq_len(nrow(s))) {
  x <- read_raw_h5(s$file[i])
  if (is.null(G)) {
    G <- x$features
    mp <- map_genes(core, G$name, G$id)
    ens_ok <- all(mp$row == cand$ENSEMBL[match(mp$gene, cand$symbol)])
    add("10 panel genes uniquely mapped by symbol", all(mp$status == "mapped"), sum(mp$status == "mapped"))
    add("panel gene Ensembl ids equal the frozen annotation snapshot", ens_ok)
    mm <- map_genes(mk, G$name, G$id)
    add("all marker genes uniquely mapped", all(mm$status == "mapped"), paste(mm$gene[mm$status != "mapped"], collapse = ";"))
    mrow <- stats::setNames(match(mm$row, G$id), mm$gene)
  }
  cc <- call_cells(x$counts, x$features$name)
  keep <- which(cc$called)
  sub <- x$counts[mrow, keep, drop = FALSE]
  rownames(sub) <- names(mrow)
  S <- marker_scores(sub, cc$total[keep])
  a <- assign_labels(S, sub["CD14",], Matrix::colSums(sub[T_GENES, , drop = FALSE]))
  cells[[i]] <- data.frame(gsm = s$gsm[i], donor = s$donor[i], time = s$time[i], library_rep = s$library_rep[i],
                           barcode = cc$barcode[keep], n_genes = cc$n_genes[keep], total = cc$total[keep],
                           mt_frac = round(cc$mt_frac[keep], 5), label = a$label, margin = round(a$margin, 5),
                           stringsAsFactors = FALSE)
  med <- vapply(names(mrow), function(g) tapply(as.numeric(log1p(1e4 * sub[g,] / cc$total[keep])), a$label, stats::median), numeric(length(unique(a$label))))
  attrs[[i]] <- data.frame(gsm = s$gsm[i], t(x$attrs), check.names = FALSE)
  files[[i]] <- data.frame(gsm = s$gsm[i], title = s$title[i], n_barcodes = ncol(x$counts), n_called = length(keep),
                           median_genes_called = stats::median(cc$n_genes[keep]), median_umis_called = stats::median(cc$total[keep]),
                           counts_integer = counts_valid(x$counts), features_identical = identical(x$features$id, G$id),
                           barcodes_unique = !anyDuplicated(colnames(x$counts)),
                           gene_expression_only = all(x$features$feature_type == "Gene Expression"))
  rm(x, sub); invisible(gc())
  message(s$gsm[i], " ", s$title[i], ": ", length(keep), " cells")
}
fl <- do.call(rbind, files)
at <- do.call(rbind, lapply(attrs, function(a) { a[setdiff(unique(unlist(lapply(attrs, names))), names(a))] <- NA; a }))
add("integer counts in every file", all(fl$counts_integer))
add("feature ids identical across files", all(fl$features_identical))
add("barcodes unique within every file", all(fl$barcodes_unique))
add("only Gene Expression features", all(fl$gene_expression_only))
chem_col <- intersect(c("chemistry_description", "chemistry"), names(at))
if (length(chem_col)) {
  ch <- tapply(at[[chem_col[1]]], s$donor[match(at$gsm, s$gsm)], function(v) length(unique(v)))
  add("one chemistry within each donor", all(ch == 1), paste(unique(at[[chem_col[1]]]), collapse = ";"))
} else add("chemistry attribute present in H5", FALSE, paste(names(at), collapse = ";"))
add("each donor's libraries have the same SRA run layout", all(tapply(s$n_runs, s$donor, function(v) length(unique(v))) == 1),
    paste(tapply(s$n_runs, s$donor, function(v) paste(unique(v), collapse = "/")), collapse = ";"))

cl <- do.call(rbind, cells)
cl$lineage <- ifelse(cl$label %in% names(LINEAGE_OF), LINEAGE_OF[cl$label], NA)
gz <- gzfile(LABELS_FILE, "w"); utils::write.table(cl, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

lab_tab <- as.data.frame(table(gsm = cl$gsm, label = cl$label), stringsAsFactors = FALSE)
lab_tab <- merge(lab_tab, s[, c("gsm", "title", "donor", "time", "time_h", "library_rep")], by = "gsm")
lab_tab$frac <- lab_tab$Freq / ave(lab_tab$Freq, lab_tab$gsm, FUN = sum)
lab_tab <- lab_tab[order(lab_tab$donor, lab_tab$time_h, lab_tab$library_rep, lab_tab$label),]
out(lab_tab, "annotation_label_fractions.csv")

grp <- aggregate(list(n_cells = cl$lineage), list(donor = cl$donor, lineage = cl$lineage, time = cl$time), length)
grp$time_h <- as.integer(sub("h", "", grp$time))
grp <- grp[order(grp$lineage, grp$donor, grp$time_h),]
grp$primary <- grp$time %in% PRIMARY
grp$below_flag <- grp$n_cells < FLAG_CELLS
out(grp, "cell_counts.csv")
lib0 <- aggregate(list(n_cells = cl$lineage[cl$time == "0h"]), list(donor = cl$donor[cl$time == "0h"], lineage = cl$lineage[cl$time == "0h"],
                                                                    library_rep = cl$library_rep[cl$time == "0h"]), length)
out(lib0, "cell_counts_0h_by_library.csv")
pg <- expand.grid(donor = sort(unique(s$donor)), lineage = c("mono", "nk"), time = PRIMARY, stringsAsFactors = FALSE)
pg$n <- grp$n_cells[match(paste(pg$donor, pg$lineage, pg$time), paste(grp$donor, grp$lineage, grp$time))]
pg$n[is.na(pg$n)] <- 0L
add(sprintf("every primary donor x lineage x time group has >= %d cells", MIN_CELLS), all(pg$n >= MIN_CELLS), min(pg$n))
add(sprintf("every 0h library has >= %d cells per lineage", MIN_CELLS), all(lib0$n_cells >= MIN_CELLS) && nrow(lib0) == 12, min(lib0$n_cells))

fr <- lab_tab[lab_tab$label %in% names(LINEAGE_OF),]
fr$lineage <- LINEAGE_OF[fr$label]
f0 <- aggregate(list(f0 = fr$frac[fr$time == "0h"]), list(donor = fr$donor[fr$time == "0h"], lineage = fr$lineage[fr$time == "0h"]), mean)
fr <- merge(fr, f0, by = c("donor", "lineage"))
fr$fold_vs_0h <- fr$frac / fr$f0
fr$flag_fold <- fr$fold_vs_0h > FOLD_FLAG | fr$fold_vs_0h < 1 / FOLD_FLAG
out(fr[order(fr$lineage, fr$donor, fr$time_h, fr$library_rep), c("donor", "lineage", "title", "time", "library_rep", "Freq", "frac", "f0", "fold_vs_0h", "flag_fold")],
    "lineage_fraction_stability.csv")
add(sprintf("no primary-time lineage fraction differs > %d-fold from 0h", FOLD_FLAG), !any(fr$flag_fold[fr$time %in% PRIMARY]),
    sum(fr$flag_fold[fr$time %in% PRIMARY]))

mq <- do.call(rbind, lapply(split(cl, list(cl$label, cl$time), drop = TRUE), function(z)
  data.frame(label = z$label[1], time = z$time[1], n = nrow(z), margin_q10 = qtl(z$margin, .1), margin_median = stats::median(z$margin),
             median_genes = stats::median(z$n_genes), median_umis = stats::median(z$total), median_mt = stats::median(z$mt_frac))))
out(mq[order(mq$label, as.integer(sub("h", "", mq$time))),], "annotation_margins.csv")

mexp <- list()
for (i in seq_len(nrow(s))) {
  x <- read_raw_h5(s$file[i])
  z <- cl[cl$gsm == s$gsm[i] & !is.na(cl$lineage),]
  sub <- x$counts[mrow, z$barcode, drop = FALSE]
  ln <- log1p(1e4 * t(t(as.matrix(sub)) / z$total)); rownames(ln) <- names(mrow)
  for (L in c("mono", "nk")) {
    j <- z$lineage == L
    mexp[[length(mexp) + 1]] <- data.frame(gsm = s$gsm[i], donor = s$donor[i], time = s$time[i], library_rep = s$library_rep[i], lineage = L,
                                           marker = names(mrow), set = rep(names(MARKERS), lengths(MARKERS))[match(names(mrow), unlist(MARKERS))],
                                           mean_lognorm = rowMeans(ln[, j, drop = FALSE]), frac_detected = rowMeans(as.matrix(sub[, j, drop = FALSE]) > 0))
  }
  rm(x); invisible(gc())
}
mexp <- do.call(rbind, mexp)
out(mexp, "marker_expression_by_lineage.csv")
add("NK cells have zero CD3 counts (rule)", all(mexp$frac_detected[mexp$lineage == "nk" & mexp$marker %in% T_GENES] == 0))
add("every CD14 monocyte has CD14 detected (rule)", all(mexp$frac_detected[mexp$lineage == "mono" & mexp$marker == "CD14"] == 1))

out(fl, "files.csv")
utils::write.csv(at, file.path(OUT, "h5_attributes.csv"), row.names = FALSE)
cc <- do.call(rbind, chk)
cc$gate <- !grepl("rule\\)$", cc$check)
out(cc, "feasibility_checks.csv")
utils::write.table(data.frame(item = c("run_id", "category", "labels_file", "labels_sha256", "panel_values_computed"),
                              value = c(RUN_ID, "feasibility and annotation gate (no panel-gene values read beyond identifier mapping)",
                                        LABELS_FILE, unname(tools::sha256sum(LABELS_FILE)), "no")),
                   file.path(OUT, "annotation_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "annotation_sessionInfo.txt"))
print(cc); print(pg)
if (!all(cc$pass)) message("FEASIBILITY GATE FAILED: ", paste(cc$check[!cc$pass], collapse = "; ")) else message("FEASIBILITY GATE PASSED")
