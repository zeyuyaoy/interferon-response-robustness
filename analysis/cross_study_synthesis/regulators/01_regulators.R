for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/independent_cohort/scripts/utils_hirisa.R")) source(f)
D <- "analysis/cross_study_synthesis/regulators"; OUT <- "results/cross_study_synthesis/regulators"
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", D), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PRIOR <- 2; MIN_TARGETS <- 5; TOL <- 1e-8; GAP <- 2
ISGF3 <- c("STAT1", "STAT2", "IRF9"); NFKB <- c("RELA", "NFKB1", "REL")
FOCUS <- c(ISGF3, "IRF1", "IRF7", "IRF3", NFKB, "NFKB")
COH <- c(k = "GSE96583", h = "HIRISA GSE306664", g = "GSE226572 (0h)")
chk <- list()
add <- function(check, pass, value = NA) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))
src <- utils::read.delim(file.path(D, "network_source.tsv"), stringsAsFactors = FALSE)
spec <- utils::read.delim(file.path(D, "specification.tsv"), stringsAsFactors = FALSE)
add("specification status is LOCKED", startsWith(spec$value[spec$key == "status"], "LOCKED"), spec$value[spec$key == "specified_utc"])
net_file <- src$value[src$key == "file"]
add("network md5 equals network_source.tsv", unname(tools::md5sum(net_file)) == src$value[src$key == "md5"], unname(tools::md5sum(net_file)))
core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene

ct <- utils::read.delim(net_file, stringsAsFactors = FALSE)
lg <- function(x) x %in% c(TRUE, "True")
ct$tf <- ct$source_genesymbol
cx <- startsWith(ct$source, "COMPLEX:")
ct$tf[cx & grepl("^(JUN|FOS)", ct$tf)] <- "AP1"
ct$tf[cx & grepl("^(REL|NFKB)", ct$tf)] <- "NFKB"
ct$w <- ifelse(lg(ct$consensus_inhibition) & !lg(ct$consensus_stimulation), -1, 1)
net <- stats::aggregate(w ~ tf + target_genesymbol, ct, sum)
net$w <- ifelse(net$w < 0, -1, 1)
names(net)[2] <- "target"
add("all non-complex sources are single gene symbols", !any(grepl("_", ct$tf[!cx])))
net <- net[!net$target %in% core,]
add("panel genes absent from all regulons used", !any(core %in% net$target), sum(core %in% net$target))

by_symbol <- function(counts, sym) {
  ok <- !is.na(sym) & nzchar(sym) & !sym %in% sym[duplicated(sym)]
  m <- counts[ok, , drop = FALSE]; rownames(m) <- sym[ok]; m
}

kang <- function() {
  sce <- load_kang("cache/Kang18_8vs8.Rda")
  cd <- as.data.frame(SummarizedExperiment::colData(sce)); cnt <- SummarizedExperiment::assay(sce, "counts")
  lin <- c("CD14+ Monocytes" = MONO, "NK cells" = NK)[as.character(cd$cell)]
  k <- which(cd$multiplets == "singlet" & !is.na(cd$ind) & cd$stim == "ctrl" & !is.na(lin) &
               Matrix::colSums(cnt) > 0 & Matrix::colSums(cnt > 0) >= 200)
  list(counts = by_symbol(cnt[, k], as.character(SummarizedExperiment::rowData(sce)$SYMBOL)),
       meta = data.frame(donor = as.character(cd$ind[k]), lineage = unname(lin[k]), cond = "ctrl"))
}

hirisa <- function() {
  s <- merge(utils::read.delim("analysis/independent_cohort/samples.tsv"), utils::read.delim("analysis/independent_cohort/downloads.tsv")[, c("gsm", "file", "md5")], by = "gsm")
  s <- s[s$cond == "ctrl",]
  add("HIRISA control files md5 equal downloads.tsv", nrow(s) == 10 && all(unname(tools::md5sum(s$file)) == s$md5), nrow(s))
  x <- lapply(seq_len(nrow(s)), function(i) {
    r <- read_hirisa_h5(s$file[i]); m <- r$meta
    k <- (if (s$lineage[i] == MONO) m$celltype.l2 == "CD14 Mono" else m$celltype.l1 == "NK") & m$n_genes >= 200
    list(counts = by_symbol(r$counts[, k, drop = FALSE], r$features$name), meta = data.frame(donor = s$donor[i], lineage = s$lineage[i], cond = "ctrl", n = sum(k)))
  })
  add("HIRISA feature symbols identical across control files", length(unique(lapply(x, function(z) rownames(z$counts)))) == 1)
  list(counts = do.call(cbind, lapply(x, `[[`, "counts")),
       meta = do.call(rbind, lapply(x, function(z) z$meta[rep(1, z$meta$n), c("donor", "lineage", "cond")])))
}

gse226572 <- function() {
  sel <- readRDS("cache/gse226572/selected_cells.rds"); cm <- sel$meta
  lin <- c("CD14 Mono" = MONO, "NK" = NK)[cm$label]
  k <- which(cm$time == "0h" & !is.na(lin))
  list(counts = by_symbol(sel$counts[, k], sel$features$name), meta = data.frame(donor = cm$donor[k], lineage = unname(lin[k]), cond = "ctrl"))
}

ulm <- function(S, W) {
  Sc <- sweep(S, 2, colMeans(S)); Wc <- sweep(W, 2, colMeans(W)); n <- nrow(S)
  sxx <- colSums(Wc^2); b <- crossprod(Wc, Sc) / sxx
  rss <- matrix(colSums(Sc^2), ncol(W), ncol(S), byrow = TRUE) - b^2 * sxx
  structure(b / sqrt(rss / (n - 2) / sxx), slope = b)
}

regulon <- function(g) {
  n_t <- table(net$tf[net$target %in% g])
  tfs <- sort(names(n_t)[n_t >= MIN_TARGETS])
  e <- net[net$tf %in% tfs & net$target %in% g,]
  W <- matrix(0, length(g), length(tfs), dimnames = list(g, tfs))
  W[cbind(match(e$target, g), match(e$tf, tfs))] <- e$w
  W
}

expected <- list(k = c("101", "1015", "1016", "1039", "107", "1244", "1256", "1488"),
                 h = c("2616BW", "3283BW", "3491BW", "3955BW", "6811BW"), g = c("D34", "D38", "D39"))
res <- list(); stats_all <- list(); SL <- list(); PS <- list()
for (co in names(COH)) {
  x <- switch(co, k = kang(), h = hirisa(), g = gse226572())
  add(paste0(COH[co], ": counts integer"), counts_valid(x$counts))
  pb <- pseudobulk(x$counts, x$meta)
  add(paste0(COH[co], ": expected donors present with both lineages"),
      identical(complete_donors(group_table(pb$meta, expected[[co]], c(MONO, NK), "ctrl")), sort(expected[[co]])) &&
        setequal(pb$meta$donor, expected[[co]]), paste(sort(unique(pb$meta$donor)), collapse = ";"))
  pb$meta$cond <- pb$meta$lineage; pb$meta$lineage <- "ctrl"
  nm <- logcpm_by_lineage(pb, PRIOR, "TMM")
  g <- nm$keep$ctrl
  pm <- pb$meta[order(pb$meta$donor),]
  S <- nm$Y[g, pm$profile[pm$cond == MONO]] - nm$Y[g, pm$profile[pm$cond == NK]]
  colnames(S) <- pm$donor[pm$cond == MONO]
  W <- regulon(g); tfs <- colnames(W); n_t <- colSums(W != 0)
  A <- ulm(S, W)
  SL[[co]] <- S
  PS[[co]] <- data.frame(universe = "primary (cohort-specific)", cohort = COH[co], tf = rep(tfs, ncol(A)), donor = rep(colnames(A), each = length(tfs)),
                         n_genes = length(g), n_targets = as.integer(n_t[tfs]), t = as.vector(A), slope = as.vector(attr(A, "slope")))
  if (co == "k") {
    f <- stats::lm(S[, 1] ~ W[, "STAT1"])
    d <- abs(summary(f)$coefficients[2, "t value"] - A["STAT1", 1])
    add("ULM t equals direct lm() fit (GSE96583, STAT1, first donor)", d < TOL, d)
    d <- abs(stats::coef(f)[2] - attr(A, "slope")["STAT1", 1])
    add("ULM slope equals direct lm() coefficient (GSE96583, STAT1, first donor)", d < TOL, d)
  }
  res[[co]] <- data.frame(cohort = COH[co], tf = rep(tfs, ncol(A)), donor = rep(colnames(A), each = length(tfs)),
                          n_targets = as.integer(n_t[tfs]), t = as.vector(A))
  stats_all[[co]] <- data.frame(cohort = COH[co], n_genes_expressed = length(g), n_tfs_scored = length(tfs),
                                n_donors = ncol(S), panel_in_expressed = sum(core %in% g))
  rm(x, pb); invisible(gc())
}
by_donor <- do.call(rbind, res); rownames(by_donor) <- NULL
summ <- stats::aggregate(t ~ cohort + tf + n_targets, by_donor, function(v) c(mean = mean(v), sd = stats::sd(v), n = length(v)))
summ <- data.frame(summ[, 1:3], mean_t = summ$t[, "mean"], sd_t = summ$t[, "sd"], n = summ$t[, "n"])
summ <- summ[order(summ$cohort, summ$tf),]
wide <- function(cn) { z <- summ[summ$cohort == cn,]; stats::setNames(z$mean_t, z$tf) }
both <- intersect(names(wide(COH["k"])), names(wide(COH["h"])))
diff <- wide(COH["h"])[both] - wide(COH["k"])[both]
top <- names(sort(abs(diff), decreasing = TRUE))[1:10]
pick <- function(tf, cn, v) { z <- summ[summ$cohort == cn & summ$tf == tf, v]; if (length(z)) z else NA }
focus <- do.call(rbind, lapply(unique(c(FOCUS, top)), function(tf) data.frame(
  tf = tf, set = if (tf %in% FOCUS) "pre-stated" else "post hoc ranking (top 10 |HIRISA - GSE96583|)",
  in_top10_post_hoc = tf %in% top,
  mean_t_GSE96583 = pick(tf, COH["k"], "mean_t"), sd_t_GSE96583 = pick(tf, COH["k"], "sd_t"), n_GSE96583 = pick(tf, COH["k"], "n"),
  mean_t_HIRISA = pick(tf, COH["h"], "mean_t"), sd_t_HIRISA = pick(tf, COH["h"], "sd_t"), n_HIRISA = pick(tf, COH["h"], "n"),
  mean_t_GSE226572_0h = pick(tf, COH["g"], "mean_t"), sd_t_GSE226572_0h = pick(tf, COH["g"], "sd_t"), n_GSE226572_0h = pick(tf, COH["g"], "n"),
  diff_HIRISA_minus_GSE96583 = if (tf %in% both) unname(diff[tf]) else NA)))
up <- function(s) all(s %in% both) && all(diff[s] > GAP)
none <- function(s) !any(diff[intersect(s, both)] > GAP)
verdict <- if (up(ISGF3) && up(NFKB)) "mixed" else if (up(ISGF3) && none(NFKB)) "ISGF3 supported" else
  if (up(NFKB) && none(ISGF3)) "NF-kB supported" else "not resolved"
add("all pre-stated ISGF3 and NF-kB TFs scored in both GSE96583 and HIRISA", all(c(ISGF3, NFKB) %in% both), paste(setdiff(c(ISGF3, NFKB), both), collapse = ";"))
gs <- sort(intersect(rownames(SL$k), rownames(SL$h)))
Ws <- regulon(gs)
add("shared universe: identical gene set scored in both cohorts", identical(rownames(SL$k[gs,]), rownames(SL$h[gs,])) && identical(rownames(Ws), gs), length(gs))
add("shared universe: 0 panel genes in regulons used", !any(core %in% rownames(Ws)[rowSums(Ws != 0) > 0]), sum(core %in% rownames(Ws)[rowSums(Ws != 0) > 0]))
sh <- do.call(rbind, lapply(c("k", "h"), function(co) {
  A <- ulm(SL[[co]][gs,], Ws)
  data.frame(universe = "shared (expressed in GSE96583 and HIRISA)", cohort = COH[co], tf = rep(colnames(Ws), ncol(A)), donor = rep(colnames(A), each = ncol(Ws)),
             n_genes = length(gs), n_targets = as.integer(colSums(Ws != 0)), t = as.vector(A), slope = as.vector(attr(A, "slope")))
}))
add("shared universe: all pre-stated ISGF3 and NF-kB TFs scored", all(c(ISGF3, NFKB) %in% colnames(Ws)), paste(setdiff(c(ISGF3, NFKB), colnames(Ws)), collapse = ";"))
sd_long <- rbind(do.call(rbind, PS[c("k", "h")]), sh)
sd_long <- sd_long[sd_long$tf %in% FOCUS,]
ag <- function(v) stats::aggregate(sd_long[[v]], sd_long[c("universe", "cohort", "tf")], function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x)))
mt <- ag("t"); ms <- ag("slope"); nt <- stats::aggregate(n_targets ~ universe + cohort + tf + n_genes, sd_long, `[`, 1)
sens <- merge(merge(data.frame(mt[1:3], mean_t = mt$x[, "mean"], sd_t = mt$x[, "sd"], n = mt$x[, "n"]),
                    data.frame(ms[1:3], mean_slope = ms$x[, "mean"], sd_slope = ms$x[, "sd"])), nt)
dk <- sens[sens$cohort == COH["k"],]; dh <- sens[sens$cohort == COH["h"],]
dd <- merge(dk[c("universe", "tf", "mean_t", "mean_slope")], dh[c("universe", "tf", "mean_t", "mean_slope")], by = c("universe", "tf"), suffixes = c("_k", "_h"))
dd$diff_t_HIRISA_minus_GSE96583 <- dd$mean_t_h - dd$mean_t_k; dd$diff_slope_HIRISA_minus_GSE96583 <- dd$mean_slope_h - dd$mean_slope_k
sens <- merge(sens, dd[c("universe", "tf", "diff_t_HIRISA_minus_GSE96583", "diff_slope_HIRISA_minus_GSE96583")], all.x = TRUE)
sens <- sens[order(sens$universe, match(sens$tf, FOCUS), sens$cohort),]
ds <- stats::setNames(dd$diff_t_HIRISA_minus_GSE96583, dd$tf)[dd$universe == unique(sh$universe)]
up_s <- function(s) all(s %in% names(ds)) && all(ds[s] > GAP)
none_s <- function(s) !any(ds[intersect(s, names(ds))] > GAP)
verdict_s <- if (up_s(ISGF3) && up_s(NFKB)) "mixed" else if (up_s(ISGF3) && none_s(NFKB)) "ISGF3 supported" else
  if (up_s(NFKB) && none_s(ISGF3)) "NF-kB supported" else "not resolved"
cc <- rbind(do.call(rbind, chk), data.frame(check = "interpretation rule verdict (descriptive)", pass = TRUE, value = verdict),
            data.frame(check = "interpretation rule verdict, shared-universe sensitivity (descriptive)", pass = TRUE, value = verdict_s),
            data.frame(check = "shared universe: genes / TFs scored", pass = TRUE, value = paste(length(gs), ncol(Ws), sep = " / ")))
st <- do.call(rbind, stats_all); st$n_tfs_both_k_h <- length(both)
wr <- function(df, f) { df$run_id <- RUN_ID; utils::write.csv(df, file.path(OUT, f), row.names = FALSE) }
wr(by_donor, "tf_activity_by_donor.csv"); wr(summ, "tf_activity_summary.csv"); wr(focus, "focus_tfs.csv"); wr(sens, "sensitivity_shared_universe.csv")
wr(rbind(cc, data.frame(check = paste0(st$cohort, ": genes / TFs scored / donors"), pass = TRUE,
                        value = paste(st$n_genes_expressed, st$n_tfs_scored, st$n_donors, sep = " / "))), "checks.csv")
writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(st); print(cc[, 1:3]); print(focus[, c(1, 4, 7, 10, 13)], digits = 3)
if (!all(cc$pass)) stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; "))
