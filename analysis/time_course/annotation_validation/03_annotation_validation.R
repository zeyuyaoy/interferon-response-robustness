for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R",
            "analysis/time_course/scripts/utils_gse226572.R")) source(f)
D <- "analysis/time_course"; A <- file.path(D, "annotation_validation"); OUT <- "results/time_course/annotation_validation"
SPEC <- file.path(A, "specification.tsv"); HIST <- file.path(A, "spec_history.tsv")
if (!spec_gate(SPEC, HIST)) stop("specification.tsv differs from the last spec_history.tsv entry; analysis refused", call. = FALSE)
STARTED <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", A, "analysis/time_course/scripts", "analysis/shared", "analysis/inputs"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(OUT, c("checks.csv", "cell_counts_by_annotation.csv", "gamma_by_annotation.csv", "gamma_by_donor_annotation.csv", "run_manifest.tsv", "sessionInfo.txt")))
spec <- utils::read.delim(SPEC, stringsAsFactors = FALSE)
AL_FILE <- "cache/gse226572/author_labels.tsv.gz"; SEL_FILE <- "cache/gse226572/selected_cells.rds"
AL_SHA <- regmatches(spec$value[spec$key == "author_labels"], regexpr("[0-9a-f]{64}", spec$value[spec$key == "author_labels"]))
tspec <- utils::read.delim(file.path(D, "specification.tsv"), stringsAsFactors = FALSE)
MARKER_SHA <- regmatches(tspec$value[tspec$key == "cell_labels"], regexpr("[0-9a-f]{64}", tspec$value[tspec$key == "cell_labels"]))
TOL <- 1e-8; PRIMARY <- c("0h", "4h", "24h"); EARLY <- "4h"; LATE <- "24h"
AUTHOR_MONO <- c("4. Classical Monocytes", "16. Classical monocytes 2"); AUTHOR_NC <- "18. Non-classical Monocytes"
CD14 <- c("CD14 Mono" = "mono", "NK" = "nk"); MONO_ALL <- c("CD14 Mono" = "mono", "CD16 Mono" = "mono", "NK" = "nk")
TR <- c(p2 = "edgeR logCPM, prior count 2", p05 = "edgeR logCPM, prior count 0.5", fixed = "log2(CPM + 1)")

chk <- list()
add <- function(check, pass, value = NA) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))
out <- function(df, name) write_result(df, name, RUN_ID, dir = OUT)
finish <- function() {
  cc <- do.call(rbind, chk); out(cc, "checks.csv")
  writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  if (!all(cc$pass)) stop("annotation validation checks failed: ", paste(cc$check[!cc$pass], collapse = "; "), call. = FALSE)
}

h <- utils::read.delim(HIST, header = FALSE, colClasses = "character")
add("specification MD5 equals the last spec_history entry", TRUE, h[nrow(h), 2])
al_ok <- file.exists(AL_FILE) && identical(unname(tools::sha256sum(AL_FILE)), AL_SHA)
add("author_labels sha256 equals the specification", al_ok, AL_SHA)
sel <- readRDS(SEL_FILE)
add("marker-rule selected cells built from the locked cell_labels", identical(sel$labels_sha256, MARKER_SHA), MARKER_SHA)
if (!al_ok) finish()
s <- merge(utils::read.delim(file.path(D, "samples.tsv"), stringsAsFactors = FALSE),
           utils::read.delim(file.path(D, "downloads.tsv"), colClasses = "character")[, c("gsm", "file", "sha256")], by = "gsm")
s <- s[s$time %in% PRIMARY,]
add("primary-time H5 files have sha256 equal to downloads.tsv", nrow(s) == 12 && all(unname(tools::sha256sum(s$file)) == s$sha256), nrow(s))

al <- utils::read.delim(gzfile(AL_FILE), stringsAsFactors = FALSE)
al$lin_primary <- ifelse(al$author_celltype %in% AUTHOR_MONO, MONO, ifelse(al$author_population == "NK", NK, NA))
al$lin_all <- ifelse(al$author_celltype %in% c(AUTHOR_MONO, AUTHOR_NC), MONO, al$lin_primary)
al <- al[!is.na(al$lin_all) & al$time %in% PRIMARY,]
mats <- list(); feats <- NULL; same <- logical(); miss <- integer()
for (i in seq_len(nrow(s))) {
  x <- read_raw_h5(s$file[i])
  if (is.null(feats)) feats <- x$features
  same[i] <- identical(x$features$id, feats$id)
  l <- al[al$gsm == s$gsm[i],]
  j <- match(l$barcode, colnames(x$counts)); miss[i] <- sum(is.na(j))
  mats[[i]] <- x$counts[, j[!is.na(j)], drop = FALSE]
  colnames(mats[[i]]) <- paste(s$gsm[i], l$barcode[!is.na(j)], sep = ":")
  rm(x); invisible(gc())
}
XA <- do.call(cbind, mats)
am <- al[match(colnames(XA), paste(al$gsm, al$barcode, sep = ":")),]
add("feature ids identical across files", all(same) && identical(feats$id, sel$features$id))
add("every author-labelled barcode found in its H5", sum(miss) == 0, sum(miss))
add("integer non-negative counts", counts_valid(XA))

core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
cand <- utils::read.delim("analysis/inputs/annotation/candidate_rows.tsv", stringsAsFactors = FALSE)
cand <- cand[cand$set == "core",]
mp <- map_genes(core, feats$name, feats$id)
add("10 panel genes uniquely mapped with Ensembl ids equal to the frozen snapshot",
    all(mp$status == "mapped") && all(mp$row == cand$ENSEMBL[match(mp$gene, cand$symbol)]))
genes <- stats::setNames(mp$row, mp$gene)

estimate <- function(X, donor, lineage, time) {
  k <- which(!is.na(lineage) & time %in% PRIMARY)
  pb <- pseudobulk(X[, k, drop = FALSE], data.frame(donor = donor[k], lineage = lineage[k], cond = time[k]))
  n2 <- logcpm_by_lineage(pb, 2); n05 <- logcpm_by_lineage(pb, 0.5)
  eff <- stats::setNames(n2$norm$eff_lib_size, n2$norm$profile)[colnames(pb$counts)]
  Ys <- list(p2 = n2$Y, p05 = n05$Y, fixed = logcpm_fixed(pb$counts, eff, 1))
  list(pb = pb, tc = lapply(Ys, function(Y) time_contrasts(Y, pb$meta, genes, c(EARLY, LATE))))
}
cmk <- sel$meta
ANN <- list(
  marker_cd14 = list(label = "marker rule: CD14 Mono vs NK (frozen primary)", e = estimate(sel$counts, cmk$donor, unname(CD14[cmk$label]), cmk$time)),
  marker_cd14_or_cd16 = list(label = "marker rule: CD14 or CD16 Mono vs NK", e = estimate(sel$counts, cmk$donor, unname(MONO_ALL[cmk$label]), cmk$time)),
  author_classical = list(label = "depositors: classical monocytes vs NK (primary of this analysis)", e = estimate(XA, am$donor, am$lin_primary, am$time)),
  author_all_monocytes = list(label = "depositors: classical or non-classical monocytes vs NK", e = estimate(XA, am$donor, am$lin_all, am$time)))
for (nm in names(ANN)) add(sprintf("18 primary profiles (%s)", nm), ncol(ANN[[nm]]$e$pb$counts) == 18, ncol(ANN[[nm]]$e$pb$counts))

donor_rows <- do.call(rbind, lapply(names(ANN), function(nm) do.call(rbind, lapply(names(TR), function(tr) {
  g <- gamma_table(ANN[[nm]]$e$tc[[tr]], EARLY, LATE)
  data.frame(annotation = nm, transformation = TR[[tr]], g)
}))))
summ <- do.call(rbind, lapply(names(ANN), function(nm) do.call(rbind, lapply(names(TR), function(tr) {
  tc <- ANN[[nm]]$e$tc[[tr]]; g <- gamma_table(tc, EARLY, LATE)
  ci <- ci_t(g$gamma)
  data.frame(annotation = nm, description = ANN[[nm]]$label, transformation = TR[[tr]],
             theta_4h = mean(tc$delta[tc$time == EARLY]), theta_24h = mean(tc$delta[tc$time == LATE]), Gamma = mean(g$gamma),
             n_donors = nrow(g), n_donor_gamma_pos = sum(g$gamma > 0), gamma_min = min(g$gamma), gamma_max = max(g$gamma),
             t_ci95_lo_df2 = ci[1], t_ci95_hi_df2 = ci[2], signflip_p_one_sided_increase = signflip_p(g$gamma),
             R_mono_4h = mean(tc$R_mono[tc$time == EARLY]), R_mono_24h = mean(tc$R_mono[tc$time == LATE]),
             R_nk_4h = mean(tc$R_nk[tc$time == EARLY]), R_nk_24h = mean(tc$R_nk[tc$time == LATE]))
}))))
counts <- do.call(rbind, lapply(names(ANN), function(nm) data.frame(annotation = nm, ANN[[nm]]$e$pb$meta[, c("donor", "lineage", "cond", "n_cells")])))
names(counts)[names(counts) == "cond"] <- "time"

g3 <- utils::read.csv("results/time_course/g3_summary.csv", stringsAsFactors = FALSE)
g3s <- utils::read.csv("results/time_course/g3_sensitivity.csv", stringsAsFactors = FALSE)
p2 <- summ$transformation == TR[["p2"]]
add("marker_cd14 reproduces frozen Gamma", abs(summ$Gamma[p2 & summ$annotation == "marker_cd14"] - as.numeric(g3$value[g3$quantity == "Gamma"])) < TOL,
    summ$Gamma[p2 & summ$annotation == "marker_cd14"])
add("marker_cd14_or_cd16 reproduces frozen sensitivity Gamma",
    abs(summ$Gamma[p2 & summ$annotation == "marker_cd14_or_cd16"] - g3s$Gamma[g3s$setting == "mono_cd14_or_cd16"]) < TOL,
    summ$Gamma[p2 & summ$annotation == "marker_cd14_or_cd16"])
add("marker_cd14 log2(CPM+1) and prior 0.5 reproduce frozen sensitivities",
    abs(summ$Gamma[summ$annotation == "marker_cd14" & summ$transformation == TR[["p05"]]] - g3s$Gamma[g3s$setting == "prior_0.5"]) < TOL &&
      abs(summ$Gamma[summ$annotation == "marker_cd14" & summ$transformation == TR[["fixed"]]] - g3s$Gamma[g3s$setting == "log2_cpm_plus_1"]) < TOL)
add("Gamma equals theta(24h) - theta(4h) in every row", max(abs(summ$Gamma - (summ$theta_24h - summ$theta_4h))) < TOL)
pr <- summ[p2 & summ$annotation == "author_classical",]
decision <- if (pr$Gamma > 0 && pr$n_donor_gamma_pos == 3) "temporal result preserved under the depositors' annotation" else "annotation-sensitive"
summ$decision <- ifelse(p2 & summ$annotation == "author_classical", decision, "")

out(counts, "cell_counts_by_annotation.csv"); out(summ, "gamma_by_annotation.csv"); out(donor_rows, "gamma_by_donor_annotation.csv")
finish()
inputs <- c(SPEC, HIST, AL_FILE, file.path(D, "specification.tsv"), "analysis/inputs/core_isgs.tsv", "analysis/inputs/annotation/candidate_rows.tsv",
            file.path(D, "samples.tsv"), file.path(D, "downloads.tsv"), "results/time_course/g3_summary.csv", "results/time_course/g3_sensitivity.csv")
man <- rbind(data.frame(item = c("run_id", "started_utc", "finished_utc", "git_sha", "spec_md5", "R"),
                        value = c(RUN_ID, STARTED, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), sha[1], h[nrow(h), 2],
                                  R.version.string)),
             data.frame(item = paste0("sha256:", inputs), value = unname(tools::sha256sum(inputs))))
utils::write.table(man, file.path(OUT, "run_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ[, c("annotation", "transformation", "theta_4h", "theta_24h", "Gamma", "n_donor_gamma_pos", "gamma_min", "gamma_max")], digits = 3)
print(donor_rows[donor_rows$transformation == TR[["p2"]],], digits = 3); cat("decision:", decision, "\n")
