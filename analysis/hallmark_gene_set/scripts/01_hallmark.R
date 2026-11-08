for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R",
            "analysis/independent_cohort/scripts/utils_hirisa.R", "analysis/time_course/scripts/utils_gse226572.R")) source(f)
D <- "analysis/hallmark_gene_set"; OUT <- "results/hallmark_gene_set"
SPEC <- file.path(D, "specification.tsv"); HIST <- file.path(D, "spec_history.tsv")
if (!spec_gate(SPEC, HIST)) stop("specification.tsv differs from the last spec_history.tsv entry; analysis refused", call. = FALSE)
STARTED <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", D, "analysis/shared", "analysis/inputs", "analysis/independent_cohort/scripts/utils_hirisa.R", "analysis/time_course/scripts"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); unlink(Sys.glob(file.path(OUT, "*")))
TOL <- 1e-8; PRIOR <- 2; OFFSET <- 1
spec <- utils::read.delim(SPEC, stringsAsFactors = FALSE)
ELIG_SHA <- regmatches(spec$value[spec$key == "eligible_genes"], regexpr("[0-9a-f]{64}", spec$value[spec$key == "eligible_genes"]))
chk <- list()
add <- function(check, pass, value = NA) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))
out <- function(df, name) write_result(df, paste0(name, ".csv"), RUN_ID, dir = OUT)

finish <- function() {
  cc <- do.call(rbind, chk); out(cc, "checks")
  writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  if (!all(cc$pass)) stop(sprintf("checks failed: %s", paste(cc$check[!cc$pass], collapse = "; ")), call. = FALSE)
}

EF <- file.path(D, "eligible_genes.tsv")
add("eligible_genes.tsv sha256 equals the specification", identical(unname(tools::sha256sum(EF)), ELIG_SHA), ELIG_SHA)
if (!all(do.call(rbind, chk)$pass)) finish()
el <- utils::read.delim(EF, stringsAsFactors = FALSE)
el <- el[el$eligible,]
add("93 eligible genes, 88 without panel genes", nrow(el) == 93 && sum(!el$in_core) == 88, nrow(el))
core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene

fixed_Y <- function(pb, nm) logcpm_fixed(pb$counts, stats::setNames(nm$norm$eff_lib_size, nm$norm$profile)[colnames(pb$counts)], OFFSET)

score <- function(pb, rows, ctrl, stim) {
  nm <- logcpm_by_lineage(pb, PRIOR)
  lapply(list(p2 = nm$Y, fixed = fixed_Y(pb, nm)), function(Y) {
    if (length(stim) == 1) return(donor_contrasts(induction(Y, pb$meta, rows, ctrl, stim), MONO, NK))
    gamma_table(time_contrasts(Y, pb$meta, rows, stim, ctrl), stim[1], stim[2])
  })
}

sce <- load_kang("cache/Kang18_8vs8.Rda")
cd <- as.data.frame(SummarizedExperiment::colData(sce))
cnt <- SingleCellExperiment::counts(sce)
lin <- c("CD14+ Monocytes" = MONO, "NK cells" = NK)[as.character(cd$cell)]
keep <- as.character(cd$multiplets) %in% "singlet" &
  !is.na(cd$ind) &
  as.character(cd$stim) %in% c("ctrl", "stim") &
  !is.na(cd$cell) &
  !is.na(lin) &
  Matrix::colSums(cnt) > 0 &
  Matrix::colSums(cnt > 0) >= 200
donors_k <- strsplit(utils::read.csv("results/primary_contrast/02_primary_estimate/analysis_set.csv")$donors, ";")[[1]]
keep <- keep & as.character(cd$ind) %in% donors_k
PBk <- pseudobulk(cnt[, keep], data.frame(donor = as.character(cd$ind[keep]), lineage = unname(lin[keep]), cond = as.character(cd$stim[keep])))
rm(sce, cnt, e); invisible(gc())
sym <- rownames(PBk$counts)
mk <- map_genes(core, sym, sym)
rows_k <- list(panel_10 = stats::setNames(mk$row, mk$gene), hallmark_93 = stats::setNames(el$gse96583_row, el$symbol),
               hallmark_minus_panel_88 = stats::setNames(el$gse96583_row[!el$in_core], el$symbol[!el$in_core]))
add("GSE96583: 32 profiles", ncol(PBk$counts) == 32, ncol(PBk$counts))
add("GSE96583: all scored rows present", all(unlist(rows_k) %in% sym))
ref_k <- utils::read.csv("results/primary_contrast/11_hallmark_gene_sets/hallmark_resolution_symbols.csv")
ref_k <- ref_k$row[ref_k$build == "grch37" & ref_k$status == "resolved"]
ref_h <- utils::read.csv("results/primary_contrast/11_hallmark_gene_sets/hallmark_repaired_primary.csv")
tk96 <- mean(score(PBk, ref_k, "ctrl", "stim")$p2$delta)
add("GSE96583: 96-gene Hallmark reproduces frozen theta", abs(tk96 - ref_h$theta[ref_h$build == "grch37"]) < TOL, tk96)

hs <- merge(utils::read.delim("analysis/independent_cohort/samples.tsv"), utils::read.delim("analysis/independent_cohort/downloads.tsv")[, c("gsm", "file", "md5")], by = "gsm")
add("HIRISA: 20 files with md5 equal to downloads.tsv", nrow(hs) == 20 && all(unname(tools::md5sum(hs$file)) == hs$md5))
hs$profile <- paste(hs$donor, hs$lineage, hs$cond, sep = "|"); hs <- hs[order(hs$profile),]
pbc <- list(); Gh <- NULL
for (i in seq_len(nrow(hs))) {
  x <- read_hirisa_h5(hs$file[i])
  if (is.null(Gh)) Gh <- x$features$id
  stopifnot(identical(x$features$id, Gh))
  k <- (if (hs$lineage[i] == MONO) x$meta$celltype.l2 == "CD14 Mono" else x$meta$celltype.l1 == "NK") & x$meta$n_genes >= 200
  pbc[[i]] <- Matrix::rowSums(x$counts[, k, drop = FALSE])
  rm(x); invisible(gc())
}
PBh <- list(counts = do.call(cbind, pbc), meta = data.frame(profile = hs$profile, donor = hs$donor, lineage = hs$lineage, cond = hs$cond))
dimnames(PBh$counts) <- list(Gh, hs$profile)
cand <- utils::read.delim("analysis/inputs/annotation/candidate_rows.tsv"); cand <- cand[cand$set == "core",]
rows_h <- list(panel_10 = stats::setNames(cand$ENSEMBL[match(core, cand$symbol)], core), hallmark_93 = stats::setNames(el$hgnc_ensembl, el$symbol),
               hallmark_minus_panel_88 = stats::setNames(el$hgnc_ensembl[!el$in_core], el$symbol[!el$in_core]))
add("HIRISA: all scored rows present", all(unlist(rows_h) %in% Gh))

sel <- readRDS("cache/gse226572/selected_cells.rds")
cm <- sel$meta; pk <- which(cm$time %in% c("0h", "4h", "24h") & cm$label %in% c("CD14 Mono", "NK"))
PBg <- pseudobulk(sel$counts[, pk], data.frame(donor = cm$donor[pk], lineage = c("CD14 Mono" = MONO, "NK" = NK)[cm$label[pk]], cond = cm$time[pk]))
rm(sel); invisible(gc())
rows_g <- list(panel_10 = stats::setNames(cand$ENSEMBL[match(core, cand$symbol)], core), hallmark_93 = rows_h$hallmark_93,
               hallmark_minus_panel_88 = rows_h$hallmark_minus_panel_88)
add("GSE226572: all scored rows present", all(unlist(rows_g) %in% rownames(PBg$counts)))

ref <- c(gse96583 = utils::read.csv("results/primary_contrast/02_primary_estimate/primary.csv")$theta,
         hirisa = with(utils::read.csv("results/independent_cohort/h2_theta.csv"), theta[donor_set == "all_donors" & grepl("prior count 2", transformation)]),
         gse226572 = as.numeric(with(utils::read.csv("results/time_course/g3_summary.csv"), value[quantity == "Gamma"])))
res <- list(); donors <- list()
for (gs in names(rows_k)) {
  fk <- score(PBk, rows_k[[gs]], "ctrl", "stim"); fh <- score(PBh, rows_h[[gs]], "ctrl", "stim")
  fg <- score(PBg, rows_g[[gs]], "0h", c("4h", "24h"))
  for (tr in c("p2", "fixed")) {
    res[[length(res) + 1]] <- data.frame(cohort = c("GSE96583", "HIRISA GSE306664", "GSE226572", "GSE226572", "GSE226572"), gene_set = gs,
                                         n_genes = length(rows_k[[gs]]), transformation = c(p2 = "edgeR logCPM, prior count 2", fixed = "log2(CPM + 1)")[[tr]],
                                         estimand = c("theta", "theta", "theta_4h", "theta_24h", "Gamma"),
                                         estimate = c(mean(fk[[tr]]$delta), mean(fh[[tr]]$delta), mean(fg[[tr]]$delta_early), mean(fg[[tr]]$delta_late), mean(fg[[tr]]$gamma)),
                                         n_donors = c(nrow(fk[[tr]]), nrow(fh[[tr]]), rep(nrow(fg[[tr]]), 3)),
                                         n_donors_positive = c(sum(fk[[tr]]$delta > 0), sum(fh[[tr]]$delta > 0), sum(fg[[tr]]$delta_early > 0), sum(fg[[tr]]$delta_late > 0), sum(fg[[tr]]$gamma > 0)),
                                         ci_lower = c(ci_t(fk[[tr]]$delta)[1], ci_t(fh[[tr]]$delta)[1], NA, NA, ci_t(fg[[tr]]$gamma)[1]),
                                         ci_upper = c(ci_t(fk[[tr]]$delta)[2], ci_t(fh[[tr]]$delta)[2], NA, NA, ci_t(fg[[tr]]$gamma)[2]))
    donors[[length(donors) + 1]] <- rbind(
      data.frame(cohort = "GSE96583", gene_set = gs, transformation = tr, donor = fk[[tr]]$donor, quantity = "delta", value = fk[[tr]]$delta),
      data.frame(cohort = "HIRISA GSE306664", gene_set = gs, transformation = tr, donor = fh[[tr]]$donor, quantity = "delta", value = fh[[tr]]$delta),
      data.frame(cohort = "GSE226572", gene_set = gs, transformation = tr, donor = fg[[tr]]$donor, quantity = "gamma", value = fg[[tr]]$gamma))
  }
  if (gs == "panel_10") {
    add("GSE96583: panel_10 reproduces frozen theta", abs(mean(fk$p2$delta) - ref[["gse96583"]]) < TOL, mean(fk$p2$delta))
    add("HIRISA: panel_10 reproduces frozen theta", abs(mean(fh$p2$delta) - ref[["hirisa"]]) < TOL, mean(fh$p2$delta))
    add("GSE226572: panel_10 reproduces frozen Gamma", abs(mean(fg$p2$gamma) - ref[["gse226572"]]) < 1e-6, mean(fg$p2$gamma))
    if (!all(do.call(rbind, chk)$pass)) finish()
  }
}
res <- do.call(rbind, res); donors <- do.call(rbind, donors)
pr <- res[res$transformation == "edgeR logCPM, prior count 2" & res$estimand %in% c("theta", "Gamma"),]
dir_tab <- do.call(rbind, lapply(c("hallmark_93", "hallmark_minus_panel_88"), function(gs) {
  p <- pr[pr$gene_set == gs,]; q <- pr[pr$gene_set == "panel_10",]
  data.frame(gene_set = gs, cohort = p$cohort, estimand = p$estimand, panel_10 = q$estimate, hallmark = p$estimate,
             persists = c(p$estimate[1] > 0, p$estimate[2] < 0, p$estimate[3] > 0), donors_same_direction_as_panel =
               c(p$n_donors_positive[1], p$n_donors[2] - p$n_donors_positive[2], p$n_donors_positive[3]), n_donors = p$n_donors)
}))
d93 <- dir_tab[dir_tab$gene_set == "hallmark_93",]
hp <- d93$persists[d93$cohort == "HIRISA GSE306664"]; gp <- d93$persists[d93$cohort == "GSE226572"]
outcome <- if (hp && gp) "both persist: findings extend beyond the ten-gene signature" else if (hp) "HIRISA persists, GSE226572 changes: the temporal result may be specific to particular ISGs" else if (gp) "HIRISA changes, GSE226572 persists: the HIRISA reversal may depend on the gene program" else "neither persists: restrict conclusions to the original panel"
out(res, "estimates"); out(donors, "donor_values"); out(dir_tab, "direction"); out(data.frame(outcome = outcome), "outcome")
out(data.frame(cohort = c("GSE96583", "HIRISA GSE306664", "GSE226572"), n_profiles = c(ncol(PBk$counts), ncol(PBh$counts), ncol(PBg$counts))), "profiles")
finish()
man <- data.frame(item = c("run_id", "started_utc", "finished_utc", "git_sha", "spec_md5", "eligible_sha256"),
                  value = c(RUN_ID, STARTED, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), sha[1], unname(tools::md5sum(SPEC)), ELIG_SHA))
utils::write.table(man, file.path(OUT, "run_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(do.call(rbind, chk)[, 1:3]); print(res[, c(1, 2, 4, 5, 6, 7, 8)]); print(dir_tab); cat(outcome, "\n")
