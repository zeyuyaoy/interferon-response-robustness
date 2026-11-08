PRIOR <- 2; TOL <- 1e-8
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "analysis/cross_study_synthesis/baseline_harmonization.R"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE)
level_p2 <- function(cpm, mean_eff) { k <- PRIOR / mean_eff; log2(cpm + k * 1e6) - log2(1 + 2 * k) }
dir.create("results/cross_study_synthesis", recursive = TRUE, showWarnings = FALSE)
chk <- list()
add <- function(check, pass, value) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = value)

k <- rd("results/primary_contrast/09_expression_levels/core_baseline_stim_induction.csv")
k <- k[k$prior_count == PRIOR,]
kn <- rd("results/primary_contrast/02_primary_estimate/normalization.csv")
me_k <- tapply(kn$eff_lib_size, kn$lineage, mean)
k$level_p2 <- level_p2(k$cpm_ctrl, me_k[k$lineage])
add("GSE96583: prior-2 level formula reproduces edgeR logcpm_ctrl", max(abs(k$level_p2 - k$logcpm_ctrl)) < TOL, max(abs(k$level_p2 - k$logcpm_ctrl)))
kang <- aggregate(list(mean_logcpm_p2 = k$logcpm_ctrl, mean_log2_cpm_plus_1 = log2(k$cpm_ctrl + 1)), k[, c("donor", "lineage")], mean)
dec <- rd("results/primary_contrast/09_expression_levels/core_level_decomposition_by_donor.csv")
dec <- dec[dec$prior_count == PRIOR,]
km <- kang[kang$lineage == "mono",]
add("GSE96583: mono prior-2 means equal ctrl_mono of the frozen decomposition",
    max(abs(km$mean_logcpm_p2 - dec$ctrl_mono[match(km$donor, dec$donor)])) < TOL, max(abs(km$mean_logcpm_p2 - dec$ctrl_mono[match(km$donor, dec$donor)])))

h <- rd("results/independent_cohort/h2_baseline.csv")
h <- h[h$cond == "ctrl",]
hl <- rd("results/independent_cohort/h2_library.csv")
h$level_p2 <- level_p2(h$cpm, hl$mean_eff_lib[match(h$lineage, hl$lineage)])
hir <- aggregate(list(mean_logcpm_p2 = h$level_p2, mean_log2_cpm_plus_1 = log2(h$cpm + 1)), h[, c("donor", "lineage")], mean)
add("HIRISA: 10 panel genes per donor x lineage", all(table(h$donor, h$lineage) == 10), 10)

g <- rd("results/time_course/g3_baseline_summary.csv")
g <- g[, c("donor", "lineage", "mean_logcpm_p2", "mean_log2_cpm_plus_1")]

out <- rbind(data.frame(cohort = "GSE96583", kang), data.frame(cohort = "HIRISA GSE306664", hir), data.frame(cohort = "GSE226572 (0h)", g))
out$donor <- as.character(out$donor)
out <- out[order(out$cohort, out$lineage, out$mean_logcpm_p2),]
out$run_id <- RUN_ID
cc <- do.call(rbind, chk); cc$run_id <- RUN_ID
utils::write.csv(out, "results/cross_study_synthesis/baseline_harmonization.csv", row.names = FALSE)
utils::write.csv(cc, "results/cross_study_synthesis/baseline_harmonization_checks.csv", row.names = FALSE)
if (!all(cc$pass)) stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; "))
print(cc[, 1:3]); print(out[out$lineage == "mono", 1:5])
