GAP <- 5.55; MIN_CELLS <- 20
SCOPE <- "analysis/cross_study_synthesis/sle_bound"
OUT <- "results/cross_study_synthesis/sle_bound"
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", SCOPE), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
chk <- list()
add <- function(check, pass, value) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))

spec <- utils::read.delim(file.path(SCOPE, "specification.tsv"), stringsAsFactors = FALSE)
add("specification locked", grepl("^LOCKED", spec$value[spec$key == "status"]), spec$value[spec$key == "specified_utc"])
genes <- utils::read.delim("analysis/inputs/core_isgs.tsv", stringsAsFactors = FALSE)$gene
x <- utils::read.csv(gzfile("cache/perez2022/perez2022_cM_NK_panel_counts.csv.gz"), stringsAsFactors = FALSE)
x$is_primary_data <- as.logical(x$is_primary_data)
add("10 panel genes found", all(genes %in% names(x)) && length(genes) == 10, sum(genes %in% names(x)))
cnt <- as.matrix(x[, genes])
add("panel counts integer and non-negative", all(cnt == round(cnt)) && all(cnt >= 0), max(abs(cnt - round(cnt))))
add("raw_sum integer and >= panel sum", all(x$raw_sum == round(x$raw_sum)) && all(x$raw_sum >= rowSums(cnt)), min(x$raw_sum - rowSums(cnt)))
add("all cells primary data", all(x$is_primary_data), mean(x$is_primary_data))
add("assay 10x 3' v2", all(x$assay == "10x 3' v2"), paste(unique(x$assay), collapse = ";"))
add("census lineage matches author label", all((x$census_lineage == "cM") == (x$author_cell_type == "cM")) && all(x$author_cell_type %in% c("cM", "NK")),
    paste(names(table(x$author_cell_type)), table(x$author_cell_type), collapse = ";"))
ac <- utils::read.csv("cache/perez2022/perez2022_author_cell_type_counts.csv", stringsAsFactors = FALSE)
add("all author cM and NK cells fetched", all(ac$n_cells_h5ad[match(c("cM", "NK"), ac$author_cell_type)] == as.vector(table(factor(x$author_cell_type, c("cM", "NK"))))),
    paste(ac$n_cells_h5ad[match(c("cM", "NK"), ac$author_cell_type)], collapse = ";"))
keep <- x$author_cell_type %in% c("cM", "NK") & x$is_primary_data
x <- x[keep,]; cnt <- cnt[keep, , drop = FALSE]

key <- paste(x$donor_id, x$author_cell_type)
pb <- rowsum(cbind(cnt, raw_sum = x$raw_sum, n_cells = 1), key)
lev <- data.frame(donor_id = sub(" [^ ]+$", "", rownames(pb)), lineage = sub("^.* ", "", rownames(pb)), n_cells = pb[, "n_cells"], lib = pb[, "raw_sum"])
lev$level <- rowMeans(log2(1e6 * pb[, genes] / pb[, "raw_sum"] + 1))
w <- reshape(lev, idvar = "donor_id", timevar = "lineage", direction = "wide")
w[is.na(w)] <- 0
ann <- function(f) tapply(as.character(x[[f]]), x$donor_id, function(v) paste(sort(unique(v)), collapse = "|"))
d <- data.frame(donor_id = w$donor_id, n_cells_cM = w$n_cells.cM, n_cells_NK = w$n_cells.NK, lib_cM = w$lib.cM, lib_NK = w$lib.NK,
                mean_log2cpm1_cM = w$level.cM, mean_log2cpm1_NK = w$level.NK, stringsAsFactors = FALSE)
for (f in c("disease", "disease_state", "self_reported_ethnicity", "Processing_Cohort", "sample_uuid")) d[[f]] <- ann(f)[d$donor_id]
d$n_samples <- lengths(strsplit(d$sample_uuid, "|", fixed = TRUE)); d$sample_uuid <- NULL
d$group <- ifelse(d$disease == "systemic lupus erythematosus", "SLE", ifelse(d$disease == "normal", "healthy", NA))
d$included <- d$n_cells_cM >= MIN_CELLS & d$n_cells_NK >= MIN_CELLS & !is.na(d$group)
d$delta_ctrl <- ifelse(d$included, d$mean_log2cpm1_cM - d$mean_log2cpm1_NK, NA)
d$run_id <- RUN_ID
d <- d[order(d$group, d$delta_ctrl),]
add("one disease label per donor", !any(grepl("|", d$disease, fixed = TRUE)), sum(grepl("|", d$disease, fixed = TRUE)))
add("one ethnicity per donor", !any(grepl("|", d$self_reported_ethnicity, fixed = TRUE)), sum(grepl("|", d$self_reported_ethnicity, fixed = TRUE)))
add("donors with cells in >1 Processing_Cohort (amendment 1: cohort indicators)", TRUE, sum(grepl("|", d$Processing_Cohort[d$included], fixed = TRUE)))
for (g in c("SLE", "healthy")) {
  add(paste0("donors total ", g), sum(d$group == g, na.rm = TRUE) > 0, sum(d$group == g, na.rm = TRUE))
  add(paste0("donors included (>= ", MIN_CELLS, " cells both lineages) ", g), sum(d$included & d$group == g) >= 10, sum(d$included & d$group == g))
}
add("min cells per lineage among included donors", min(d$n_cells_cM[d$included], d$n_cells_NK[d$included]) >= MIN_CELLS, min(d$n_cells_cM[d$included], d$n_cells_NK[d$included]))
add("donors with >1 sample pooled", TRUE, sum(d$n_samples[d$included] > 1))

s <- d[d$included,]
s$group <- factor(s$group, c("healthy", "SLE"))
tt <- t.test(delta_ctrl ~ group, data = s)
eff <- data.frame(model = "unadjusted Welch", estimate = unname(tt$estimate[2] - tt$estimate[1]), ci_low = -tt$conf.int[2], ci_high = -tt$conf.int[1],
                  p = tt$p.value, n_healthy = sum(s$group == "healthy"), n_sle = sum(s$group == "SLE"), terms = "", dropped_aliased = "")
et <- table(s$self_reported_ethnicity)
s$ancestry <- ifelse(s$self_reported_ethnicity %in% names(et)[et >= 3], s$self_reported_ethnicity, "Other")
for (k in 1:4) s[[paste0("cohort", k)]] <- as.integer(sapply(strsplit(s$Processing_Cohort, "|", fixed = TRUE), function(v) any(as.numeric(v) == k)))
fit <- lm(delta_ctrl ~ group + ancestry + cohort1 + cohort2 + cohort3 + cohort4, data = s)
al <- names(coef(fit))[is.na(coef(fit))]
ci <- confint(fit)["groupSLE",]
eff <- rbind(eff, data.frame(model = "lm adjusted ancestry + Processing_Cohort membership indicators", estimate = unname(coef(fit)["groupSLE"]), ci_low = ci[1], ci_high = ci[2],
                             p = summary(fit)$coefficients["groupSLE", 4], n_healthy = sum(s$group == "healthy"), n_sle = sum(s$group == "SLE"),
                             terms = paste(names(coef(fit)), collapse = ";"), dropped_aliased = paste(al, collapse = ";")))
add("adjusted SLE term estimable", !"groupSLE" %in% al, paste(al, collapse = ";"))
add("Welch estimate equals group mean difference", abs(eff$estimate[1] - diff(tapply(s$delta_ctrl, s$group, mean))) < 1e-10, eff$estimate[1])

q <- do.call(rbind, lapply(levels(s$group), function(g) { v <- s$delta_ctrl[s$group == g]
  data.frame(group = g, n = length(v), mean = mean(v), sd = sd(v), q05 = quantile(v, 0.05, names = FALSE), q50 = quantile(v, 0.5, names = FALSE), q95 = quantile(v, 0.95, names = FALSE)) }))
q$run_id <- RUN_ID
qh <- q[q$group == "healthy",]; qs <- q[q$group == "SLE",]
ci_max <- max(abs(c(eff$ci_low, eff$ci_high)))
pct <- c(qs$q95 - qh$q05, qs$q05 - qh$q95)
pct_max <- pct[which.max(abs(pct))]
bound <- data.frame(model = c("bound: max |CI limit| over unadjusted and adjusted", "bound: extreme group percentile contrast (SLE q95 - healthy q05 or SLE q05 - healthy q95)"),
                    estimate = c(ci_max, pct_max), ci_low = NA, ci_high = NA, p = NA, n_healthy = qh$n, n_sle = qs$n, terms = "", dropped_aliased = "")
eff <- rbind(eff, bound)
eff$fraction_of_gap <- abs(eff$estimate) / GAP
eff$gap <- GAP
eff$run_id <- RUN_ID
add("bound fractions finite", all(is.finite(eff$fraction_of_gap)), paste(round(eff$fraction_of_gap, 3), collapse = ";"))

cc <- do.call(rbind, chk); cc$run_id <- RUN_ID
utils::write.csv(d, file.path(OUT, "donor_delta_ctrl.csv"), row.names = FALSE)
utils::write.csv(eff, file.path(OUT, "sle_effect.csv"), row.names = FALSE)
utils::write.csv(q, file.path(OUT, "group_quantiles.csv"), row.names = FALSE)
utils::write.csv(cc, file.path(OUT, "checks.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
if (!all(cc$pass)) stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; "))
print(cc[, 1:3]); print(eff[, c("model", "estimate", "ci_low", "ci_high", "p", "fraction_of_gap")]); print(q)
