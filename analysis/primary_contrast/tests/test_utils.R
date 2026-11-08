source("analysis/primary_contrast/scripts/config.R")
for (f in Sys.glob("analysis/shared/utils_*.R")) source(f)
suppressPackageStartupMessages(library(Matrix))

ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE); cat("ok -", msg, "\n") }

err_class <- function(expr) tryCatch({ expr; "none" }, not_estimable = function(e) "not_estimable",
                                     error = function(e) "error")

set.seed(1)
G <- 300
meta <- expand.grid(cell = 1:5, cond = c("ctrl", "stim"), lineage = c("mono", "nk"),
                    donor = c("d1", "d2", "d3"), stringsAsFactors = FALSE)
meta$cell <- paste0("c", seq_len(nrow(meta)))
counts <- Matrix(rpois(G * nrow(meta), 20), G, nrow(meta), sparse = TRUE,
                 dimnames = list(paste0("g", 1:G), meta$cell))

pb <- pseudobulk(counts, meta)
ok(ncol(pb$counts) == 12 && identical(colnames(pb$counts), pb$meta$profile), "pseudobulk columns align with metadata")
for (p in pb$meta$profile) {
  m <- pb$meta[pb$meta$profile == p,]
  cells <- meta$cell[meta$donor == m$donor &
                       meta$lineage == m$lineage &
                       meta$cond == m$cond]
  stopifnot(all(pb$counts[, p] == rowSums(as.matrix(counts[, cells]))))
}
ok(TRUE, "group sums equal per-group column sums")
ok(sum(pb$counts) == sum(counts) && all(pb$meta$n_cells == 5), "total counts conserved, n_cells correct")

perm <- sample(ncol(counts))
pb2 <- pseudobulk(counts[, perm], meta[perm,])
ok(identical(pb$counts, pb2$counts) && identical(pb$meta, pb2$meta), "pseudobulk invariant to cell order")
f1 <- theta_pipeline(counts, meta, paste0("g", 1:10), 2, "mono", "nk")
f2 <- theta_pipeline(counts[, perm], meta[perm,], paste0("g", 1:10), 2, "mono", "nk")
ok(isTRUE(all.equal(f1$theta, f2$theta)), "theta invariant to cell order")

drop <- !(meta$donor == "d2" &
  meta$lineage == "nk" &
  meta$cond == "stim")
gt <- group_table(meta[drop,], c("d1", "d2", "d3"), c("mono", "nk"), c("ctrl", "stim"))
ok(nrow(gt) == 12 && sum(gt$n_cells == 0) == 1, "group_table reports empty group as 0")
ok(identical(complete_donors(gt), c("d1", "d3")), "complete_donors drops donor with empty group")
f3 <- theta_pipeline(counts[, drop], meta[drop,], paste0("g", 1:10), 2, "mono", "nk")
ok(nrow(f3$delta) == 2 && !"d2" %in% f3$delta$donor, "missing group: donor contrast is dropped, not imputed")

z <- counts; z[, meta$donor == "d1" &
                   meta$lineage == "mono" &
                   meta$cond == "ctrl"] <- 0
ok(err_class(theta_pipeline(z, meta, paste0("g", 1:10), 2, "mono", "nk")) == "not_estimable",
   "zero pseudobulk library -> not_estimable")

gm <- map_genes(c("A", "B", "C"), c("A", "B", "B", NA), c("r1", "r2", "r3", "r4"))
ok(identical(gm$status, c("mapped", "ambiguous", "missing")) && identical(gm$row, c("r1", NA, NA)),
   "map_genes flags mapped/ambiguous/missing")

ok(counts_valid(counts) &&
     !counts_valid(Matrix(c(1, 2.5), 1)) &&
     !counts_valid(Matrix(c(1, -1), 1)),
   "counts_valid rejects non-integer and negative")

ok(signflip_p(1:8) == 1 / 256, "sign-flip all positive = 1/256")
ok(signflip_p(-(1:8)) == 1, "sign-flip all negative = 1")
ok(signflip_p(rep(0, 8)) == 1, "sign-flip all zero = 1")
pm <- signflip_p(c(3, -1, 2, 1, -2, 4, 1, 0.5)); ok(pm > 1 / 256 && pm < 1, "sign-flip mixed in (1/256, 1)")
ok(length(ci_t(1)) == 2 && all(is.na(ci_t(1))), "t-interval NA for n = 1")

d <- c(1, 2, 3, 4, 5, 6, 7, 8)
l <- loo(d, letters[1:8])
ok(isTRUE(all.equal(l$theta_minus[1], mean(d[-1]))) && isTRUE(all.equal(l$influence, l$theta_minus - 4.5)),
   "LODO arithmetic")

grp <- rep(c("a", "b"), c(10, 4))
i <- subsample_idx(grp, 4, 1)
ok(length(i) == 8 &&
     !anyDuplicated(i) &&
     all(table(grp[i]) == 4), "subsample balanced without replacement")
ok(identical(i, subsample_idx(grp, 4, 1)), "subsample reproducible by seed")
ok(err_class(subsample_idx(grp, 5, 1)) == "not_estimable", "subsample target > group -> not_estimable")

base <- cbind(mono = 1:20 / 2, nk = 20:1 / 2); rownames(base) <- paste0("g", 1:20)
core <- c("g1", "g2"); elig <- setdiff(rownames(base), c(core, "g3"))
pools <- ref_pools(base, core, elig, 5)
ok(all(lengths(lapply(pools, `[[`, "gene")) == 5) && !any(c(core, "g3") %in% unlist(lapply(pools, `[[`, "gene"))),
   "reference pools exclude core/excluded genes")
set.seed(2); ps <- replicate(200, sample_panel(pools))
ok(all(apply(ps, 2, function(p) !anyDuplicated(p))), "reference panel: no duplicates within panel")
tiny <- ref_pools(base, core, c("g10"), 50)
ok(nrow(tiny$g1) == 1 && err_class(sample_panel(tiny)) == "not_estimable", "exhausted pool -> not_estimable")

R <- data.frame(donor = rep(letters[1:8], 2), lineage = rep(c("mono", "nk"), each = 8), R = 1)
lz <- loo(rep(0, 8), letters[1:8])
cz <- robustness_criteria(0, lz, rep(0.1, 100), R, "mono", "nk", CFG, TRUE)
ok(all(is.na(cz$criterion_met[cz$criterion %in% c("donor_robustness_lodo", "subsampling_robustness")])),
   "theta == 0 -> directional criteria NA")
sub <- c(rep(1, 95), rep(NA, 5))
c1 <- robustness_criteria(1, loo(rep(1, 8), letters[1:8]), sub, R, "mono", "nk", CFG, TRUE)
ok(c1$criterion_met[c1$criterion == "subsampling_robustness"], "95/100 agreeing -> met")
c2 <- robustness_criteria(1, loo(rep(1, 8), letters[1:8]), c(rep(1, 94), rep(NA, 6)), R, "mono", "nk", CFG, TRUE)
ok(!c2$criterion_met[c2$criterion == "subsampling_robustness"] &&
     grepl("of 100 attempted", c2$threshold[c2$criterion == "subsampling_robustness"]),
   "failed draws count as disagreeing; denominator stays 100")
c3 <- robustness_criteria(1, NULL, NULL, R, "mono", "nk", CFG, FALSE)
ok(all(is.na(c3$criterion_met)), "complete-donor set -> 8-donor criteria NA")

.status$rows <- list()
run_component("a", function() 1)
run_component("b", function() not_estimable("data"))
run_component("c", function() stop("bug"))
run_component("d", function() { component_warning("w"); 1 })
ok(identical(status_table()$status, c("completed", "not_estimable", "failed", "completed_with_warning")),
   "run_component maps statuses")


hg <- data.frame(symbol = c("A", "B", "C", "D", "E"), hgnc_id = paste0("HGNC:", 1:5),
                 hgnc_ensembl = paste0("ENS", 1:5))
cand <- data.frame(symbol = c("A", "B", "B", "C", "C", "D", "D"),
                   row = c("A", "oldB", "X", "C_ENS3", "C_ENS3alt", "D_ENS4", "D_ENS9"),
                   SYMBOL = c("A", "oldB", "X", "C", "C", "D", "D"),
                   ENSEMBL = c("ENS1", "ENS2", "ENS8", "ENS3", "ENS3alt", "ENS4", "ENS9"))
ens <- data.frame(ENSEMBL = c("ENS1", "ENS2", "ENS8", "ENS3", "ENS3alt", "ENS4", "ENS9"),
                  in_current = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
                  seq_region = c("1", "2", "3", "6", "HSCHR6_MHC_COX", "7", NA),
                  ens_hgnc_id = c("HGNC:1", "HGNC:2", "HGNC:9", "HGNC:3", "HGNC:3", "HGNC:4", NA))
rs <- resolve_symbols(cand, hg, ens, primary = as.character(1:22))
st <- setNames(rs$symbols$status, rs$symbols$symbol)
ok(identical(unname(st[c("A", "B", "C", "D", "E")]), c("resolved", "resolved", "resolved", "unresolved", "unresolved")),
   "resolve_symbols: renamed and alt-haplotype resolve; unannotated duplicate and absent gene do not")
cl <- setNames(rs$rows$class, rs$rows$row)
ok(cl[["X"]] == "distinct_gene" &&
     cl[["C_ENS3alt"]] == "alt_annotation" &&
     cl[["D_ENS9"]] == "unresolved_not_in_current_ensembl" &&
     rs$symbols$row[rs$symbols$symbol == "B"] == "oldB",
   "resolve_symbols: row classes (distinct gene, alt annotation, unannotated)")
ok(setequal(linked_rows(rs), c("A", "oldB", "C_ENS3", "C_ENS3alt", "D_ENS4")),
   "linked_rows: selected, HGNC-ID and alt-annotation rows; not distinct or unannotated")

f <- theta_pipeline(counts, meta, paste0("g", 1:10), 2, "mono", "nk")
pm <- f$pb$meta; Y <- f$norm$Y
tg <- sapply(paste0("g", 1:10), function(g) mean(sapply(c("d1", "d2", "d3"), function(d) {
  p <- function(l, c) pm$profile[pm$donor == d & pm$lineage == l & pm$cond == c]
  (Y[g, p("mono", "stim")] - Y[g, p("mono", "ctrl")]) - (Y[g, p("nk", "stim")] - Y[g, p("nk", "ctrl")])
})))
ok(isTRUE(all.equal(mean(tg), f$theta, tolerance = 1e-12)), "theta = mean of gene-level contrasts")

eff <- f$norm$norm$eff_lib_size[match(colnames(f$pb$counts), f$norm$norm$profile)]
Yf <- logcpm_fixed(f$pb$counts, eff, 1)
ok(max(abs(Yf - log2(sweep(f$pb$counts, 2, eff, "/") * 1e6 + 1))) < 1e-12, "logcpm_fixed = log2(1e6 C / E + 1)")
cst <- sapply(c("mono", "nk"), function(l) {
  j <- pm$lineage == l
  cL <- 1e6 * 2 / mean(eff[j])
  D <- Y[, j] - logcpm_fixed(f$pb$counts[, j], eff[j], cL)
  max(D) - min(D)
})
ok(all(cst < 1e-10), "edgeR logCPM = log2(CPM + c_L) + lineage constant, c_L = 1e6 p / mean_L(E)")

cat("\nALL TESTS PASSED\n")
