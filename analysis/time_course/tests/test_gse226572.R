for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/time_course/scripts/utils_gse226572.R")) source(f)
suppressPackageStartupMessages(library(Matrix))
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE); cat("ok -", msg, "\n") }

S <- rbind(c(2, 1, 0, 0, 0, 0, 0, 0), c(0, 0, 1.5, 0.2, 0, 0, 0, 0), c(0, 0, 1.5, 0, 0, 0, 0, 0),
           c(1, 1, 0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0, 0, 0), c(1.2, 0, 0, 0, 0, 0, 0, 0))
colnames(S) <- names(MARKERS)
a <- assign_labels(S, cd14_counts = c(3, 0, 0, 2, 0, 0), t_counts = c(0, 0, 2, 0, 0, 0))
ok(identical(a$label, c("CD14 Mono", "NK", "unassigned", "unassigned", "unassigned", "unassigned")),
   "labels: argmax, T-marker exclusion for NK, tie and all-zero unassigned, CD14 count required")
ok(isTRUE(all.equal(a$margin, c(1, 1.3, 1.5, 0, 0, 1.2))), "label margin is best minus second-best score")

mc <- Matrix(c(10, 0, 0, 5), 2, 2, sparse = TRUE, dimnames = list(c("g1", "g2"), c("c1", "c2")))
sc <- marker_scores(mc, total = c(100, 50), sets = list(A = "g1", B = c("g1", "g2")))
ok(isTRUE(all.equal(unname(sc[, "A"]), c(log1p(1000), 0))) && isTRUE(all.equal(unname(sc[1, "B"]), log1p(1000) / 2)),
   "marker score is mean log1p(1e4 * count / total)")

cnt <- Matrix(0, 4, 3, sparse = TRUE, dimnames = list(c("A", "B", "MT-1", "C"), c("x", "y", "z")))
cnt[, "x"] <- c(5, 5, 1, 1); cnt[, "y"] <- c(1, 0, 9, 0); cnt[, "z"] <- c(1, 1, 0, 1)
cc <- call_cells(cnt, rownames(cnt), min_genes = 3, max_mt = 0.2)
ok(identical(cc$called, c(TRUE, FALSE, TRUE)) && isTRUE(all.equal(cc$mt_frac, c(1 / 12, 0.9, 0))), "cell calling uses gene count and mito fraction")

set.seed(2)
G <- 60; genes <- paste0("g", 1:G)
meta <- expand.grid(cell = 1:4, rep = 1:2, time = c("0h", "4h", "24h"), lineage = c("mono", "nk"), donor = c("D1", "D2", "D3"), stringsAsFactors = FALSE)
meta <- meta[!(meta$rep == 2 & meta$time != "0h"),]
meta$cond <- meta$time
cells <- Matrix(rpois(G * nrow(meta), 30), G, nrow(meta), sparse = TRUE, dimnames = list(genes, paste0("c", seq_len(nrow(meta)))))
pb <- pseudobulk(cells, meta)
ok(length(unique(pb$meta$donor)) == 3 && ncol(pb$counts) == 18, "two 0h libraries pool into one profile; donors stay at 3")
j <- meta$donor == "D2" & meta$lineage == "nk" & meta$time == "0h"
ok(all(pb$counts[, "D2|nk|0h"] == rowSums(as.matrix(cells[, j]))) && pb$meta$n_cells[pb$meta$profile == "D2|nk|0h"] == 8,
   "pooled 0h profile is the sum of both libraries")

Y <- matrix(0, 3, 18, dimnames = list(c("p1", "p2", "p3"), pb$meta$profile))
for (k in seq_len(nrow(pb$meta))) {
  m <- pb$meta[k,]
  base <- if (m$lineage == "mono") 1 else 2
  add <- c("0h" = 0, "4h" = if (m$lineage == "mono") 3 else 1, "24h" = if (m$lineage == "mono") 2 else 1.5)[[m$cond]]
  Y[, k] <- base + add + (m$donor == "D3") * c(0, 0, 3) * (m$cond == "24h") * (m$lineage == "nk")
}
tc <- time_contrasts(Y, pb$meta, c("p1", "p2", "p3"), c("24h", "4h"))
ok(identical(sort(unique(tc$time)), c("24h", "4h")) && all(table(tc$time) == 3), "one contrast per donor per time")
g <- gamma_table(tc)
ok(isTRUE(all.equal(g$gamma[g$donor == "D1"], (2 - 1.5) - (3 - 1))) && isTRUE(all.equal(g$gamma[g$donor == "D3"], (2 - (1.5 + 3 / 3)) - (3 - 1))),
   "gamma = delta(24h) - delta(4h), matched by donor")
ok(isTRUE(all.equal(mean(g$gamma), (-1.5 - 1.5 - 2.5) / 3)), "Gamma weights donors equally")
ok(identical(g$sign_change, c(FALSE, FALSE, TRUE)), "sign change flagged only when delta changes sign")

Yr <- matrix(rnorm(3 * 18), 3, 18, dimnames = dimnames(Y))
gr <- gamma_table(time_contrasts(Yr, pb$meta, c("p1", "p2", "p3"), c("4h", "24h")))
pf <- function(d, l, t) pb$meta$profile[pb$meta$donor == d & pb$meta$lineage == l & pb$meta$cond == t]
direct <- sapply(gr$donor, function(d) mean((Yr[, pf(d, "mono", "24h")] - Yr[, pf(d, "mono", "4h")]) -
                                              (Yr[, pf(d, "nk", "24h")] - Yr[, pf(d, "nk", "4h")])))
ok(isTRUE(all.equal(unname(direct), gr$gamma)), "gamma equals [M(24h) - M(4h)] - [N(24h) - N(4h)] without 0h terms")
Y0 <- Yr; z <- pb$meta$cond == "0h"; Y0[, z] <- Y0[, z] + matrix(rnorm(3 * sum(z), sd = 5), 3)
g0 <- gamma_table(time_contrasts(Y0, pb$meta, c("p1", "p2", "p3"), c("4h", "24h")))
ok(isTRUE(all.equal(g0$gamma, gr$gamma)) && !isTRUE(all.equal(g0$delta_late, gr$delta_late)),
   "changing the 0h reference values changes delta but not gamma")

ok(isTRUE(all.equal(signflip_p(c(-1, -2, -3)), 1)) && isTRUE(all.equal(signflip_p(c(1, 2, 3)), 1 / 8)),
   "one-sided sign-flip p with 3 donors has minimum 1/8 (two-sided 1/4)")
eff <- c(a = 2e5, b = 1e6)
cm <- matrix(c(0, 10), 1, 2, dimnames = list("x", c("a", "b")))
ok(isTRUE(all.equal(as.vector(logcpm_fixed(cm, eff, 1)), unname(log2(c(0, 10) / eff * 1e6 + 1)))) && isTRUE(all.equal(1e6 * 2 / mean(eff), 1e6 * 2 / 6e5)),
   "fixed offset and CPM added by prior count")

tmp <- tempfile(); hist <- tempfile()
writeLines(c("key\tvalue", "a\t1"), tmp)
writeLines(paste("2026-01-01T00:00:00Z", unname(tools::md5sum(tmp)), "locked", sep = "\t"), hist)
ok(spec_gate(tmp, hist), "spec gate accepts the locked specification")
cat("b\t2\n", file = tmp, append = TRUE)
ok(!spec_gate(tmp, hist), "spec gate refuses a modified specification")
cat("\nALL TESTS PASSED\n")
