TOL <- 1e-12
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "analysis/time_course/crossing"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
OUT <- "results/time_course/crossing"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE)
chk <- list()
add <- function(check, pass, value) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))

spec <- utils::read.delim("analysis/time_course/crossing/specification.tsv", stringsAsFactors = FALSE)
add("specification is locked", grepl("^LOCKED", spec$value[spec$key == "status"]), spec$value[spec$key == "status"])

tr <- rd("results/time_course/g3_trajectory.csv")
s <- utils::read.delim("analysis/time_course/samples.tsv", stringsAsFactors = FALSE)
s <- s[s$time_h > 0,]
grid <- lapply(split(s$time_h, s$donor), function(x) sort(unique(x)))
obs <- lapply(split(tr$time_h, tr$donor), sort)
add("donors are D34, D38, D39", identical(sort(names(obs)), c("D34", "D38", "D39")), paste(names(obs), collapse = ";"))
add("trajectory times per donor equal the samples.tsv grid excluding 0h", identical(obs[names(grid)], grid),
    paste(sapply(names(obs), function(d) paste0(d, ":", paste(obs[[d]], collapse = ","))), collapse = " "))
add("one trajectory row per donor x time", !anyDuplicated(tr[, c("donor", "time_h")]), nrow(tr))
add("time label equals time_h", all(tr$time == paste0(tr$time_h, "h")), "")
add("delta equals R_mono - R_nk", max(abs(tr$delta - (tr$R_mono - tr$R_nk))) < TOL, max(abs(tr$delta - (tr$R_mono - tr$R_nk))))
add("no missing or exactly zero delta", all(is.finite(tr$delta) & tr$delta != 0), min(abs(tr$delta)))
add("single upstream run_id", length(unique(tr$run_id)) == 1, paste(unique(tr$run_id), collapse = ";"))
dbt <- rd("results/time_course/g3_delta_by_time.csv")
sens_times <- tapply(dbt$time, dbt$setting, function(x) paste(sort(unique(x)), collapse = ","))
add("sensitivity settings lack full trajectories (primary only)", all(sens_times[names(sens_times) != "primary"] == "24h,4h"),
    paste(names(sens_times), sens_times, sep = ":", collapse = " "))

cr <- do.call(rbind, lapply(split(tr, tr$donor), function(z) {
  z <- z[order(z$time_h),]
  i <- which(diff(sign(z$delta)) != 0)
  if (!length(i)) return(NULL)
  j <- i + 1
  data.frame(setting = "primary", donor = z$donor[1], crossing = seq_along(i), n_crossings_donor = length(i),
             direction = ifelse(z$delta[i] < 0, "neg_to_pos", "pos_to_neg"),
             t_lo_h = z$time_h[i], t_hi_h = z$time_h[j], delta_lo = z$delta[i], delta_hi = z$delta[j],
             t_linear_interp_h = z$time_h[i] + (z$time_h[j] - z$time_h[i]) * z$delta[i] / (z$delta[i] - z$delta[j]))
}))
rownames(cr) <- NULL
add("interpolation points lie inside their intervals", all(cr$t_linear_interp_h > cr$t_lo_h & cr$t_linear_interp_h < cr$t_hi_h), "")
add("crossing intervals join adjacent observed times",
    all(mapply(function(d, a, b) match(b, grid[[d]]) - match(a, grid[[d]]) == 1, cr$donor, cr$t_lo_h, cr$t_hi_h)), "")

merge_iv <- function(lo, hi) {
  if (!length(lo)) return(matrix(numeric(0), ncol = 2))
  o <- order(lo); lo <- lo[o]; hi <- hi[o]; m <- c(lo[1], hi[1]); res <- NULL
  for (k in seq_along(lo)[-1]) if (lo[k] <= m[2]) m[2] <- max(m[2], hi[k]) else { res <- rbind(res, m); m <- c(lo[k], hi[k]) }
  unname(rbind(res, m))
}
intersect_iv <- function(a, b) {
  res <- NULL
  for (i in seq_len(nrow(a))) for (j in seq_len(nrow(b))) {
    lo <- max(a[i, 1], b[j, 1]); hi <- min(a[i, 2], b[j, 2])
    if (lo <= hi) res <- rbind(res, c(lo, hi))
  }
  if (is.null(res)) matrix(numeric(0), ncol = 2) else merge_iv(res[, 1], res[, 2])
}
fmt <- function(m) if (!nrow(m)) "empty" else paste(sprintf("[%g, %g] h", m[, 1], m[, 2]), collapse = " U ")
donors <- sort(names(obs))
rng <- do.call(rbind, lapply(c("all", "neg_to_pos", "pos_to_neg"), function(dir) {
  x <- if (dir == "all") cr else cr[cr$direction == dir,]
  per <- lapply(donors, function(d) merge_iv(x$t_lo_h[x$donor == d], x$t_hi_h[x$donor == d]))
  u <- merge_iv(x$t_lo_h, x$t_hi_h)
  it <- Reduce(intersect_iv, per)
  data.frame(setting = "primary", direction = dir, n_crossings = nrow(x), n_donors_with_crossing = length(unique(x$donor)),
             donor_intervals = paste(donors, sapply(per, fmt), sep = ": ", collapse = "; "),
             union = fmt(u), union_min_h = if (nrow(u)) min(u) else NA, union_max_h = if (nrow(u)) max(u) else NA,
             intersection = fmt(it), intersection_is_single_time = nrow(it) > 0 && all(it[, 1] == it[, 2]),
             interp_min_h = if (nrow(x)) min(x$t_linear_interp_h) else NA, interp_max_h = if (nrow(x)) max(x$t_linear_interp_h) else NA)
}))
add("every donor has at least one sign change", rng$n_donors_with_crossing[rng$direction == "all"] == 3, rng$n_donors_with_crossing[rng$direction == "all"])

cr$run_id <- RUN_ID; rng$run_id <- RUN_ID
cc <- do.call(rbind, chk); cc$run_id <- RUN_ID
utils::write.csv(cr, file.path(OUT, "crossings_by_donor.csv"), row.names = FALSE)
utils::write.csv(rng, file.path(OUT, "crossing_range.csv"), row.names = FALSE)
utils::write.csv(cc, file.path(OUT, "checks.csv"), row.names = FALSE)
if (!all(cc$pass)) stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; "))
print(cc[, 1:3]); print(cr[, 1:10]); print(rng[, -ncol(rng)])
