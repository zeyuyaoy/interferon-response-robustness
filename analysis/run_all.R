source("analysis/config.R")
for (f in Sys.glob("analysis/utils_*.R")) source(f)

STEPS <- c("01_quality_control", "02_primary_estimate", "03_donor_influence", "04_cell_subsampling",
           "05_sensitivity_analyses", "06_reference_panels", "07_robustness_criteria", "08_pseudocount_sensitivity",
           "09_expression_levels", "10_normalization_diagnostics", "11_hallmark_gene_sets", "12_verification",
           "13_summary_figures", "14_fixed_offset_measurement")
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- tryCatch(length(system2("git", c("status", "--porcelain", "--", ".", "':(exclude)results'", "':(exclude)figures'"),
                                 stdout = TRUE)) > 0, error = function(e) FALSE)
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
unlink(c(CFG$results_dir, CFG$figures_dir), recursive = TRUE)
dir.create(CFG$figures_dir)
ST <- new.env()
ST$figs <- list(); ST$dims <- list()
message("run_id: ", RUN_ID)

for (s in STEPS) {
  CFG$step <- s
  CFG$step_dir <- file.path(CFG$results_dir, s)
  dir.create(CFG$step_dir, recursive = TRUE)
  tryCatch(source(file.path("analysis", paste0(s, ".R"))), error = function(e) {
    .status$rows[[length(.status$rows) + 1]] <-
      data.frame(component = paste0("script:", s), status = "failed", reason = conditionMessage(e))
    message("[failed] ", s, ": ", conditionMessage(e))
  })
}

write_result(status_table(), "status.csv", RUN_ID, dir = CFG$results_dir)
writeLines(utils::capture.output(sessionInfo()), file.path(CFG$results_dir, "sessionInfo.txt"))
