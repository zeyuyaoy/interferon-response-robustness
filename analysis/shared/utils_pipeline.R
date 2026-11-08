load_kang <- function(path, fetch = muscData::Kang18_8vs8) {
  if (!file.exists(path)) return(fetch())
  e <- new.env()
  nm <- load(path, envir = e)
  e[[nm[1]]]
}

MONO <- "mono"; NK <- "nk"
LIN <- c(mono = "CD14+ monocytes", nk = "NK cells")

.status <- new.env()
.status$rows <- list()
.status$warn <- character()

not_estimable <- function(reason) {
  stop(structure(class = c("not_estimable", "error", "condition"),
                 list(message = reason, call = NULL)))
}

component_warning <- function(msg) .status$warn <- c(.status$warn, msg)

run_component <- function(name, fn) {
  .status$warn <- character()
  res <- tryCatch(
    list(ok = TRUE, value = fn()),
    not_estimable = function(e) list(ok = FALSE, status = "not_estimable",
                                     reason = conditionMessage(e)),
    error = function(e) list(ok = FALSE, status = "failed",
                             reason = paste0(conditionMessage(e), " [call: ",
                                             paste(deparse(conditionCall(e)), collapse = " "), "]"))
  )
  if (isTRUE(res$ok)) {
    status <- if (length(.status$warn)) "completed_with_warning" else "completed"
    reason <- paste(.status$warn, collapse = "; ")
  } else {
    status <- res$status
    reason <- res$reason
  }
  .status$rows[[length(.status$rows) + 1]] <-
    data.frame(component = name, status = status, reason = reason)
  message(sprintf("[%s] %s %s", status, name, if (nzchar(reason)) paste("-", reason) else ""))
  invisible(if (isTRUE(res$ok)) res$value else NULL)
}

status_table <- function() do.call(rbind, .status$rows)

write_result <- function(df, file, run_id, dir = CFG$step_dir) {
  df <- as.data.frame(df)
  df$run_id <- if (nrow(df)) run_id else character()
  utils::write.csv(df, file.path(dir, file), row.names = FALSE)
  invisible(df)
}

reldiff <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300), 0)
qtl <- function(x, p) stats::quantile(x, p, names = FALSE, type = 7)

write_gz <- function(df, file, dir = CFG$step_dir) {
  gz <- gzfile(file.path(dir, file), "w"); utils::write.csv(df, gz, row.names = FALSE); close(gz)
}


analysis_cells <- function(keep = NULL) {
  m <- ST$meta
  sel <- m$donor %in% ST$analysis$donors
  if (!is.null(keep)) sel <- sel & keep
  list(counts = ST$counts[, sel, drop = FALSE], meta = m[sel,])
}

need_core <- function() {
  if (is.null(ST$analysis)) not_estimable("no analysis set")
  if (is.null(ST$core_rows)) not_estimable("core genes not uniquely mapped")
}

need_primary <- function() {
  need_core()
  if (is.null(ST$primary)) not_estimable("primary estimate unavailable")
}

need_fits <- function() {
  need_primary()
  if (is.null(ST$fit_sens)) not_estimable("prior.count sensitivity fit unavailable")
}

hallmark_symbols <- function(file) unique(utils::read.delim(file)$gene_symbol)
