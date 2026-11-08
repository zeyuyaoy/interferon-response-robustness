rds_slot_reader <- function(path, target, max_len = 5e6, chunk = 2^26) {
  con <- gzcon(pipe(paste("gzip -dc", shQuote(path)), "rb"))
  on.exit(close(con))
  st <- new.env()
  st$refs <- list(); st$skipped <- 0; st$n_items <- 0; st$found <- NULL
  int <- function(n = 1) readBin(con, "integer", n, size = 4, endian = "big")
  skip <- function(bytes) {
    st$skipped <- st$skipped + bytes
    while (bytes > 0) { k <- min(bytes, chunk); r <- readBin(con, "raw", k); if (length(r) != k) stop("unexpected end of stream"); bytes <- bytes - k }
  }
  len <- function() { n <- int(); if (n == -1L) { u <- int(); l <- int(); n <- u * 2^32 + l }; n }
  addref <- function(x) { st$refs[[length(st$refs) + 1]] <- x; x }
  chars <- function() {
    n <- int()
    if (n == -1L) return(NA_character_)
    if (n == 0L) return("")
    rawToChar(readBin(con, "raw", n))
  }
  strvec <- function() { int(); n <- int(); vapply(seq_len(n), function(i) item(), character(1)) }
  attach_attr <- function(v, a) {
    if (is.null(a)) return(v)
    if (is.list(a)) for (nm in names(a)) if (!is.null(a[[nm]]) && nzchar(nm) && !inherits(a[[nm]], "rds_skipped"))
      v <- tryCatch({ attr(v, nm) <- a[[nm]]; v }, error = function(e) v)
    v
  }
  item <- function() {
    st$n_items <- st$n_items + 1
    flags <- int()
    type <- bitwAnd(flags, 255L)
    has_attr <- bitwAnd(flags, 512L) != 0
    has_tag <- bitwAnd(flags, 1024L) != 0
    if (type %in% c(254L, 253L, 252L, 251L, 250L, 242L, 241L)) return(NULL)
    if (type == 255L) { i <- bitwShiftR(flags, 8L); if (i == 0L) i <- int(); return(st$refs[[i]]) }
    if (type %in% c(247L, 248L, 249L)) return(addref(structure(list(strvec()), class = "rds_ns")))
    if (type == 9L) return(chars())
    if (type == 1L) return(addref(as.character(item())))
    if (type == 4L) {
      i <- length(st$refs) + 1; st$refs[[i]] <- structure(list(), class = "rds_env")
      int(); item(); item(); item(); item()
      return(st$refs[[i]])
    }
    if (type %in% c(2L, 6L, 3L, 5L, 17L, 239L, 240L)) {
      out <- list(); nms <- character()
      repeat {
        if (has_attr) item()
        tag <- if (has_tag) item() else ""
        val <- item()
        if (identical(tag, target) && is.null(st$found) && is.list(val) && any(c("orig.ident", "nCount_RNA") %in% attr(val, "names"))) { st$found <- val; stop(structure(class = c("rds_found", "condition"), list(message = "found", call = NULL))) }
        out[length(out) + 1] <- list(val); nms <- c(nms, if (is.null(tag)) "" else as.character(tag))
        flags <- int()
        type2 <- bitwAnd(flags, 255L)
        if (!type2 %in% c(2L, 6L, 3L, 5L, 17L, 239L, 240L)) {
          if (type2 == 254L) break
          rest <- item_after_flags(flags)
          out[length(out) + 1] <- list(rest); nms <- c(nms, "")
          break
        }
        has_attr <- bitwAnd(flags, 512L) != 0
        has_tag <- bitwAnd(flags, 1024L) != 0
      }
      names(out) <- nms
      return(out)
    }
    item_after_flags(flags)
  }
  item_after_flags <- function(flags) {
    type <- bitwAnd(flags, 255L)
    has_attr <- bitwAnd(flags, 512L) != 0
    if (type %in% c(254L, 253L, 252L, 251L, 250L, 242L, 241L)) return(NULL)
    if (type == 255L) { i <- bitwShiftR(flags, 8L); if (i == 0L) i <- int(); return(st$refs[[i]]) }
    if (!type %in% c(10L, 13L, 14L, 15L, 16L, 19L, 20L, 24L, 25L, 22L, 23L, 7L, 8L, 238L, 9L, 1L, 4L, 247L, 248L, 249L)) {
      if (type %in% c(2L, 6L, 3L, 5L, 17L, 239L, 240L)) stop("pairlist continuation not expected here")
      stop(sprintf("unsupported SEXP type %d", type))
    }
    if (type %in% c(9L, 1L, 4L, 247L, 248L, 249L)) stop("unexpected type in continuation")
    v <- switch(as.character(type),
      "10" = , "13" = { n <- len(); if (n > max_len) { skip(4 * n); structure(list(), class = "rds_skipped") } else { x <- int(n); if (type == 10L) as.logical(x) else x } },
      "14" = { n <- len(); if (n > max_len) { skip(8 * n); structure(list(), class = "rds_skipped") } else readBin(con, "double", n, size = 8, endian = "big") },
      "15" = { n <- len(); skip(16 * n); structure(list(), class = "rds_skipped") },
      "24" = { n <- len(); skip(n); structure(list(), class = "rds_skipped") },
      "16" = { n <- len(); if (n > max_len) { for (i in seq_len(n)) item(); structure(list(), class = "rds_skipped") } else vapply(seq_len(n), function(i) { x <- item(); if (is.null(x)) NA_character_ else x }, character(1)) },
      "19" = , "20" = { n <- len(); if (n > max_len) { for (i in seq_len(n)) item(); structure(list(), class = "rds_skipped") } else { l <- vector("list", n); for (i in seq_len(n)) { x <- item(); if (!is.null(x)) l[[i]] <- x }; l } },
      "25" = structure(list(), class = "rds_s4"),
      "22" = { addref(NULL); item(); item(); structure(list(), class = "rds_extptr") },
      "23" = { addref(NULL); structure(list(), class = "rds_weakref") },
      "7" = , "8" = { n <- int(); readBin(con, "raw", n); structure(list(), class = "rds_builtin") },
      "238" = { info <- item(); state <- item(); a <- item(); structure(list(info = info, state = state, attr = a), class = "rds_altrep") })
    if (type == 238L) {
      a_alt <- v$attr
      cls <- v$info[[1]]
      if (identical(cls, "compact_intseq") && is.numeric(v$state) && length(v$state) == 3) v <- as.integer(seq(v$state[2], by = v$state[3], length.out = v$state[1]))
      else if (identical(cls, "compact_realseq") && is.numeric(v$state) && length(v$state) == 3) v <- seq(v$state[2], by = v$state[3], length.out = v$state[1])
      else if (identical(cls, "wrap_integer") || identical(cls, "wrap_real") || identical(cls, "wrap_string") || identical(cls, "wrap_logical")) v <- v$state[[1]]
      else if (identical(cls, "deferred_string")) v <- as.character(v$state[[1]])
    }
    if (type == 238L) return(if (inherits(v, "rds_altrep")) v else attach_attr(v, a_alt))
    if (has_attr) v <- attach_attr(v, item())
    v
  }
  con_header <- readBin(con, "raw", 2)
  if (!identical(rawToChar(con_header), "X\n")) stop("not an XDR serialization")
  ver <- int(3)
  if (ver[1] == 3L) { n <- int(); readBin(con, "raw", n) }
  res <- tryCatch({ item(); NULL }, rds_found = function(e) st$found)
  list(value = res, found = !is.null(res), n_items = st$n_items, skipped_bytes = st$skipped, version = ver)
}

as_metadata_frame <- function(v) {
  if (!is.data.frame(v)) stop("meta.data was not reconstructed as a data.frame")
  data.frame(cell = attr(v, "row.names"), v, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
}
