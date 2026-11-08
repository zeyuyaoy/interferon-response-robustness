f <- "release/checksums.sha256"
p <- sub("^[0-9a-f]+  ", "", readLines(f))
tracked <- c("analysis", "results")
p <- p[file.exists(p) & !grepl(paste0("^(", paste(tracked, collapse = "|"), ")/"), p)]
p <- sort(unique(c(p, unlist(lapply(tracked, list.files, recursive = TRUE, full.names = TRUE)), "release/update_checksums.R")))
h <- unname(tools::sha256sum(p))
stopifnot(!anyNA(h))
writeLines(paste0(h, "  ", p), f)
cat(length(p), "files;", sum(grepl(paste0("^(", paste(tracked, collapse = "|"), ")/"), p)), "under", paste(tracked, collapse = ", "), "\n")
