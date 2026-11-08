f <- "release/checksums.sha256"
p <- sort(unique(c(".gitignore", "LICENSE", "README.md",
                   unlist(lapply(c("analysis", "results", "references", "release"), list.files,
                                 recursive = TRUE, full.names = TRUE)))))
p <- setdiff(p, f)
h <- unname(tools::sha256sum(p))
stopifnot(!anyNA(h))
writeLines(paste0(h, "  ", p), f)
cat(length(p), "files checksummed\n")
