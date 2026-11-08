f <- "release/checksums.sha256"
p <- sub("^[0-9a-f]+  ", "", readLines(f))
p <- sort(unique(c(p[!startsWith(p, "hirisa/")], list.files("hirisa", recursive = TRUE, full.names = TRUE),
                   "release/update_checksums.R")))
h <- unname(tools::sha256sum(p))
stopifnot(!anyNA(h))
writeLines(paste0(h, "  ", p), f)
cat(length(p), "files,", sum(startsWith(p, "hirisa/")), "under hirisa/\n")
