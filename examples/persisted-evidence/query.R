#!/usr/bin/env Rscript
# Dependencies: base R and the matching installed Rosalind native executable.
# Read counters as strings: base R's numeric type cannot preserve every uint64.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript query.R evidence-dataset.manifest.json subset.bed output.tsv")
}
binary <- Sys.getenv("ROSALIND_BINARY", unset = "rosalind")
command <- c("dataset", "extract", "--dataset", args[[1L]], "--regions", args[[2L]],
             "--fields", "depths,alleles", "--format", "tsv", "-o", args[[3L]])
status <- system2(binary, shQuote(command))
if (status != 0L) stop(sprintf("Native evidence query failed (exit %d)", status))
evidence <- read.delim(args[[3L]], colClasses = "character", check.names = FALSE,
                       comment.char = "", quote = "", na.strings = NULL)
stopifnot(is.character(evidence$callable_depth))
print(head(evidence[c("#contig", "pos", "callable_depth", "a", "c", "g", "t")]))
cat(sprintf("Verified query: %d positions. Receipt: %s.manifest.json\n",
            nrow(evidence), args[[3L]]))
