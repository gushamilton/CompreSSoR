#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
path <- Sys.getenv("COMPRESSOR_SCALING_SOURCE")
x <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(path))), showProgress = FALSE)
setnames(x, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "p"))
debug_rows <- as.integer(Sys.getenv("COMPRESSOR_DEBUG_ROWS", unset = "1000000"))
x <- x[seq_len(min(debug_rows, nrow(x)))]
lengths <- c(`1`=248956422,`2`=242193529,`3`=198295559,`4`=190214555,`5`=181538259,`6`=170805979,`7`=159345973,`8`=145138636,`9`=138394717,`10`=133797422,`11`=135086622,`12`=133275309,`13`=114364328,`14`=107043718,`15`=101991189,`16`=90338345,`17`=83257441,`18`=80373285,`19`=58617616,`20`=64444167,`21`=46709983,`22`=50818468)
off <- c(0,cumsum(as.numeric(lengths)))[seq_along(lengths)]; names(off) <- names(lengths)
base <- c(A=0L,C=1L,G=2L,T=3L)
pos <- unname(off[as.character(x$chrom)]) + as.numeric(x$pos) - 1
sub <- as.integer(base[x$ref])*4L + as.integer(base[x$alt])
key <- paste(pos, sub, sep=":")
cat("rows", nrow(x), "unique_full", length(unique(key)), "duplicate_count", sum(duplicated(key)), "first_duplicate", anyDuplicated(key), "\n")
if (anyDuplicated(key)) {
  duplicate_key <- key[which(duplicated(key))[1L]]
  print(x[key == duplicate_key])
  cat("source_key values:\n")
  print(paste(x$chrom[key == duplicate_key], x$pos[key == duplicate_key],
              x$ref[key == duplicate_key], x$alt[key == duplicate_key], sep = ":"))
}
