library(Biostrings)
library(stringr) 
library(seqinr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!= 3) 
   stop("Three arguments must be supplied: Fprobe, Rprobe and the sequence", call=FALSE)[1]

seqfile <- args[3]
sequence <- read.fasta(seqfile, as.string = TRUE, seqonly = TRUE)[[1]]
Fprobe <- args[1]
Rprobe <- args[2]
sequence <- DNAString(str_replace_all(sequence, "[\r\n]", ""))
rg <- ranges(matchProbePair(Fprobe, Rprobe, sequence, algorithm="auto", logfile=NULL, verbose=FALSE))
result_list <- list(start(rg), end(rg))
result_list