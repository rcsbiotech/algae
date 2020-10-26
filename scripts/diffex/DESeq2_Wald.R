# Public data algae pipeline: Wald test #
## Given a set of salmon quant transcriptomes, calculate
## differentially expressed genes between pairs of conditions,
## passed by a "C_condition" on a "metadata.txt" tsv file.

#!/usr/bin/Rscript
# /usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

# Help section ----
if("--help" %in% args) {
  cat("
      deparse_taxa.R            Split taxonomy table with uniques and collapsed
      
      Arguments:
      
      --input=Path             - Path to salmon quant files
      --metadata=metadata.tsv  - Path to metadata.txt file, with header and one line
                                 per sample
      --outDir=folderPath      - Name to write all differential expression output
      --threads=int            - Number of threads to run DESeq2 parallel
      --help                   - Print this help
 
      Example:
      ./deparse_taxa.R --x=\"input1.txt\" --out=\"output.txt\" 
      
      Rafael Correia da Silva
      GCCRC - Genomics for Climate Change Research Center
      www.gccrc.unicamp.br
      rcs.biotec@gmail.com
      \n\n")
  
  q(save="no")
}
# Help section end ----

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1


## Missing for X
if(is.null(argsL[['input']])) {
  sink(stderr())
  cat("\nERROR: Missing input file 'input' !\n\n")
  sink()
  q(save="no")
}

## Missing for X
if(is.null(argsL[['metadata']])) {
  sink(stderr())
  cat("\nERROR: Missing metadata file 'metadata' !\n\n")
  sink()
  q(save="no")
}

## Missing output
if(is.null(argsL[['outDir']])) {
  sink(stderr())
  cat("\nERROR: Missing output directory !\n\n")
  sink()
  q(save="no")
}

### Dummy examples (to test implementation)

### Libraries
require("DESeq2")
require("tximport")
require("gtools")
require("BiocParallel")