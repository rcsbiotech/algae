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
      DESeq2_Wald.R            Split taxonomy table with uniques and collapsed
      
      Arguments:
      
      --input=Path             - Path to salmon quant files directory
      --metadata=metadata.tsv  - Path to metadata.txt file, with header and one line
                                 per sample
      --map=gene_trans_map     - Path to gene to transcript map
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

### [r2c] Dummy missing libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tximport")

### [r2c] Dummy vars
dummy.input <- "test/01_deseq2_Wald/results/06_diffex/salmon_quant_PRJNA609760"
dummy.metadata <- read.delim("E:/Workspace/GIT/algae2/test/01_deseq2_Wald/data/intel/PRJNA609760/metadata.txt")
dummy.map <- read.delim("E:/Workspace/GIT/algae2/test/01_deseq2_Wald/results/04_trinity_assembly/trinity_PRJNA609760/gene_trans_map.transposed", stringsAsFactors=FALSE)
dummy.threads <- 6

### Register thread number
# register(MulticoreParam(i.threads))
# [r2c] SnowParam
# register(SnowParam(6))

### Quant files path
files <- file.path(dummy.input, dummy.metadata$Run, "quant.sf")

### Extract all treatments
treatmentCols <- grep(pattern = "^C_", x = colnames(dummy.metadata))
timeCourseCols <- grep(pattern = "^TC_", x = colnames(dummy.metadata))

### Import with tximport
txi <- tximport(files, type="salmon", tx2gene = dummy.map)

## Wald Test Function per treatment
DESeq2_Wald <- function(txiData, metadata, design)

# head(txi$counts, 4)
dds <- DESeqDataSetFromTximport(txi, colData = dummy.metadata, design =~ C_Temperature)
deseq.dds <- DESeq(dds, parallel=T, fitType='local')

## Generate combinations of treatments
# Treatment combinations
treat.combs <- as.data.frame(permutations(n = length(levels(dummy.metadata$C_Temperature)),
                            r = 2,
                            v = levels(dummy.metadata$C_Temperature),
                            repeats.allowed = F))

## Paste as text to generate outputs
treat.combs$text.combs = paste(treat.combs[,1], treat.combs[,2], sep='_')














