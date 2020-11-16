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


# Parsing user args ----
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1


# Missing inputs/outputs ----
## Missing input salmon
if(is.null(argsL[['input']])) {
  sink(stderr())
  cat("\nERROR: Missing input salmon_dir 'input' !\n\n")
  sink()
  q(save="no")
}

## Missing for metadata
if(is.null(argsL[['metadata']])) {
  sink(stderr())
  cat("\nERROR: Missing metadata file 'metadata' !\n\n")
  sink()
  q(save="no")
}

## Missing map
if(is.null(argsL[['map']])) {
  sink(stderr())
  cat("\nERROR: Missing map file 'map' !\n\n")
  sink()
  q(save="no")
}

## Missing threads
if(is.null(argsL[['threads']])) {
  sink(stderr())
  cat("\nERROR: Missing threads number !\n\n")
  sink()
  q(save="no")
}

## Missing outdir
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

### Input vars renaming
i.salmon = argsL[['input']]
i.metadata = argsL[['metadata']]
i.map = argsL[['map']]
i.outdir = argsL[['outDir']]
i.threads = argsL[['threads']]


### [r2c] Dummy vars
# dummy.input <- "test/01_deseq2_Wald/results/06_diffex/salmon_quant_PRJNA609760"
# dummy.metadata <- read.delim("E:/Workspace/GIT/algae2/test/01_deseq2_Wald/data/intel/PRJNA609760/metadata.txt", stringsAsFactors=T)
# dummy.map <- read.delim("E:/Workspace/GIT/algae2/test/01_deseq2_Wald/results/04_trinity_assembly/trinity_PRJNA609760/gene_trans_map.transposed", stringsAsFactors=FALSE)
# dummy.threads <- 6
# dummy.outdir <- "test/01_deseq2_Wald/results/06_diffex/output"

### Register thread number
register(MulticoreParam(i.threads))
# [r2c] SnowParam
# register(SnowParam(1))

### Quant files path
files <- file.path(i.salmon, i.metadata$Run, "quant.sf")

### Import with tximport
txi <- tximport(files, type="salmon", i.map = dummy.map)

## Wald Test Function per treatment (Wald test and output)
  ## 1. DESeq2 analysis over DESeq2 object for a given set of metadata
    ## a. Generates combinatorics over treatment combinations
  ## 2. Stores all outputs for every metadata treatment control pair
  ## 3. Stores the matrix in a (nxnxp) data.frame() structure

DESeq2_Wald <- function(InTxiData, InMetadata, InInterestColumn) {

  ## Default formulae
  InFormula = paste("~", InInterestColumn, sep="")
  
  ## Create DDS object with given metadata
  ddsWald <- DESeqDataSetFromTximport(InTxiData,
                                      colData = InMetadata,
                                      design = as.formula(InFormula))
  
  ## Runs DESeq2
  out.ddsWald <- DESeq(ddsWald, parallel=T, fitType = "local")
  
  ## Extracts outputs
  # Treatment combinations
  treat.combs <- as.data.frame(permutations(n = length(levels(InMetadata[,InInterestColumn])),
                                            r = 2,
                                            v = levels(InMetadata[,InInterestColumn]),
                                            repeats.allowed = F))
  
  # Names of treatments
  treat.combs$text.combs = paste(treat.combs[,1], treat.combs[,2], sep='_')
  
  # DESeq out list
  deseq.out <- list()
  
  # Loops over treatment
  for (trat in 1:nrow(treat.combs)) {
    ## Gets first trat name
    trat1 <- as.character(treat.combs[trat,1])
    ## Gets second trat name
    trat2 <- as.character(treat.combs[trat,2])
    ## Gets paired name
    trat.name <- as.character(treat.combs[trat,3])
    ## Get numbers from DESeq2
    deseq.tmp.results <- results(out.ddsWald, contrast=c(InInterestColumn, trat1, trat2))
    ## Filter out bad stuff
    deseq.tmp.results <- deseq.tmp.results[!is.na(deseq.tmp.results$padj),]
    deseq.tmp.results <- deseq.tmp.results[abs(deseq.tmp.results$log2FoldChange) > 1, ]
    ## Saves in a list, treatment by treatment
    deseq.out[trat.name] <- deseq.tmp.results
  }
  
  ## Final output tables
  DiffExGenesLists <- deseq.out
  return(DiffExGenesLists)
}
  

## DESeq2_Wald F'n Testing ----
# txiData <- txi
# metadata <- dummy.metadata
# interestColumn = "C_Dummy"
## Works:
# Trat_Temperature <- DESeq2_Wald(txiData, metadata, interestColumn)

### Extract all treatments
treatmentCols <- grep(pattern = "^C_", x = colnames(i.metadata))
timeCourseCols <- grep(pattern = "^TC_", x = colnames(i.metadata))

## Loop: generates one DESeq2 object per column ----
allDiffExWald <- list()

for (cols in treatmentCols){
  ## Generates colNames
  ## Grabs the interestColumn
  interestColumn <- colnames(i.metadata[cols])
  
  ## Runs diffEx
  tmp.diffex <- DESeq2_Wald(txi, i.metadata, interestColumn)
  
  ## Stores output
  ### Parses list into objects
  all.df <- unlist(tmp.diffex)
  
  for (listElement in 1:length(all.df)) {
    ## Saves DF inside another list
    allDiffExWald[names(all.df[listElement])] <- all.df[listElement]
    }
  
  
}


# Writing outputs ----
# 1. One table per condition, with condition as a name
# 2. Saves to target output directory
for (df in 1:length(allDiffExWald)) {
  curTable <- as.data.frame(allDiffExWald[df])
  curTableName <- names(allDiffExWald[df])
  
  ## Only L2FC, p-val, p-adj
  curTable <- curTable[,c(2,5,6)]
  colnames(curTable) <- c("L2FC", "p-val", "p-adj")
  curTable$gene_id <- row.names(curTable)
  
  ## Merge with Isoform
  curTable <- merge(curTable, i.map, by="gene_id", all.x=T)
  
  ## OutTable
  outTable <- curTable[,c(5,1,2,3,4)]
  
  ## Generates output path
  outpath = paste(dummy.outdir, "DiffEx_", curTableName, ".tsv", sep = "")
  
  ## Escreve a tabela
  write.table(x = outTable, file = outpath, quote = F,
    sep = "\t", row.names = F, col.names = T)
}

# Testing zone ----

# Test: Not taking a specific, new column
# Fixed: mismatched object #
# test <- DESeq2_Wald(txiData, metadata, interestColumn)





