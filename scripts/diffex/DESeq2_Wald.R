#!/usr/bin/env Rscript

# Public data algae pipeline: Wald test #
## Given a set of salmon quant transcriptomes, calculate
## differentially expressed genes between pairs of conditions,
## passed by a "C_condition" on a "metadata.txt" tsv file.

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

# Help section ----
if("--help" %in% args) {
  cat("
      DESeq2_Wald.R            Simple pairwise comparison for transcripts
      
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
      
      Example run:
      Rscript ./scripts/diffex/DESeq2_Wald.R \
	    --input='results/06_diffex/salmon_quant_PRJNA401507' \
	    --metadata='data/intel/PRJNA401507/metadata.txt' \
	    --map='results/04_trinity_assembly/trinity_PRJNA401507/Trinity.fasta.gene_trans_map.transposed' \
	    --outDir='results/06_diffex/diffex_PRJNA401507' \
	    --threads=40
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

### Libraries
require("DESeq2", quietly = T)
require("tximport", quietly = T)
require("gtools", quietly = T)
require("BiocParallel", quietly = T)

### Input vars renaming
i.salmon = argsL[['input']]
i.threads = argsL[['threads']]
i.outdir = argsL[['outDir']]

### Inputs
i.metadata <- read.delim(file = argsL[['metadata']], stringsAsFactors=TRUE)
i.map <- read.delim(argsL[['map']], stringsAsFactors=FALSE, header = F)

# Fold change filter.
## (2021-02) Maybe will change to user input at a later time.
i.l2fc <- 02

#### Add colnames to tx2gene map
colnames(i.map) <- c("transcript_id", "gene_id")

### File validation tests
if (! is.data.frame(i.metadata)) {
  stop("Metadata file doesn't appear to be a table")
}

### Register thread number
register(MulticoreParam(i.threads))
# [r2c] SnowParam
# register(SnowParam(1))

# Dummy vars ----
### [r2c] Dummy vars
# i.salmon <- "test/01_deseq2_Wald/results/06_diffex/salmon_quant_PRJNA609760"
# i.threads <- 6
# i.outdir <- "test/01_deseq2_Wald/results/06_diffex/output"
# i.metadata <- read.delim("test/01_deseq2_Wald/data/intel/PRJNA609760/metadata.txt")
# i.map <- read.delim("test/01_deseq2_Wald/results/04_trinity_assembly/trinity_PRJNA609760/gene_trans_map.transposed", stringsAsFactors = F, header = F)

# tximport ----
files <- file.path(i.salmon, i.metadata$Run, "quant.sf")

### Import with tximport
txi <- tximport(files, type="salmon", tx2gene = i.map)

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
    # deseq.tmp.results <- deseq.tmp.results[!is.na(deseq.tmp.results$padj),]
    # deseq.tmp.results <- deseq.tmp.results[abs(deseq.tmp.results$log2FoldChange) > 1, ]
    ## Round numbers (less disk space)
    deseq.tmp.results$pvalue <- round(deseq.tmp.results$pvalue, 4)
    deseq.tmp.results$log2FoldChange <- round(deseq.tmp.results$log2FoldChange, 4)
    deseq.tmp.results$padj <- round(deseq.tmp.results$padj, 4)
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

## [r2c]
# backup <- allDiffExWald
# allDiffExWald <- backup

# Generate master table ----
## New addition. Generates one large table with all differential expression.
## Must harness allDiffExWald list object.

## Start with empty df
# all.df <- data.frame(matrix(nrow=dim(allDiffExWald[[1]])[1],ncol=0))
count_flag = 0

for (df in names(allDiffExWald)){
    
    ## Adds gene names
    # allDiffExWald[df][[1]]$gene <- row.names(allDiffExWald[df][[1]])
    
    ## Captures each diffex.df
    current.df <- as.data.frame(allDiffExWald[df])
    
    ## Add gene name
    current.df$gene <- row.names(allDiffExWald[df][[1]])
    
    ## Keeps only useful info
    current.df <- current.df[,c(2,5,6,7)]
    
    ## Stores the information.
    ## if 0: new table
    ## if =! 0: merge tables
    
    if (count_flag == 0) {
        keep.df <- current.df
    } else {
        keep.df <- merge(keep.df, current.df, by = "gene")
        
    }
    
    
    ## Adds to counter
    count_flag = count_flag + 1
    
}




# Writing outputs ----

## Output: master table (keep.df)
### Merge gene and isoform
i.map.tmp <- i.map
colnames(i.map.tmp)[2] <- "gene"
keep.df <- merge(keep.df, i.map.tmp, by = "gene")

### Put gene as col2, isoform as col1
allColsLength <- dim(keep.df)[2]
keep.df <- keep.df[,c(1,allColsLength,2:(allColsLength-1))]

## Generates output path
outpath_master = paste(i.outdir, "MasterTable", ".tsv", sep = "")

## Escreve a tabela
write.table(x = keep.df, file = outpath_master, quote = F,
            sep = "\t", row.names = F, col.names = T)


## Output: Per condition table ----
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
  
  ## OutTable: pick only useful columns
  outTable <- curTable[,c(5,1,2,3,4)]
  
  ## Filter clean outs: L2FC
  outTable <- outTable[abs(outTable$L2FC) > 1,]
  
  ## Filter by: pval
  outTable <- outTable[abs(outTable$`p-val`) < 0.12,]
  
  ## Filter out: NA p-adj
  outTable <- outTable[! is.na(outTable$`p-adj`),]
  
  ## Significant at 0.05
  outTable$sigpadj005 <- ifelse(outTable$`p-adj` < 0.05, "yes_padj005", "no_padj005")
  outTable$sigpadj010 <- ifelse(outTable$`p-adj` < 0.10, "yes_padj010", "no_padj010")
  outTable$sigpval005 <- ifelse(outTable$`p-val` < 0.05, "yes_pval005", "no_pval005")
  outTable$sigpval010 <- ifelse(outTable$`p-val` < 0.10, "yes_pval010", "no_pval010")
  
  ## Generates output path
  outpath = paste(i.outdir, "DiffEx_", curTableName, ".tsv", sep = "")
  
  ## Escreve a tabela
  write.table(x = outTable, file = outpath, quote = F,
    sep = "\t", row.names = F, col.names = T)
}

# Testing zone ----

# Test: Not taking a specific, new column
# Fixed: mismatched object #
# test <- DESeq2_Wald(txiData, metadata, interestColumn)





