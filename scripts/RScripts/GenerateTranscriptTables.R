#!/usr/bin/env Rscript

# GenerateTranscriptTables.R

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
      GenerateTranscriptTables.R
      
      ---
      From differentially expression tables and some annotation,
      parses and outputs with the information:
      
      1. Nยบ of treatments in which the gene was DE (column Multisignificance)
      2. What are those treatments (column Stress)
      
      ---
      
      Arguments:
      
      --master=path.tsv        - Path to DESeq2 MasterTable
      --samba=samba.tsv        - Path to samba coding/noncoding table
      --emapper=emapper.tsv    - Path to eggnog-mapper clean table
      --out=filename.tsv       - Path to output (TranscriptTable)
      --help                   - Print this help

      Example run:
      Rscript ./scripts/diffex/GenerateTranscriptTables.R \
	    --master='results/06_diffex/salmon_quant_PRJNA401507' \
	    --samba='data/intel/PRJNA401507/metadata.txt' \
	    --emapper='results/04_trinity_assembly/trinity_PRJNA401507/Trinity.fasta.gene_trans_map.transposed' \
	    --out='results/06_diffex/diffex_PRJNA401507'
	    
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
if(is.null(argsL[['master']])) {
    sink(stderr())
    cat("\nERROR: Missing input masterTable 'master' !\n\n")
    sink()
    q(save="no")
}

## Missing for metadata
if(is.null(argsL[['samba']])) {
    sink(stderr())
    cat("\nERROR: Missing samba tsv file 'samba' !\n\n")
    sink()
    q(save="no")
}

## Missing map
if(is.null(argsL[['emapper']])) {
    sink(stderr())
    cat("\nERROR: Missing emapper file 'emapper' !\n\n")
    sink()
    q(save="no")
}


## Missing outdir
if(is.null(argsL[['out']])) {
    sink(stderr())
    cat("\nERROR: Missing output filename!\n\n")
    sink()
    q(save="no")
}

### Libraries
require("stringi", quietly = T)
require("gtools", quietly = T)

# Description: Stores information prior to sequence clustering
# Inputs:
# 1. Samba.tsv coding and noncoding annotation
# 2. MasterTable.tsv diffex conditions table
# 3. Emapper.tsv crude sequence annotation
# Outputs:
# One TSV file, with headers:
# [1] Transcripts
# [2] Genes
# [3] Peptides
# [4] StressActivity
# [5] Annotation

#### [ Test files ] ####
## 1. Samba coding/noncoding TSV
# df.samba <- read.delim(
#     file = "./test/02_transcriptTable/samba.tsv",
#     stringsAsFactors = F
# )
# 
# ## 2. RNA-Seq MasterTable
# df.masterTable <- read.delim(
#     file = "./test/02_transcriptTable/MasterTable.tsv",
#     stringsAsFactors = F
# )
# 
# ## 3. Emapper table
# df.emapper <- read.delim(
#     file = "./test/02_transcriptTable/test.emapper.annotations",
#     stringsAsFactors = F
# )

#### Actual input ####

# Paths
path.samba.df <- argsL[['samba']]
path.master.df <- argsL[['master']]
path.emapper.df <- argsL[['emapper']]
path.outfile <- argsL[['out']]

# Read files
df.samba <- read.delim(path.samba.df, stringsAsFactors = F)
df.masterTable <- read.delim(path.master.df, stringsAsFactors = F)
df.emapper <- read.delim(path.emapper.df, stringsAsFactors = F)

# Read tables
colnames(df.emapper)[1] <- "transcript_id"
colnames(df.samba)[1] <- "transcript_id"

# [ Annotate stress activity ] ####
## 1. Get all columns with sig padj
## 2. Split per category
## 3. Bind back, split with ";"
padj.cols <- grep("padj", colnames(df.masterTable))

# [[ Case 1. Generic A vs B ]] ####
## Split first 6 columns, since output is always doubled
valid.cols.limit <- (length(padj.cols)/2)
padj.cols <- padj.cols[1:valid.cols.limit]

# Generates a new data frame with transcripts
df.to.transcriptTable <- df.masterTable[,c("transcript_id", "gene")]

## Counts how many times it was significant in stresses
df.to.transcriptTable$Multisignificance <- 0

## Test value per value, for each valid column
for (column_index in 1:length(padj.cols)) {
    ## Makes a tmp replacement
    df.tmp <- df.masterTable
    ## Get one sig column
    column <- padj.cols[column_index]
    ## Get column names
    colname <- colnames(df.masterTable)[column]
    ## Removes ".padj" from end of string
    sigcolname <- strsplit(colname, ".padj")[[1]]
    ## Split treatments:
    ### [1] T1
    ### [2] T2.padj
    ## Replace NAs with 1
    df.masterTable[,column] <- na.replace(df.tmp[,column], 1)
    ## Get significant values
    sigvals <- ifelse(df.tmp[,column] <= 0.05,
                      yes = sigcolname,
                      no = "NS")
    ## Based on which transcript, paste treatment name for
    ## significant ones
    df.to.transcriptTable$Stress <- paste(
        df.to.transcriptTable$Stress,
        sigvals, sep=";")
    
    ## Sum one per significant treatment
    ### Get all significants as 1
    sumvec <- na.replace(as.integer(df.masterTable[,column] <= 0.05), 0)
    ### Sum to table
    df.to.transcriptTable$Multisignificance <- 
        df.to.transcriptTable$Multisignificance + sumvec
    
    
    
}

## Remove bad stuff from start of transcript table
df.to.transcriptTable$Stress <- gsub("^;", "", df.to.transcriptTable$Stress)

## Add KO and annotation
df.emapper.sub <- df.emapper[,c("eggNOG.annot", "KEGG_KOs", "transcript_id")]

## Merge out
df.out <- merge(y = df.emapper.sub,
                x = df.to.transcriptTable,
                by = "transcript_id",
                all.x = T)

## Merge again, with samba
df.out <- merge(x = df.out,
                y = df.samba[,c(1,3)],
                by = "transcript_id",
                all.x = T)


## Output table
## Escreve a tabela
write.table(x = df.out, file = path.outfile, quote = F,
            sep = "\t", row.names = F, col.names = T)







































