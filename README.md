# Algae Prospection Project

My algae gene prospection project.

## Current protocol

01. Download SRA from prefetch
02. Check MD5 with vdb-validate
03. Convert to FASTQ with fastq-dump
04. Download some genome and annotation
05. Map reads to genome and annotation
06. De novo assembly (Trinity, if unmappable)
07. TransDecoder to get peptides
08. Eggnog-mapper annotation
09. Salmon quantification
10. [wip] DESeq2 Wald test (RScript)
11. [to do] DESeq2 LRT test (RScript)

## To do (2020-10-26)

1. Fork decision: if alignment < 60%, do denovo, else, do genome-guided
2. Add conda specification to each pipe step/environments
3. Format to recommended Snakemake structure
4. Add logs
5. Add salmon fork: if genome-guided, build a decoy transcriptome for additional mapping precision
6. Google Cloud auth//integration in a private (.gitignore) file, and safe auth options
7. DiffEx plots (RScript)
9. Merge annotated genes with sequences
