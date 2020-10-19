# Algae Prospection Project

My algae gene prospection project.

## Current protocol

01. Download SRA from prefetch
02. Check MD5 with vdb-validate
03. Convert to FASTQ with fastq-dump
04. Download some genome and annotation
05. Map reads to genome and annotation
06. De novo assembly
07. (WIP) TransDecoder
08. Eggnog-mapper annotation
09. Salmon quantification
10. DESeq2 Wald test
11. DESeq2 LRT test

## To do (2020-10-19)

1. Fork decision: if alignment < 60%, do denovo, else, do genome-guided
2. Add conda specification to each pipe step


