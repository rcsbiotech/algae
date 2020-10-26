# Algae Prospection Project

My algae gene prospection project.

## Current protocol

01. Download SRA from prefetch
02. Check MD5 with vdb-validate
03. Convert to FASTQ with fastq-dump
04. Download some genome and annotation
05. Map reads to genome and annotation
06. De novo assembly
07. TransDecoder
08. [to do] Eggnog-mapper annotation
09. [wip] Salmon quantification
10. DESeq2 Wald test (RScript)
11. DESeq2 LRT test (RScript)

## To do (2020-10-26)

1. Fork decision: if alignment < 60%, do denovo, else, do genome-guided
2. Add conda specification to each pipe step/environments
3. Format to recommended Snakemake structure
4. Add logs
5. Add salmon fork: if genome-guided, build a decoy transcriptome for additional mapping precision
6. Google Cloud auth//integration in a private (.gitignore) file, and safe auth options

## Current results tree

results
├── 00_sra_files
│   └── PRJNA609760
├── 01_vdb_logs
│   └── PRJNA609760
├── 02_fastq_dump
│   └── PRJNA609760
├── 03_STAR_alignment
│   └── PRJNA609760
├── 04_trinity_assembly
│   └── trinity_PRJNA609760
└── 05_annotation
    ├── 01_transdecoder_PRJNA609760
    └── 01_transdecoder_PRJNA609760.__checkpoints_longorfs