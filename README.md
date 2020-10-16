# Algae Prospection Project

My algae gene prospection project.

## Current protocol

1. Download SRA from prefetch
2. Check MD5 with vdb-validate
3. Convert to FASTQ with fastq-dump
4. Download some genome and annotation
5. Map reads to genome and annotation
6. NEXT: decision fork (de novo/genome-guided)

## To do (2020-10-16)

1. Fork decision: if alignment < 60%, do denovo, else, do genome-guided

