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
10. DESeq2 Wald test (RScript)
11. Single-end branch
12. Diversified, fast annotation


## About the metadata file

- Must have 4 columns:
    1. "Bioproject" column with the Bioproject ID for all files
    2. "Run" column - SRR Run - one ID per file
    3. "Layout" column - either PAIRED or SINGLE
    4. Any number of treatment columns starting with "C_"
    5. Any number of dependent treatments starting with "TC_"
- No spaces or dashes are allowed in any header or data;

## Releases

- 0.1: From reads to assemble
- 0.2: (Current) Differential expression with DESeq2 (Wald/LRT test)
- 0.3: (wip) Multiple annotation sources.

## About the pipeline - Single-end, Paired-end and mixed reads.

- Paired-end and Single-end: fully implemented and functional. R2 read is passed as parameter, not input, throughout the entire pipeline.
- Mixed samples: since it's the rarest case, still to be implemented.

## To do (2020-10-26)

1. Fork decision: if alignment < 60%, do denovo, else, do genome-guided
2. Add conda specification to each pipe step/environments
3. Format to recommended Snakemake structure
4. (done) Add logs
5. Add salmon fork: if genome-guided, build a decoy transcriptome for additional mapping precision
6. Google Cloud auth//integration in a private (.gitignore) file, and safe auth options
7. DiffEx plots (RScript)
9. Merge annotated genes with sequences
10. Add "snakemake" as base environment
11. Add an external reference to databases