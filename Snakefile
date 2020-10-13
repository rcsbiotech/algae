## YAML configs
# configfile:"config/config.yaml"

## SRA Accessions
ACCESSION=["SRR11215921","SRR11215922","SRR11215923"]
BIOPROJECT="PRJNA609760"

## RUN Setup
# Choose between { Paired | Single | Mixed }
VAR_LIB_LAYOUT="PAIRED"
# Choose between { Genome-Guided/GG | de novo }
VAR_ASSEMBLY_STRATEGY="GG"

## LINKS/REFS/NETWORK Setup
### 1. Search for genome in Taxonomy, or nearest
### 2. Click the latest/most interesting
### 3. On top left, right click copy link "Download sequences in FASTA for genome"
REF_URL_NCBI_GENOME=["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/754/935/GCF_002754935.1_ASM275493v1/GCF_002754935.1_ASM275493v1_genomic.fna.gz",]
REF_URL_NCBI_ANNOTATION=["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/754/935/GCF_002754935.1_ASM275493v1/GCF_002754935.1_ASM275493v1_genomic.gff.gz",]
REF_NCBI_GENOME_CODE=["GCF_002754935.1_ASM275493v1",]

## All outputs rule
### Check for commas (,) when adding new rules or it'll fizzle!
rule all:
    input:
        sra=expand("results/00_sra_files/{Bioproject}/{accession}/{accession}.sra",
            accession=ACCESSION, Bioproject=BIOPROJECT),
        vdb=expand("results/01_vdb_logs/{Bioproject}/{accession}.vdb.txt",
            accession=ACCESSION, Bioproject=BIOPROJECT),
        fqf=expand("results/02_fastq_dump/{Bioproject}/{accession}/{accession}_1.fastq",
            accession=ACCESSION, Bioproject=BIOPROJECT),
        ref_genome=expand("data/genomes/{Bioproject}/{GenomeCode}.fasta",
            Bioproject=BIOPROJECT, GenomeCode=REF_NCBI_GENOME_CODE),
        ref_annot=expand("data/genomes/{Bioproject}/{GenomeCode}.gff",
            Bioproject=BIOPROJECT, GenomeCode=REF_NCBI_GENOME_CODE)

## Baixa os SRA
rule prefetch:
    params:
        outdir="results/00_sra_files/{Bioproject}/",
        force="all"
    output:
        "results/00_sra_files/{Bioproject}/{accession}/{accession}.sra"
    shell:
        "prefetch {ACCESSION} --force {params.force} "
        "--output-directory {params.outdir}"

## Valida o md5hash
rule md5validate:
    input:
        "results/00_sra_files/{Bioproject}/{accession}/{accession}.sra"
    output:
        "results/01_vdb_logs/{Bioproject}/{accession}.vdb.txt"
    shell:
        "vdb-validate {input} 2> {output}"
        
## Fastq-dump for paired end
rule fastqdump:
    input:
        "results/00_sra_files/{Bioproject}/{accession}/{accession}.sra"
    output:
        fqdir=directory("results/02_fastq_dump/{Bioproject}/{accession}"),
        fq1="results/02_fastq_dump/{Bioproject}/{accession}/{accession}_1.fastq",
        fq2="results/02_fastq_dump/{Bioproject}/{accession}/{accession}_2.fastq"
    shell:
        "fastq-dump --split-e {input} "
        "--outdir {output.fqdir}"
        
## Imports reference genome
rule import_reference_genome:
    output:
        reference_genome="data/genome/{REF_NCBI_GENOME_CODE}.fna",
        reference_annotation="data/genome/{REF_NCBI_GENOME_CODE}.gff"
    shell:
        "wget -qO {REF_URL_NCBI_GENOME} >> data/genome/{REF_NCBI_GENOME_CODE}.fna.gz && "
        "wget -qO {REF_URL_NCBI_ANNOTATION} >> data/genome/{REF_NCBI_GENOME_CODE}.gff.gz && "
        "gzip -d data/genome/*.gz"

## Genome download
### IF, genome-guided, downloads
### ELSE, run nothing, just print a message




