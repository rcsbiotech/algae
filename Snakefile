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
# Reads length (on average)
VAR_READS_LENGTH=150

## LINKS/REFS/NETWORK Setup
### 1. Search for genome in Taxonomy, or nearest
### 2. Click the latest/most interesting
### 3. On top left, right click copy link "Download sequences in FASTA for genome"
REF_URL_NCBI_GENOME=["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/754/935/GCF_002754935.1_ASM275493v1/GCF_002754935.1_ASM275493v1_genomic.fna.gz",]
REF_URL_NCBI_ANNOTATION=["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/754/935/GCF_002754935.1_ASM275493v1/GCF_002754935.1_ASM275493v1_genomic.gtf.gz",]
REF_NCBI_GENOME_CODE=["GCF_002754935.1_ASM275493v1",]

## THREADS per process
THREADS_STAR_INDEX=10

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
        ref_genome=expand("data/genomes/{Bioproject}/{GenomeCode}.fna",
            Bioproject=BIOPROJECT, GenomeCode=REF_NCBI_GENOME_CODE),
        ref_annot=expand("data/genomes/{Bioproject}/{GenomeCode}.gtf",
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
        reference_genome="data/genomes/{BIOPROJECT}/{REF_NCBI_GENOME_CODE}.fna",
        reference_annotation="data/genomes/{BIOPROJECT}/{REF_NCBI_GENOME_CODE}.gtf"
    shell:
        "wget -q -O data/genomes/{BIOPROJECT}/{REF_NCBI_GENOME_CODE}.fna.gz {REF_URL_NCBI_GENOME} && "
        "wget -q -O data/genomes/{BIOPROJECT}/{REF_NCBI_GENOME_CODE}.gtf.gz {REF_URL_NCBI_ANNOTATION} && "
        "gzip -d data/genomes/{BIOPROJECT}/*.gz"


## STAR Alignment
## 1. Generate genome index
## 2. Align to genome
rule star_genome_generate:
    input:
        ref_genome="data/genomes/{BIOPROJECT}/{REF_NCBI_GENOME_CODE}.fna",
        ref_annotation="data/genomes/{BIOPROJECT}/{REF_NCBI_GENOME_CODE}.gtf"
    output:
        star_index=directory("data/genome/{BIOPROJECT}/star_index")
    params:
        star_overhang={VAR_READS_LENGTH}-1
    threads: {THREADS_STAR_INDEX}
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output.star_index} "
        "--genomeFastaFiles {input.ref_genome} "
        "--sjdbGTFfile {input.ref_annotation} "
        "--genomeSAindexNbases 10 "
        "--sjdbOverhang {params.star_overhang}"
        
## Genome download
### IF, genome-guided, downloads
### ELSE, run nothing, just print a message




