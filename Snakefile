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
# STAR Overhang
VAR_OVERHANG=VAR_READS_LENGTH-1

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
            Bioproject=BIOPROJECT, GenomeCode=REF_NCBI_GENOME_CODE),
	star_index=expand("data/genomes/{Bioproject}/star_index",
            Bioproject=BIOPROJECT),
        bam=expand("results/03_STAR_alignment/{Bioproject}/{accession}.Aligned.out.bam",
            Bioproject=BIOPROJECT, accession=ACCESSION)

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
        genome=expand("data/genomes/{code}/{ref}.fna",
            code=BIOPROJECT, ref=REF_NCBI_GENOME_CODE),
        annotf=expand("data/genomes/{code}/{ref}.gtf",
            code=BIOPROJECT, ref=REF_NCBI_GENOME_CODE)
    output:
        directory("data/genomes/{BIOPROJECT}/star_index")
    params:
        star_overhang={VAR_OVERHANG}
    threads: 10
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.annotf} "
        "--genomeSAindexNbases 10 "
        "--sjdbOverhang {params.star_overhang}"

## Star align only for paired
## (to-do) OBS: there's a shell method for acquiring basenames (fqb), consider
## Swapping to python lambda functions if possible

## EXPAND adds all files
## W/OUT EXPAND runs one per file!
rule star_align:
    input:
        index=expand("data/genomes/{code}/star_index", code=BIOPROJECT),
        annotf=expand("data/genomes/{code}/{ref}.gtf", code=BIOPROJECT, ref=REF_NCBI_GENOME_CODE),
        fq1="results/02_fastq_dump/{code}/{acc}/{acc}_1.fastq", 
        fq2="results/02_fastq_dump/{code}/{acc}/{acc}_2.fastq"
    output:
        "results/03_STAR_alignment/{code}/{acc}.Aligned.out.bam"
    params:
        multimapNmax=20,
        outSAMtype="BAM Unsorted",
        outdir="results/03_STAR_alignment/{code}",
        unmapped="Within"
    threads: 30
    shell:
        "fqb=$(basename {input.fq1} _1.fastq) && "
        "STAR "
        "--runThreadN {threads} "
        "--genomeDir {input.index} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--sjdbGTFfile {input.annotf} "
        "--outSAMstrandField intronMotif "
        "--outFilterIntronMotifs RemoveNoncanonical "
        "--outFilterMultimapNmax {params.multimapNmax} "
        "--outFileNamePrefix {params.outdir}/${{fqb}}. "
        "--outSAMunmapped {params.unmapped} "
        "--outSAMtype {params.outSAMtype}"
        
## Genome download
### IF, genome-guided, downloads
### ELSE, run nothing, just print a message




