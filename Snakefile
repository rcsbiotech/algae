## [.. Instructions (wip) ..] ##
# 1. Create a metadata.txt file inside data/intel/{BIOPROJECT}
# 2. Create sample blank files {ACCESSION}.blank inside data/intel/{BIOPROJECT/samples
# 3. Run the pipeline


## YAML configs
# configfile:"config/config.yaml"

# Libs
import pandas as pd

# [.. General vars ..] #
VAR_ANALYSIS_ID="TEST001"

# [ .. Public data vars .. ] #
## SRA Accessions
### Download method: {PREFETCH | WGET}
VAR_DOWNLOAD_METHOD="WGET"
### To-do: read from metadata - S4/json file
BIOPROJECT="PRJNA609760"

# Metadata
## Parses metadata from "metadata.txt"
metadata_path = ("data/intel/", BIOPROJECT, "/metadata.txt")
metadata = pd.read_table("".join(metadata_path)).set_index("Run", drop=False)
ACCESSION = list(metadata.index)
## Works @ 2020-10-20
# print(ACCESSION)

## RUN Setup
# Choose between { Paired | Single | Mixed }
VAR_LIB_LAYOUT="PAIRED"
# Reads length (on average)
VAR_READS_LENGTH=150
# STAR Overhang
VAR_OVERHANG=VAR_READS_LENGTH-1

## Assembly variables
### "de_novo" or "genome_guided"
### TO DO:"de_novo" if no genome or poor alignment representation in gene space;
### "genome_guided" otherwise
VAR_ASSEMBLY_STRATEGY="de_novo"

## LINKS/REFS/NETWORK Setup
### 1. Search for genome in Taxonomy, or nearest
### 2. Click the latest/most interesting
### 3. On top left, right click copy link "Download sequences in FASTA for genome"
REF_URL_NCBI_GENOME=["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/754/935/GCF_002754935.1_ASM275493v1/GCF_002754935.1_ASM275493v1_genomic.fna.gz",]
REF_URL_NCBI_ANNOTATION=["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/754/935/GCF_002754935.1_ASM275493v1/GCF_002754935.1_ASM275493v1_genomic.gtf.gz",]
REF_NCBI_GENOME_CODE=["GCF_002754935.1_ASM275493v1",]

## THREADS per process
THREADS_STAR_INDEX=10

## "All outputs" rule
### For the pipeline to run back to back, every output must be passed as an input
### to all, otherwise one must explicit tell the shell which output must be generated.
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
            Bioproject=BIOPROJECT, accession=ACCESSION),
        trinity_asm=expand("results/04_trinity_assembly/trinity_{Bioproject}/Trinity.fasta",
            Bioproject=BIOPROJECT),
        trinity_renamed=expand("results/04_trinity_assembly/trinity_{Bioproject}/Trinity.renamed.fasta",
            Bioproject=BIOPROJECT),
        longest_orfs=expand("results/05_annotation/01_transdecoder_{Bioproject}/longest_orfs.pep",
            Bioproject=BIOPROJECT),
        salmon_index=expand("results/06_diffex/salmon_index_{Bioproject}",
            Bioproject=BIOPROJECT),
        salmon_quant=expand("results/06_diffex/salmon_quant_{Bioproject}/{acc}",
            Bioproject=BIOPROJECT, acc=ACCESSION)

## Baixa os SRA
rule prefetch:
    input:
        "data/intel/{Bioproject}/samples/{accession}.blank"
    params:
        outdir="results/00_sra_files/{Bioproject}/{accession}",
        force="all"
    output:
        "results/00_sra_files/{Bioproject}/{accession}/{accession}.sra"
    run:
        if VAR_DOWNLOAD_METHOD == "PREFETCH":
            shell("prefetch {wildcards.accession} --force {params.force} "
            "--output-directory {params.outdir}")
        if VAR_DOWNLOAD_METHOD == "WGET":
            shell("mkdir -p {params.outdir} && "
                "wget http://sra-download.ncbi.nlm.nih.gov/srapub/{wildcards.accession} "
                "-O {params.outdir}/{wildcards.accession}.sra")
        

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
        
# [ .. Trinity assembly .. ] #
# some R code join(map()) is run inside shell due to 
# trinity needing input fqfiles as a comma-separated objects,
# not spaces, which is default snakemake behaviour.
# to-do: implement genome-guided bam assembly.
# also to-do: recheck previous alignments.
rule trinity_asm:
    input:
        fq_read1=expand("results/02_fastq_dump/{code}/{acc}/{acc}_1.fastq",
            code=BIOPROJECT, acc=ACCESSION),
        fq_read2=expand("results/02_fastq_dump/{code}/{acc}/{acc}_2.fastq",
            code=BIOPROJECT, acc=ACCESSION)
    output:
        trinity_fasta="results/04_trinity_assembly/trinity_{code}/Trinity.fasta"
    params:
        trinity_dir=directory("results/04_trinity_assembly/trinity_{code}"),
        max_memory="600G",
        max_read_cov=20
    threads: 50
    run:
        if VAR_ASSEMBLY_STRATEGY == "de_novo":
            read1_parsed = ",".join(map(str, input.fq_read1))
            read2_parsed = ",".join(map(str, input.fq_read2))
            shell("Trinity "
            "--seqType fq "
            "--left {read1_parsed} "
            "--right {read2_parsed} "
            "--output {params.trinity_dir} "
            "--max_memory {params.max_memory} "
            "--normalize_max_read_cov {params.max_read_cov} "
            "--CPU {threads}")
        # Placeholder: still gotta add bam > sam > sort sam.    
        if VAR_ASSEMBLY_STRATEGY == "genome_guided":
            shell("Trinity --version")

# rename fasta headers by analysis code
rule after_asm_rename:
    input:
        trinity_fasta="results/04_trinity_assembly/trinity_{Bioproject}/Trinity.fasta"
    output:
        ren_out="results/04_trinity_assembly/trinity_{Bioproject}/Trinity.renamed.fasta"
    params:
        ID=expand("{id}", id=VAR_ANALYSIS_ID),
        trinity_dir="results/04_trinity_assembly/trinity_{Bioproject}"
    run:
        shell("sed 's/TRINITY/{params.ID}_TRINITY/' {input} > {output.ren_out} && "
            "sed -i 's/ .*//' {output.ren_out} && "
            "map=$(ls {params.trinity_dir}/*gene_trans_map) && "
            "sed -i 's/TRINITY/{params.ID}_TRINITY/g' ${{map}}")
            
# --- Quantification --- #

## Salmon step 1. Transcriptome index
### The index is needed to map reads to a transcriptome, and count.
rule salmon_index:
    input:
        trinity_fasta="results/04_trinity_assembly/trinity_{Bioproject}/Trinity.renamed.fasta"
    output:
        salmon_index=directory("results/06_diffex/salmon_index_{Bioproject}")
    threads: 50
    run:
        shell("salmon index "
        "--transcripts {input.trinity_fasta} "
        "--index {output.salmon_index} "
        "--threads {threads}")

## Salmon step 2. Actual quantification
## About parameters (--validateMappings)
### Enables selective alignment of the sequencing reads when mapping them to the 
### transcriptome. This can improve both the sensitivity and specificity of mapping 
### and, as a result, can improve quantification accuracy (salmon docs).

rule salmon_quant:
    input:
        fq1="results/02_fastq_dump/{Bioproject}/{acc}/{acc}_1.fastq", 
        fq2="results/02_fastq_dump/{Bioproject}/{acc}/{acc}_2.fastq",
        salmon_index="results/06_diffex/salmon_index_{Bioproject}"
    output:
        salmon_quant=directory("results/06_diffex/salmon_quant_{Bioproject}/{acc}")
    params:
        mapdir="results/04_trinity_assembly/trinity_{Bioproject}",
        libType="A",
        scoreFraction="0.70"
    threads: 8
    run:
        shell(
            "map=$(ls {params.mapdir}/*gene_trans_map) && "
            "salmon quant "
            "--index {input.salmon_index} "
            "--libType {params.libType} "
            "--mates1 {input.fq1} "
            "--mates2 {input.fq2} "
            "-o {output.salmon_quant} "
            "--geneMap ${{map}} "
            "--validateMappings "
            "--gcBias --seqBias "
            "--threads {threads}"
        )







# --- DiffEx ---- #
## Inputs:
### - A set of salmon directories
### - A metadata.txt file ([to do] to be specified on instructions.md)
## Outputs:
### - All sets of diffex genes with values (p, p-adj, deltaFC, deltaquant)
           



            
# --- Annotation --- #
## TransDecoder - pull peptides
rule transdecoder_get_peptides:
    input:
        "results/04_trinity_assembly/trinity_{Bioproject}/Trinity.renamed.fasta"
    output:
        "results/05_annotation/01_transdecoder_{Bioproject}/longest_orfs.pep"
    params:
        inpdir="results/04_trinity_assembly/trinity_{Bioproject}",
        outdir="results/05_annotation/01_transdecoder_{Bioproject}",
        minsize=70
    run:
        shell("map=$(ls {params.inpdir}/*gene_trans_map) && "
	"TransDecoder.LongOrfs -t {input} "
        "-m {params.minsize} "
        "--output_dir {params.outdir} "
        "--gene_trans_map ${{map}}")



