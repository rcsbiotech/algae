## YAML configs
# configfile:"config/config.yaml"

## SRA Accessions
ACCESSION=["SRR11215921","SRR11215922","SRR11215923"]
BIOPROJECT="PRJNA609760"

## Regra geral
rule all:
    input:
        sra=expand("results/00_sra_files/{Bioproject}/{accession}/{accession}.sra",
            accession=ACCESSION, Bioproject=BIOPROJECT),
        vdb=expand("results/01_vdb_logs/{Bioproject}/{accession}.vdb.txt",
            accession=ACCESSION, Bioproject=BIOPROJECT)

## Baixa os SRA
rule prefetch:
    params:
        outdir="results/00_sra_files/{Bioproject}/"
    output:
        "results/00_sra_files/{Bioproject}/{accession}/{accession}.sra"
    shell:
        "prefetch {ACCESSION} "
        "--output-directory {params.outdir}"

## Valida o md5hash
rule md5validate:
    input:
        raw_sra=expand("results/00_sra_files/{Bioproject}/{accession}/{accession}.sra",
            accession=ACCESSION, Bioproject=BIOPROJECT)
    output:
        vdb_log=expand("results/01_vdb_logs/{Bioproject}/{accession}.vdb.txt",
            accession=ACCESSION, Bioproject=BIOPROJECT)
    shell:
        "vdb-validate {input.raw_sra} >> "
        "{output.vdb_log}"

