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

