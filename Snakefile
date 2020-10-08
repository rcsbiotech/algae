## YAML configs
# configfile:"config/config.yaml"

## SRA Accessions
ACCESSION=["SRR11215921","SRR11215922","SRR11215923"]
Bioproject="PRJNA609760"

## Regra geral
rule all:
    input:
        sra=expand("results/00_sra_files/{accession}.sra", accession=ACCESSION)

## Script do tutorial
rule prefetch:
    params:
        outdir="./results/00_sra_files/{Bioproject}/"
    output:
        "results/00_sra_files/{Bioproject}/{accession}/{accession}.sra"
    shell:
        "prefetch {ACCESSION} "
        "--output-directory {params.outdir}"

