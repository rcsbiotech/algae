# Diretório
base_dir <- "."

# Nomes das amostras a serem analisadas

sample_id = c("SRR5684403","SRR5684404","SRR5684405","SRR5684406","SRR5684407","SRR5684408","SRR5684409","SRR5684410","SRR5684411","SRR5684412","SRR5684413","SRR5684414","SRR5684415","SRR5684416","SRR5684417")

# Lista dos diretórios (paste)
# lista dos diretorios criados pelo kallisto (output). Lembrar de manter a msm ordem q no sample_id

kallisto_out = c(file.path(base_dir,"tempo1_kallisto_transcripts_index.out"),file.path(base_dir,"tempo4_kallisto_transcripts_index.out"))