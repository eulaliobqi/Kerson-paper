#!/bin/bash
# Setup do ambiente mamba para o projeto Kerson-paper
# Executar UMA VEZ no servidor: eulalio@200.235.143.10
# Tempo estimado: 10-15 minutos
#
# Uso:
#   bash analyses/00_setup_environment.sh

set -euo pipefail

ENV_NAME="kerson-paper"

echo "================================================================="
echo "Criando ambiente mamba: ${ENV_NAME}"
echo "================================================================="

# Verificar se mamba está disponível
command -v mamba &>/dev/null || { echo "ERRO: mamba não encontrado. Instale miniforge3."; exit 1; }

# Criar ambiente (ignorar se já existir)
if mamba env list | grep -q "^${ENV_NAME}"; then
    echo "Ambiente '${ENV_NAME}' já existe. Atualizando pacotes..."
else
    mamba create -n "$ENV_NAME" python=3.11 -y
fi

# Instalar todas as dependências
mamba install -n "$ENV_NAME" -c bioconda -c conda-forge -y \
    bedtools \
    samtools \
    hmmer \
    mafft \
    kaks-calculator \
    last \
    mcscanx \
    biopython \
    wget \
    pymol-open-source \
    r-base \
    r-tidyverse \
    r-pheatmap \
    r-ggplot2 \
    r-rcolorbrewer \
    r-gggenes \
    r-scales \
    r-gridextra

echo ""
echo "================================================================="
echo "Ambiente '${ENV_NAME}' pronto!"
echo "================================================================="
echo ""
echo "Ativar com:"
echo "  mamba activate ${ENV_NAME}"
echo ""
echo "Verificar instalação:"
echo "  bedtools --version"
echo "  hmmscan -h | head -2"
echo "  mafft --version"
echo "  Rscript -e 'library(pheatmap); library(gggenes); cat(\"R OK\n\")'"
