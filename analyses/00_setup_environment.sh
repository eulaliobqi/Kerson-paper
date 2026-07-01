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

command -v mamba &>/dev/null || { echo "ERRO: mamba não encontrado. Instale miniforge3."; exit 1; }

if mamba env list | grep -q "^${ENV_NAME}"; then
    echo "Ambiente '${ENV_NAME}' já existe. Atualizando pacotes..."
else
    mamba create -n "$ENV_NAME" python=3.11 -y
fi

# Bioinformática (bioconda + conda-forge)
mamba install -n "$ENV_NAME" -c bioconda -c conda-forge -y \
    bedtools \
    samtools \
    hmmer \
    mafft \
    last \
    mcscanx \
    biopython \
    wget

# R e pacotes de visualização
mamba install -n "$ENV_NAME" -c conda-forge -c bioconda -y \
    r-base \
    r-tidyverse \
    r-pheatmap \
    r-ggplot2 \
    r-rcolorbrewer \
    r-gggenes \
    r-scales \
    r-gridextra \
    r-ggnewscale \
    r-ggrepel

# PyMOL (separado — às vezes tem conflito de canal)
mamba install -n "$ENV_NAME" -c conda-forge -y pymol-open-source || \
    echo "AVISO: pymol-open-source não instalado (opcional para Dia 5; pode instalar depois)"

# Ka/Ks — implementado via BioPython (NG86); nenhum executável externo necessário
# O módulo Bio.codonalign.codonseq.cal_dn_ds já está incluído no biopython instalado acima
echo ""
echo "=== Ka/Ks (Dia 4): BioPython NG86 já disponível no ambiente ==="
mamba run -n "${ENV_NAME}" python3 -c "from Bio.codonalign import CodonAlignment; print('BioPython Ka/Ks OK')"

echo ""
echo "================================================================="
echo "Ambiente '${ENV_NAME}' pronto!"
echo "================================================================="
echo ""
echo "Ativar com:"
echo "  mamba activate ${ENV_NAME}"
echo ""
echo "Verificar instalação:"
echo "  mamba run -n ${ENV_NAME} bedtools --version"
echo "  mamba run -n ${ENV_NAME} hmmscan -h | head -2"
echo "  mamba run -n ${ENV_NAME} mafft --version"
echo "  mamba run -n ${ENV_NAME} Rscript -e 'library(pheatmap); library(gggenes); cat(\"R OK\n\")'"
