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
    r-gridextra

# PyMOL (separado — às vezes tem conflito de canal)
mamba install -n "$ENV_NAME" -c conda-forge -y pymol-open-source || \
    echo "AVISO: pymol-open-source não instalado (opcional para Dia 5; pode instalar depois)"

# KaKs_Calculator 2.0 — binário não disponível no conda; instalar manualmente
echo ""
echo "=== KaKs_Calculator 2.0 (Dia 4) ==="
KAKS_BIN="$(conda run -n ${ENV_NAME} python3 -c "import sys; print(sys.prefix)")/bin/KaKs_Calculator"
if [ ! -f "$KAKS_BIN" ]; then
    echo "Baixando KaKs_Calculator 2.0..."
    TMPDIR_KAKS=$(mktemp -d)
    # Binário Linux 64-bit compilado (GitHub mirror)
    wget -q -O "${TMPDIR_KAKS}/KaKs_Calculator2.0.tar.gz" \
        "https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download" \
        --timeout=60 || true

    if [ -f "${TMPDIR_KAKS}/KaKs_Calculator2.0.tar.gz" ]; then
        tar -xzf "${TMPDIR_KAKS}/KaKs_Calculator2.0.tar.gz" -C "$TMPDIR_KAKS"
        # Tentar compilar
        SRC_DIR=$(find "$TMPDIR_KAKS" -maxdepth 3 -name "Makefile" | head -1 | xargs -r dirname)
        if [ -n "$SRC_DIR" ]; then
            make -C "$SRC_DIR" 2>/dev/null && \
            cp "${SRC_DIR}/KaKs_Calculator" "$KAKS_BIN" && \
            chmod +x "$KAKS_BIN" && \
            echo "KaKs_Calculator instalado em: $KAKS_BIN" || true
        fi
        rm -rf "$TMPDIR_KAKS"
    fi

    # Fallback: yn00 do PAML como alternativa
    if [ ! -f "$KAKS_BIN" ]; then
        echo "KaKs_Calculator não compilado — instalando yn00 (PAML) como alternativa..."
        mamba install -n "$ENV_NAME" -c bioconda -y paml 2>/dev/null && \
            echo "yn00 disponível como alternativa a KaKs_Calculator" || \
            echo "AVISO: instale KaKs_Calculator 2.0 manualmente (sourceforge.net/projects/kakscalculator2)"
    fi
else
    echo "KaKs_Calculator já instalado."
fi

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
