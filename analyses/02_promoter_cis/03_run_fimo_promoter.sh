#!/usr/bin/env bash
# 03_run_fimo_promoter.sh
# Análise de elementos cis via FIMO (MEME Suite) + JASPAR2024 CORE Plants
# Alternativa server-side ao PlantCARE — reprodutível, scriptável, citável.
#
# Referência: Grant CE, Bailey TL, Noble WS (2011) FIMO: scanning for occurrences
#   of a given motif. Bioinformatics 27(7):1017-1018.
# Base de motifs: JASPAR2024 CORE Plantae non-redundant (Rauluseviciute et al. 2024)
#
# Pré-requisito no repositório:
#   databases/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt
#   (baixar no Windows em https://jaspar.elixir.no/download/data/2024/CORE/ e commitar)
#
# Uso:
#   bash 03_run_fimo_promoter.sh [upstream_fasta] [output_dir]
#   bash 03_run_fimo_promoter.sh rlp_upstream_2kb.fa fimo_out_49genes

set -euo pipefail
cd "$(dirname "$0")"

ENV="kerson-paper"
UPSTREAM_FA="${1:-rlp_upstream_2kb.fa}"
OUT_DIR="${2:-fimo_out_49genes}"
DB_DIR="databases"
JASPAR_DB="${DB_DIR}/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
PVAL_THRESH="1e-4"    # p-value cutoff para hits (equivale a PlantCARE 'High' confidence)
QVAL_THRESH="0.05"    # FDR cutoff

echo "================================================"
echo "FIMO — Cis-element scan em LRR-RLP promotores"
echo "================================================"
echo "  Upstream FASTA : $UPSTREAM_FA"
echo "  Banco de motifs: $JASPAR_DB"
echo "  p-value cutoff : $PVAL_THRESH"
echo "  Saída          : $OUT_DIR"
echo ""

# ── Verificar dependências ────────────────────────────────────────────────────
if ! mamba run -n "$ENV" fimo --version &>/dev/null 2>&1; then
    echo "Instalando MEME Suite via mamba..."
    mamba install -n "$ENV" -c bioconda -c conda-forge -y meme
fi

# ── Verificar arquivos de entrada ─────────────────────────────────────────────
if [ ! -f "$UPSTREAM_FA" ]; then
    echo "ERRO: '$UPSTREAM_FA' não encontrado."
    echo "Rode primeiro: python3 fetch_upstream_local.py (Windows) e commitar no repo."
    exit 1
fi

if [ ! -f "$JASPAR_DB" ]; then
    echo "ERRO: banco JASPAR não encontrado em '$JASPAR_DB'."
    echo ""
    echo "Para baixar no Windows:"
    echo "  Invoke-WebRequest -Uri 'https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt' \\"
    echo "    -OutFile 'analyses/02_promoter_cis/databases/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt'"
    echo ""
    echo "Depois: git add analyses/02_promoter_cis/databases/ && git commit -m 'Add JASPAR2024 plant motifs'"
    echo "No servidor: git pull"
    exit 1
fi

# ── Contar sequências ─────────────────────────────────────────────────────────
N_SEQ=$(grep -c "^>" "$UPSTREAM_FA" || true)
N_MOT=$(grep -c "^MOTIF" "$JASPAR_DB" || true)
echo "Sequências no FASTA : $N_SEQ"
echo "Motifs no JASPAR2024: $N_MOT"
echo ""

# ── Rodar FIMO ───────────────────────────────────────────────────────────────
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

mamba run -n "$ENV" fimo \
    --o "$OUT_DIR" \
    --thresh "$PVAL_THRESH" \
    --qv-thresh "$QVAL_THRESH" \
    --max-stored-scores 2000000 \
    --parse-genomic-coord \
    "$JASPAR_DB" \
    "$UPSTREAM_FA"

echo ""
echo "FIMO concluído!"
N_HITS=$(tail -n +2 "$OUT_DIR/fimo.tsv" | grep -v "^#" | wc -l || echo "?")
echo "  Total de hits: $N_HITS"
echo "  Arquivo principal: $OUT_DIR/fimo.tsv"
echo ""
echo "Próximo passo:"
echo "  python3 parse_fimo_results.py $OUT_DIR/fimo.tsv"
echo "  Rscript 02_plot_plantcare_heatmap.R fimo_parsed_counts.csv"
