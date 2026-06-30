#!/bin/bash
# Pipeline Ka/Ks para pares de genes LRR-RLP duplicados em tomate
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Pipeline: CDS FASTA (SGN ITAG4.0) → MAFFT → KaKs_Calculator 2.0 → sumário
#
# Instalar KaKs_Calculator (se ausente):
#   mamba install -n kerson-paper -c bioconda kaks-calculator

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAIRS_FILE="${SCRIPT_DIR}/gene_pairs.tsv"
OUTDIR="${SCRIPT_DIR}/results"
DB_DIR="/home/eulalio/databases/itag4.0"
CDS_FA="${DB_DIR}/ITAG4.0_CDS.fa"

mkdir -p "$OUTDIR"

# ── Exportar variáveis antes dos heredocs Python que as usam ──────────────────
export OUTDIR_PY="$OUTDIR"
export CDS_FA_PY="$CDS_FA"

# ── CDS do ITAG4.0 ────────────────────────────────────────────────────────────
# Fonte 1: arquivo local no repo com os 21 genes necessários para Ka/Ks
REPO_CDS="${SCRIPT_DIR}/ITAG4.0_CDS_kaks_genes.fasta"

if [ ! -f "$CDS_FA" ]; then
    if [ -f "$REPO_CDS" ]; then
        echo "Usando CDS local do repo (21 genes pré-extraídos)..."
        cp "$REPO_CDS" "$CDS_FA"
    else
        echo "Baixando CDS ITAG4.0 do SGN (35 MB — pode falhar no servidor UFV)..."
        wget -c -q --show-progress -O "${DB_DIR}/ITAG4.0_CDS.fasta" \
            "https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_CDS.fasta"
        mv "${DB_DIR}/ITAG4.0_CDS.fasta" "$CDS_FA"
    fi
fi

# ── Processar cada par ────────────────────────────────────────────────────────
while IFS=$'\t' read -r gene_a gene_b tipo evidencia; do
    [[ "$gene_a" =~ ^# ]] && continue
    [[ -z "$gene_a" ]] && continue

    echo ""
    echo "============================================================"
    echo "Par: ${gene_a} × ${gene_b} | ${tipo}"
    echo "============================================================"

    PAIR_DIR="${OUTDIR}/${gene_a}_vs_${gene_b}"
    mkdir -p "$PAIR_DIR"

    # Exportar variáveis do par para o Python heredoc
    export GENE_A="$gene_a"
    export GENE_B="$gene_b"
    export PAIR_DIR_PY="$PAIR_DIR"

    # 1. Extrair CDS (pega a primeira isoforma de cada gene)
    python3 << 'PYEOF'
import os, sys
from Bio import SeqIO

genes     = {os.environ["GENE_A"]: None, os.environ["GENE_B"]: None}
cds_file  = os.environ["CDS_FA_PY"]
pair_dir  = os.environ["PAIR_DIR_PY"]

for rec in SeqIO.parse(cds_file, "fasta"):
    for g in list(genes):
        if g in rec.id and genes[g] is None:
            genes[g] = rec
            break
    if all(v is not None for v in genes.values()):
        break

for gene, rec in genes.items():
    if rec is None:
        print(f"AVISO: CDS não encontrado para {gene}", file=sys.stderr)
        continue
    outpath = os.path.join(pair_dir, f"{gene}_CDS.fa")
    with open(outpath, "w") as f:
        SeqIO.write(rec, f, "fasta")
    print(f"  Extraído: {rec.id} ({len(rec.seq)} nt) → {outpath}")
PYEOF

    CDS_A="${PAIR_DIR}/${gene_a}_CDS.fa"
    CDS_B="${PAIR_DIR}/${gene_b}_CDS.fa"

    [ ! -f "$CDS_A" ] && { echo "SKIP: CDS ausente para ${gene_a}"; continue; }
    [ ! -f "$CDS_B" ] && { echo "SKIP: CDS ausente para ${gene_b}"; continue; }

    # 2. Concatenar e alinhar com MAFFT
    cat "$CDS_A" "$CDS_B" > "${PAIR_DIR}/combined.fa"
    echo "  Alinhando com MAFFT..."
    mafft --auto --quiet "${PAIR_DIR}/combined.fa" > "${PAIR_DIR}/aligned.fa"

    # 3. Converter para AXT (formato do KaKs_Calculator)
    export PAIR_NAME_PY="${gene_a}-${gene_b}"
    export ALIGNED_FA_PY="${PAIR_DIR}/aligned.fa"
    export AXT_OUT_PY="${PAIR_DIR}/pair.axt"

    python3 << 'PYEOF'
import os
from Bio import SeqIO

pair_name   = os.environ["PAIR_NAME_PY"]
aligned_fa  = os.environ["ALIGNED_FA_PY"]
axt_out     = os.environ["AXT_OUT_PY"]

aln = list(SeqIO.parse(aligned_fa, "fasta"))
if len(aln) < 2:
    print("ERRO: alinhamento com < 2 sequências")
    exit(1)

with open(axt_out, "w") as f:
    f.write(f"{pair_name}\n")
    f.write(str(aln[0].seq).replace("-", "") + "\n")
    f.write(str(aln[1].seq).replace("-", "") + "\n")
    f.write("\n")
print(f"  AXT gerado: {axt_out}")
PYEOF

    # 4. KaKs_Calculator 2.0
    echo "  Calculando Ka/Ks (método NG)..."
    if command -v KaKs_Calculator &>/dev/null; then
        KaKs_Calculator \
            -i "${PAIR_DIR}/pair.axt" \
            -o "${PAIR_DIR}/kaks_results.txt" \
            -m NG 2>/dev/null
        echo "  Resultado:"
        cat "${PAIR_DIR}/kaks_results.txt"
    else
        echo "  AVISO: KaKs_Calculator não encontrado."
        echo "  Instale: mamba install -n kerson-paper -c bioconda kaks-calculator"
    fi

done < "$PAIRS_FILE"

# ── Compilar resultados ───────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "COMPILANDO RESULTADOS..."
echo "============================================================"

export RESULTS_DIR_PY="$OUTDIR"
export SUMMARY_OUT_PY="${SCRIPT_DIR}/kaks_summary.tsv"

python3 << 'PYEOF'
import os, glob

results_dir = os.environ["RESULTS_DIR_PY"]
out_summary = os.environ["SUMMARY_OUT_PY"]

header = "pair\tsequence\tmethod\tka\tks\tka_ks\tpvalue\tlen"
rows   = []

for f in sorted(glob.glob(os.path.join(results_dir, "**/kaks_results.txt"), recursive=True)):
    pair = os.path.basename(os.path.dirname(f))
    with open(f) as fp:
        lines = fp.readlines()
    for line in lines[1:]:
        parts = line.strip().split("\t")
        if len(parts) >= 7:
            rows.append(f"{pair}\t" + "\t".join(parts[:7]))

if rows:
    with open(out_summary, "w") as f:
        f.write(header + "\n")
        f.write("\n".join(rows) + "\n")
    print(f"Sumário salvo: {out_summary}")
    print(f"Total de pares analisados: {len(rows)}")
else:
    print("Nenhum resultado para compilar.")
PYEOF

echo ""
echo "Pipeline Ka/Ks concluído."
echo "Interpretação: Ka/Ks < 1 = purificante | ≈1 = neutro | > 1 = positiva"
