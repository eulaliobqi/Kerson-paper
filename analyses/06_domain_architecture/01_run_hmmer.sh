#!/bin/bash
# Anotação de domínios Pfam para os 49 LRR-RLPs
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Entrada: proteins_49rlp.fa  (gerado por 00_fetch_proteins.sh)
# Saída:   hmmer_out/hmmer_domains.tsv  → input para plot_domain_architecture.R
#
# Pré-requisito: bash 00_fetch_proteins.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFAM_DB="/home/eulalio/databases/pfam/Pfam-A.hmm"
PROTEINS_FA="${SCRIPT_DIR}/proteins_49rlp.fa"
OUTDIR="${SCRIPT_DIR}/hmmer_out"

mkdir -p "$OUTDIR"

# ── Download Pfam-A (se ausente) ──────────────────────────────────────────────
if [ ! -f "$PFAM_DB" ]; then
    echo "Baixando Pfam-A.hmm (~1.5 GB)..."
    mkdir -p /home/eulalio/databases/pfam
    wget -c -P /home/eulalio/databases/pfam \
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip /home/eulalio/databases/pfam/Pfam-A.hmm.gz
    echo "Indexando Pfam-A..."
    hmmpress "$PFAM_DB"
fi

# ── Verificar sequências de entrada ──────────────────────────────────────────
if [ ! -f "$PROTEINS_FA" ]; then
    echo "Arquivo de proteínas não encontrado: ${PROTEINS_FA}"
    echo "Execute primeiro: bash 00_fetch_proteins.sh"
    exit 1
fi

N_SEQS=$(grep -c "^>" "$PROTEINS_FA")
echo "Sequências encontradas: ${N_SEQS}"

# ── HMMER hmmscan vs Pfam-A ───────────────────────────────────────────────────
echo "Rodando hmmscan (${N_SEQS} sequências vs Pfam-A, 16 CPUs)..."
hmmscan \
    --domtblout "${OUTDIR}/rlp_pfam_domtbl.txt" \
    --cpu 16 \
    -E 0.001 \
    --domE 0.001 \
    --noali \
    "$PFAM_DB" \
    "$PROTEINS_FA" \
    > "${OUTDIR}/hmmscan.log" 2>&1

echo "hmmscan concluído. Log: ${OUTDIR}/hmmscan.log"

# ── Parsear domtblout → TSV limpo ─────────────────────────────────────────────
# FIX: exportar OUTDIR antes do heredoc Python que o usa
export OUTDIR_PY="$OUTDIR"

python3 << 'PYEOF'
import os
from collections import Counter

outdir  = os.environ["OUTDIR_PY"]
domtbl  = os.path.join(outdir, "rlp_pfam_domtbl.txt")
out_tsv = os.path.join(outdir, "hmmer_domains.tsv")

rows = []
with open(domtbl) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split()
        if len(cols) < 23:
            continue
        target_name = cols[0]   # domínio Pfam
        query_name  = cols[3]   # proteína
        query_len   = int(cols[5])
        evalue      = float(cols[12])
        score       = float(cols[13])
        ali_from    = int(cols[17])
        ali_to      = int(cols[18])
        description = " ".join(cols[22:])
        rows.append((query_name, query_len, target_name, ali_from, ali_to,
                     evalue, score, description))

rows.sort(key=lambda x: (x[0], x[3]))

header = "gene_id\tprotein_len\tdomain\tstart\tend\tevalue\tscore\tdescription"
with open(out_tsv, "w") as f:
    f.write(header + "\n")
    for r in rows:
        f.write("\t".join(str(v) for v in r) + "\n")

print(f"Domínios encontrados: {len(rows)}")
print(f"TSV salvo: {out_tsv}")

domain_counts = Counter(r[2] for r in rows)
print("\nTop 15 domínios Pfam mais frequentes:")
for dom, cnt in domain_counts.most_common(15):
    print(f"  {dom}: {cnt}")
PYEOF

echo ""
echo "PRÓXIMO PASSO:"
echo "  Rscript ${SCRIPT_DIR}/02_plot_domain_architecture.R ${OUTDIR}/hmmer_domains.tsv"
