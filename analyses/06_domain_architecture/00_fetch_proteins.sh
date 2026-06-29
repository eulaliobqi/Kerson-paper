#!/bin/bash
# Extrai sequências proteicas dos 49 LRR-RLPs do ITAG4.0
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Entrada:  ../../ids_49_rlp_tomato.txt  (49 IDs)
# Saída:    proteins_49rlp.fa            (input para 01_run_hmmer.sh)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
IDS_FILE="${REPO_ROOT}/ids_49_rlp_tomato.txt"
DB_DIR="/home/eulalio/databases/itag4.0"
PROTEOME="${DB_DIR}/ITAG4.0_proteins.fa"
OUT_FA="${SCRIPT_DIR}/proteins_49rlp.fa"

mkdir -p "$DB_DIR"

# ── Download proteoma ITAG4.0 (se ausente) ────────────────────────────────────
if [ ! -f "$PROTEOME" ]; then
    echo "Baixando proteoma ITAG4.0..."
    wget -c -P "$DB_DIR" \
        "https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_proteins.fasta.gz"
    gunzip "${DB_DIR}/ITAG4.0_proteins.fasta.gz"
    mv "${DB_DIR}/ITAG4.0_proteins.fasta" "$PROTEOME"
fi

N_TOTAL=$(grep -c "^>" "$PROTEOME")
echo "Proteoma carregado: ${N_TOTAL} sequências"
echo "Filtrando pelos 49 IDs em: ${IDS_FILE}"

# ── Filtrar sequências pelos 49 IDs ──────────────────────────────────────────
export IDS_FILE_PY="$IDS_FILE"
export PROTEOME_PY="$PROTEOME"
export OUT_FA_PY="$OUT_FA"

python3 << 'PYEOF'
import os, sys
from Bio import SeqIO

ids_file   = os.environ["IDS_FILE_PY"]
proteome   = os.environ["PROTEOME_PY"]
out_fa     = os.environ["OUT_FA_PY"]

# Carregar IDs dos 49 RLPs
with open(ids_file) as f:
    target_ids = {line.strip() for line in f if line.strip()}

print(f"IDs alvo: {len(target_ids)}")

# Buscar no proteoma — pegar PRIMEIRA isoforma de cada gene
found     = {}
not_found = set(target_ids)

for rec in SeqIO.parse(proteome, "fasta"):
    for gene_id in list(not_found):
        if gene_id in rec.id:
            found[gene_id] = rec
            not_found.discard(gene_id)
            break

print(f"Encontrados: {len(found)} / {len(target_ids)}")

if not_found:
    print(f"AVISO: {len(not_found)} IDs não encontrados:", file=sys.stderr)
    for gid in sorted(not_found):
        print(f"  {gid}", file=sys.stderr)

# Escrever FASTA filtrado
records = [found[gid] for gid in sorted(found)]
with open(out_fa, "w") as f:
    SeqIO.write(records, f, "fasta")

print(f"FASTA salvo: {out_fa} ({len(records)} sequências)")
PYEOF

echo ""
echo "Próximo passo: bash 01_run_hmmer.sh"
