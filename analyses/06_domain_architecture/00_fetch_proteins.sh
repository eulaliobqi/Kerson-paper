#!/bin/bash
# 00_fetch_proteins.sh
# Extrai sequências proteicas dos 49 LRR-RLPs do ITAG4.0
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Estratégia (em ordem de preferência):
#   1. Arquivo local no repo (analyses/06_domain_architecture/ITAG4.0_proteins.fasta)
#      — commitado no git, obtido via `git pull` sem depender de rede externa
#   2. SGN FTP (ftp.solgenomics.net — inacessível no servidor UFV)
#   3. Cópia manual via SCP do notebook Windows
#
# Entrada:  ../../ids_49_rlp_tomato.txt  (49 IDs)
# Saída:    proteins_49rlp.fa            (input para 01_run_hmmer.sh)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
IDS_FILE="${REPO_ROOT}/ids_49_rlp_tomato.txt"
DB_DIR="/home/eulalio/databases/itag4.0"
PROTEOME="${DB_DIR}/ITAG4.0_proteins.fa"
REPO_PROTEOME="${SCRIPT_DIR}/ITAG4.0_proteins.fasta"
OUT_FA="${SCRIPT_DIR}/proteins_49rlp.fa"

mkdir -p "$DB_DIR"

echo "================================================================="
echo " 00_fetch_proteins.sh — Proteoma ITAG4.0 → 49 LRR-RLPs"
echo " IDs: ${IDS_FILE}"
echo " Saída: ${OUT_FA}"
echo "================================================================="
echo ""

# ── [1/3] Proteoma ITAG4.0 ────────────────────────────────────────────────────
if [ ! -f "$PROTEOME" ]; then
    echo "[1/3] Configurando proteoma ITAG4.0..."

    # Fonte 1: arquivo local no repositório (commitado no Windows, obtido via git pull)
    if [ -f "$REPO_PROTEOME" ]; then
        echo "  Usando cópia local do repo (git pull incluiu ITAG4.0_proteins.fasta)..."
        cp "$REPO_PROTEOME" "$PROTEOME"
        echo "  Proteoma copiado de ${REPO_PROTEOME}"

    # Fonte 2: SGN FTP — URL sem .gz (arquivo descomprimido, 12 MB)
    else
        SGN_URL="https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_proteins.fasta"
        echo "  Arquivo local ausente — tentando SGN FTP (pode falhar no servidor UFV)..."
        if wget -q --spider --timeout=15 "$SGN_URL" 2>/dev/null; then
            echo "  SGN acessível — baixando (~12 MB)..."
            wget -c -q --show-progress -O "$PROTEOME" "$SGN_URL"
            echo "  Proteoma obtido via SGN FTP"
        else
            echo ""
            echo "  ERRO: Proteoma indisponível."
            echo "  No NOTEBOOK WINDOWS, já está em:"
            echo "    C:\\Users\\eulal\\kerson-paper\\analyses\\06_domain_architecture\\ITAG4.0_proteins.fasta"
            echo "  Copie via SCP:"
            echo "    scp C:\\Users\\eulal\\kerson-paper\\analyses\\06_domain_architecture\\ITAG4.0_proteins.fasta \\"
            echo "        eulalio@200.235.143.10:${PROTEOME}"
            exit 1
        fi
    fi
else
    echo "[1/3] Proteoma já disponível: ${PROTEOME}"
fi

# Verificar integridade
N_TOTAL=$(grep -c "^>" "$PROTEOME" 2>/dev/null || echo 0)
if [ "$N_TOTAL" -eq 0 ]; then
    echo "ERRO: proteoma vazio ou corrompido: ${PROTEOME}"
    echo "Remova e re-execute: rm -f ${PROTEOME}"
    exit 1
fi
echo "  Proteoma carregado: ${N_TOTAL} sequências"

# ── Filtrar os 49 IDs ─────────────────────────────────────────────────────────
echo ""
echo "[2/3] Filtrando pelos 49 IDs de ${IDS_FILE}..."
echo "       (ITAG4.0 usa sufixo .1 ou .2 nas isoformas — pegando a 1ª por gene)"

export IDS_FILE_PY="$IDS_FILE"
export PROTEOME_PY="$PROTEOME"
export OUT_FA_PY="$OUT_FA"

python3 - << 'PYEOF'
import os, sys
from Bio import SeqIO

ids_file = os.environ["IDS_FILE_PY"]
proteome = os.environ["PROTEOME_PY"]
out_fa   = os.environ["OUT_FA_PY"]

# Carregar IDs alvo
with open(ids_file) as f:
    target_ids = [l.strip() for l in f if l.strip()]
target_set = set(target_ids)

print(f"IDs alvo: {len(target_ids)}")

# Buscar no proteoma: pegar a PRIMEIRA isoforma de cada gene
# ITAG4.0 nomeia isoformas como Solyc01g005730.3.1, Solyc01g005730.3.2 ...
found     = {}
not_found = set(target_ids)

for rec in SeqIO.parse(proteome, "fasta"):
    if not not_found:
        break
    for gid in list(not_found):
        # Busca exata: gene_id deve aparecer como prefixo do ID da sequência
        # Ex: rec.id = "Solyc01g005730.3.1" → contém "Solyc01g005730"
        if rec.id.startswith(gid) or f"|{gid}" in rec.id or gid in rec.id:
            found[gid] = rec
            not_found.discard(gid)
            break

print(f"Encontrados: {len(found)} / {len(target_ids)}")

if not_found:
    print(f"AVISO: {len(not_found)} genes sem sequência no proteoma:", file=sys.stderr)
    for gid in sorted(not_found):
        print(f"  {gid}", file=sys.stderr)

if not found:
    print("ERRO CRITICO: nenhuma proteína encontrada.", file=sys.stderr)
    sys.exit(1)

# Escrever FASTA na ordem do IDs_FILE (mantém ordem reproducível)
records = [found[gid] for gid in target_ids if gid in found]
with open(out_fa, "w") as fh:
    SeqIO.write(records, fh, "fasta")

print(f"FASTA salvo: {out_fa} ({len(records)} sequencias)")
PYEOF

# Verificação final
N_OUT=$(grep -c "^>" "$OUT_FA" 2>/dev/null || echo 0)

echo ""
echo "[3/3] Verificacao da saida..."
echo "  Proteinas extraidas: ${N_OUT} / 49"

if [ "$N_OUT" -lt 45 ]; then
    echo "AVISO: menos de 45 proteinas extraidas. Verifique o proteoma e os IDs."
fi

echo ""
echo "================================================================="
echo " CONCLUIDO: ${OUT_FA}"
echo " ${N_OUT} sequencias prontas para HMMER"
echo "================================================================="
echo ""
echo "Proximo passo:"
echo "  bash ${SCRIPT_DIR}/01_run_hmmer.sh"
