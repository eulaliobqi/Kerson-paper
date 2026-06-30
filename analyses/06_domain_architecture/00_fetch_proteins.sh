#!/bin/bash
# 00_fetch_proteins.sh
# Extrai sequências proteicas dos 49 LRR-RLPs do ITAG4.0
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# NOTA (2026-06-30): SGN FTP e EnsemblPlants FTP são INACESSÍVEIS no servidor UFV.
# Estratégia atualizada (em ordem de preferência):
#   1. NCBI FTP  — ftp.ncbi.nlm.nih.gov (GCF_000188115.5 = SL4.0, proteínas RefSeq)
#   2. SGN FTP   — ftp.solgenomics.net (ITAG4.0, fonte oficial — inacessível)
#   3. EnsemblPlants FTP — ftp.ensemblgenomes.ebi.ac.uk (inacessível)
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

echo "================================================================="
echo " 00_fetch_proteins.sh — Proteoma ITAG4.0 → 49 LRR-RLPs"
echo " IDs: ${IDS_FILE}"
echo " Saída: ${OUT_FA}"
echo "================================================================="
echo ""

# ── Download do proteoma ITAG4.0 ──────────────────────────────────────────────
if [ ! -f "$PROTEOME" ]; then
    echo "[1/3] Baixando proteoma ITAG4.0..."

    # URL 1: NCBI FTP (GCF_000188115.5 = SL4.0, RefSeq proteins — acessível no servidor UFV)
    NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL4.0/GCF_000188115.5_SL4.0_protein.faa.gz"

    # URL 2: SGN FTP (ITAG4.0 canônico — inacessível no servidor UFV mas mantido como fallback)
    SGN_URL="https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_proteins.fasta.gz"

    # URL 3: EnsemblPlants FTP (inacessível no servidor UFV mas mantido como fallback)
    ENSEMBL_BASE="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/solanum_lycopersicum/pep/"

    echo "  Testando URL primária: NCBI FTP (GCF_000188115.5)..."
    if wget -q --spider --timeout=20 "$NCBI_URL" 2>/dev/null; then
        echo "  NCBI FTP acessível — baixando (~40 MB)..."
        wget -c -q --show-progress -O "${DB_DIR}/ncbi_proteins.faa.gz" "$NCBI_URL"
        gunzip -c "${DB_DIR}/ncbi_proteins.faa.gz" > "$PROTEOME"
        rm -f "${DB_DIR}/ncbi_proteins.faa.gz"
        echo "  Proteoma obtido via NCBI FTP (RefSeq GCF_000188115.5)"
    else
        echo "  NCBI FTP inacessível — tentando SGN FTP..."
        if wget -q --spider --timeout=15 "$SGN_URL" 2>/dev/null; then
            echo "  SGN FTP acessível — baixando (~50 MB)..."
            wget -c -q --show-progress -O "${DB_DIR}/ITAG4.0_proteins.fasta.gz" "$SGN_URL"
            gunzip -c "${DB_DIR}/ITAG4.0_proteins.fasta.gz" > "$PROTEOME"
            rm -f "${DB_DIR}/ITAG4.0_proteins.fasta.gz"
            echo "  Proteoma obtido via SGN"
        else
            echo "  SGN FTP inacessível — tentando EnsemblPlants FTP..."
            LISTING=$(wget -q --timeout=20 -O- "$ENSEMBL_BASE" 2>/dev/null || true)
            ENSEMBL_FILE=""
            if [ -n "$LISTING" ]; then
                ENSEMBL_FILE=$(echo "$LISTING" | grep -oP 'Solanum_lycopersicum\.SL4\.0\.\d+\.pep\.all\.fa\.gz' | head -1 || true)
            fi

            if [ -n "$ENSEMBL_FILE" ]; then
                ENSEMBL_URL="${ENSEMBL_BASE}${ENSEMBL_FILE}"
                echo "  Baixando ${ENSEMBL_URL}..."
                wget -c -q --show-progress -O "${DB_DIR}/ensembl_proteins.fa.gz" "$ENSEMBL_URL"
                gunzip -c "${DB_DIR}/ensembl_proteins.fa.gz" > "$PROTEOME"
                rm -f "${DB_DIR}/ensembl_proteins.fa.gz"
                echo "  Proteoma obtido via EnsemblPlants"
            else
                echo ""
                echo "  ERRO: Nenhuma fonte de proteoma acessível."
                echo "  Baixe manualmente no NOTEBOOK WINDOWS e copie via SCP:"
                echo "    wget -O proteins.faa.gz '${NCBI_URL}'"
                echo "    gunzip proteins.faa.gz"
                echo "    scp proteins.faa eulalio@200.235.143.10:${PROTEOME}"
                exit 1
            fi
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
