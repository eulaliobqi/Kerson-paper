#!/bin/bash
# 01_run_hmmer.sh
# Anotação de domínios Pfam para os 49 LRR-RLPs
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Entrada: proteins_49rlp.fa        (gerado por 00_fetch_proteins.sh)
# Saída:   hmmer_out/hmmer_domains.tsv  → input para 02_plot_domain_architecture.R
#
# Tempo estimado: 10–30 min (16 CPUs, Pfam-A completo ~21.000 domínios)
# Pré-requisito: bash 00_fetch_proteins.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFAM_DIR="/home/eulalio/databases/pfam"
PFAM_DB="${PFAM_DIR}/Pfam-A.hmm"
PROTEINS_FA="${SCRIPT_DIR}/proteins_49rlp.fa"
OUTDIR="${SCRIPT_DIR}/hmmer_out"
DOMTBL="${OUTDIR}/rlp_pfam_domtbl.txt"
DOMAINS_TSV="${OUTDIR}/hmmer_domains.tsv"

mkdir -p "$OUTDIR" "$PFAM_DIR"

echo "================================================================="
echo " 01_run_hmmer.sh — HMMER hmmscan vs Pfam-A"
echo " Proteinas: ${PROTEINS_FA}"
echo " Saida:     ${DOMAINS_TSV}"
echo "================================================================="
echo ""

# ── [1/4] Verificar sequências de entrada ────────────────────────────────────
if [ ! -f "$PROTEINS_FA" ]; then
    echo "ERRO: arquivo de proteinas nao encontrado: ${PROTEINS_FA}"
    echo "Execute primeiro:"
    echo "  bash ${SCRIPT_DIR}/00_fetch_proteins.sh"
    exit 1
fi

N_SEQS=$(grep -c "^>" "$PROTEINS_FA")
echo "[1/4] Sequencias de entrada: ${N_SEQS}"
if [ "$N_SEQS" -eq 0 ]; then
    echo "ERRO: ${PROTEINS_FA} esta vazio."
    exit 1
fi

# ── [2/4] Download + press de Pfam-A (se necessário) ─────────────────────────
# Verificar se Pfam-A.hmm existe E se os índices do hmmpress existem
PFAM_NEEDS_PRESS=0

if [ ! -f "$PFAM_DB" ]; then
    echo "[2/4] Baixando Pfam-A.hmm (~1.5 GB)..."

    # URL primária: EBI FTP (mais estável)
    PFAM_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    # URL de fallback: InterPro FTP (mesma release)
    PFAM_URL2="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/Pfam-A.hmm.gz"

    if wget -q --spider --timeout=15 "$PFAM_URL" 2>/dev/null; then
        wget -c -q --show-progress -O "${PFAM_DIR}/Pfam-A.hmm.gz" "$PFAM_URL"
    elif wget -q --spider --timeout=15 "$PFAM_URL2" 2>/dev/null; then
        wget -c -q --show-progress -O "${PFAM_DIR}/Pfam-A.hmm.gz" "$PFAM_URL2"
    else
        echo "ERRO: Pfam-A.hmm inacessivel via FTP."
        echo "Baixe manualmente:"
        echo "  wget -c -O ${PFAM_DIR}/Pfam-A.hmm.gz '${PFAM_URL}'"
        exit 1
    fi

    gunzip "${PFAM_DIR}/Pfam-A.hmm.gz"
    PFAM_NEEDS_PRESS=1
fi

# Verificar se índices do hmmpress existem (.h3i, .h3m, .h3f, .h3p)
# Se Pfam-A.hmm existe mas índices não, re-rodar hmmpress
if [ -f "$PFAM_DB" ]; then
    if [ ! -f "${PFAM_DB}.h3i" ] || [ ! -f "${PFAM_DB}.h3m" ]; then
        echo "  Indices hmmpress ausentes — executando hmmpress..."
        PFAM_NEEDS_PRESS=1
    fi
fi

if [ "$PFAM_NEEDS_PRESS" -eq 1 ]; then
    echo "[2/4] Indexando Pfam-A com hmmpress (pode levar 2-5 min)..."
    hmmpress "$PFAM_DB"
    echo "  hmmpress concluido"
else
    N_MODELS=$(grep -c "^NAME" "$PFAM_DB" 2>/dev/null || echo "desconhecido")
    echo "[2/4] Pfam-A disponivel: ~${N_MODELS} modelos HMM"
fi

# ── [3/4] HMMER hmmscan ───────────────────────────────────────────────────────
if [ -f "$DOMTBL" ] && [ -s "$DOMTBL" ]; then
    echo ""
    echo "[3/4] hmmscan já executado: ${DOMTBL}"
    echo "      Para re-executar: rm ${DOMTBL} && bash $0"
else
    echo ""
    echo "[3/4] Executando hmmscan (${N_SEQS} sequencias vs Pfam-A, 16 CPUs)..."
    echo "      Tempo estimado: 10–30 min dependendo do hardware"
    echo ""

    hmmscan \
        --domtblout "$DOMTBL" \
        --cpu 16 \
        -E 0.001 \
        --domE 0.001 \
        --noali \
        "$PFAM_DB" \
        "$PROTEINS_FA" \
        > "${OUTDIR}/hmmscan.log" 2>&1

    echo "hmmscan concluido."
    echo "Log: ${OUTDIR}/hmmscan.log"
fi

# Verificar saída do hmmscan
if [ ! -s "$DOMTBL" ]; then
    echo "ERRO: domtblout vazio — verifique ${OUTDIR}/hmmscan.log"
    exit 1
fi

N_HITS=$(grep -vc "^#" "$DOMTBL" 2>/dev/null || echo 0)
echo "  Hits encontrados no domtblout: ${N_HITS}"

# ── [4/4] Parsear domtblout → TSV limpo ──────────────────────────────────────
echo ""
echo "[4/4] Convertendo domtblout para TSV limpo..."

export DOMTBL_PY="$DOMTBL"
export DOMAINS_TSV_PY="$DOMAINS_TSV"
export PROTEINS_FA_PY="$PROTEINS_FA"

python3 - << 'PYEOF'
import os, sys
from collections import Counter

domtbl    = os.environ["DOMTBL_PY"]
out_tsv   = os.environ["DOMAINS_TSV_PY"]
prot_fa   = os.environ["PROTEINS_FA_PY"]

# Mapear ID da isoforma (Solyc01g005730.3.1) → gene_id (Solyc01g005730)
from Bio import SeqIO
import re
isoform_to_gene = {}
for rec in SeqIO.parse(prot_fa, "fasta"):
    # Extrair gene_id removendo sufixos .N.N do final do ID ITAG4.0
    m = re.match(r"(Solyc\d+g\d+)", rec.id)
    if m:
        isoform_to_gene[rec.id] = m.group(1)

rows = []
with open(domtbl) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split()
        if len(cols) < 23:
            continue
        # Formato domtblout HMMER3:
        # col 0: target name (domínio Pfam)
        # col 3: query name  (proteína)
        # col 5: query len
        # col 11: i-evalue (domain)
        # col 12: score (domain)
        # col 17: hmm from, 18: hmm to, 19: ali from, 20: ali to
        target_name = cols[0]    # Pfam domain name (ex: LRR_1)
        query_name  = cols[3]    # proteína (ex: Solyc01g005730.3.1)
        query_len   = int(cols[5])
        ievalue     = float(cols[11])
        score       = float(cols[12])
        hmm_from    = int(cols[15])
        hmm_to      = int(cols[16])
        ali_from    = int(cols[17])
        ali_to      = int(cols[18])
        env_from    = int(cols[19])
        env_to      = int(cols[20])
        acc         = cols[1]    # Pfam accession (ex: PF00560.33)
        description = " ".join(cols[22:])

        # Mapear isoforma para gene_id
        gene_id = isoform_to_gene.get(query_name, query_name)

        rows.append((
            gene_id, query_name, query_len, target_name, acc,
            ali_from, ali_to, env_from, env_to,
            ievalue, score, description
        ))

# Ordenar por gene_id depois por posição no alinhamento
rows.sort(key=lambda x: (x[0], x[5]))

header = (
    "gene_id\tisoform_id\tprotein_len\tdomain\tpfam_acc"
    "\tali_start\tali_end\tenv_start\tenv_end\tievalue\tscore\tdescription"
)
with open(out_tsv, "w") as fh:
    fh.write(header + "\n")
    for r in rows:
        fh.write("\t".join(str(v) for v in r) + "\n")

n_genes   = len(set(r[0] for r in rows))
n_domains = len(rows)
print(f"Dominios encontrados: {n_domains}")
print(f"Genes com dominios: {n_genes}")
print(f"TSV salvo: {out_tsv}")

domain_counts = Counter(r[3] for r in rows)
print("\nTop 15 dominios Pfam mais frequentes:")
for dom, cnt in domain_counts.most_common(15):
    print(f"  {dom}: {cnt}")
PYEOF

echo ""
echo "================================================================="
echo " CONCLUIDO: ${DOMAINS_TSV}"
echo "================================================================="
echo ""
echo "PROXIMO PASSO:"
echo "  Rscript ${SCRIPT_DIR}/02_plot_domain_architecture.R ${DOMAINS_TSV}"
