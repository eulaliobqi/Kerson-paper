#!/bin/bash
# Extrai sequências 2 kb upstream dos 7 genes LRR-RLP do tomate (ITAG4.0)
# Executar no servidor: eulalio@200.235.143.10
# Dependências: wget, samtools, bedtools, python3 (biopython)
# Ambiente: mamba activate kerson-paper

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DB_DIR="/home/eulalio/databases/itag4.0"
GENOME="${DB_DIR}/S_lycopersicum_chromosomes.4.00.fa"
GFF3="${DB_DIR}/ITAG4.0_gene_models.gff3"
UPSTREAM=2000
OUT_FA="${SCRIPT_DIR}/rlp_upstream_2kb.fa"
GENES_FILE="${REPO_ROOT}/analyses/genes_7rlp.txt"

mkdir -p "$DB_DIR"

# ── Download genoma ITAG4.0 (se ausente) ──────────────────────────────────────
if [ ! -f "$GENOME" ]; then
    echo "[1/5] Baixando genoma ITAG4.0..."
    wget -c -P "$DB_DIR" \
        "https://ftp.solgenomics.net/tomato_genome/assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa.gz"
    gunzip "${DB_DIR}/S_lycopersicum_chromosomes.4.00.fa.gz"
fi

if [ ! -f "$GFF3" ]; then
    echo "[2/5] Baixando anotação GFF3 (ITAG4.1 → fallback REST API)..."
    GFF3_GZ="${DB_DIR}/gene_models.gff3.gz"

    # Tentar ITAG4.1 (mesmo IDs Solyc, anotação mais recente)
    ITAG41="https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff3.gz"
    if wget -q --spider "$ITAG41" 2>/dev/null; then
        wget -c -O "$GFF3_GZ" "$ITAG41"
        gunzip -c "$GFF3_GZ" > "$GFF3"
        echo "GFF3 obtido via ITAG4.1 SGN"
    else
        echo "ITAG4.1 indisponível — usando Ensembl Plants REST API para os 7 genes..."
        export GENES_PATH_PRE="$GENES_FILE"
        export GENOME_FAI_PRE="${GENOME}.fai"
        # Criar FAI antes de precisar do GFF3
        [ -f "${GENOME}.fai" ] || samtools faidx "$GENOME"
        python3 << 'PYEOF' > "$GFF3"
import sys, os, json, urllib.request as req

genes_file = os.environ["GENES_PATH_PRE"]
fai        = os.environ["GENOME_FAI_PRE"]

# Detectar nomenclatura de cromossomos no FASTA
chrom_names = []
with open(fai) as f:
    for line in f:
        chrom_names.append(line.split("\t")[0])

# Mapear inteiro Ensembl → nome real no FASTA
def map_chrom(ec):
    for c in chrom_names:
        if ec == c or c.endswith(f"ch{int(ec):02d}") or c == f"ch{ec}" or c == f"Chr{ec}":
            return c
    return ec  # retornar como está se não encontrar

with open(genes_file) as f:
    genes = [l.strip() for l in f if l.strip()]

print("##gff-version 3")
for gene in genes:
    url = f"https://rest.ensembl.org/lookup/id/{gene}?content-type=application/json;expand=0"
    try:
        with req.urlopen(url, timeout=15) as r:
            d = json.loads(r.read())
        chrom  = map_chrom(str(d["seq_region_name"]))
        start  = d["start"]
        end    = d["end"]
        strand = "+" if d["strand"] == 1 else "-"
        attrs  = f"ID={gene};Name={gene}"
        print(f"{chrom}\tEnsemblPlants\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")
        print(f"  {gene}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
    except Exception as e:
        print(f"  AVISO: {gene} — {e}", file=sys.stderr)
PYEOF
        echo "GFF3 mínimo gerado via Ensembl Plants REST API"
    fi
fi

# ── Indexar genoma ────────────────────────────────────────────────────────────
[ -f "${GENOME}.fai" ] || samtools faidx "$GENOME"
cut -f1,2 "${GENOME}.fai" > /tmp/chrom.sizes

echo "[3/5] Extraindo coordenadas dos genes do GFF3..."

# FIX: exportar variáveis ANTES do heredoc Python que as usa
export GFF3_PATH="$GFF3"
export GENES_PATH="$GENES_FILE"

python3 << 'PYEOF' > /tmp/rlp_genes.bed
import sys, os

genes_file = os.environ["GENES_PATH"]
gff3_path  = os.environ["GFF3_PATH"]

with open(genes_file) as f:
    genes = [line.strip() for line in f if line.strip()]

found = set()
with open(gff3_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 9 or cols[2] != "gene":
            continue
        attrs = cols[8]
        for gene in genes:
            if gene in attrs and gene not in found:
                chrom  = cols[0]
                start  = int(cols[3]) - 1  # BED é 0-based
                end    = int(cols[4])
                strand = cols[6]
                print(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}")
                found.add(gene)
                break

missing = set(genes) - found
if missing:
    print(f"AVISO: genes não encontrados: {missing}", file=sys.stderr)
PYEOF

echo "Genes encontrados:"
awk '{print $4, $1, $2, $3, $6}' /tmp/rlp_genes.bed

echo ""
echo "[4/5] Calculando regiões upstream com bedtools flank..."
bedtools flank \
    -i /tmp/rlp_genes.bed \
    -g /tmp/chrom.sizes \
    -l ${UPSTREAM} \
    -r 0 \
    -s \
    > /tmp/rlp_upstream.bed

echo "[5/5] Extraindo sequências FASTA..."
bedtools getfasta \
    -fi "$GENOME" \
    -bed /tmp/rlp_upstream.bed \
    -name \
    -s \
    > "$OUT_FA"

echo ""
echo "=========================================================="
echo "CONCLUÍDO: ${OUT_FA}"
echo "=========================================================="
echo ""
echo "PRÓXIMOS PASSOS:"
echo "1. Acesse: https://bioinformatics.psb.ugent.be/webtools/plantcare/html/"
echo "2. Selecione 'Search Promoter' → cole as sequências de ${OUT_FA}"
echo "3. Baixe os resultados como TXT (opção 'Download results')"
echo "4. Salve como: ${SCRIPT_DIR}/plantcare_results.txt"
echo "5. Execute: Rscript ${SCRIPT_DIR}/02_plot_plantcare_heatmap.R"
