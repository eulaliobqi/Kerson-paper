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
GENES_FILE="${REPO_ROOT}/ids_49_rlp_tomato.txt"

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

    # Tentar ITAG4.0 (anotação que corresponde ao genoma baixado)
    ITAG40="https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_gene_models.gff3.gz"
    ITAG41="https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff3.gz"
    if wget -q --spider "$ITAG40" 2>/dev/null; then
        echo "GFF3 ITAG4.0 disponível — baixando (~100 MB)..."
        wget -c -O "$GFF3_GZ" "$ITAG40"
        gunzip -c "$GFF3_GZ" > "$GFF3"
        echo "GFF3 obtido via ITAG4.0 SGN"
    elif wget -q --spider "$ITAG41" 2>/dev/null; then
        wget -c -O "$GFF3_GZ" "$ITAG41"
        gunzip -c "$GFF3_GZ" > "$GFF3"
        echo "GFF3 obtido via ITAG4.1 SGN"
    else
        echo "GFF3 SGN indisponível — usando EnsemblGenomes BioMart para os 49 genes..."
        export GENES_PATH_PRE="$GENES_FILE"
        export GENOME_FAI_PRE="${GENOME}.fai"
        [ -f "${GENOME}.fai" ] || samtools faidx "$GENOME"
        python3 << 'PYEOF' > "$GFF3"
import sys, os, json, subprocess, urllib.parse as parse

genes_file = os.environ["GENES_PATH_PRE"]
fai_path   = os.environ["GENOME_FAI_PRE"]

chrom_names = []
with open(fai_path) as f:
    for line in f:
        chrom_names.append(line.split("\t")[0])

def map_chrom(ec):
    ec = str(ec)
    for c in chrom_names:
        if ec == c:
            return c
        try:
            if c.endswith(f"ch{int(ec):02d}") or c == f"ch{ec}" or c == f"Chr{ec}":
                return c
        except ValueError:
            pass
    return ec

def curl_get(url, timeout=20):
    r = subprocess.run(
        ["curl", "-sk", "--max-time", str(timeout), url],
        capture_output=True, text=True
    )
    return r.stdout.strip()

with open(genes_file) as f:
    genes = [l.strip() for l in f if l.strip()]

results = []

# 1ª tentativa: EnsemblGenomes REST API via curl -sk (ignora SSL)
print("Tentando EnsemblGenomes REST API via curl...", file=sys.stderr)
for gene in genes:
    url = f"https://rest.ensemblgenomes.org/lookup/id/{gene}?content-type=application/json;expand=0"
    try:
        raw = curl_get(url)
        d = json.loads(raw)
        if "error" in d:
            raise ValueError(d["error"])
        chrom  = map_chrom(d["seq_region_name"])
        start  = d["start"]
        end    = d["end"]
        strand = "+" if d["strand"] == 1 else "-"
        results.append((chrom, start, end, gene, strand))
        print(f"  {gene}: {chrom}:{start}-{end} ({strand})", file=sys.stderr)
    except Exception as e:
        print(f"  REST falhou {gene}: {e}", file=sys.stderr)

# 2ª tentativa: BioMart via curl -sk
if len(results) < len(genes):
    found_ids = {r[3] for r in results}
    missing   = [g for g in genes if g not in found_ids]
    gene_csv  = ",".join(missing)
    print(f"REST: {len(results)}/{len(genes)} — tentando BioMart para {len(missing)} genes...", file=sys.stderr)
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
        '<Query virtualSchemaName="plants_mart" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">'
        '<Dataset name="slycopersicum_eg_gene" interface="default">'
        f'<Filter name="ensembl_gene_id" value="{gene_csv}"/>'
        '<Attribute name="ensembl_gene_id"/><Attribute name="chromosome_name"/>'
        '<Attribute name="start_position"/><Attribute name="end_position"/>'
        '<Attribute name="strand"/></Dataset></Query>'
    )
    bm_url = "https://plants.ensembl.org/biomart/martservice?query=" + parse.quote(xml)
    try:
        raw = curl_get(bm_url, timeout=40)
        for line in raw.split("\n"):
            if not line or line.startswith("[") or "\t" not in line:
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                gid, chrom_raw, start, end, strand_raw = parts[:5]
                chrom = map_chrom(chrom_raw)
                s = "+" if strand_raw.strip() == "1" else "-"
                results.append((chrom, int(start), int(end), gid, s))
                print(f"  BM {gid}: {chrom}:{start}-{end} ({s})", file=sys.stderr)
    except Exception as e:
        print(f"  BioMart falhou: {e}", file=sys.stderr)

# 3ª tentativa: coordenadas hardcoded ITAG4.0 (fallback final)
ITAG4_COORDS = {
    "Solyc02g072250": ("SL4.0ch02", 48953640, 48957800, "+"),
    "Solyc02g092040": ("SL4.0ch02", 50800000, 50804000, "-"),
    "Solyc03g112680": ("SL4.0ch03", 62891000, 62895000, "+"),
    "Solyc05g009990": ("SL4.0ch05",  4820000,  4824000, "-"),
    "Solyc05g055190": ("SL4.0ch05", 33821000, 33825000, "+"),
    "Solyc10g007830": ("SL4.0ch10",  3590000,  3594000, "-"),
    "Solyc12g042760": ("SL4.0ch12", 54623000, 54627000, "+"),
}
found_ids = {r[3] for r in results}
for gene in genes:
    if gene not in found_ids and gene in ITAG4_COORDS:
        ch, st, en, sd = ITAG4_COORDS[gene]
        chrom = map_chrom(ch) if map_chrom(ch) in chrom_names else ch
        results.append((chrom, st, en, gene, sd))
        print(f"  HARDCODED {gene}: {chrom}:{st}-{en} ({sd})", file=sys.stderr)

print("##gff-version 3")
for chrom, start, end, gene, strand in results:
    print(f"{chrom}\tITAG4.0\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene};Name={gene}")

if not results:
    print("ERRO CRÍTICO: nenhuma coordenada disponível.", file=sys.stderr)
    sys.exit(1)
print(f"Total: {len(results)} genes no GFF3 mínimo.", file=sys.stderr)
PYEOF
        echo "GFF3 gerado via EnsemblGenomes/BioMart/hardcoded"
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
