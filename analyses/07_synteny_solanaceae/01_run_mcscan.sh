#!/bin/bash
# Sintenia interespécie: tomate × batata × pimenta (LAST + MCScanX)
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Dependências: wget, last, mcscanx, python3
#   mamba install -n kerson-paper -c bioconda last mcscanx -y

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DB_DIR="/home/eulalio/databases/solanaceae_synteny"
OUTDIR="${SCRIPT_DIR}/mcscan_results"
GENES_FILE="${REPO_ROOT}/analyses/genes_7rlp.txt"

mkdir -p "$DB_DIR" "$OUTDIR"

# ── URLs dos genomas Solanaceae ───────────────────────────────────────────────
TOMATO_GFF="https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_gene_models.gff3.gz"
TOMATO_PEP="https://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_proteins.fasta.gz"
POTATO_GFF="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/226/075/GCA_000226075.1_DM_v6.1/GCA_000226075.1_DM_v6.1_genomic.gff.gz"
POTATO_PEP="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/226/075/GCA_000226075.1_DM_v6.1/GCA_000226075.1_DM_v6.1_protein.faa.gz"
PEPPER_GFF="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/710/875/GCA_000710875.1_Pepper_1.55/GCA_000710875.1_Pepper_1.55_genomic.gff.gz"
PEPPER_PEP="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/710/875/GCA_000710875.1_Pepper_1.55/GCA_000710875.1_Pepper_1.55_protein.faa.gz"

# ── 1. Download anotações ──────────────────────────────────────────────────────
for pair in "tomato|${TOMATO_GFF}" "tomato_pep|${TOMATO_PEP}" \
            "potato|${POTATO_GFF}" "potato_pep|${POTATO_PEP}" \
            "pepper|${PEPPER_GFF}" "pepper_pep|${PEPPER_PEP}"; do
    name="${pair%%|*}"
    url="${pair##*|}"
    ext="${url##*.gz}"
    out="${DB_DIR}/${name}.${url##*.}"
    out="${DB_DIR}/${name}.gz"
    final="${DB_DIR}/${name}"
    [ -f "$final" ] && continue
    echo "Baixando ${name}..."
    wget -q -c -O "$out" "$url"
    gunzip -c "$out" > "$final"
done

# ── 2. GFF → BED para MCScanX ─────────────────────────────────────────────────
export DB_DIR_PY="$DB_DIR"

python3 << 'PYEOF'
import os, re

db = os.environ["DB_DIR_PY"]
species_gff = {
    "tomato": os.path.join(db, "tomato"),
    "potato": os.path.join(db, "potato"),
    "pepper": os.path.join(db, "pepper"),
}

for sp, gff in species_gff.items():
    if not os.path.exists(gff):
        print(f"AVISO: {gff} não encontrado — pulando {sp}")
        continue
    bed_out = gff + ".bed"
    rows = []
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            chrom  = cols[0]
            start  = int(cols[3]) - 1
            end    = int(cols[4])
            strand = cols[6]
            m = re.search(r"ID=([^;]+)", cols[8])
            gene_id = m.group(1) if m else "."
            rows.append(f"{sp}_{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}")
    with open(bed_out, "w") as f:
        f.write("\n".join(rows) + "\n")
    print(f"{sp}: {len(rows)} genes → {bed_out}")
PYEOF

# ── 3. Comparação proteínas com LAST ─────────────────────────────────────────
for sp1 in tomato potato pepper; do
    for sp2 in tomato potato pepper; do
        [ "$sp1" = "$sp2" ] && continue
        BLAST_OUT="${OUTDIR}/${sp1}_vs_${sp2}.blast"
        [ -f "$BLAST_OUT" ] && continue
        echo "Comparando proteínas: ${sp1} × ${sp2}..."
        lastdb -p "${OUTDIR}/${sp2}_db" "${DB_DIR}/${sp2}_pep"
        lastal -f BlastTab "${OUTDIR}/${sp2}_db" "${DB_DIR}/${sp1}_pep" > "$BLAST_OUT"
    done
done

# ── 4. Concatenar e rodar MCScanX ─────────────────────────────────────────────
cat "${DB_DIR}/tomato.bed" "${DB_DIR}/potato.bed" "${DB_DIR}/pepper.bed" \
    > "${OUTDIR}/all_species.gff" 2>/dev/null || true

cat "${OUTDIR}/tomato_vs_potato.blast" \
    "${OUTDIR}/tomato_vs_pepper.blast" \
    "${OUTDIR}/potato_vs_pepper.blast" \
    > "${OUTDIR}/all_species.blast" 2>/dev/null || true

if command -v MCScanX &>/dev/null; then
    echo "Rodando MCScanX..."
    MCScanX "${OUTDIR}/all_species" -a -e 1e-10 -s 5 -m 25 -w 5
    echo "MCScanX concluído."
else
    echo ""
    echo "MCScanX não encontrado. Instale: mamba install -c bioconda mcscanx"
    echo "Ou use TBtools-II (GUI): Synteny → MCScanX Wrapper"
fi

# ── 5. Extrair blocos com os 7 RLPs focais ────────────────────────────────────
export COLLINEARITY_PY="${OUTDIR}/all_species.collinearity"
export OUTFILE_PY="${OUTDIR}/rlp_synteny_blocks.txt"
export GENES_FILE_PY="$GENES_FILE"

python3 << 'PYEOF'
import os, sys

collinearity_file = os.environ["COLLINEARITY_PY"]
out_file          = os.environ["OUTFILE_PY"]
genes_file        = os.environ["GENES_FILE_PY"]

with open(genes_file) as f:
    genes_of_interest = [line.strip() for line in f if line.strip()]

if not os.path.exists(collinearity_file):
    print(f"Arquivo collinearity ainda não gerado: {collinearity_file}")
    sys.exit(0)

with open(collinearity_file) as f_in, open(out_file, "w") as f_out:
    block = []
    keep  = False
    for line in f_in:
        if line.startswith("#"):
            if block and keep:
                f_out.writelines(block)
            block = [line]
            keep  = False
        else:
            block.append(line)
            if any(g in line for g in genes_of_interest):
                keep = True
    if block and keep:
        f_out.writelines(block)

print(f"Blocos de sintenia com os 7 RLPs salvos: {out_file}")
PYEOF

echo ""
echo "VISUALIZAÇÃO FINAL: TBtools-II (GUI)"
echo "  1. Advanced Circos Plot → ${OUTDIR}/all_species.collinearity"
echo "  2. Filtrar por cromossomos dos 7 RLPs"
