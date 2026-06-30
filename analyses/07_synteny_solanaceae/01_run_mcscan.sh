#!/bin/bash
# Sintenia interespécie: tomate × batata × pimenta (LAST + MCScanX)
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Dependências: wget, last, mcscanx, python3
#   mamba install -n kerson-paper -c bioconda last mcscanx -y
#
# Tomato: GFF e proteínas vêm do repositório (git pull já os incluiu).
# Batata/pimenta: baixados do NCBI (acessível no servidor UFV).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DB_DIR="/home/eulalio/databases/solanaceae_synteny"
OUTDIR="${SCRIPT_DIR}/mcscan_results"
GENES_FILE="${REPO_ROOT}/analyses/genes_7rlp.txt"

mkdir -p "$DB_DIR" "$OUTDIR"

# ── Todos os arquivos vêm do repositório (git pull inclui tudo) ───────────────
# Tomato: ITAG4.0 (SGN, 2019)
REPO_TOMATO_GFF="${SCRIPT_DIR}/ITAG4.0_gene_models.gff.gz"
REPO_TOMATO_PEP="${REPO_ROOT}/analyses/06_domain_architecture/ITAG4.0_proteins.fasta"
# Potato: EnsemblPlants release-61 (Solanum tuberosum SolTub_3.0)
REPO_POTATO_GFF="${SCRIPT_DIR}/potato_gff.gff3.gz"
REPO_POTATO_PEP="${SCRIPT_DIR}/potato_pep.fa.gz"
# Pepper: EnsemblPlants release-61 (Capsicum annuum ASM51225v2)
REPO_PEPPER_GFF="${SCRIPT_DIR}/pepper_gff.gff3.gz"
REPO_PEPPER_PEP="${SCRIPT_DIR}/pepper_pep.fa.gz"

echo "================================================================="
echo " 01_run_mcscan.sh — Sintenia Solanaceae (MCScanX)"
echo "================================================================="
echo ""

# ── 1. Preparar tomate (cópias locais do repo) ────────────────────────────────
echo "[1/5] Preparando arquivos de tomate..."
TOMATO_GFF_FINAL="${DB_DIR}/tomato.gff"
TOMATO_PEP_FINAL="${DB_DIR}/tomato_pep"

if [ ! -f "$TOMATO_GFF_FINAL" ]; then
    if [ -f "$REPO_TOMATO_GFF" ]; then
        echo "  Descomprimindo GFF tomato do repo..."
        gunzip -c "$REPO_TOMATO_GFF" > "$TOMATO_GFF_FINAL"
    else
        echo "  ERRO: ${REPO_TOMATO_GFF} não encontrado."
        echo "  Certifique-se de ter rodado: git pull"
        exit 1
    fi
else
    echo "  GFF tomato já disponível."
fi

if [ ! -f "$TOMATO_PEP_FINAL" ]; then
    if [ -f "$REPO_TOMATO_PEP" ]; then
        echo "  Copiando proteínas tomato do repo..."
        cp "$REPO_TOMATO_PEP" "$TOMATO_PEP_FINAL"
    else
        echo "  ERRO: ${REPO_TOMATO_PEP} não encontrado."
        exit 1
    fi
else
    echo "  Proteínas tomato já disponíveis."
fi

# ── 2. Download batata e pimenta (NCBI — acessível no servidor) ───────────────
echo "[2/5] Preparando anotações batata/pimenta (do repositório)..."

copy_repo_file() {
    local repo_src="$1"
    local dest_name="$2"
    local dest="${DB_DIR}/${dest_name}"
    if [ -f "$dest" ]; then
        echo "  ${dest_name} já disponível."
        return
    fi
    if [ ! -f "$repo_src" ]; then
        echo "  ERRO: ${repo_src} não encontrado — rode: git pull"
        exit 1
    fi
    echo "  Descomprimindo ${dest_name}..."
    gunzip -c "$repo_src" > "$dest"
    echo "  OK: ${dest_name}"
}

copy_repo_file "$REPO_POTATO_GFF" "potato.gff"
copy_repo_file "$REPO_POTATO_PEP" "potato_pep"
copy_repo_file "$REPO_PEPPER_GFF" "pepper.gff"
copy_repo_file "$REPO_PEPPER_PEP" "pepper_pep"

# ── 3. GFF → BED para MCScanX ─────────────────────────────────────────────────
echo "[3/5] Convertendo GFF → BED..."
export DB_DIR_PY="$DB_DIR"

python3 << 'PYEOF'
import os, re, gzip

db = os.environ["DB_DIR_PY"]
species_gff = {
    "tomato": os.path.join(db, "tomato.gff"),
    "potato": os.path.join(db, "potato.gff"),
    "pepper": os.path.join(db, "pepper.gff"),
}

for sp, gff in species_gff.items():
    if not os.path.exists(gff):
        print(f"  AVISO: {gff} ausente — pulando {sp}")
        continue
    bed_out = os.path.join(db, f"{sp}.bed")
    rows = []
    open_fn = gzip.open if gff.endswith(".gz") else open
    with open_fn(gff, "rt", errors="ignore") as f:
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
    print(f"  {sp}: {len(rows)} genes -> {bed_out}")
PYEOF

# ── 4. LAST: comparação proteína × proteína ───────────────────────────────────
echo "[4/5] Comparando proteínas com LAST..."

if ! command -v lastdb &>/dev/null; then
    echo "  AVISO: LAST não encontrado."
    echo "  Instale: mamba install -n kerson-paper -c bioconda last -y"
else
    SPECIES=("tomato" "potato" "pepper")
    for sp2 in "${SPECIES[@]}"; do
        pep2="${DB_DIR}/${sp2}_pep"
        [ -f "$pep2" ] || continue
        DB_PREFIX="${OUTDIR}/${sp2}_db"
        [ -f "${DB_PREFIX}.prj" ] || lastdb -p "$DB_PREFIX" "$pep2"
        for sp1 in "${SPECIES[@]}"; do
            [ "$sp1" = "$sp2" ] && continue
            pep1="${DB_DIR}/${sp1}_pep"
            [ -f "$pep1" ] || continue
            BLAST_OUT="${OUTDIR}/${sp1}_vs_${sp2}.blast"
            [ -f "$BLAST_OUT" ] && continue
            echo "  ${sp1} × ${sp2}..."
            lastal -f BlastTab "$DB_PREFIX" "$pep1" > "$BLAST_OUT"
        done
    done
fi

# ── 5. MCScanX ────────────────────────────────────────────────────────────────
echo "[5/5] Rodando MCScanX..."

# Combinar BEDs e BLASTs disponíveis
cat "${DB_DIR}/tomato.bed" "${DB_DIR}/potato.bed" "${DB_DIR}/pepper.bed" \
    > "${OUTDIR}/all_species.gff" 2>/dev/null || true

for pair in "tomato_vs_potato" "tomato_vs_pepper" "potato_vs_pepper"; do
    [ -f "${OUTDIR}/${pair}.blast" ] && cat "${OUTDIR}/${pair}.blast"
done > "${OUTDIR}/all_species.blast" 2>/dev/null || true

if ! command -v MCScanX &>/dev/null; then
    echo "  AVISO: MCScanX não encontrado."
    echo "  Instale: mamba install -n kerson-paper -c bioconda mcscanx -y"
else
    MCScanX "${OUTDIR}/all_species" -a -e 1e-10 -s 5 -m 25 -w 5
    echo "  MCScanX concluído."

    # Extrair blocos com os 7 RLPs focais
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
    print(f"Arquivo collinearity não gerado: {collinearity_file}")
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

print(f"Blocos de sintenia com os 7 RLPs: {out_file}")
PYEOF
fi

echo ""
echo "================================================================="
echo " CONCLUIDO"
echo " Resultados em: ${OUTDIR}/"
echo "================================================================="
echo ""
echo "VISUALIZACAO: TBtools-II"
echo "  Advanced Circos Plot -> ${OUTDIR}/all_species.collinearity"
