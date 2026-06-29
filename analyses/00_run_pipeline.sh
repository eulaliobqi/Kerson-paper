#!/bin/bash
# Pipeline Kerson-paper — análises bioinformáticas (Dias 2–7)
# Executar no servidor: eulalio@200.235.143.10
# Após: git pull origin main
#
# Ordem de execução:
#   Dia 2: PlantCARE (cis-elements)
#   Dia 3: Atlas de expressão (TFGD)
#   Dia 4: Ka/Ks (duplicações)
#   Dia 5: RMSD + MolProbity (qualidade 3D)
#   Dia 6: Arquitetura de domínios Pfam
#   Dia 7: Sintenia Solanaceae (MCScanX)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SERVER="eulalio@200.235.143.10"

echo "================================================================="
echo "Kerson-paper — Pipeline Bioinformática"
echo "Dir: ${SCRIPT_DIR}"
echo "================================================================="
echo ""

# ── Verificações iniciais ─────────────────────────────────────────────────────
check_tool() {
    command -v "$1" &>/dev/null || echo "  AUSENTE: $1 (instalar antes de prosseguir)"
}

echo "Verificando dependências..."
check_tool bedtools
check_tool samtools
check_tool mafft
check_tool hmmscan
check_tool python3
check_tool Rscript
echo ""

# ── DIA 2: Elementos cis nos promotores ──────────────────────────────────────
run_dia2() {
    echo "=== DIA 2: Elementos cis — PlantCARE ==="
    cd "${SCRIPT_DIR}/02_promoter_cis"
    bash 01_fetch_upstream.sh
    echo ""
    echo "PAUSA MANUAL: enviar rlp_upstream_2kb.fa ao PlantCARE"
    echo "  https://bioinformatics.psb.ugent.be/webtools/plantcare/html/"
    echo "  Salvar resultado como: ${SCRIPT_DIR}/02_promoter_cis/plantcare_results.txt"
    echo ""
    echo "Após receber resultado do PlantCARE, executar:"
    echo "  Rscript 02_plot_plantcare_heatmap.R"
    cd "$SCRIPT_DIR"
}

# ── DIA 3: Atlas de expressão (TFGD) ─────────────────────────────────────────
run_dia3() {
    echo "=== DIA 3: Atlas de Expressão — TFGD ==="
    cd "${SCRIPT_DIR}/03_expression_atlas"
    python3 01_fetch_expression.py
    Rscript 02_plot_expression_heatmap.R
    echo "Figuras geradas: expression_atlas_heatmap.pdf + expression_atlas_absolute.pdf"
    cd "$SCRIPT_DIR"
}

# ── DIA 4: Ka/Ks ─────────────────────────────────────────────────────────────
run_dia4() {
    echo "=== DIA 4: Ka/Ks — KaKs_Calculator 2.0 ==="
    cd "${SCRIPT_DIR}/04_kaks"
    bash 01_run_kaks_pipeline.sh
    echo "Resultado: ${SCRIPT_DIR}/04_kaks/kaks_summary.tsv"
    cd "$SCRIPT_DIR"
}

# ── DIA 5: RMSD + MolProbity ─────────────────────────────────────────────────
run_dia5() {
    echo "=== DIA 5: RMSD + MolProbity ==="
    cd "${SCRIPT_DIR}/05_rmsd_quality"
    echo "Colocar PDBs em: ${SCRIPT_DIR}/05_rmsd_quality/pdb_models/"
    echo "  Fontes: Swiss-Model (https://swissmodel.expasy.org/)"
    echo "          AlphaFold3 (https://alphafoldserver.com/)"
    echo ""
    python3 01_calc_rmsd_molprobity.py --pdb-dir ./pdb_models/ --mode both
    echo "Script PyMOL gerado. Executar:"
    echo "  pymol -c rmsd_calc.py ./pdb_models/"
    cd "$SCRIPT_DIR"
}

# ── DIA 6: Arquitetura de domínios ───────────────────────────────────────────
run_dia6() {
    echo "=== DIA 6: Arquitetura de Domínios — Pfam/HMMER ==="
    cd "${SCRIPT_DIR}/06_domain_architecture"
    echo "REQUERIDO: colocar proteínas FASTA em proteins_49rlp.fa"
    echo "  (extrair dos 49 RLPs identificados pelo RLPredictOme)"
    echo ""
    if [ -f "proteins_49rlp.fa" ]; then
        bash 01_run_hmmer.sh
        Rscript 02_plot_domain_architecture.R hmmer_out/hmmer_domains.tsv
    else
        echo "Testando modo DEMO (sem proteins_49rlp.fa)..."
        Rscript 02_plot_domain_architecture.R  # modo demo
    fi
    cd "$SCRIPT_DIR"
}

# ── DIA 7: Sintenia Solanaceae ────────────────────────────────────────────────
run_dia7() {
    echo "=== DIA 7: Sintenia Solanaceae ==="
    cd "${SCRIPT_DIR}/07_synteny_solanaceae"
    bash 01_run_mcscan.sh
    echo "Blocos de sintenia: mcscan_results/rlp_synteny_blocks.txt"
    echo "Visualização final: TBtools-II → Advanced Circos Plot"
    cd "$SCRIPT_DIR"
}

# ── Seleção de análise ────────────────────────────────────────────────────────
DIA="${1:-all}"

case "$DIA" in
    dia2|2) run_dia2 ;;
    dia3|3) run_dia3 ;;
    dia4|4) run_dia4 ;;
    dia5|5) run_dia5 ;;
    dia6|6) run_dia6 ;;
    dia7|7) run_dia7 ;;
    all)
        run_dia2
        echo ""
        run_dia3
        echo ""
        run_dia4
        echo ""
        run_dia5
        echo ""
        run_dia6
        echo ""
        run_dia7
        ;;
    *)
        echo "Uso: bash 00_run_pipeline.sh [2|3|4|5|6|7|all]"
        echo "  2 = PlantCARE     3 = TFGD expressão"
        echo "  4 = Ka/Ks         5 = RMSD + MolProbity"
        echo "  6 = Domínios      7 = Sintenia"
        ;;
esac

echo ""
echo "Pipeline concluído. Sincronizar resultados:"
echo "  git add analyses/ && git commit -m 'Add bioinformatics results' && git push"
