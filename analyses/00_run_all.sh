#!/bin/bash
# 00_run_all.sh — Orquestrador do pipeline bioinformático Kerson-paper
# Análises genome-wide de 49 LRR-RLPs em Solanum lycopersicum (ITAG4.0)
#
# Executar no servidor: eulalio@200.235.143.10
# Ambiente: mamba activate kerson-paper
#
# Uso:
#   bash 00_run_all.sh --all              # Executa Dias 2-7 em sequência
#   bash 00_run_all.sh --step 6           # Executa apenas o Dia 6
#   bash 00_run_all.sh --step dia6        # Idem (aceita "dia6" ou "6")
#   bash 00_run_all.sh --step 2,3,6       # Múltiplos steps (separados por vírgula)
#   bash 00_run_all.sh --list             # Lista todos os steps disponíveis
#   bash 00_run_all.sh --dry-run --all    # Mostra o que seria executado sem rodar
#
# Steps:
#   2 = Elementos cis nos promotores (PlantCARE, requer interação manual)
#   3 = Atlas de expressão (TFGD / RNA-Seq)
#   4 = Ka/Ks — seleção purificadora (PAML/yn00)
#   5 = Qualidade dos modelos 3D (RMSD + MolProbity)
#   6 = Arquitetura de domínios Pfam (HMMER) — pode rodar AGORA
#   7 = Sintenia Solanaceae (MCScanX)

set -uo pipefail
# Nota: sem -e global; cada step tem seu próprio tratamento de erro

# ── Caminhos ──────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_FILE="${SCRIPT_DIR}/pipeline_run.log"
GENES_FILE="${REPO_ROOT}/ids_49_rlp_tomato.txt"

# ── Helpers de log ────────────────────────────────────────────────────────────
log() {
    local ts
    ts="$(date '+%Y-%m-%d %H:%M:%S')"
    local msg="[${ts}] $*"
    echo "$msg"
    echo "$msg" >> "$LOG_FILE"
}

log_err() {
    local ts
    ts="$(date '+%Y-%m-%d %H:%M:%S')"
    local msg="[${ts}] ERRO: $*"
    echo "$msg" >&2
    echo "$msg" >> "$LOG_FILE"
}

log_sep() {
    local line="================================================================="
    echo "$line"
    echo "$line" >> "$LOG_FILE"
}

# ── Flags globais ─────────────────────────────────────────────────────────────
DRY_RUN=0
STEPS_TO_RUN=()

# ── Verificar se está dentro de um screen ─────────────────────────────────────
check_screen() {
    if [ -z "${STY:-}" ] && [ -z "${TMUX:-}" ]; then
        log_err "Processo NÃO está dentro de screen ou tmux!"
        log_err "Processos longos podem ser mortos por SIGTTOU ao desconectar SSH."
        log_err "Execute:"
        log_err "  screen -S kerson-paper"
        log_err "  mamba activate kerson-paper"
        log_err "  bash ${SCRIPT_DIR}/00_run_all.sh $*"
        echo ""
        read -r -p "Continuar mesmo assim? [s/N]: " REPLY
        if [[ ! "$REPLY" =~ ^[Ss]$ ]]; then
            exit 1
        fi
    else
        log "Sessao screen/tmux detectada: ${STY:-${TMUX}} — OK"
    fi
}

# ── Verificar ferramenta disponível ──────────────────────────────────────────
require_tool() {
    local tool="$1"
    if ! command -v "$tool" &>/dev/null; then
        log_err "Ferramenta ausente: ${tool}"
        log_err "Instale com: mamba install -n kerson-paper -c bioconda ${tool}"
        return 1
    fi
    return 0
}

check_deps_base() {
    log "Verificando dependencias basicas..."
    local ok=1
    for tool in python3 Rscript samtools bedtools; do
        if command -v "$tool" &>/dev/null; then
            log "  [OK] $tool"
        else
            log_err "  [AUSENTE] $tool"
            ok=0
        fi
    done
    [ "$ok" -eq 1 ] || { log_err "Instale as dependencias antes de continuar."; exit 1; }
}

# ── git pull automático ──────────────────────────────────────────────────────
git_pull() {
    log "Sincronizando repositorio (git pull)..."
    if git -C "$REPO_ROOT" pull --ff-only 2>>"$LOG_FILE"; then
        log "  git pull OK"
    else
        log_err "  git pull falhou — continuando com versao local"
    fi
}

# ══════════════════════════════════════════════════════════════════════════════
# DIA 2 — Elementos cis nos promotores (PlantCARE)
# ══════════════════════════════════════════════════════════════════════════════
run_step2() {
    log_sep
    log "DIA 2: Elementos cis — promotores 2 kb upstream + PlantCARE"

    local DIR="${SCRIPT_DIR}/02_promoter_cis"
    local UPSTREAM_FA="${DIR}/rlp_upstream_2kb.fa"
    local PLANTCARE="${DIR}/plantcare_results.txt"
    local COORDS_TSV="${DIR}/coords_49genes.tsv"

    # Verificar/gerar coordenadas
    if [ ! -f "$COORDS_TSV" ]; then
        log "  coords_49genes.tsv ausente — executando 00_get_coords_all49.py..."
        if [ "$DRY_RUN" -eq 0 ]; then
            python3 "${DIR}/00_get_coords_all49.py" 2>>"$LOG_FILE" || {
                log_err "00_get_coords_all49.py falhou — abortando Dia 2"
                return 1
            }
        else
            log "  [DRY-RUN] python3 ${DIR}/00_get_coords_all49.py"
        fi
    else
        log "  coords_49genes.tsv presente ($(tail -n +2 "$COORDS_TSV" | wc -l) genes)"
    fi

    # Gerar upstream FASTA
    if [ ! -f "$UPSTREAM_FA" ]; then
        require_tool bedtools || return 1
        require_tool samtools || return 1
        log "  Extraindo regioes upstream com 01_fetch_upstream.sh..."
        if [ "$DRY_RUN" -eq 0 ]; then
            bash "${DIR}/01_fetch_upstream.sh" 2>>"$LOG_FILE" || {
                log_err "01_fetch_upstream.sh falhou — abortando Dia 2"
                return 1
            }
        else
            log "  [DRY-RUN] bash ${DIR}/01_fetch_upstream.sh"
        fi
    else
        log "  upstream FASTA ja existe: ${UPSTREAM_FA}"
        log "  Sequencias: $(grep -c '^>' "$UPSTREAM_FA") genes"
    fi

    # Pausa manual — PlantCARE requer submissao web
    log ""
    log "  *** PAUSA MANUAL NECESSARIA ***"
    log "  O PlantCARE NAO tem API automatica — submissao manual obrigatoria."
    log ""
    log "  1. Acesse: https://bioinformatics.psb.ugent.be/webtools/plantcare/html/"
    log "  2. Clique em 'Search Promoter'"
    log "  3. Cole o conteudo de: ${UPSTREAM_FA}"
    log "  4. Execute a analise e clique em 'Download results'"
    log "  5. Salve o arquivo TXT como: ${PLANTCARE}"
    log "  6. Depois execute: bash ${SCRIPT_DIR}/00_run_all.sh --step 2b"
    log ""

    if [ "$DRY_RUN" -eq 0 ]; then
        read -r -p "Pressione ENTER apos salvar plantcare_results.txt (ou Ctrl+C para sair): "
    fi

    # Gerar heatmap PlantCARE
    if [ -f "$PLANTCARE" ]; then
        log "  Gerando heatmap de cis-elements..."
        if [ "$DRY_RUN" -eq 0 ]; then
            Rscript "${DIR}/02_plot_plantcare_heatmap.R" 2>>"$LOG_FILE" || {
                log_err "02_plot_plantcare_heatmap.R falhou"
                return 1
            }
        else
            log "  [DRY-RUN] Rscript ${DIR}/02_plot_plantcare_heatmap.R"
        fi
        log "  Figuras geradas em: ${DIR}/"
    else
        log_err "  plantcare_results.txt ausente — execute o PlantCARE e repita --step 2"
        return 1
    fi

    log "DIA 2 CONCLUIDO"
}

# ══════════════════════════════════════════════════════════════════════════════
# DIA 3 — Atlas de expressão (TFGD)
# ══════════════════════════════════════════════════════════════════════════════
run_step3() {
    log_sep
    log "DIA 3: Atlas de Expressao — TFGD"

    local DIR="${SCRIPT_DIR}/03_expression_atlas"

    require_tool python3 || return 1

    log "  Baixando dados de expressao do TFGD..."
    if [ "$DRY_RUN" -eq 0 ]; then
        python3 "${DIR}/01_fetch_expression.py" 2>>"$LOG_FILE" || {
            log_err "01_fetch_expression.py falhou"
            return 1
        }
    else
        log "  [DRY-RUN] python3 ${DIR}/01_fetch_expression.py"
    fi

    log "  Gerando heatmaps de expressao..."
    if [ "$DRY_RUN" -eq 0 ]; then
        Rscript "${DIR}/02_plot_expression_heatmap.R" 2>>"$LOG_FILE" || {
            log_err "02_plot_expression_heatmap.R falhou"
            return 1
        }
    else
        log "  [DRY-RUN] Rscript ${DIR}/02_plot_expression_heatmap.R"
    fi

    log "  Verificar saida: ${DIR}/expression_atlas_heatmap.pdf"
    log "DIA 3 CONCLUIDO"
}

# ══════════════════════════════════════════════════════════════════════════════
# DIA 4 — Ka/Ks (duplicações, seleção purificadora)
# ══════════════════════════════════════════════════════════════════════════════
run_step4() {
    log_sep
    log "DIA 4: Ka/Ks — KaKs_Calculator / PAML yn00"

    local DIR="${SCRIPT_DIR}/04_kaks"
    local PAIRS="${DIR}/gene_pairs.tsv"

    require_tool mafft || return 1

    if [ ! -f "$PAIRS" ]; then
        log_err "  gene_pairs.tsv ausente: ${PAIRS}"
        log_err "  Gere o arquivo com os pares de genes paralogos antes de continuar."
        return 1
    fi

    log "  Pares de genes: $(wc -l < "$PAIRS")"
    if [ "$DRY_RUN" -eq 0 ]; then
        bash "${DIR}/01_run_kaks_pipeline.sh" 2>>"$LOG_FILE" || {
            log_err "01_run_kaks_pipeline.sh falhou"
            return 1
        }
    else
        log "  [DRY-RUN] bash ${DIR}/01_run_kaks_pipeline.sh"
    fi

    log "  Verificar saida: ${DIR}/kaks_summary.tsv"
    log "DIA 4 CONCLUIDO"
}

# ══════════════════════════════════════════════════════════════════════════════
# DIA 5 — Qualidade dos modelos 3D (RMSD + MolProbity)
# ══════════════════════════════════════════════════════════════════════════════
run_step5() {
    log_sep
    log "DIA 5: Qualidade dos modelos 3D (RMSD + MolProbity)"

    local DIR="${SCRIPT_DIR}/05_rmsd_quality"
    local PDB_DIR="${DIR}/pdb_models"

    if [ ! -d "$PDB_DIR" ] || [ -z "$(ls -A "$PDB_DIR" 2>/dev/null)" ]; then
        log_err "  Diretorio de PDBs vazio ou ausente: ${PDB_DIR}"
        log_err "  Coloque os arquivos .pdb em ${PDB_DIR}/"
        log_err "  Fontes: Swiss-Model (https://swissmodel.expasy.org/)"
        log_err "          AlphaFold3  (https://alphafoldserver.com/)"
        return 1
    fi

    local N_PDBS
    N_PDBS=$(ls "${PDB_DIR}"/*.pdb 2>/dev/null | wc -l)
    log "  PDBs encontrados: ${N_PDBS}"

    if [ "$DRY_RUN" -eq 0 ]; then
        python3 "${DIR}/01_calc_rmsd_molprobity.py" \
            --pdb-dir "$PDB_DIR" \
            --mode both 2>>"$LOG_FILE" || {
            log_err "01_calc_rmsd_molprobity.py falhou"
            return 1
        }
    else
        log "  [DRY-RUN] python3 ${DIR}/01_calc_rmsd_molprobity.py --pdb-dir ${PDB_DIR} --mode both"
    fi

    log "  Verificar saida: ${DIR}/rmsd_results.tsv e rmsd_calc.py (para PyMOL)"
    log "DIA 5 CONCLUIDO"
}

# ══════════════════════════════════════════════════════════════════════════════
# DIA 6 — Arquitetura de domínios Pfam (HMMER)
# ══════════════════════════════════════════════════════════════════════════════
run_step6() {
    log_sep
    log "DIA 6: Arquitetura de Dominios Pfam — HMMER"

    local DIR="${SCRIPT_DIR}/06_domain_architecture"
    local PROTEINS="${DIR}/proteins_49rlp.fa"
    local DOMAINS_TSV="${DIR}/hmmer_out/hmmer_domains.tsv"

    require_tool hmmscan || return 1

    # Etapa 6a: Buscar proteínas
    if [ ! -f "$PROTEINS" ]; then
        log "  proteins_49rlp.fa ausente — executando 00_fetch_proteins.sh..."
        if [ "$DRY_RUN" -eq 0 ]; then
            bash "${DIR}/00_fetch_proteins.sh" 2>>"$LOG_FILE" || {
                log_err "00_fetch_proteins.sh falhou"
                return 1
            }
        else
            log "  [DRY-RUN] bash ${DIR}/00_fetch_proteins.sh"
        fi
    else
        log "  proteins_49rlp.fa presente ($(grep -c '^>' "$PROTEINS") sequencias)"
    fi

    # Etapa 6b: HMMER
    if [ ! -f "$DOMAINS_TSV" ]; then
        log "  Executando hmmscan vs Pfam-A (pode levar 10-30 min)..."
        if [ "$DRY_RUN" -eq 0 ]; then
            bash "${DIR}/01_run_hmmer.sh" 2>>"$LOG_FILE" || {
                log_err "01_run_hmmer.sh falhou"
                return 1
            }
        else
            log "  [DRY-RUN] bash ${DIR}/01_run_hmmer.sh"
        fi
    else
        log "  hmmer_domains.tsv presente ($(tail -n +2 "$DOMAINS_TSV" | wc -l) dominios)"
    fi

    # Etapa 6c: Plot
    log "  Gerando plot de arquitetura de dominios..."
    if [ "$DRY_RUN" -eq 0 ]; then
        Rscript "${DIR}/02_plot_domain_architecture.R" "$DOMAINS_TSV" 2>>"$LOG_FILE" || {
            log_err "02_plot_domain_architecture.R falhou"
            return 1
        }
    else
        log "  [DRY-RUN] Rscript ${DIR}/02_plot_domain_architecture.R ${DOMAINS_TSV}"
    fi

    log "  Verificar saida: ${DIR}/domain_architecture.pdf"
    log "DIA 6 CONCLUIDO"
}

# ══════════════════════════════════════════════════════════════════════════════
# DIA 7 — Sintenia Solanaceae (MCScanX)
# ══════════════════════════════════════════════════════════════════════════════
run_step7() {
    log_sep
    log "DIA 7: Sintenia Solanaceae — MCScanX"

    local DIR="${SCRIPT_DIR}/07_synteny_solanaceae"

    if [ "$DRY_RUN" -eq 0 ]; then
        bash "${DIR}/01_run_mcscan.sh" 2>>"$LOG_FILE" || {
            log_err "01_run_mcscan.sh falhou"
            return 1
        }
    else
        log "  [DRY-RUN] bash ${DIR}/01_run_mcscan.sh"
    fi

    log "  Verificar saida: ${DIR}/mcscan_results/rlp_synteny_blocks.txt"
    log "  Visualizacao final: TBtools-II -> Advanced Circos Plot"
    log "DIA 7 CONCLUIDO"
}

# ══════════════════════════════════════════════════════════════════════════════
# Dispatch de steps
# ══════════════════════════════════════════════════════════════════════════════
run_step() {
    local step="$1"
    # Normalizar: "dia6", "6", "step6" → "6"
    step="${step//[Dd][Ii][Aa]/}"
    step="${step//[Ss][Tt][Ee][Pp]/}"
    step="${step// /}"

    case "$step" in
        2) run_step2 ;;
        3) run_step3 ;;
        4) run_step4 ;;
        5) run_step5 ;;
        6) run_step6 ;;
        7) run_step7 ;;
        *)
            log_err "Step desconhecido: '$step' (validos: 2 3 4 5 6 7)"
            return 1
            ;;
    esac
}

list_steps() {
    echo ""
    echo "Steps disponíveis:"
    echo "  2  — PlantCARE: elementos cis nos promotores 2 kb upstream"
    echo "       (requer interação manual para submissão ao PlantCARE)"
    echo "  3  — Atlas de expressão via TFGD"
    echo "  4  — Ka/Ks: seleção purificadora em pares de genes paralogos"
    echo "       (requer gene_pairs.tsv preenchido)"
    echo "  5  — Qualidade de modelos 3D: RMSD + MolProbity"
    echo "       (requer PDBs em 05_rmsd_quality/pdb_models/)"
    echo "  6  — Arquitetura de dominios Pfam via HMMER  *** PODE RODAR AGORA ***"
    echo "  7  — Sintenia Solanaceae via MCScanX"
    echo ""
    echo "Exemplos:"
    echo "  bash 00_run_all.sh --step 6          # rodar apenas Dia 6 (HMMER)"
    echo "  bash 00_run_all.sh --step 6,7        # rodar Dia 6 e 7"
    echo "  bash 00_run_all.sh --all             # rodar todos"
    echo "  bash 00_run_all.sh --dry-run --all   # simular sem executar"
    echo ""
}

# ── Parsing de argumentos ─────────────────────────────────────────────────────
if [ $# -eq 0 ]; then
    echo "Uso: bash 00_run_all.sh [--step N] [--all] [--list] [--dry-run]"
    list_steps
    exit 0
fi

PARSED_STEPS=""
RUN_ALL=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --step|-s)
            shift
            PARSED_STEPS="${PARSED_STEPS},${1}"
            shift
            ;;
        --step=*)
            PARSED_STEPS="${PARSED_STEPS},${1#--step=}"
            shift
            ;;
        --all|-a)
            RUN_ALL=1
            shift
            ;;
        --dry-run|-n)
            DRY_RUN=1
            shift
            ;;
        --list|-l)
            list_steps
            exit 0
            ;;
        --help|-h)
            head -30 "$0"
            list_steps
            exit 0
            ;;
        *)
            echo "Argumento desconhecido: $1"
            echo "Use: bash 00_run_all.sh --help"
            exit 1
            ;;
    esac
done

# ── Inicialização ─────────────────────────────────────────────────────────────
{
    log_sep
    log "Kerson-paper — Pipeline Bioinformatica"
    log "Repositorio: ${REPO_ROOT}"
    log "Log: ${LOG_FILE}"
    [ "$DRY_RUN" -eq 1 ] && log "[DRY-RUN ATIVO — nenhum script sera executado]"

    # Verificar screen (apenas para --all ou steps longos)
    if [ "$RUN_ALL" -eq 1 ] || echo "$PARSED_STEPS" | grep -qE "6|7|4"; then
        check_screen "$@" || true
    fi

    # git pull antes de qualquer coisa
    git_pull

    # Verificar dependências básicas
    check_deps_base

    # ── Executar steps ────────────────────────────────────────────────────────
    FAILED_STEPS=()

    if [ "$RUN_ALL" -eq 1 ]; then
        for s in 2 3 4 5 6 7; do
            run_step "$s" || FAILED_STEPS+=("$s")
        done
    else
        # Processar steps individuais: "2,3,6" ou "6"
        IFS=',' read -ra STEP_LIST <<< "${PARSED_STEPS#,}"
        for s in "${STEP_LIST[@]}"; do
            s="${s// /}"
            [ -z "$s" ] && continue
            run_step "$s" || FAILED_STEPS+=("$s")
        done
    fi

    # ── Relatório final ───────────────────────────────────────────────────────
    log_sep
    if [ ${#FAILED_STEPS[@]} -eq 0 ]; then
        log "PIPELINE CONCLUIDO COM SUCESSO"
        log ""
        log "Sincronizar resultados:"
        log "  cd ${REPO_ROOT}"
        log "  git add analyses/ results/"
        log "  git commit -m 'Add bioinformatics results'"
        log "  git push"
    else
        log_err "PIPELINE CONCLUIDO COM ERROS"
        log_err "Steps com falha: ${FAILED_STEPS[*]}"
        log_err "Verificar log completo: ${LOG_FILE}"
        exit 1
    fi
}
