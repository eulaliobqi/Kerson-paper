# Kerson-paper — Análise Genome-wide de LRR-RLPs em *Solanum lycopersicum*

Pipeline bioinformático completo para o artigo científico sobre identificação e caracterização de Leucine-Rich Repeat Receptor-Like Proteins (LRR-RLPs) em tomate.

**Aluno:** Fólkerson Marinho Ferreira (doutorado, UFV)  
**Orientadora:** Profª Elizabeth Pacheco Batista Fontes  
**Co-orientadores:** Prof. Pedro Augusto Braga Reis · Dr. Eulalio Gutemberg Dias dos Santos

---

## Sobre o projeto

49 LRR-RLPs foram identificados no genoma de *Solanum lycopersicum* (ITAG4.0) via RLPredictOme. Destes, 7 foram selecionados com base em ortologia com *Arabidopsis thaliana* para análise funcional aprofundada (qRT-PCR em seca e infecção por begomovírus ToYSV).

Este repositório contém os scripts para 6 análises bioinformáticas complementares que compõem o artigo final.

---

## Requisitos (servidor)

- Servidor Linux com GPU (testado em Debian, RTX 5070 Ti, 32 cores)
- [Miniforge3](https://github.com/conda-forge/miniforge) com `mamba`
- Git, wget

---

## Setup inicial (executar uma única vez)

```bash
# 1. Clonar o repositório
git clone https://github.com/eulaliobqi/Kerson-paper ~/kerson-paper
cd ~/kerson-paper

# 2. Criar ambiente conda com todas as dependências
bash analyses/00_setup_environment.sh

# 3. Ativar o ambiente (necessário antes de cada sessão de análise)
mamba activate kerson-paper
```

**Dependências instaladas automaticamente:**  
`bedtools`, `samtools`, `hmmer`, `mafft`, `kaks-calculator`, `last`, `mcscanx`,  
`biopython`, `pymol-open-source`, `R` + `tidyverse`, `pheatmap`, `ggplot2`, `gggenes`

---

## Dados de entrada

| Arquivo | Conteúdo |
|---------|----------|
| `ids_49_rlp_tomato.txt` | 49 IDs dos RLPs identificados por RLPredictOme (ITAG4.0) |
| `analyses/genes_7rlp.txt` | 7 genes focais para análise funcional |
| `analyses/gene_metadata.tsv` | Metadados dos 7 genes (ortólogos, função, n° LRR, cromossomo) |
| `analyses/04_kaks/gene_pairs.tsv` | 13 pares de genes para análise Ka/Ks |

---

## Análises — Passo a Passo

### Dia 2 — Elementos Cis nos Promotores (PlantCARE)

**Objetivo:** Identificar elementos cis-regulatórios nos 2 kb upstream dos 7 RLPs focais.  
**Figura gerada:** `analyses/02_promoter_cis/plantcare_heatmap.pdf`

```bash
mamba activate kerson-paper
cd ~/kerson-paper

# Passo 1: Extrair sequências upstream do genoma ITAG4.0
bash analyses/02_promoter_cis/01_fetch_upstream.sh
# → Gera: analyses/02_promoter_cis/rlp_upstream_2kb.fa
```

**Passo 2 (manual):**
1. Acesse: https://bioinformatics.psb.ugent.be/webtools/plantcare/html/
2. Selecione "Search Promoter" → cole o conteúdo de `rlp_upstream_2kb.fa`
3. Baixe os resultados como TXT → salve como `analyses/02_promoter_cis/plantcare_results.txt`

```bash
# Passo 3: Gerar heatmap (funciona em modo DEMO mesmo sem plantcare_results.txt)
Rscript analyses/02_promoter_cis/02_plot_plantcare_heatmap.R
# → Gera: plantcare_heatmap.pdf + plantcare_counts.csv
```

**Saída esperada:** Heatmap com 7 genes × categorias de elementos cis  
(ABA/seca, JA, SA/defesa, GA, luz, desenvolvimento, circadiano)

---

### Dia 3 — Atlas de Expressão Pública (TFGD)

**Objetivo:** Perfil de expressão dos 7 RLPs em tecidos, seca e infecção viral.  
**Figura gerada:** `analyses/03_expression_atlas/expression_atlas_heatmap.pdf`

```bash
mamba activate kerson-paper
cd ~/kerson-paper

# Buscar dados do TFGD (Cornell) — cai para DEMO se indisponível
python3 analyses/03_expression_atlas/01_fetch_expression.py
# → Gera: analyses/03_expression_atlas/expression_matrix.csv

# Gerar heatmaps (z-score + absoluto)
Rscript analyses/03_expression_atlas/02_plot_expression_heatmap.R
# → Gera: expression_atlas_heatmap.pdf + expression_atlas_absolute.pdf
```

**Dica:** Se os dados online estiverem indisponíveis, edite `expression_matrix.csv`  
com valores de FPKM de estudos RNA-Seq relevantes do NCBI GEO.

---

### Dia 4 — Pressão Seletiva Ka/Ks (KaKs_Calculator 2.0)

**Objetivo:** Calcular Ka/Ks para 13 pares de RLPs duplicados, inferindo modo de evolução.  
**Saída:** `analyses/04_kaks/kaks_summary.tsv`

```bash
mamba activate kerson-paper
cd ~/kerson-paper

# Rodar pipeline completo (ITAG4.0 CDS → MAFFT → KaKs_Calculator)
bash analyses/04_kaks/01_run_kaks_pipeline.sh
```

**Pares analisados (13 total):**
- 11 pares de **duplicação em tandem** (clusters no chr01, chr06, chr07, chr12)
- 2 pares de **duplicação segmental** (chr02 e chr05/chr03)

**Interpretação:**  
`Ka/Ks < 1` = seleção purificadora (função conservada)  
`Ka/Ks ≈ 1` = evolução neutra  
`Ka/Ks > 1` = seleção positiva (neofuncionalização)

---

### Dia 5 — Qualidade dos Modelos 3D (RMSD + MolProbity)

**Objetivo:** Quantificar RMSD entre os 7 modelos 3D e validar qualidade estrutural.  
**Saídas:** `rmsd_matrix.csv` + template `molprobity_template.tsv`

```bash
mamba activate kerson-paper
cd ~/kerson-paper/analyses/05_rmsd_quality

# Passo 1: Colocar os 7 PDBs no diretório pdb_models/
# Fontes recomendadas:
#   AlphaFold DB: https://alphafold.ebi.ac.uk/ (buscar por Solyc ID)
#   Swiss-Model: https://swissmodel.expasy.org/
# Nomear como: Solyc05g055190.pdb, Solyc03g112680.pdb, etc.

# Passo 2: Gerar script PyMOL e template MolProbity
python3 01_calc_rmsd_molprobity.py --pdb-dir ./pdb_models/ --mode both

# Passo 3: Calcular RMSD
pymol -c rmsd_calc.py ./pdb_models/
# → Gera: rmsd_matrix.csv
```

**Passo 4 (MolProbity — manual):**
1. Acesse: https://molprobity.biochem.duke.edu/
2. Submeta cada PDB → anote Ramachandran favored (%), clashscore
3. Preencha `molprobity_template.tsv`

**Métricas alvo:** Ramachandran favored > 95%, Clashscore < 20

---

### Dia 6 — Arquitetura de Domínios (49 RLPs via Pfam)

**Objetivo:** Mapear domínios proteicos (LRR, sinal, TM) em todos os 49 RLPs.  
**Figura gerada:** `analyses/06_domain_architecture/domain_architecture.pdf`

```bash
mamba activate kerson-paper
cd ~/kerson-paper

# Passo 1: Extrair sequências dos 49 RLPs do ITAG4.0
bash analyses/06_domain_architecture/00_fetch_proteins.sh
# → Gera: analyses/06_domain_architecture/proteins_49rlp.fa

# Passo 2: Anotar domínios Pfam (hmmscan, ~1.5 GB Pfam-A, 16 CPUs, ~30 min)
bash analyses/06_domain_architecture/01_run_hmmer.sh
# → Gera: analyses/06_domain_architecture/hmmer_out/hmmer_domains.tsv

# Passo 3: Gerar figura de arquitetura de domínios
Rscript analyses/06_domain_architecture/02_plot_domain_architecture.R \
    analyses/06_domain_architecture/hmmer_out/hmmer_domains.tsv
# → Gera: domain_architecture.pdf + domain_architecture.png
```

**Nota:** O script de plot funciona em modo DEMO (sem `hmmer_domains.tsv`) para  
testar a visualização antes de rodar o HMMER.

---

### Dia 7 — Sintenia Solanaceae (LAST + MCScanX)

**Objetivo:** Identificar blocos de sintenia envolvendo os RLPs em tomate × batata × pimenta.  
**Saídas:** `mcscan_results/all_species.collinearity` + `rlp_synteny_blocks.txt`

```bash
mamba activate kerson-paper
cd ~/kerson-paper

# Pipeline completo: download genomas → LAST → MCScanX → filtrar RLPs focais
bash analyses/07_synteny_solanaceae/01_run_mcscan.sh
```

**Visualização final (TBtools-II):**
1. Abra TBtools-II → Synteny → Advanced Circos Plot
2. Carregue `mcscan_results/all_species.collinearity`
3. Filtre por blocos em `rlp_synteny_blocks.txt`

---

## Saídas esperadas (figuras para o artigo)

| Figura | Arquivo | Análise |
|--------|---------|---------|
| Fig. nova — Elementos cis | `plantcare_heatmap.pdf` | Dia 2 |
| Fig. nova — Atlas expressão | `expression_atlas_heatmap.pdf` | Dia 3 |
| Fig. existente (expandir) | `domain_architecture.pdf` | Dia 6 |
| Fig. existente (expandir) | Via TBtools-II | Dia 7 |
| Supl. nova — Ka/Ks | `kaks_summary.tsv` | Dia 4 |
| Supl. nova — Qualidade 3D | `rmsd_matrix.csv` + `molprobity_template.tsv` | Dia 5 |

---

## Executar tudo de uma vez

```bash
mamba activate kerson-paper
cd ~/kerson-paper
bash analyses/00_run_pipeline.sh all
```

Ou por análise individual:

```bash
bash analyses/00_run_pipeline.sh 2   # PlantCARE
bash analyses/00_run_pipeline.sh 3   # Expressão
bash analyses/00_run_pipeline.sh 4   # Ka/Ks
bash analyses/00_run_pipeline.sh 5   # RMSD/MolProbity
bash analyses/00_run_pipeline.sh 6   # Domínios
bash analyses/00_run_pipeline.sh 7   # Sintenia
```

---

## Citação

> Ferreira, F.M. *et al.* (2025). Genome-Wide Identification of Leucine-Rich Repeat  
> Receptor-Like Proteins in Tomato and Implications in Begomoviral Infection and  
> Drought Stress. *Frontiers in Plant Science* (em preparação).

---

## Contato

eulalio.santos@ufv.br · github.com/eulaliobqi
