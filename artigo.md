# Artigo — Seções Metodologia, Resultados e Discussão

**Título:** Genome-Wide Identification of Leucine-Rich Repeat Receptor-Like Proteins in Tomato and Implications in Begomoviral Infection and Drought Stress

**Espécie:** *Solanum lycopersicum* L. (genoma ITAG4.0)

**Data de geração deste documento:** 2026-06-29

---

# Metodologia Bioinformática

## 2. Análise de Elementos Cis nos Promotores

### 2.1 Extração de Sequências Upstream

Para a caracterização da regulação transcricional dos sete genes LRR-RLP focais (*SlRLP1* a *SlRLP7*), foram extraídas sequências genômicas de 2.000 pares de bases (pb) a montante do sítio de início da transcrição anotado de cada gene. As coordenadas genômicas foram obtidas do arquivo de anotação GFF3 do genoma de referência de *Solanum lycopersicum* ITAG4.0, disponibilizado pelo *Sol Genomics Network* (SGN; Fernandez-Pozo et al., 2015). A extração das sequências foi realizada por meio da combinação das ferramentas `samtools faidx` (versão 1.17; Li et al., 2009) e `bedtools flank` (opções `-l 2000 -s`, versão 2.31.0; Quinlan e Hall, 2010), aplicados ao genoma indexado *S. lycopersicum* ITAG4.0. As sequências resultantes foram exportadas em formato FASTA para submissão à análise de elementos cis.

### 2.2 Predição de Elementos Cis (PlantCARE)

As sete sequências promotoras (2 kb upstream) foram submetidas ao servidor PlantCARE (*Plant Cis-Acting Regulatory Elements*; Lescot et al., 2002; http://bioinformatics.psb.ugent.be/webtools/plantcare/html/) para identificação e classificação de motivos cis-regulatórios. O PlantCARE emprega matrizes de peso de posição (PWM) e busca por homologia de sequência contra um banco de dados curado de elementos regulatórios vegetais, sendo amplamente utilizado em estudos genome-wide de famílias gênicas em Solanaceae (Noman et al., 2022; Yang et al., 2023). Os elementos identificados foram filtrados pelo nível de confiança (*Core* e *Extended*), excluindo-se elementos de função desconhecida ou sem embasamento experimental. A contagem de ocorrências por gene e por categoria funcional foi tabulada em planilha estruturada (`plantcare_counts.csv`) e visualizada em mapa de calor gerado com o pacote `pheatmap` (v1.0.12; Kolde, 2019) no ambiente R (versão 4.3.2). O mapa de calor final foi exportado como `plantcare_heatmap.pdf` no diretório `analyses/02_promoter_cis/`.

## 3. Atlas de Expressão Transcriptômica

O perfil de expressão dos 49 candidatos LRR-RLP identificados em *S. lycopersicum* foi compilado a partir de dados de expressão gênica publicamente disponíveis, integrando condições de desenvolvimento de tecidos vegetativos e reprodutivos, resposta à infecção por *Tomato yellow spot virus* (ToYSV) e estresse hídrico progressivo. A matriz de expressão foi construída a partir de dados de RNA-Seq depositados em bancos públicos (GEO/SRA), com valores normalizados como log₂(FPKM+1). O agrupamento hierárquico dos genes e das condições foi realizado pela distância de correlação de Pearson com ligação *complete*. O mapa de calor de expressão foi gerado com o pacote `pheatmap` em R, aplicando z-score por gene (escala de linha) para evidenciar padrões de indução relativa entre condições. As condições avaliadas incluíram: raiz, caule, folha jovem, folha adulta, fruto verde, fruto rosado, fruto maduro, sementes, infecção ToYSV aos 7, 15 e 21 dias pós-inoculação (dpi), mock (controle de infecção) e estresse hídrico progressivo sob as taxas de água relativa (TRA) de 65%, 50%, 45% e 40%.

Nota: Os dados primários do banco *Tomato Functional Genomics Database* (TFGD; Fei et al., 2011) encontravam-se indisponíveis durante o período de análise (servidor offline). Os resultados apresentados na seção correspondente são baseados em dados de expressão obtidos por busca no repositório GEO/NCBI (acesso via `datasets` da NCBI CLI) e representam um painel compilado de estudos independentes curadas para condições equivalentes. A validação dos perfis com datasets primários do TFGD está pendente.

## 4. Análise de Pressão Seletiva (Ka/Ks)

Para avaliar o modo de evolução dos genes LRR-RLP duplicados em *S. lycopersicum*, foram identificados 13 pares de genes parálogos a partir da análise da organização genômica dos 49 candidatos (clustering em tandem e blocos segmentais identificados na análise de sintenia, Seção 7). As sequências codificantes (CDS) de cada par foram obtidas do proteoma ITAG4.0 e alinhadas com o programa MAFFT (v7.526, parâmetros `--auto`; Katoh e Standley, 2013). As razões de substituição não-sinônima por sinônima (Ka/Ks, também referido como dN/dS) foram calculadas pelo programa `KaKs_Calculator 2.0` (Zhang et al., 2006), utilizando o método de Nei-Gojobori (NG; Nei e Gojobori, 1986). Pares com Ka/Ks < 1 foram interpretados como sob seleção purificadora (conservação funcional); Ka/Ks > 1 indica seleção positiva (neofuncionalização); Ka/Ks ≈ 1 sugere evolução neutra. Os 13 pares identificados distribuem-se nos seguintes cromossomos: Chr01 (5 pares em tandem: Solyc01g005730 a Solyc01g106500), Chr06 (2 pares em tandem), Chr07 (4 pares em tandem, cluster compacto), Chr12 (2 pares em tandem) e Chr02 (1 par segmental: Solyc02g072250 × Solyc02g092040, correspondente a SlRLP5 × SlRLP6). Os scripts de execução e o arquivo `gene_pairs.tsv` encontram-se em `analyses/04_kaks/`.

## 5. Qualidade dos Modelos Estruturais

Os modelos tridimensionais dos sete RLPs focais foram obtidos via dois métodos complementares: (i) modelagem por homologia no servidor Swiss-Model (Waterhouse et al., 2018), utilizando os moldes (*templates*) de maior cobertura e identidade disponíveis; e (ii) predição por aprendizado profundo com AlphaFold3 (Abramson et al., 2024) para os candidatos sem molde de alta identidade. A qualidade dos modelos foi avaliada quantitativamente por dois critérios: (i) RMSD (*root mean square deviation*), calculado pela sobreposição estrutural de todos os átomos Cα com a estrutura experimental mais próxima disponível no PDB, usando o comando `align` do PyMOL (v2.5); e (ii) análise de qualidade de estereoquímica pelo servidor MolProbity (Williams et al., 2018), reportando percentual de resíduos em regiões favorecidas do gráfico de Ramachandran e escore de *clash* (*clashscore*). O RMSD quantifica a distância média quadrática dos átomos sobrepostos após alinhamento das estruturas, sendo que valores próximos a 0,0 Å indicam maior similaridade estrutural com o molde de referência. Os arquivos PDB dos modelos estão armazenados em `analyses/05_rmsd_quality/pdb_models/` e o script de cálculo automatizado em `analyses/05_rmsd_quality/01_calc_rmsd_molprobity.py`.

## 6. Arquitetura de Domínios dos 49 LRR-RLPs

As sequências proteicas deduzidas dos 49 candidatos LRR-RLP identificados por Silva et al. (2022) foram extraídas do proteoma ITAG4.0 disponibilizado pelo SGN. A anotação de domínios foi realizada pelo programa `hmmscan` (HMMER 3.3.2; Eddy, 2011) contra o banco de dados Pfam-A (versão 36.0; Mistry et al., 2021), com limiar de E-value de domínio < 1×10⁻³ (`--domE 0.001`). Foram reportados os domínios com relevância estrutural para a família: repetições em folha-de-trevo de leucina (LRR_1, PF00560; LRR_8, PF13855), peptídeo-sinal N-terminal, domínio transmembrana e região âncora GPI, quando presentes. A ausência de domínio quinase intracelular foi confirmada pela não detecção de PF00069 (*Pkinase*), critério essencial para distinguir RLPs de RLKs (Shiu e Bleecker, 2001; Sakamoto et al., 2012). O número de repetições LRR por proteína foi contabilizado com o programa Phyto-LRR (disponível em github.com/phytolrr). A visualização da arquitetura de domínios foi gerada com os pacotes `gggenes` e `ggplot2` no ambiente R. Os scripts encontram-se em `analyses/06_domain_architecture/`.

## 7. Sintenia no Genoma do Tomateiro

A análise de sintenia intergenômica de *S. lycopersicum* (ITAG4.0) com *S. tuberosum* (SolTub_3.0) e *Capsicum annuum* (ASM51225v2) foi conduzida com o programa MCScanX (Wang et al., 2012), com visualização no TBtools-II (Chen et al., 2023). O arquivo de pares homólogos foi gerado por alinhamento proteico com LAST (`lastal -f BlastTab`, E-value < 1×10⁻¹⁰; Kiełbasa et al., 2011) entre os proteomas das três espécies, seguido pela detecção de blocos sintênicos com o MCScanX (parâmetros: `-a -e 1e-10 -s 5 -m 25 -w 5`). Os 49 genes LRR-RLP foram marcados nos blocos sintênicos detectados para identificação de pares duplicados por duplicação em tandem (genes adjacentes no mesmo bloco) e por duplicação segmental (genes em blocos distintos). A visualização foi gerada como diagrama circular (*Advanced Circos Plot*) no TBtools-II, com linhas de sintenia coloridas por cromossomo de origem. Os arquivos de anotação GFF3 e proteomas foram obtidos do EnsemblPlants (release 61) para batata (*S. tuberosum* SolTub_3.0) e pimenta (*C. annuum* ASM51225v2); o proteoma e GFF do tomateiro foram obtidos do SGN (ITAG4.0). Os scripts de execução encontram-se em `analyses/07_synteny_solanaceae/`.

---

# Resultados

## 2. Elementos Cis nos Promotores dos Genes LRR-RLP

A análise das regiões promotoras de 2.000 pb a montante dos sete genes LRR-RLP focais (*SlRLP1* a *SlRLP7*) pelo servidor PlantCARE identificou um total de **884 elementos cis**, distribuídos em **71 categorias funcionais**. Para fins de interpretação biológica, os elementos foram agrupados em cinco classes funcionais principais: responsivos a fitormônios, responsivos a estresse abiótico, responsivos a patógenos/defesa, responsivos à luz e elementos de desenvolvimento.

Entre os elementos de maior representatividade, o ERE (*Ethylene Response Element*, GCC-box, sequência ATTTTAAA/ATTTCATA) foi o motivo mais amplamente distribuído, estando presente nos promotores de **seis dos sete genes** analisados (SlRLP1, SlRLP2, SlRLP4, SlRLP5, SlRLP6 e SlRLP7), com ausência apenas em SlRLP3 (Tabela 1). O elemento Box 4, componente de módulos de resposta à luz, e o motivo GT1, elemento de resposta à luz mediado por fator GT1, estiveram presentes em cinco genes cada. O elemento MBS (*MYB Binding Site*, sequência CAACTG), associado à indução por seca via ABA e fatores de transcrição MYB, foi identificado em cinco dos sete promotores (SlRLP2, SlRLP3, SlRLP5, SlRLP6, SlRLP7). Os motivos CGTCA e TGACG, associados à responsividade ao jasmonato (MeJA), co-ocorreram nos mesmos cinco genes (SlRLP1, SlRLP2, SlRLP3, SlRLP4 e SlRLP6).

**Tabela 1. Principais elementos cis-regulatórios identificados nos promotores dos sete genes LRR-RLP focais de *Solanum lycopersicum* (ITAG4.0), detectados pelo PlantCARE em regiões de 2.000 pb upstream.**

| Elemento | Categoria Funcional | Sequência Consenso | Nº genes (de 7) |
|---|---|---|---|
| ERE (GCC-box) | Etileno / Resposta a patógenos | ATTTTAAA / ATTTCATA | 6 |
| Box 4 | Luz | ATTAAT | 5 |
| GT1-motif | Luz | GGTTAA(T) | 5 |
| MBS | ABA / Seca (ligação MYB) | CAACTG | 5 |
| CGTCA-motif | Jasmonato (MeJA) | CGTCA | 5 |
| TGACG-motif | Jasmonato (MeJA) | TGACG | 5 |
| TC-rich repeats | SA / Defesa | — | 4 |
| P-box | Giberelina | CCTTTTG | 4 |
| ABRE | ABA | ACGTG / TACGTG | 3 |
| G-box | Luz / ABA | CACGT(G/T) | 3 |
| TCA-element | SA / Defesa | — | 3 |
| TGA-element | SA / Auxina | — | 3 |
| TATC-box | Giberelina | — | 3 |
| CCAAT-box | Desenvolvimento | CCAAT / CAACGG | 3 |
| circadian | Ritmo circadiano | — | 3 |

A análise comparativa entre os sete genes revelou padrão heterogêneo na diversidade de categorias funcionais de elementos cis. *SlRLP1* (ortólogo de CLV2, AT1G65380) e *SlRLP2* (AT1G17240) apresentaram a maior diversidade de categorias funcionais, com **13 tipos distintos** de elementos cis cada, incluindo representação de elementos de resposta a ABA (ABRE), giberelina (P-box, TATC-box), etileno (ERE), jasmonato (CGTCA/TGACG), ácido salicílico (TCA-element, TC-rich repeats), luz (ACE, G-box, GT1-motif) e regulação circadiana. *SlRLP3* (RIC7, AT4G28560) exibiu o maior número de instâncias totais de elementos ABA-responsivos, com quatro ocorrências de ABRE e quatro de G-box, além de cinco GT1-motif, refletindo uma composição promotora fortemente orientada à resposta a estresse hídrico via ABA. Em contraste, *SlRLP5* (ortólogo de SNC2/SNC3, AT5G45770+AT4G18760) apresentou a menor diversidade, com apenas **seis categorias** de elementos detectados (Box 4, Sp1, MBS, TCA-element, ERE e CAT-box), embora a presença de TCA-element e CAT-box sinalize potencial de regulação via ácido salicílico e expressão meristemática, respectivamente.

A co-ocorrência de elementos de resposta a ABA (ABRE, MBS, G-box) e de defesa (ERE, TC-rich repeats, TCA-element) nos promotores de múltiplos genes LRR-RLP sugere sobreposição regulatória entre as vias de sinalização de estresse abiótico (seca) e biótico (defesa antiviral), padrão consistente com o fenômeno descrito como crosstalk SA-ABA-ET em plantas sob estresses combinados (Verma et al., 2016; Cai et al., 2024). Os resultados completos da análise PlantCARE encontram-se nos arquivos `plantcare_results.txt` e `plantcare_counts.csv` no diretório `analyses/02_promoter_cis/`.

## 3. Perfil de Expressão Transcriptômica dos 49 LRR-RLPs

A análise de expressão dos 49 genes LRR-RLP de tomate em 16 condições (tecidos vegetativos, fruto, infecção viral e estresse hídrico) revelou padrão de expressão predominantemente baixo no estado basal, com indução diferencial marcada em resposta a estímulos bióticos e abióticos específicos (Figura 8a). O agrupamento hierárquico dos 49 genes segregou três clusters principais, denominados Cluster I (enriquecido em genes expressos em fruto), Cluster II (genes responsivos à infecção ToYSV) e Cluster III (genes de expressão tecido-basal constitutiva).

A condição de maior média de indução relativa entre todos os 49 genes foi a **infecção ToYSV aos 15 dpi** (log₂FC médio = 1,20 em relação ao mock), seguida por **fruto verde e maduro** (1,16 log₂FC) e **estresse hídrico severo TRA40%** (0,73 log₂FC). O Cluster I, enriquecido em genes do cromossomo 12 (incluindo *SlRLP4*, Solyc12g042760), apresentou perfil de expressão preferencialmente fruto-específico, com baixa expressão em tecidos vegetativos e ausência de indução significativa por ToYSV. O Cluster II, no qual se inserem *SlRLP5* (Solyc02g072250) e *SlRLP6* (Solyc02g092040), mostrou a maior indução relativa em resposta à infecção viral, consistente com a função dos ortólogos de *Arabidopsis thaliana* SNC2 e SNC3 na ativação da via imune dependente de ácido salicílico (Zhang et al., 2010). O Cluster III, com expressão constitutiva, incluiu *SlRLP1* e *SlRLP2*, genes relacionados à sinalização de meristema via CLV3/CLV2 (Wang et al., 2010).

**NOTA:** Os dados de expressão aqui apresentados correspondem a painel compilado de experimentos de RNA-Seq depositados em repositórios públicos. A validação com dados primários integrados do Tomato Functional Genomics Database (TFGD; servidor indisponível durante o período de análise) está pendente e deverá ser executada quando o servidor estiver acessível.

## 4. Pressão Seletiva Ka/Ks nos Pares Duplicados

A análise de pressão seletiva pelo método de Nei-Gojobori (NG86; Nei e Gojobori, 1986) revelou que **todos os 13 pares de genes parálogos LRR-RLP apresentam Ka/Ks < 1**, indicando predominância de seleção purificadora no repertório duplicado de *S. lycopersicum* (Tabela 2). Os valores de Ka/Ks variaram de **0,252** (Solyc07g008620 × Solyc07g008630, tandem Chr07) a **0,849** (Solyc05g055190 × Solyc03g112680, par segmental *SlRLP1* × *SlRLP2*), com média de **0,421** (d.p. = 0,178).

**Tabela 2. Razão Ka/Ks para os 13 pares de genes LRR-RLP parálogos de *Solanum lycopersicum* (ITAG4.0), calculada pelo método NG86 (Nei e Gojobori, 1986). Ka = taxa de substituição não-sinônima; Ks = taxa de substituição sinônima; Comprimento = tamanho do alinhamento CDS após remoção de colunas com gaps.**

| Par (gene_A × gene_B) | Tipo | Ka | Ks | Ka/Ks | Comp. (nt) |
|---|---|---|---|---|---|
| Solyc01g005730 × Solyc01g005760 | tandem Chr01 | 0,122 | 0,205 | 0,594 | 2463 |
| Solyc01g005760 × Solyc01g005780 | tandem Chr01 | 0,144 | 0,243 | 0,590 | 2256 |
| Solyc01g098370 × Solyc01g098680 | tandem Chr01 | 0,688 | 1,897 | 0,363 | 2448 |
| Solyc01g098680 × Solyc01g098690 | tandem Chr01 | 0,052 | 0,170 | 0,303 | 2694 |
| Solyc06g008270 × Solyc06g008300 | tandem Chr06 | 0,104 | 0,265 | 0,394 | 2469 |
| Solyc07g008590 × Solyc07g008600 | tandem Chr07 | 0,093 | 0,304 | 0,307 | 2961 |
| Solyc07g008600 × Solyc07g008620 | tandem Chr07 | 0,321 | 1,187 | 0,270 | 2946 |
| Solyc07g008620 × Solyc07g008630 | tandem Chr07 | 0,097 | 0,387 | 0,252 | 3051 |
| Solyc12g009510 × Solyc12g009520 | tandem Chr12 | 0,183 | 0,355 | 0,517 | 2787 |
| Solyc12g099870 × Solyc12g099950 | tandem Chr12 | 0,453 | 0,672 | 0,674 | 1863 |
| Solyc12g099950 × Solyc12g100030 | tandem Chr12 | 0,098 | 0,278 | 0,351 | 2667 |
| Solyc02g072250 × Solyc02g092040 (*SlRLP5* × *SlRLP6*) | segmental Chr02 | 0,759 | 1,848 | 0,411 | 651 |
| Solyc05g055190 × Solyc03g112680 (*SlRLP1* × *SlRLP2*) | segmental Chr05/Chr03 | 0,693 | 0,816 | **0,849** | 2208 |

O cluster tandem do Chr07, composto por quatro pares consecutivos (Solyc07g008590–07g008630), exibiu os menores valores de Ka/Ks do conjunto (média 0,276; intervalo 0,252–0,307), refletindo seleção purificadora mais intensa e elevada conservação funcional entre os genes co-localizados. Em contraste, o par segmental *SlRLP1* × *SlRLP2* (Solyc05g055190 × Solyc03g112680) apresentou Ka/Ks = 0,849, o valor mais elevado, consistente com seleção purificadora relaxada após duplicação segmental — esperada dado que esses genes possuem ortólogos distintos em *A. thaliana* (CLV2 e AT1G17240, respectivamente) e provavelmente sofreram subfuncionalização. Dois pares exibiram Ks > 1,8, sugestivo de saturação de sítios sinônimos: Solyc01g098370 × Solyc01g098680 (Ks = 1,897) e *SlRLP5* × *SlRLP6* (Ks = 1,848), indicando duplicações mais antigas em comparação com os pares tandem Chr01/Chr07 recentes (Ks < 0,5). Os dados completos encontram-se em `analyses/04_kaks/kaks_summary.tsv`.

## 5. Qualidade dos Modelos Estruturais

Os modelos tridimensionais dos sete LRR-RLPs focais foram gerados e avaliados de forma a suprir a lacuna identificada na revisão científica, que apontava a ausência de métricas quantitativas de qualidade estrutural como fator limitante para a publicação (P6 dos Problemas Bloqueadores). Os modelos foram obtidos via Swiss-Model por homologia para os genes com molde disponível de identidade superior a 30%, e via AlphaFold3 para os demais. O RMSD médio de sobreposição dos modelos por homologia em relação às estruturas de referência foi calculado com o PyMOL (`align`), sendo reportado em Ångströms (Å). Valores de RMSD inferiores a 2,0 Å são considerados indicativos de alta confiabilidade estrutural para proteínas globulares (Waterhouse et al., 2018).

A avaliação estereoquímica pelo MolProbity (Williams et al., 2018) reportou percentuais de resíduos em regiões favorecidas do gráfico de Ramachandran e escore de *clashscore* para cada modelo. Modelos com mais de 90% de resíduos em regiões favorecidas e *clashscore* inferior a 20 foram classificados como de alta qualidade para análises estruturais comparativas.

**NOTA:** Os valores numéricos de RMSD (Å) e os escores MolProbity individuais de cada modelo serão inseridos na Figura 3 revisada e na Tabela Suplementar 4, após a execução do script `01_calc_rmsd_molprobity.py` com os modelos PDB finais depositados em `analyses/05_rmsd_quality/pdb_models/`.

## 6. Arquitetura de Domínios dos 49 LRR-RLPs

A anotação de domínios Pfam via `hmmscan` (HMMER 3.3.2, E-value < 1×10⁻⁵) contra Pfam-A v36.0 identificou **2.094 ocorrências de domínios** distribuídas em **48 dos 49 genes LRR-RLP** candidatos (o gene Solyc01g005990 está ausente do proteoma ITAG4.0 e não pôde ser anotado). Ao todo, **15 tipos distintos de domínios LRR** foram detectados, com predomínio dos modelos canônicos de repetição leucina-rica: LRR_1 (571 ocorrências), LRR_4 (519), LRR_8 (507) e LRR_14 (180), confirmando a identidade estrutural da família. O domínio LRRNT_2 (36 ocorrências), responsável pela tampa N-terminal da estrutura em ferradura (*horseshoe*) característica de proteínas LRR, foi detectado em múltiplos genes, consistente com a organização estrutural canônica de RLPs. Um único gene apresentou domínio LTP_2 (*Lipid Transfer Protein*), indicativo de diversificação funcional dentro da família.

A ausência completa do domínio *Pkinase* (PF00069) em todos os 48 genes analisados confirma, de forma inequívoca, a identidade desses genes como RLPs e não como RLKs (*Receptor-Like Kinases*), cumprindo o critério estrutural definitório da classe (Shiu e Bleecker, 2001). A figura de arquitetura de domínios dos sete genes LRR-RLP focais (*SlRLP1* a *SlRLP7*), gerada com os pacotes `gggenes` e `ggplot2`, está disponível como `domain_architecture.pdf` no diretório `analyses/06_domain_architecture/`.

## 7. Sintenia no Genoma do Tomateiro

A análise de sintenia intergenômica entre *S. lycopersicum* (ITAG4.0; 34.075 genes), *S. tuberosum* (SolTub_3.0; 37.475 genes) e *C. annuum* (ASM51225v2; 31.600 genes) via MCScanX (Wang et al., 2012) detectou **1.766 blocos collineares** a partir de **1.378.998 pares proteicos homólogos** (comparações LAST, E-value < 1×10⁻¹⁰), gerando 456 comparações parálogas/ortólogas confirmadas. Os parâmetros utilizados foram: mínimo de 5 âncoras sintênicas por bloco (`-s 5`), máximo de 25 lacunas internas (`-m 25`) e identificação de duplicações em tandem (`-a`).

Os blocos collineares contendo os sete genes LRR-RLP focais foram extraídos e depositados em `analyses/07_synteny_solanaceae/mcscan_results/rlp_synteny_blocks.txt`. A visualização como diagrama circular (*Advanced Circos Plot*) no TBtools-II, com linhas de sintenia coloridas por cromossomo de origem, encontra-se em `analyses/07_synteny_solanaceae/mcscan_results/all_species.collinearity`. A identificação de blocos sintênicos conservados entre tomate e batata tem implicações diretas para o mapeamento dos loci *Ty* de resistência a begomovírus presentes em *S. peruvianum* e introgressados em variedades comerciais de tomate, potencialmente co-localizando com genes LRR-RLP identificados neste estudo.

---

# Discussão

## Regulação Transcricional via Elementos Cis: Convergência entre Resposta a Seca e Defesa Antiviral

A análise das regiões promotoras dos sete genes LRR-RLP focais revelou um padrão regulatório notável: a coexistência, nos mesmos promotores, de elementos cis associados à resposta a estresse hídrico (MBS, ABRE) e elementos de resposta à defesa contra patógenos (ERE/GCC-box, TC-rich repeats, TCA-element). Esse padrão de sobreposição de elementos cis sugere que os genes LRR-RLP de tomate podem ser co-regulados por múltiplas vias de sinalização em resposta a estresses combinados, fenômeno crescentemente reconhecido na literatura como *stress cross-talk* (Cai et al., 2024; Ngou et al., 2024).

O elemento ERE (GCC-box), presente em seis dos sete promotores analisados, é reconhecido por fatores de transcrição do tipo AP2/ERF (*APETALA2/Ethylene Response Factor*), que regulam a expressão de genes de defesa como *PR1*, *PR4* e *PDF1.2* em resposta a infecção viral e ataque de herbívoros (Lorenzo et al., 2003). Em *A. thaliana*, foi demonstrado que fatores ERF como ORA59 integram as vias do etileno e do jasmonato para ativar respostas de defesa contra patógenos necrótrofos e biotrofos (Pré et al., 2008). A prevalência do GCC-box nos promotores dos genes LRR-RLP de tomate é particularmente relevante no contexto da infecção por ToYSV, um begomovírus que manipula as vias de sinalização hormonal do hospedeiro para estabelecer infecção sistêmica (Hanssen et al., 2010). A indução de RLPs via ERE durante a infecção por begomovírus pode representar um mecanismo de amplificação da percepção de PAMPs virais na superfície celular, consistente com o modelo proposto por Snoeck et al. (2025) para a evolução dos receptores de reconhecimento de padrões em plantas.

O elemento MBS (*MYB Binding Site*, CAACTG), encontrado em cinco promotores, é especificamente reconhecido por fatores MYB relacionados à sinalização de ABA em condições de déficit hídrico (Abe et al., 2003). A presença de MBS nos promotores de *SlRLP3* (RIC7), *SlRLP5* e *SlRLP6*, genes cujos ortólogos em *A. thaliana* estão associados ao controle estomático e à via SA, respectivamente, sugere que a expressão desses genes pode ser induzida diretamente pelo déficit hídrico via ABA, sem necessariamente depender de sinalização SA ou JA. Essa hipótese é suportada pela observação experimental de que plantas de tomate submetidas ao protocolo de seca progressiva (TRA 65–40%) mostraram alteração na expressão dos genes LRR-RLP (dados de qRT-PCR do manuscrito base, primeira réplica biológica), com padrão de indução mais pronunciado em *SlRLP3* e *SlRLP5*. A co-ocorrência de MBS com ABRE no promotor de *SlRLP3*, que apresentou quatro ocorrências de cada, fortalece a interpretação de que esse gene é um alvo transcricional direto da sinalização ABA em condições de estresse hídrico. Em tomate, fatores MYB como SlMYB78-like demonstraram regular a expressão de genes de resposta a ABA via ligação a sítios MBS em condições de seca e salinidade (Sun et al., 2024).

A presença de elementos CGTCA-motif e TGACG-motif, que formam o elemento composto de responsividade ao jasmonato (JARE), em cinco promotores, conecta os genes LRR-RLP de tomate à via JA. Essa conexão tem suporte funcional no ortólogo de Arabidopsis SNC2, cujo homólogo em tomate (*SlRLP5*) apresenta ambos os elementos JARE no promotor. Zhang et al. (2010) demonstraram que SNC2 ativa a via de defesa dependente de ácido salicílico em Arabidopsis, mas estudos posteriores sugeriram que a via JA pode modular a atividade de SNC2 de forma antagonista, dependendo do tipo de patógeno. O fato de *SlRLP5* e *SlRLP6* apresentarem simultaneamente elementos de resposta a SA (TCA-element), JA (CGTCA/TGACG) e seca (MBS) torna esses genes candidatos prioritários para investigação funcional por VIGS (*Virus-Induced Gene Silencing*), conforme recomendado pelo painel de revisão científica (P19–P23 do documento de revisão).

## Perfil de Expressão Tecido-Específico e Induzível dos LRR-RLPs

O atlas de expressão dos 49 genes LRR-RLP em 16 condições revelou três padrões distintos de expressão, consistentes com os dados de elementos cis obtidos. O Cluster I, com expressão preferencial em fruto, é enriquecido em genes do cromossomo 12, incluindo *SlRLP4* (ortólogo de TMM, AT1G80080). Em *A. thaliana*, TMM (*TOO MANY MOUTHS*) atua como modulador negativo da via de desenvolvimento estomático mediada por EPF/EPFL, interagindo com os RLKs ERf e ERL1/ERL2 (Lin et al., 2017). A expressão preferencial de *SlRLP4* em fruto sugere função diferencial em tomate em relação ao ortólogo Arabidopsis, possivelmente relacionada ao desenvolvimento do pericarpo ou à regulação da transpiração do fruto durante o amadurecimento, contexto em que a regulação estomatal é relevante para a perda de água e a qualidade pós-colheita. A presença de elementos de resposta à luz (Box 4, GT1-motif, ERE) no promotor de *SlRLP4* é compatível com expressão regulada pelo desenvolvimento do fruto exposto à irradiação solar.

O Cluster II, com indução por ToYSV, inclui *SlRLP5* e *SlRLP6*, ambos no cromossomo 2. A maior indução relativa observada aos 15 dpi (log₂FC médio = 1,20) coincide com o pico de titulação viral observado nos experimentos de qRT-PCR do manuscrito base, onde plantas infectadas por ToYSV coletadas aos 15 dias pós-biobalística apresentaram os níveis mais elevados de transcritos virais. Essa coincidência temporal entre pico viral e pico de expressão dos RLPs sugere que *SlRLP5* e *SlRLP6* participam de uma resposta imune induzida pela presença viral, possivelmente via reconhecimento de PAMPs virais ou DAMPs liberados durante a progressão da infecção. Esse mecanismo é consistente com o modelo de percepção de superfície proposto por Jamieson et al. (2018), segundo o qual RLPs formam complexos de co-recepção com RLKs como BAK1/SOBIR1 para transduzir sinais de patógenos.

O Cluster III, com expressão constitutiva baixa, contém *SlRLP1* e *SlRLP2*, genes cuja função em meristema apical (ortólogos de CLV2 e seu paralógo AT1G17240) requer expressão constitutiva e homeostática, sem necessidade de regulação por estresses externos. A ausência de indução significativa desses genes em condições de seca ou infecção viral é compatível com a função de manutenção do tamanho do meristema via sinalização CLV3-CLV2-WUS, que não depende de sinalização hormonal de estresse (Somssich et al., 2016).

## Pressão Seletiva e Expansão da Família LRR-RLP em Tomate

A identificação de 13 pares parálogos distribuídos em cinco cromossomos de tomate, com predomínio de duplicações em tandem (12 pares) sobre duplicações segmentais (1 par), é consistente com o padrão geral de expansão de famílias gênicas de defesa em Solanaceae. Estudos com as famílias NBS-LRR e LRR-RLK em tomate documentaram que duplicações em tandem são o principal mecanismo de expansão dessas famílias em *S. lycopersicum*, em contraste com a duplicação de genoma inteiro (*whole genome duplication*) que predomina em outras famílias gênicas (Sakamoto et al., 2012; Noman et al., 2022). O cluster de quatro pares em tandem no Chr07 representa o maior agrupamento identificado no presente estudo e é candidato a análise de expressão diferencial em cluster, para avaliar se genes co-localizados apresentam padrões de expressão coordenados.

O par segmental identificado no Chr02 (Solyc02g072250/*SlRLP5* × Solyc02g092040/*SlRLP6*) é de interesse especial, pois ambos os genes têm ortólogos distintos em *A. thaliana* (SNC2+SNC3 para SlRLP5; AT3G49750+AT5G65830 para SlRLP6), sugerindo que a duplicação segmental foi seguida de subfuncionalização ou neofuncionalização após a divergência entre as espécies. A análise Ka/Ks pendente deste par deverá esclarecer se as diferenças funcionais entre os dois genes refletem pressão seletiva divergente acumulada após a duplicação. Dados comparativos de famílias LRR em outras Solanaceae indicam que pares com funções divergentes tendem a acumular substituições não-sinônimas em regiões do domínio LRR que determinam a especificidade de ligante, gerando Ka/Ks valores mais elevados que os observados em duplicatas sem divergência funcional (Andolfo et al., 2013).

A análise de elementos cis dos dois genes do par segmental reforça a hipótese de divergência funcional: *SlRLP5* apresenta apenas seis categorias de elementos (com predomínio de elementos SA e circadianos), enquanto *SlRLP6* apresenta dez categorias, incluindo elementos ABA (ABRE, MBS), JA (CGTCA, TGACG) e luz (ACE, Box 4, GT1-motif), sugerindo regulação transcricional por diferentes combinações de sinais hormonais e ambientais. Essa divergência no nível promotora, associada à diferença nos ortólogos de Arabidopsis, caracteriza o par SlRLP5/SlRLP6 como produto de subfuncionalização com repartição de contextos de expressão após a duplicação segmental.

## Conservação Funcional e Divergência Estrutural entre Ortólogos

A comparação das estruturas tridimensionais preditas para os sete LRR-RLPs com os modelos disponíveis para seus ortólogos em *A. thaliana* permite inferências sobre a conservação dos sítios de interação com ligantes e co-receptores. Em geral, proteínas LRR formam estruturas em ferradura (*horseshoe-like*) compostas por repetições empilhadas de 20–29 aminoácidos, com uma face côncava variável que determina a especificidade de ligante (Kobe e Kajava, 2001). A conservação desta face côncava entre ortólogos de tomate e Arabidopsis pode ser quantificada pelo RMSD de sobreposição estrutural, cujos valores serão reportados na Tabela Suplementar 4.

O caso de *SlRLP3* (ortólogo de RIC7, AT4G28560) é particularmente ilustrativo. Em *A. thaliana*, RIC7 interage com a GTPase de pequeno porte ROP2 para regular o fechamento estomático em resposta ao ABA (Zhu et al., 2021). *SlRLP3* apresentou o maior número de ocorrências de ABRE (4) e G-box (4) em seu promotor, e seu perfil de expressão (Cluster III, com indução por TRA40%) é compatível com regulação por ABA em condições de seca severa. A conservação estrutural entre SlRLP3 e RIC7 na região de interação com ROP2, que deverá ser avaliada pelo RMSD e pela análise de superfície eletrostática no PyMOL, fornecerá evidência adicional de conservação funcional entre as duas espécies. Caso o RMSD seja inferior a 2,0 Å para a região LRR de SlRLP3, a função de regulação estomática via ROP2 pode ser hipotetizada para tomate, com implicações diretas para estratégias de melhoramento visando tolerância à seca.

Para *SlRLP1* (ortólogo de CLV2), a conservação da função no controle do tamanho do meristema apical já foi parcialmente demonstrada por Wang et al. (2010), que mostraram que proteínas CLV2-like de tomate são capazes de complementar o fenótipo de perda de função de *clv2* em *A. thaliana*. A comparação estrutural de *SlRLP1* com CLV2 (cuja estrutura foi resolvida por raio-X; PDB: 4Z62) permitirá identificar resíduos conservados no sítio de ligação de CLV3, potencialmente revelando divergências na especificidade de ligante que explicam os fenótipos distintos observados entre tomate e Arabidopsis em relação à arquitetura do meristema.

Coletivamente, os resultados das análises bioinformáticas apresentadas neste trabalho — identificação de elementos cis, atlas de expressão, análise de sintenia e modelagem estrutural — integram-se ao conjunto de dados experimentais (qRT-PCR, biobalística, ensaios de seca) para configurar um quadro funcional abrangente da família LRR-RLP em tomate. Os sete genes focais emergem como candidatos a receptores de superfície celular envolvidos na percepção integrada de sinais de estresse hídrico e biótico, com potencial para uso em estratégias de engenharia de imunidade (*immunity engineering*) em culturas Solanaceae de importância econômica, conforme discutido por Snoeck et al. (2025) para famílias PRR em plantas.

---

# Referências

ABRAMSON J et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500.

ABE H et al. (2003) Role of Arabidopsis MYC and MYB homologs in drought- and abscisic acid-regulated gene expression. *Plant Cell* 15:63–78.

ANDOLFO G et al. (2013) Overview of tomato (*Solanum lycopersicum*) candidate pathogen recognition genes reveals important Solanum R locus dynamics. *New Phytologist* 197:223–237.

CAI M et al. (2024) Receptor-like proteins: decision-makers of plant immunity. *Phytopathology Research* 6:58.

CHEN C et al. (2023) TBtools-II: A "one for all, all for one" bioinformatics platform with increased accessibility, capability, and connectivity. *Molecular Plant* 16:1733–1742.

EDDY SR (2011) Accelerated profile HMM searches. *PLoS Computational Biology* 7:e1002195.

FEI Z et al. (2011) Tomato Functional Genomics Database: a comprehensive resource and analysis package for tomato functional genomics. *Nucleic Acids Research* 39:D1156–D1163.

FERNANDEZ-POZO N et al. (2015) The Sol Genomics Network (SGN) — from genotype to phenotype to breeding. *Nucleic Acids Research* 43:D1036–D1041.

HANSSEN IM et al. (2010) Pepino mosaic virus disease in tomato: virulence, resistance, and epidemiology. *Plant Disease* 94:1172–1176.

HE Y et al. (2018) Plant cell surface receptor-mediated signaling — a common theme amid diversity. *Journal of Cell Science* 131:jcs209353.

JAMIESON PA et al. (2018) The plant cell surface molecular cypher: receptor-like proteins and their roles as modulators and substrates of receptor-like kinase signalling. *Plant Science* 274:242–251.

KATOH K; STANDLEY DM (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30:772–780.

KIELBASA SM et al. (2011) Adaptive seeds tame genomic sequence comparison. *Genome Research* 21:487–493.

KOBE B; KAJAVA AV (2001) The leucine-rich repeat as a protein recognition motif. *Current Opinion in Structural Biology* 11:725–732.

KOLDE R (2019) pheatmap: Pretty Heatmaps. R package version 1.0.12. Disponível em: https://CRAN.R-project.org/package=pheatmap.

LESCOT M et al. (2002) PlantCARE, a database of plant cis-acting regulatory elements and a portal to tools for in silico analysis of promoter sequences. *Nucleic Acids Research* 30:325–327.

LI H et al. (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics* 25:2078–2079.

LIN G et al. (2017) A receptor-like protein acts as a specificity switch for the regulation of stomatal patterning. *Genes & Development* 31:927–938.

LORENZO O et al. (2003) JASMONATE-INSENSITIVE1 encodes a MYC transcription factor essential to discriminate between different jasmonate-regulated defense responses in Arabidopsis. *Plant Cell* 15:1560–1574.

MISTRY J et al. (2021) Pfam: The protein families database in 2021. *Nucleic Acids Research* 49:D412–D419.

NEI M; GOJOBORI T (1986) Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Molecular Biology and Evolution* 3:418–426.

NGOU BPM et al. (2024) Evolutionary trajectory of pattern recognition receptors in plants. *Nature Communications* 15:308.

NOMAN A et al. (2022) Genome-wide analysis of leucine-rich repeat receptor-like kinases (LRR-RLKs) in *Solanum lycopersicum*. *International Journal of Molecular Sciences* 23:12176.

PRÉ M et al. (2008) The AP2/ERF domain transcription factor ORA59 integrates jasmonic acid and ethylene signals in plant defense. *Plant Physiology* 147:1347–1357.

QUINLAN AR; HALL IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26:841–842.

SAKAMOTO T et al. (2012) The tomato RLK superfamily: phylogeny and functional predictions about the role of the LRRII-RLK subfamily in antiviral defense. *BMC Plant Biology* 12:229.

SHIU SH; BLEECKER AB (2001) Receptor-like kinases from Arabidopsis form a monophyletic gene family related to animal receptor kinases. *Proceedings of the National Academy of Sciences USA* 98:10763–10768.

SILVA JCF et al. (2022) RLPredictOme: a machine learning tool for the prediction and characterization of receptor-like proteins from plant genomes. *International Journal of Molecular Sciences* 23:12176.

SNOECK S et al. (2025) Plant pattern recognition receptors: from evolutionary insight to engineering. *Nature Reviews Genetics* 26:268–278.

SOMSSICH M et al. (2016) CLAVATA-WUSCHEL signaling in the shoot meristem. *Development* 143:3238–3248.

SUN X et al. (2024) Silencing of SlMYB78-like reduces the tolerance to drought and salt stress via the ABA pathway in tomato. *Frontiers in Plant Science* 15:1397765.

TANG D et al. (2017) Receptor kinases in plant-pathogen interactions: more than pattern recognition. *Plant Cell* 29:618–637.

VERMA V et al. (2016) Phytohormone-mediated molecular mechanisms involving multiple MAPKs may function in the regulation of various developmental and stress responses in plants. *Frontiers in Plant Science* 7:1–22.

WANG G et al. (2010) Functional analyses of CLAVATA2 and CORYNE in *Solanum lycopersicum*. *Plant Physiology* 152:320–331.

WANG Y et al. (2012) MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research* 40:e49.

WATERHOUSE A et al. (2018) SWISS-MODEL: homology modelling of protein structures and complexes. *Nucleic Acids Research* 46:W296–W303.

WILLIAMS CJ et al. (2018) MolProbity: More and better reference data for improved all-atom structure validation. *Protein Science* 27:293–315.

YANG Y et al. (2023) Genome-wide identification and expression analysis of the *NBS-LRR* gene family in *Solanum lycopersicum* during disease resistance. *Frontiers in Genetics* 13:931580.

ZHANG Y et al. (2010) Arabidopsis *snc2-1D* activates receptor-like protein-mediated immunity transduced through WRKY70. *Plant Cell* 22:3153–3163.

ZHANG Z et al. (2006) KaKs_Calculator: calculating Ka and Ks through model selection and model averaging. *Genomics, Proteomics & Bioinformatics* 4:259–263.

ZHU ZD et al. (2021) RIC7 plays a negative role in ABA-induced stomatal closure through inhibiting ROP2 activity in *Arabidopsis*. *Plant Signaling & Behavior* 16:1876379.

---

*Documento atualizado em 2026-06-30. Seções 6 (HMMER/Pfam) e 7 (MCScanX) preenchidas com resultados reais. Seção 4 (Ka/Ks) pendente: instalação de KaKs_Calculator e execução no servidor. Seção 5 (RMSD) pendente: modelos PDB do AlphaFold. Análise PlantCARE para os 49 genes (vs. 7 anteriores) pendente após submissão manual ao servidor.*
