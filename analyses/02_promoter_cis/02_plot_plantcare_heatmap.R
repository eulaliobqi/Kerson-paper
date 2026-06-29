#!/usr/bin/env Rscript
# Heatmap de elementos cis nos promotores (2kb) dos 7 genes LRR-RLP
# Uso: Rscript 02_plot_plantcare_heatmap.R [plantcare_results.txt]
#
# Formato do arquivo de entrada (PlantCARE download):
#   Sequence <TAB> Signal <TAB> Location <TAB> Strand <TAB> Sequence <TAB> Function
# Se o arquivo não existir, roda em modo DEMO com dados sintéticos.

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Metadados dos genes ────────────────────────────────────────────────────────
gene_meta <- tibble(
  gene_id = c("Solyc05g055190","Solyc03g112680","Solyc05g009990",
               "Solyc12g042760","Solyc02g072250","Solyc02g092040","Solyc10g007830"),
  label   = c("SlRLP1\n(CLV2)","SlRLP2","SlRLP3\n(RIC7)","SlRLP4\n(TMM)",
               "SlRLP5\n(SNC2/3)","SlRLP6","SlRLP7")
)

# ── Dicionário de elementos cis por categoria funcional ───────────────────────
cis_dict <- list(
  "Luminoso"         = c("ACE","Box 4","G-box","GT1-motif","I-box","Sp1","LAMP-element"),
  "ABA / Seca"       = c("ABRE","CACGTG-motif","DRE core","MBS","MYB recognition"),
  "JA"               = c("CGTCA-motif","TGACG-motif"),
  "SA / Defesa"      = c("SARE","TCA-element","W-box","TC-rich repeats","TGA-element"),
  "GA"               = c("GARE-motif","P-box","TATC-box"),
  "Auxina"           = c("AuxRR-core","TGA-element"),
  "Etileno"          = c("GCC-box","ERE"),
  "Desenvolvimento"  = c("CAT-box","CCAAT-box","as-2-element"),
  "Circadiano"       = c("circadian","Evening Element")
)

cis_df <- tibble(
  motif    = unlist(cis_dict),
  category = rep(names(cis_dict), lengths(cis_dict))
) %>% distinct(motif, .keep_all = TRUE)

# ── Ler resultados do PlantCARE ───────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "plantcare_results.txt"

if (!file.exists(input_file)) {
  message("[DEMO] Arquivo '", input_file, "' não encontrado — usando dados simulados")
  set.seed(2024)
  results <- expand_grid(gene_id = gene_meta$gene_id,
                         motif   = cis_df$motif) %>%
    mutate(count = sample(c(rep(0L,4), 1L, 2L, 3L, 4L), n(), replace = TRUE))
} else {
  raw <- read_tsv(input_file, col_names = TRUE, show_col_types = FALSE)

  # Renomear colunas (PlantCARE usa: Sequence, Signal, Location, Strand, Seq, Function)
  colnames(raw) <- tolower(gsub("\\s+","_", colnames(raw)))
  if (!"sequence" %in% colnames(raw)) stop("Coluna 'Sequence' não encontrada no arquivo.")

  results <- raw %>%
    rename(gene_id = sequence, motif = signal) %>%
    group_by(gene_id, motif) %>%
    summarise(count = n(), .groups = "drop")
}

# ── Construir matriz ──────────────────────────────────────────────────────────
all_motifs <- cis_df$motif

mat_long <- results %>%
  filter(motif %in% all_motifs) %>%
  complete(gene_id = gene_meta$gene_id, motif = all_motifs,
           fill = list(count = 0L)) %>%
  left_join(gene_meta, by = "gene_id") %>%
  arrange(match(gene_id, gene_meta$gene_id))

mat <- mat_long %>%
  pivot_wider(id_cols = label, names_from = motif, values_from = count) %>%
  column_to_rownames("label") %>%
  as.matrix()

# Remover motivos ausentes em todos os genes
mat <- mat[, colSums(mat) > 0, drop = FALSE]

# Ordenar colunas por categoria
col_order <- cis_df %>%
  filter(motif %in% colnames(mat)) %>%
  arrange(match(category, names(cis_dict))) %>%
  pull(motif)
mat <- mat[, col_order, drop = FALSE]

# ── Anotação de colunas (categorias) ─────────────────────────────────────────
col_anno <- data.frame(
  Categoria = cis_df$category[match(colnames(mat), cis_df$motif)],
  row.names = colnames(mat)
)

n_cats   <- length(unique(col_anno$Categoria))
cat_cols <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_cats),
  unique(col_anno$Categoria)
)

# ── Plot ──────────────────────────────────────────────────────────────────────
outfile <- "plantcare_heatmap.pdf"
pdf(outfile, width = max(14, ncol(mat) * 0.55 + 4), height = 6)

pheatmap(
  mat,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  color         = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#C0392B"))(60),
  border_color  = "grey75",
  cellwidth     = 20,
  cellheight    = 24,
  fontsize       = 9,
  fontsize_row   = 10,
  fontsize_col   = 8,
  angle_col      = 45,
  annotation_col = col_anno,
  annotation_colors = list(Categoria = cat_cols),
  main           = "Elementos cis regulatórios em promotores de LRR-RLPs (2 kb upstream)",
  display_numbers = TRUE,
  number_format   = "%d",
  number_color    = "black",
  gaps_col = cumsum(table(col_anno$Categoria)[unique(col_anno$Categoria)])[-n_cats]
)

dev.off()
message("Heatmap salvo: ", outfile)

# ── Resumo ────────────────────────────────────────────────────────────────────
cat("\n=== Top elementos cis (presentes em ≥2 genes) ===\n")
freq <- colSums(mat > 0)
sort(freq[freq >= 2], decreasing = TRUE) %>%
  as.data.frame() %>%
  setNames("n_genes") %>%
  print()

cat("\n=== Elementos por gene ===\n")
rowSums(mat > 0) %>% sort(decreasing = TRUE) %>% print()

# Salvar tabela CSV de apoio para o artigo
write_csv(as.data.frame(mat) %>% rownames_to_column("gene"), "plantcare_counts.csv")
message("Tabela salva: plantcare_counts.csv")
