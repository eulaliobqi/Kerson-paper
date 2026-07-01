#!/usr/bin/env Rscript
# Heatmap PlantCARE dos 7 genes LRR-RLP focais — lê plantcare_counts.csv diretamente
# Uso: Rscript 03_plot_focal7_heatmap.R [plantcare_counts.csv]

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

script_dir <- tryCatch(dirname(normalizePath(sys.frame(0)$ofile)), error = function(e) getwd())
setwd(script_dir)

args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "plantcare_counts.csv"

if (!file.exists(input_file)) stop("Arquivo não encontrado: ", input_file)

raw <- read_csv(input_file, show_col_types = FALSE)

# Corrigir nomes de genes com \n embutido (ex: "SlRLP1\n(CLV2)" → "SlRLP1 (CLV2)")
raw[[1]] <- str_replace_all(raw[[1]], "\\s*\\n\\s*", " ")

mat <- raw %>%
  column_to_rownames(colnames(raw)[1]) %>%
  as.matrix()

# Garantir inteiros
storage.mode(mat) <- "integer"

# Ordenar genes na ordem correta (CLV2 primeiro)
gene_order <- c(
  "SlRLP1 (CLV2)", "SlRLP2", "SlRLP3 (RIC7)",
  "SlRLP4 (TMM)", "SlRLP5 (SNC2/3)", "SlRLP6", "SlRLP7"
)
present <- intersect(gene_order, rownames(mat))
mat <- mat[present, , drop = FALSE]

# ── Anotação de colunas por categoria funcional ───────────────────────────────
cis_cat <- c(
  "ACE"              = "Luminoso",
  "Box 4"            = "Luminoso",
  "G-box"            = "Luminoso",
  "GT1-motif"        = "Luminoso",
  "I-box"            = "Luminoso",
  "Sp1"              = "Luminoso",
  "LAMP-element"     = "Luminoso",
  "ABRE"             = "ABA / Seca",
  "CACGTG-motif"     = "ABA / Seca",
  "DRE core"         = "ABA / Seca",
  "MBS"              = "ABA / Seca",
  "MYB recognition"  = "ABA / Seca",
  "CGTCA-motif"      = "JA",
  "TGACG-motif"      = "JA",
  "SARE"             = "SA / Defesa",
  "TCA-element"      = "SA / Defesa",
  "W-box"            = "SA / Defesa",
  "TC-rich repeats"  = "SA / Defesa",
  "TGA-element"      = "SA / Defesa",
  "GARE-motif"       = "GA",
  "P-box"            = "GA",
  "TATC-box"         = "GA",
  "AuxRR-core"       = "Auxina",
  "GCC-box"          = "Etileno",
  "ERE"              = "Etileno",
  "CAT-box"          = "Desenvolvimento",
  "CCAAT-box"        = "Desenvolvimento",
  "as-2-element"     = "Desenvolvimento",
  "circadian"        = "Circadiano",
  "Evening Element"  = "Circadiano"
)

present_cols <- intersect(names(cis_cat), colnames(mat))
mat <- mat[, present_cols, drop = FALSE]

col_anno <- data.frame(
  Categoria = cis_cat[present_cols],
  row.names = present_cols
)

cat_order <- c("Luminoso","ABA / Seca","JA","SA / Defesa","GA","Auxina","Etileno","Desenvolvimento","Circadiano")
present_cats <- intersect(cat_order, unique(col_anno$Categoria))
n_cats <- length(present_cats)

cat_pal <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(min(8, n_cats), "Set2"))(n_cats),
  present_cats
)

# Reordenar colunas por categoria
col_order <- col_anno %>%
  rownames_to_column("motif") %>%
  mutate(cat_f = factor(Categoria, levels = present_cats)) %>%
  arrange(cat_f) %>%
  pull(motif)
mat      <- mat[, col_order, drop = FALSE]
col_anno <- col_anno[col_order, , drop = FALSE]

gaps_col <- cumsum(table(col_anno$Categoria)[present_cats])
gaps_col <- gaps_col[-length(gaps_col)]

# ── Heatmap principal ─────────────────────────────────────────────────────────
outfile <- "plantcare_focal7_heatmap.pdf"
pdf(outfile, width = max(10, ncol(mat) * 0.65 + 4), height = 5)

pheatmap(
  mat,
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  color            = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#C0392B"))(60),
  border_color     = "grey75",
  cellwidth        = 28,
  cellheight       = 30,
  fontsize          = 10,
  fontsize_row      = 11,
  fontsize_col      = 9,
  angle_col         = 45,
  annotation_col    = col_anno,
  annotation_colors = list(Categoria = cat_pal),
  main              = "Elementos cis regulatórios em promotores dos 7 LRR-RLPs focais (2 kb upstream)",
  display_numbers   = TRUE,
  number_format     = "%d",
  number_color      = "black",
  gaps_col          = gaps_col
)

dev.off()
message("Heatmap salvo: ", outfile)

# PNG para inserção no docx/artigo
png_file <- "plantcare_focal7_heatmap.png"
png(png_file, width = max(1200, ncol(mat) * 60 + 400), height = 500, res = 120)
pheatmap(
  mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color        = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#C0392B"))(60),
  border_color = "grey75",
  cellwidth = 28, cellheight = 30,
  fontsize = 10, angle_col = 45,
  annotation_col = col_anno,
  annotation_colors = list(Categoria = cat_pal),
  main = "Elementos cis regulatórios — 7 LRR-RLPs (2 kb upstream)",
  display_numbers = TRUE, number_format = "%d", number_color = "black",
  gaps_col = gaps_col
)
dev.off()
message("PNG salvo: ", png_file)

# ── Estatísticas ───────────────────────────────────────────────────────────────
cat("\n=== Elementos presentes em ≥2 genes ===\n")
freq <- colSums(mat > 0)
sort(freq[freq >= 2], decreasing = TRUE) %>% as.data.frame() %>%
  setNames("n_genes") %>% print()

cat("\n=== Total de elementos por gene ===\n")
rowSums(mat > 0) %>% sort(decreasing = TRUE) %>% print()

cat("\n=== Elementos ABA/Seca e SA/Defesa ===\n")
cat("ABA/Seca:\n")
mat[, col_anno$Categoria == "ABA / Seca", drop=FALSE] %>% print()
cat("SA/Defesa:\n")
mat[, col_anno$Categoria == "SA / Defesa", drop=FALSE] %>% print()
