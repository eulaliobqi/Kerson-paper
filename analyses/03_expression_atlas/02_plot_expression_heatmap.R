#!/usr/bin/env Rscript
# Atlas de expressão pública dos 7 LRR-RLPs em tomate
# Uso: Rscript 02_plot_expression_heatmap.R [expression_matrix.csv]
#
# Entrada: CSV com colunas gene_id, label, <condição1>, <condição2>, ...
# Valores esperados: log2(FPKM+1) ou log2FC normalizado

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "expression_matrix.csv"

if (!file.exists(input_file)) stop("Arquivo não encontrado: ", input_file,
                                    "\nExecute primeiro: python3 01_fetch_expression.py")

raw <- read_csv(input_file, show_col_types = FALSE)

# Condições de interesse (ordem de exibição no heatmap)
priority_conds <- c(
  # Tecidos vegetativos
  "Raiz","Caule","Folha jovem","Folha adulta",
  # Tecidos reprodutivos
  "Flor","Fruto verde","Fruto maduro",
  # Seca
  "Seca TRA65%","Seca TRA50%","Seca TRA45%","Seca TRA40%",
  # Vírus
  "ToYSV 7dpi","ToYSV 15dpi","ToYSV 21dpi",
  # Mock
  "Mock vírus","Mock seca"
)

# ── Construir matriz ──────────────────────────────────────────────────────────
expr_cols <- setdiff(colnames(raw), c("gene_id","label"))

# Filtrar apenas condições prioritárias que existem nos dados
keep_conds <- intersect(priority_conds, expr_cols)
# Adicionar demais colunas não listadas ao final
extra_conds <- setdiff(expr_cols, keep_conds)
ordered_conds <- c(keep_conds, extra_conds)

mat <- raw %>%
  select(label, all_of(ordered_conds)) %>%
  column_to_rownames("label") %>%
  as.matrix()

# Escalar por gene (z-score por linha) para realçar padrões relativos
mat_z <- t(scale(t(mat)))
mat_z[is.nan(mat_z)] <- 0

# ── Anotação de colunas (grupos de condição) ──────────────────────────────────
make_group <- function(cond_names) {
  case_when(
    str_detect(cond_names, "Raiz|Caule|Folha|Flor|Fruto") ~ "Tecido",
    str_detect(cond_names, "Seca|TRA")                     ~ "Seca",
    str_detect(cond_names, "ToYSV|dpi")                    ~ "Vírus",
    str_detect(cond_names, "Mock")                         ~ "Controle",
    TRUE                                                   ~ "Outro"
  )
}

col_anno <- data.frame(
  Condição = make_group(colnames(mat_z)),
  row.names = colnames(mat_z)
)

group_cols <- c(
  "Tecido"    = "#2196F3",
  "Seca"      = "#FF9800",
  "Vírus"     = "#F44336",
  "Controle"  = "#4CAF50",
  "Outro"     = "#9E9E9E"
)

# ── Plot heatmap de expressão ─────────────────────────────────────────────────
outfile <- "expression_atlas_heatmap.pdf"
pdf(outfile, width = max(12, ncol(mat_z) * 0.5 + 4), height = 5)

pheatmap(
  mat_z,
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  color          = colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
  border_color   = "grey80",
  cellwidth      = 22,
  cellheight     = 26,
  fontsize        = 9,
  fontsize_row    = 10,
  fontsize_col    = 8,
  angle_col       = 45,
  annotation_col  = col_anno,
  annotation_colors = list(Condição = group_cols[unique(col_anno$Condição)]),
  main            = "Perfil de expressão dos LRR-RLPs em S. lycopersicum (z-score)",
  gaps_col        = {
    grp <- col_anno$Condição
    cumsum(rle(grp)$lengths)[-length(rle(grp)$lengths)]
  },
  treeheight_row  = 30
)

dev.off()
message("Heatmap salvo: ", outfile)

# ── Plot de expressão absoluta (FPKM) ─────────────────────────────────────────
outfile2 <- "expression_atlas_absolute.pdf"
pdf(outfile2, width = max(12, ncol(mat) * 0.5 + 4), height = 5)

pheatmap(
  mat,
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  color          = colorRampPalette(c("#F7F7F7","#fee090","#d73027"))(80),
  border_color   = "grey80",
  cellwidth      = 22,
  cellheight     = 26,
  fontsize        = 9,
  fontsize_row    = 10,
  fontsize_col    = 8,
  angle_col       = 45,
  annotation_col  = col_anno,
  annotation_colors = list(Condição = group_cols[unique(col_anno$Condição)]),
  main            = "Perfil de expressão dos LRR-RLPs — valores absolutos [log₂(FPKM+1)]",
  gaps_col        = {
    grp <- col_anno$Condição
    cumsum(rle(grp)$lengths)[-length(rle(grp)$lengths)]
  }
)

dev.off()
message("Heatmap absoluto salvo: ", outfile2)

# ── Top condições indutoras ───────────────────────────────────────────────────
cat("\n=== Condições de maior indução (média dos genes) ===\n")
colMeans(mat, na.rm=TRUE) %>% sort(decreasing=TRUE) %>% head(10) %>% print()
