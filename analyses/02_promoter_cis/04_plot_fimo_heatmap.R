#!/usr/bin/env Rscript
# 04_plot_fimo_heatmap.R — Heatmap TFs (JASPAR2024) × 49 LRR-RLPs
# Uso: Rscript 04_plot_fimo_heatmap.R fimo_parsed_counts.csv [fimo_parsed_tf_families.tsv]
# Saída: fimo_tf_heatmap.pdf (174 mm double-column) + PNG 300 dpi

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

script_dir <- tryCatch(dirname(normalizePath(sys.frame(0)$ofile)), error = function(e) getwd())
setwd(script_dir)

args        <- commandArgs(trailingOnly = TRUE)
counts_file <- if (length(args) > 0) args[1] else "fimo_parsed_counts.csv"
ann_file    <- if (length(args) > 1) args[2] else "fimo_parsed_tf_families.tsv"

# ── Metadados dos 49 genes ────────────────────────────────────────────────────
focal_ids <- c("Solyc05g055190","Solyc03g112680","Solyc05g009990",
               "Solyc12g042760","Solyc02g072250","Solyc02g092040","Solyc10g007830")

gene_meta <- tibble(
  gene_id = c(
    "Solyc01g005730","Solyc01g005760","Solyc01g005780","Solyc01g005990",
    "Solyc01g006550","Solyc01g008390","Solyc01g009690","Solyc01g009700",
    "Solyc01g073680","Solyc01g087510","Solyc01g098370","Solyc01g098680",
    "Solyc01g098690","Solyc01g099250","Solyc01g106500",
    "Solyc02g021770","Solyc02g072250","Solyc02g092040",
    "Solyc03g082780","Solyc03g083510","Solyc03g112680",
    "Solyc04g014400",
    "Solyc05g009990","Solyc05g054900","Solyc05g055190",
    "Solyc06g008270","Solyc06g008300","Solyc06g033920",
    "Solyc07g005150","Solyc07g008590","Solyc07g008600",
    "Solyc07g008620","Solyc07g008630","Solyc07g008640",
    "Solyc08g016270","Solyc08g077740",
    "Solyc09g005090",
    "Solyc10g007830","Solyc10g076500",
    "Solyc11g011180",
    "Solyc12g006020","Solyc12g009510","Solyc12g009520","Solyc12g013680",
    "Solyc12g042760","Solyc12g049190","Solyc12g099870",
    "Solyc12g099950","Solyc12g100030"
  ),
  label = c(
    "SlRLP01","SlRLP02","SlRLP03","SlRLP04","SlRLP05","SlRLP06","SlRLP07","SlRLP08",
    "SlRLP09","SlRLP10","SlRLP11","SlRLP12","SlRLP13","SlRLP14","SlRLP15",
    "SlRLP16","SlRLP17* (SNC2/3)","SlRLP18*",
    "SlRLP19","SlRLP20","SlRLP21*",
    "SlRLP22",
    "SlRLP23* (RIC7)","SlRLP24","SlRLP25* (CLV2)",
    "SlRLP26","SlRLP27","SlRLP28",
    "SlRLP29","SlRLP30","SlRLP31","SlRLP32","SlRLP33","SlRLP34",
    "SlRLP35","SlRLP36",
    "SlRLP37",
    "SlRLP38* (Seca/Vírus)","SlRLP39",
    "SlRLP40",
    "SlRLP41","SlRLP42","SlRLP43","SlRLP44","SlRLP45* (TMM)","SlRLP46",
    "SlRLP47","SlRLP48","SlRLP49"
  )
) %>% mutate(focal = gene_id %in% focal_ids)

# Ordem de categorias para colunas
cat_order <- c("SA / Defesa","JA / Defesa","ABA / Defesa","ABA / Seca",
               "Etileno / Defesa","Luz / Desenvolvimento","Luz / GA",
               "GA / Desenvolvimento","Brassinoesteroides","Circadiano",
               "Desenvolvimento","Meristema","Estresse térmico","Outros")

# ── Ler dados ─────────────────────────────────────────────────────────────────
if (!file.exists(counts_file)) {
  message("[DEMO] fimo_parsed_counts.csv nao encontrado — modo demonstracao")
  set.seed(42)
  demo_tfs <- c("WRKY6","WRKY18","WRKY40","WRKY70","MYC2","MYC3",
                "ABF1","ABF2","MYB2","MYB96","ERF1","ORA59","PIF4","BES1",
                "HY5","NAC16","NAC29","CCA1","HSF1","TCP4")
  mat <- matrix(
    sample(0:4, length(gene_meta$label) * length(demo_tfs), replace = TRUE,
           prob = c(0.50, 0.22, 0.14, 0.09, 0.05)),
    nrow = length(gene_meta$label),
    dimnames = list(gene_meta$label, demo_tfs)
  )
  tf_ann <- tibble(
    tf_name  = demo_tfs,
    family   = c("WRKY","WRKY","WRKY","WRKY","bHLH","bHLH",
                 "bZIP","bZIP","MYB","MYB","AP2/ERF","AP2/ERF",
                 "bHLH","BES1","HY5","NAC","NAC","MYB-related","HSF","TCP"),
    category = c("SA / Defesa","SA / Defesa","SA / Defesa","SA / Defesa",
                 "JA / Defesa","JA / Defesa",
                 "ABA / Defesa","ABA / Defesa","ABA / Seca","ABA / Seca",
                 "Etileno / Defesa","Etileno / Defesa",
                 "Luz / Desenvolvimento","Brassinoesteroides",
                 "Luz / Desenvolvimento","Desenvolvimento","Desenvolvimento",
                 "Circadiano","Estresse térmico","Desenvolvimento")
  )
} else {
  raw <- read_csv(counts_file, show_col_types = FALSE)
  mat <- raw %>%
    rename(gene_id = gene) %>%
    left_join(gene_meta %>% select(gene_id, label), by = "gene_id") %>%
    mutate(label = coalesce(label, gene_id)) %>%
    select(-gene_id) %>%
    column_to_rownames("label") %>%
    as.matrix()
  # Ordenar linhas por gene_meta
  mat <- mat[intersect(gene_meta$label, rownames(mat)), , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  if (file.exists(ann_file)) {
    tf_ann <- read_tsv(ann_file, show_col_types = FALSE)
  } else {
    tf_ann <- tibble(tf_name  = colnames(mat),
                     family   = "Desconhecida",
                     category = "Outros")
  }
}

# ── Filtro: TFs presentes em ≥ 3 dos 49 genes ────────────────────────────────
keep <- colnames(mat)[colSums(mat > 0) >= 3]
mat     <- mat[, keep, drop = FALSE]
tf_ann  <- tf_ann %>% filter(tf_name %in% keep)

if (ncol(mat) == 0) stop("Nenhum TF passou o filtro de >=3 genes. Verifique o parse.")

# Ordenar colunas por categoria → família → nome
tf_ord <- tf_ann %>%
  mutate(cat_rank = match(category, cat_order, nomatch = 99)) %>%
  arrange(cat_rank, family, tf_name) %>%
  filter(tf_name %in% colnames(mat))
mat <- mat[, tf_ord$tf_name[tf_ord$tf_name %in% colnames(mat)], drop = FALSE]

message(sprintf("Heatmap: %d genes × %d TFs", nrow(mat), ncol(mat)))

# ── Anotações de colunas e linhas ─────────────────────────────────────────────
col_anno <- data.frame(
  Família   = tf_ord$family[tf_ord$tf_name %in% colnames(mat)],
  Categoria = tf_ord$category[tf_ord$tf_name %in% colnames(mat)],
  row.names = colnames(mat)
)

n_fam  <- length(unique(col_anno$Família))
n_cat  <- length(unique(col_anno$Categoria))

fam_cols <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                     "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
                     "#8DA0CB","#E78AC3","#A6D854"))(n_fam),
  unique(col_anno$Família)
)
cat_cols <- setNames(
  colorRampPalette(brewer.pal(min(n_cat, 9), "Set1"))(n_cat),
  unique(col_anno$Categoria)
)

row_anno <- data.frame(
  Focal = ifelse(rownames(mat) %in% gene_meta$label[gene_meta$focal], "Sim", "Não"),
  row.names = rownames(mat)
)

ann_colors <- list(
  Focal     = c("Sim" = "#1A252F", "Não" = "#ECF0F1"),
  Família   = fam_cols,
  Categoria = cat_cols
)

# Gaps entre categorias
cats_present <- col_anno$Categoria[!duplicated(col_anno$Categoria)]
n_cats_p <- length(cats_present)
if (n_cats_p > 1) {
  gap_cols <- cumsum(table(col_anno$Categoria)[cats_present])[-n_cats_p]
} else {
  gap_cols <- integer(0)
}

# ── Dimensões de publicação: 174 mm double-column ────────────────────────────
W_in  <- 174 / 25.4                          # 6.85 pol
H_in  <- max(8, nrow(mat) * 0.17 + 4.5)     # altura dinâmica

plot_args <- list(
  mat,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  color             = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#7B241C"))(60),
  border_color      = "grey82",
  fontsize          = 7,
  fontsize_row      = 7.5,
  fontsize_col      = 7,
  angle_col         = 45,
  annotation_col    = col_anno,
  annotation_row    = row_anno,
  annotation_colors = ann_colors,
  annotation_names_col = TRUE,
  annotation_names_row = FALSE,
  display_numbers   = TRUE,
  number_format     = "%d",
  number_color      = "grey20",
  gaps_col          = gap_cols,
  silent            = TRUE
)

outfile <- "fimo_tf_heatmap.pdf"
cairo_pdf(outfile, width = W_in, height = H_in, family = "Helvetica")
do.call(pheatmap, c(plot_args, list(
  main = "Fatores de transcrição em promotores de LRR-RLPs (JASPAR2024; 2 kb; q < 0.05)"
)))
dev.off()
message(sprintf("PDF publicacao: %s (%.0f x %.0f mm)", outfile, W_in*25.4, H_in*25.4))

png(sub("\\.pdf$",".png", outfile), width = W_in, height = H_in,
    units = "in", res = 300, type = "cairo")
do.call(pheatmap, c(plot_args, list(main = "FIMO — LRR-RLP TF binding sites (JASPAR2024)")))
dev.off()
message("PNG 300 dpi: ", sub("\\.pdf$",".png", outfile))

# ── Resumo ────────────────────────────────────────────────────────────────────
cat("\n=== Top 20 TFs (genes com binding site) ===\n")
presence <- colSums(mat > 0)
sort(presence, decreasing = TRUE)[1:min(20, length(presence))] %>%
  as.data.frame() %>% setNames("n_genes") %>% print()

cat("\n=== Genes focais — TFs exclusivos ===\n")
focal_rows <- rownames(mat)[rownames(mat) %in% gene_meta$label[gene_meta$focal]]
for (r in focal_rows) {
  specific <- names(which(mat[r,] > 0 & colSums(mat[-which(rownames(mat)==r),] > 0) == 0))
  if (length(specific) > 0)
    cat(sprintf("  %-24s: %s\n", r, paste(specific, collapse=", ")))
}
