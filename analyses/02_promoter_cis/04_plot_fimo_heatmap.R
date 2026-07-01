#!/usr/bin/env Rscript
# 04_plot_fimo_heatmap.R
# Heatmap de fatores de transcrição (JASPAR2024) nos promotores (2kb)
# dos 49 LRR-RLPs — alternativa server-side ao PlantCARE.
#
# Uso:
#   Rscript 04_plot_fimo_heatmap.R fimo_parsed_counts.csv [fimo_parsed_tf_families.tsv]

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

script_dir <- tryCatch(dirname(normalizePath(sys.frame(0)$ofile)), error = function(e) getwd())
setwd(script_dir)

# ── Argumentos ───────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
counts_file <- if (length(args) > 0) args[1] else "fimo_parsed_counts.csv"
ann_file    <- if (length(args) > 1) args[2] else "fimo_parsed_tf_families.tsv"

# ── Metadados dos 49 genes ───────────────────────────────────────────────────
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

# ── Ler contagens FIMO ───────────────────────────────────────────────────────
if (!file.exists(counts_file)) {
  message("[DEMO] Gerando dados simulados (FIMO não rodou ainda)...")
  set.seed(42)
  demo_tfs <- c("WRKY6","WRKY18","WRKY40","WRKY70","ABF1","ABF2","ABF3","ABF4",
                "MYB2","MYB96","ERF1","ORA59","PIF1","PIF3","PIF4","BES1","HY5",
                "NAC16","NAC29","TCP4","MYC2","MYC3","EIN3","ABI3","ABI5",
                "GATA17","DOF2.1","DOF5.1","HSF1","CCA1")
  mat <- matrix(
    sample(0:5, length(gene_meta$gene_id) * length(demo_tfs), replace = TRUE,
           prob = c(0.45, 0.25, 0.14, 0.08, 0.05, 0.03)),
    nrow = length(gene_meta$gene_id),
    dimnames = list(gene_meta$label, demo_tfs)
  )
  tf_ann <- tibble(
    tf_name  = demo_tfs,
    family   = c(rep("WRKY",4), rep("bZIP",4), rep("MYB",2), rep("AP2/ERF",2),
                 rep("bHLH",5), "NAC","NAC","TCP","bHLH","bHLH","EIL","bZIP","bZIP",
                 "GATA","Dof","Dof","HSF","MYB"),
    category = c(rep("SA / Defesa",4), rep("ABA / Defesa",4), rep("ABA / Seca",2),
                 rep("Etileno / Defesa",2), "Luz / Desenvolvimento","Luz / Desenvolvimento",
                 "Luz / Desenvolvimento","Brassinoesteroides","Luz / Desenvolvimento",
                 "Desenvolvimento","Desenvolvimento","Desenvolvimento",
                 "JA / Desenvolvimento","JA / Desenvolvimento","Etileno","ABA / Defesa",
                 "ABA / Defesa","Desenvolvimento","Luz / GA","Luz / GA",
                 "Estresse térmico","Circadiano")
  )
} else {
  raw <- read_csv(counts_file, show_col_types = FALSE)
  mat_df <- raw %>%
    left_join(gene_meta %>% select(gene_id = gene_id, label), by = c("gene" = "gene_id")) %>%
    mutate(label = if_else(is.na(label), gene, label)) %>%
    column_to_rownames("label") %>%
    select(-gene) %>%
    as.matrix()
  mat <- mat_df[gene_meta$label[gene_meta$label %in% rownames(mat_df)], ]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  if (file.exists(ann_file)) {
    tf_ann <- read_tsv(ann_file, show_col_types = FALSE)
  } else {
    tf_ann <- tibble(tf_name = colnames(mat), family = "Desconhecida", category = "Outros")
  }
}

# ── Filtrar: manter apenas TFs presentes em ≥ 3 genes ───────────────────────
keep_tfs <- colnames(mat)[colSums(mat > 0) >= 3]
mat <- mat[, keep_tfs, drop = FALSE]
tf_ann <- tf_ann %>% filter(tf_name %in% keep_tfs)

# Ordenar colunas por categoria funcional
cat_order <- c("SA / Defesa","ABA / Defesa","ABA / Seca","Etileno / Defesa",
               "JA / Desenvolvimento","Luz / Desenvolvimento","Luz / GA",
               "GA / Desenvolvimento","Brassinoesteroides","Circadiano",
               "Desenvolvimento","Desenvolvimento floral","Meristema",
               "Estresse térmico","Outros")
tf_ann_ord <- tf_ann %>%
  filter(tf_name %in% colnames(mat)) %>%
  mutate(cat_rank = match(category, cat_order, nomatch = 99)) %>%
  arrange(cat_rank, family, tf_name)
mat <- mat[, tf_ann_ord$tf_name[tf_ann_ord$tf_name %in% colnames(mat)], drop = FALSE]

# ── Anotação de colunas ───────────────────────────────────────────────────────
col_anno <- data.frame(
  Família    = tf_ann_ord$family[tf_ann_ord$tf_name %in% colnames(mat)],
  Categoria  = tf_ann_ord$category[tf_ann_ord$tf_name %in% colnames(mat)],
  row.names  = colnames(mat)
)

n_fam  <- length(unique(col_anno$Família))
n_cat  <- length(unique(col_anno$Categoria))
fam_cols <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                     "#A65628","#F781BF","#999999","#66C2A5","#FC8D62"))(n_fam),
  unique(col_anno$Família)
)
cat_cols <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(n_cat),
  unique(col_anno$Categoria)
)

# Anotação de linhas: destacar genes focais
row_anno <- data.frame(
  Focal = ifelse(gene_meta$focal[gene_meta$label %in% rownames(mat)], "Sim", "Não"),
  row.names = rownames(mat)
)
row_ann_cols <- list(
  Focal = c("Sim" = "#2C3E50", "Não" = "#ECF0F1"),
  Família   = fam_cols,
  Categoria = cat_cols
)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig_w <- max(14, ncol(mat) * 0.38 + 5)
fig_h <- max(10, nrow(mat) * 0.28 + 3)

outfile <- "fimo_tf_heatmap.pdf"
cairo_pdf(outfile, width = fig_w, height = fig_h)
pheatmap(
  mat,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  color           = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#7B241C"))(50),
  border_color    = "grey85",
  cellwidth       = 18,
  cellheight      = 16,
  fontsize        = 8,
  fontsize_row    = 9,
  fontsize_col    = 7.5,
  angle_col       = 45,
  annotation_col  = col_anno,
  annotation_row  = row_anno,
  annotation_colors = row_ann_cols,
  main            = "Fatores de transcrição nos promotores de LRR-RLPs (JASPAR2024; 2 kb upstream; q < 0.05)",
  display_numbers = TRUE,
  number_format   = "%d",
  number_color    = "black"
)
dev.off()
message("Figura salva: ", outfile, " (", nrow(mat), " genes × ", ncol(mat), " TFs)")

# PNG para preview
png(sub("\\.pdf$",".png", outfile),
    width = fig_w * 96, height = fig_h * 96, res = 150)
pheatmap(
  mat,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  color           = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#7B241C"))(50),
  border_color    = "grey85",
  cellwidth       = 18,
  cellheight      = 16,
  fontsize        = 8,
  fontsize_row    = 9,
  fontsize_col    = 7.5,
  angle_col       = 45,
  annotation_col  = col_anno,
  annotation_row  = row_anno,
  annotation_colors = row_ann_cols,
  main            = "Fatores de transcrição — LRR-RLP promotores (JASPAR2024)",
  display_numbers = TRUE,
  number_format   = "%d",
  number_color    = "black"
)
dev.off()
message("PNG preview salvo: ", sub("\\.pdf$",".png", outfile))
