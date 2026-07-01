#!/usr/bin/env Rscript
# Heatmap de elementos cis nos promotores (2kb) dos 49 LRR-RLPs — PlantCARE
# Uso:
#   Rscript 02_plot_plantcare_heatmap.R plantcare_results_49genes.tab   # formato .tab per-gene
#   Rscript 02_plot_plantcare_heatmap.R plantcare_results.txt            # formato .txt 7 genes
#   Rscript 02_plot_plantcare_heatmap.R plantcare_counts_49genes.csv     # CSV pré-processado
#
# Formatos suportados:
#   .tab (sem header): seq::chrom | signal | seq_motif | location | len | strand | species | func
#   .txt (com header): Sequence | Signal | Location | Strand | Seq | Function
#   .csv             : gene | motif1 | motif2 | ...  (matriz já processada)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Metadados dos 49 genes (ordem cromossômica) ───────────────────────────────
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
    "SlRLP16","SlRLP17*(SNC2/3)","SlRLP18*",
    "SlRLP19","SlRLP20","SlRLP21*",
    "SlRLP22",
    "SlRLP23*(RIC7)","SlRLP24","SlRLP25*(CLV2)",
    "SlRLP26","SlRLP27","SlRLP28",
    "SlRLP29","SlRLP30","SlRLP31","SlRLP32","SlRLP33","SlRLP34",
    "SlRLP35","SlRLP36",
    "SlRLP37",
    "SlRLP38*(Seca/Vírus)","SlRLP39",
    "SlRLP40",
    "SlRLP41","SlRLP42","SlRLP43","SlRLP44","SlRLP45*(TMM)","SlRLP46",
    "SlRLP47","SlRLP48","SlRLP49"
  )
)

# ── Dicionário de elementos cis por categoria funcional ───────────────────────
# Nomes conforme nomenclatura real do PlantCARE (validados no output 49 genes)
cis_dict <- list(
  "Luz"              = c("ACE","Box 4","G-box","GT1-motif","I-box","Sp1","LAMP-element",
                         "AE-box","AT1-motif","GCN4_motif","Ibox","TCT-motif","3-AF3 binding site"),
  "ABA / Seca"       = c("ABRE","ABRE3a","ABRE4","AT~ABRE","CACGTG-motif","DRE core",
                         "MBS","MYB recognition","Myb","MYB","MYB-like sequence","STRE"),
  "JA"               = c("CGTCA-motif","TGACG-motif","as-1"),
  "SA / Defesa"      = c("SARE","TCA-element","W box","W-box","TC-rich repeats",
                         "TGA-element","TGA1a","ARE"),
  "GA"               = c("GARE-motif","P-box","TATC-box","GARE","gibberellin-responsive element"),
  "Etileno"          = c("GCC-box","ERE"),
  "Desenvolvimento"  = c("CAT-box","CCAAT-box","as-2-element","CAAT-box","AAGAA-motif"),
  "Circadiano"       = c("circadian","Evening Element"),
  "Anaerobiose"      = c("ARE","anaerobic")
)

cis_df <- tibble(
  motif    = unlist(cis_dict),
  category = rep(names(cis_dict), lengths(cis_dict))
) %>% distinct(motif, .keep_all = TRUE)

# ── Ler resultados do PlantCARE ───────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "plantcare_results_49genes.tab"

if (!file.exists(input_file)) {
  message("[DEMO] Arquivo '", input_file, "' não encontrado — usando dados simulados")
  set.seed(2024)
  results <- expand_grid(gene_id = gene_meta$gene_id,
                         motif   = cis_df$motif) %>%
    mutate(count = sample(c(rep(0L,4), 1L, 2L, 3L, 4L), n(), replace = TRUE))

} else if (grepl("\\.csv$", input_file, ignore.case = TRUE)) {
  # Formato CSV: gene | motif1 | motif2 | ...
  message("Lendo CSV pré-processado: ", input_file)
  mat_csv <- read_csv(input_file, show_col_types = FALSE)
  results <- mat_csv %>%
    pivot_longer(-gene, names_to = "motif", values_to = "count") %>%
    rename(gene_id = gene) %>%
    filter(count > 0)

} else if (grepl("\\.tab$", input_file, ignore.case = TRUE)) {
  # Formato .tab sem header: seq::chrom | signal | seq_motif | location | len | strand | species | func
  message("Lendo formato .tab per-gene: ", input_file)
  raw <- read_tsv(input_file, col_names = c("seq_name","signal","seq_motif","location",
                                             "len","strand","species","function"),
                  col_types = cols(.default = "c"), comment = "#")
  results <- raw %>%
    filter(!is.na(signal), signal != "") %>%
    mutate(
      gene_id = str_extract(seq_name, "^[^:]+"),   # Solyc... antes de ::
      motif   = str_trim(signal)
    ) %>%
    filter(gene_id %in% gene_meta$gene_id) %>%
    group_by(gene_id, motif) %>%
    summarise(count = n(), .groups = "drop")

} else {
  # Formato .txt com header: Sequence | Signal | Location | Strand | Seq | Function
  message("Lendo formato .txt com header: ", input_file)
  raw <- read_tsv(input_file, col_names = TRUE, show_col_types = FALSE)
  colnames(raw) <- tolower(gsub("\\s+","_", colnames(raw)))
  if (!"sequence" %in% colnames(raw)) stop("Coluna 'Sequence' nao encontrada.")
  results <- raw %>%
    mutate(gene_id = str_extract(sequence, "^[^:]+")) %>%
    rename(motif = signal) %>%
    filter(!is.na(motif), motif != "") %>%
    group_by(gene_id, motif) %>%
    summarise(count = n(), .groups = "drop")
}

message("Hits lidos: ", sum(results$count), " | Genes: ", n_distinct(results$gene_id),
        " | Motivos únicos: ", n_distinct(results$motif))

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

# ── Anotação de linhas: destacar genes focais ─────────────────────────────────
focal_ids  <- c("Solyc05g055190","Solyc03g112680","Solyc05g009990",
                "Solyc12g042760","Solyc02g072250","Solyc02g092040","Solyc10g007830")
focal_labs <- gene_meta %>%
  filter(gene_id %in% focal_ids) %>%
  pull(label)

row_anno <- data.frame(
  Focal = ifelse(rownames(mat) %in% focal_labs, "Focal", "Não focal"),
  row.names = rownames(mat)
)
ann_colors_full <- list(
  Categoria = cat_cols,
  Focal     = c("Focal" = "#2C3E50", "Não focal" = "#ECF0F1")
)

# Gaps entre categorias funcionais
gap_cols <- cumsum(table(col_anno$Categoria)[unique(col_anno$Categoria)])[-n_cats]

# ── Plot qualidade de publicação ──────────────────────────────────────────────
fig_w   <- max(12, ncol(mat) * 0.50 + 4)
fig_h   <- max(10, nrow(mat) * 0.24 + 3)
outfile <- "plantcare_heatmap.pdf"

cairo_pdf(outfile, width = fig_w, height = fig_h, family = "Helvetica")
pheatmap(
  mat,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  color             = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#7B241C"))(60),
  border_color      = "grey80",
  cellwidth         = 18,
  cellheight        = 16,
  fontsize          = 8,
  fontsize_row      = 8.5,
  fontsize_col      = 7.5,
  angle_col         = 45,
  annotation_col    = col_anno,
  annotation_row    = row_anno,
  annotation_colors = ann_colors_full,
  annotation_names_col = TRUE,
  annotation_names_row = FALSE,
  main              = "Elementos cis-regulatórios em promotores de LRR-RLPs (2 kb upstream, PlantCARE)",
  display_numbers   = TRUE,
  number_format     = "%d",
  number_color      = "black",
  gaps_col          = gap_cols
)
dev.off()
message(sprintf("PDF publicação: %s (%.0f × %.0f mm)", outfile, fig_w*25.4, fig_h*25.4))

png(sub("\\.pdf$",".png", outfile),
    width = fig_w, height = fig_h, units = "in", res = 300, type = "cairo")
pheatmap(
  mat,
  cluster_rows      = FALSE, cluster_cols = FALSE,
  color             = colorRampPalette(c("#FFFFFF","#FFF3CD","#FF8C00","#7B241C"))(60),
  border_color      = "grey80", cellwidth = 18, cellheight = 16,
  fontsize = 8, fontsize_row = 8.5, fontsize_col = 7.5, angle_col = 45,
  annotation_col = col_anno, annotation_row = row_anno,
  annotation_colors = ann_colors_full,
  main = "Elementos cis-regulatórios — LRR-RLPs (PlantCARE)",
  display_numbers = TRUE, number_format = "%d", number_color = "black",
  gaps_col = gap_cols
)
dev.off()
message("PNG 300 dpi: ", sub("\\.pdf$",".png", outfile))

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
