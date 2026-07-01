#!/usr/bin/env Rscript
# Visualização Ka/Ks dos 13 pares parálogos de LRR-RLPs
# Uso: Rscript 02_plot_kaks.R [kaks_summary.tsv]

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

script_dir <- tryCatch(dirname(normalizePath(sys.frame(0)$ofile)), error = function(e) getwd())
setwd(script_dir)

args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "kaks_summary.tsv"

if (!file.exists(input_file)) stop("Arquivo não encontrado: ", input_file)

raw <- read_tsv(input_file, show_col_types = FALSE)

# Esperado: Sequence Method Ka Ks Ka/Ks P-Value Length
# Coluna 1 tem formato "GeneA-GeneB"
colnames(raw) <- c("pair","method","Ka","Ks","KaKs","pvalue","length")

df <- raw %>%
  mutate(
    # Extrair genes A e B
    gene_A = str_extract(pair, "^[^-]+"),
    gene_B = str_extract(pair, "(?<=-).+"),
    # Classificar tipo
    tipo = if_else(
      str_extract(gene_A, "Solyc(\\d+)", group=1) == str_extract(gene_B, "Solyc(\\d+)", group=1),
      "Tandem", "Segmental"
    ),
    # Rótulo curto para o gráfico
    label = paste0(
      str_replace(gene_A, "Solyc(\\d+)g(\\d+)", "\\1g\\2"),
      "\n×\n",
      str_replace(gene_B, "Solyc(\\d+)g(\\d+)", "\\1g\\2")
    ),
    pair_order = reorder(pair, KaKs)
  )

# ── Plot principal: barras ordenadas por Ka/Ks ───────────────────────────────
pal <- c("Tandem" = "#2196F3", "Segmental" = "#FF9800")

p1 <- ggplot(df, aes(x = pair_order, y = KaKs, fill = tipo)) +
  geom_col(width = 0.7, color = "white") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 0.6, y = 1.04, label = "Ka/Ks = 1 (neutro)", color = "#E74C3C",
           hjust = 0, size = 3.2) +
  geom_text(aes(label = sprintf("%.3f", KaKs)), vjust = -0.4, size = 3.0, color = "grey30") +
  scale_fill_manual(values = pal, name = "Tipo de duplicação") +
  scale_y_continuous(limits = c(0, 1.15), expand = c(0,0)) +
  labs(
    title    = "Pressão seletiva em pares parálogos de LRR-RLPs (S. lycopersicum)",
    subtitle = sprintf("Método: NG86 (Nei & Gojobori, 1986) — todos os %d pares sob seleção purificadora (Ka/Ks < 1)", nrow(df)),
    x        = "Par de genes",
    y        = "Ka/Ks"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    legend.position  = "top",
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave("kaks_barplot.pdf", p1, width = 12, height = 6)
ggsave("kaks_barplot.png", p1, width = 12, height = 6, dpi = 300)
message("Figura salva: kaks_barplot.pdf / .png")

# ── Plot complementar: Ka vs Ks scatter ──────────────────────────────────────
p2 <- ggplot(df, aes(x = Ks, y = Ka, color = tipo, label = tipo)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(size = KaKs), alpha = 0.85) +
  scale_color_manual(values = pal, name = "Tipo") +
  scale_size_continuous(range = c(2, 6), name = "Ka/Ks", limits = c(0, 1)) +
  annotate("text", x = max(df$Ks)*0.9, y = max(df$Ks)*0.9, label = "Ka = Ks",
           angle = 42, color = "grey50", size = 3) +
  labs(
    title = "Ka vs Ks — pares parálogos LRR-RLP",
    x     = "Ks (taxa sinônima)",
    y     = "Ka (taxa não-sinônima)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right")

ggsave("kaks_scatter.pdf", p2, width = 7, height = 6)
ggsave("kaks_scatter.png", p2, width = 7, height = 6, dpi = 300)
message("Figura salva: kaks_scatter.pdf / .png")

# ── Estatísticas resumo ───────────────────────────────────────────────────────
cat("\n=== Resumo Ka/Ks por tipo de duplicação ===\n")
df %>%
  group_by(tipo) %>%
  summarise(
    n       = n(),
    Ka_mean = mean(Ka),
    Ks_mean = mean(Ks),
    KaKs_mean = mean(KaKs),
    KaKs_min  = min(KaKs),
    KaKs_max  = max(KaKs)
  ) %>%
  print()

cat("\n=== Par com menor Ka/Ks (maior conservação) ===\n")
print(df %>% arrange(KaKs) %>% select(pair, tipo, Ka, Ks, KaKs) %>% head(3))

cat("\n=== Par com maior Ka/Ks (menor conservação) ===\n")
print(df %>% arrange(desc(KaKs)) %>% select(pair, tipo, Ka, Ks, KaKs) %>% head(3))
