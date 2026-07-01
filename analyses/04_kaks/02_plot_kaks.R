#!/usr/bin/env Rscript
# 02_plot_kaks.R — Figura Ka/Ks em qualidade de publicação
# Painel A: barplot ordenado por Ka/Ks (13 pares)
# Painel B: scatter Ka vs Ks com anotação dos pares focais
#
# Uso: Rscript 02_plot_kaks.R [kaks_summary.tsv]
# Saída: kaks_figure.pdf (174 mm × 110 mm, double-column journal)
#        kaks_figure.png (300 dpi)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(patchwork)
})

script_dir <- tryCatch(dirname(normalizePath(sys.frame(0)$ofile)), error = function(e) getwd())
setwd(script_dir)

args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "kaks_summary.tsv"

if (!file.exists(input_file)) stop("Arquivo não encontrado: ", input_file)

raw <- read_tsv(input_file, show_col_types = FALSE)
# TSV tem 8 colunas: pair sequence method ka ks ka_ks pvalue len
colnames(raw) <- c("pair","sequence","method","Ka","Ks","KaKs","pvalue","length")

# Pares focais para anotação
FOCAL_PAIRS <- c(
  "Solyc05g055190-Solyc03g112680",   # SlRLP1 × SlRLP2 (CLV2-like, segmental)
  "Solyc02g072250-Solyc02g092040"    # SlRLP5 × SlRLP6 (SNC2/3, segmental)
)

df <- raw %>%
  mutate(
    gene_A = str_extract(pair, "^[^-]+"),
    gene_B = str_extract(pair, "(?<=-).+"),
    chrom_A = str_extract(gene_A, "(?<=Solyc)\\d+"),
    chrom_B = str_extract(gene_B, "(?<=Solyc)\\d+"),
    tipo = if_else(chrom_A == chrom_B, "Tandem", "Segmental"),
    focal = pair %in% FOCAL_PAIRS,
    # Rótulo compacto: cromossomo + posição reduzida
    label_short = paste0(
      "Chr", chrom_A, "\n",
      str_replace(gene_A, "Solyc\\d+g0*(\\d+)", "\\1"),
      "×",
      str_replace(gene_B, "Solyc\\d+g0*(\\d+)", "\\1")
    ),
    # Rótulo completo para genes focais
    label_focal = case_when(
      pair == "Solyc05g055190-Solyc03g112680" ~ "SlRLP1×SlRLP2\n(CLV2-like)",
      pair == "Solyc02g072250-Solyc02g092040" ~ "SlRLP5×SlRLP6\n(SNC2/3)",
      TRUE ~ NA_character_
    ),
    pair_order = reorder(pair, KaKs)
  )

# ── Paleta Okabe-Ito (colorblind-friendly) ────────────────────────────────────
pal <- c("Tandem" = "#0072B2", "Segmental" = "#E69F00")

# ── Painel A: barplot ─────────────────────────────────────────────────────────
pA <- ggplot(df, aes(x = pair_order, y = KaKs, fill = tipo)) +
  geom_col(width = 0.72, color = NA) +
  # Destaque para pares focais
  geom_col(data = filter(df, focal), aes(x = pair_order, y = KaKs),
           fill = NA, color = "#CC0000", linewidth = 0.6, width = 0.72) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#CC0000",
             linewidth = 0.5, alpha = 0.8) +
  annotate("text", x = 0.55, y = 1.06, label = "Ka/Ks = 1",
           color = "#CC0000", hjust = 0, size = 2.6, fontface = "italic") +
  geom_text(aes(label = sprintf("%.3f", KaKs), color = focal),
            vjust = -0.35, size = 2.5) +
  scale_fill_manual(values = pal, name = "Duplicação") +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#CC0000"), guide = "none") +
  scale_y_continuous(limits = c(0, 1.18), expand = c(0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  scale_x_discrete(labels = setNames(df$label_short, df$pair)) +
  labs(
    x = NULL,
    y = expression(italic("K")[a] * "/" * italic("K")[s])
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 6.5,
                                     color = ifelse(levels(df$pair_order) %in% FOCAL_PAIRS,
                                                    "#CC0000", "grey20")),
    axis.text.y      = element_text(size = 7.5),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.4),
    legend.position  = c(0.85, 0.82),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.key.size  = unit(0.35, "cm"),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "grey90"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(4, 6, 2, 4, "mm")
  )

# ── Painel B: scatter Ka vs Ks ────────────────────────────────────────────────
max_val <- max(c(df$Ka, df$Ks)) * 1.08

pB <- ggplot(df, aes(x = Ks, y = Ka, fill = tipo)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55",
              linewidth = 0.4) +
  annotate("text", x = max_val * 0.88, y = max_val * 0.94,
           label = expression(italic("K")[a] * " = " * italic("K")[s]),
           angle = 43, color = "grey55", size = 2.5) +
  geom_point(aes(size = KaKs), shape = 21, color = "white", stroke = 0.4, alpha = 0.9) +
  # Círculo extra para pares focais
  geom_point(data = filter(df, focal),
             aes(x = Ks, y = Ka, size = KaKs),
             shape = 21, fill = NA, color = "#CC0000", stroke = 1.0) +
  geom_label_repel(
    data = filter(df, focal),
    aes(label = label_focal),
    size = 2.5, color = "#CC0000",
    box.padding = 0.5, point.padding = 0.4,
    segment.color = "#CC0000", segment.size = 0.4,
    min.segment.length = 0, max.overlaps = 10,
    fill = "white", label.size = 0.2
  ) +
  scale_fill_manual(values = pal, name = "Duplicação", guide = "none") +
  scale_size_continuous(
    range  = c(2.5, 6.5),
    name   = expression(italic("K")[a] * "/" * italic("K")[s]),
    limits = c(0, 1),
    breaks = c(0.25, 0.5, 0.75, 1.0)
  ) +
  scale_x_continuous(limits = c(0, max_val), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, max_val * 0.6), expand = c(0.02, 0)) +
  labs(
    x = expression(italic("K")[s] * " (substituições sinônimas/sítio)"),
    y = expression(italic("K")[a] * " (substituições não-sinônimas/sítio)")
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text        = element_text(size = 7.5),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.4),
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.key.size  = unit(0.3, "cm"),
    panel.grid       = element_blank(),
    plot.margin = margin(4, 4, 2, 4, "mm")
  )

# ── Combinar painéis A+B ──────────────────────────────────────────────────────
fig <- pA + pB +
  plot_layout(widths = c(1.8, 1.0)) +
  plot_annotation(
    tag_levels = "A",
    caption    = sprintf(
      "Todos os %d pares apresentam Ka/Ks < 1 (método NG86; Nei & Gojobori, 1986 via BioPython). Barras vermelhas = pares segmentais focais (SlRLP1×SlRLP2; SlRLP5×SlRLP6).",
      nrow(df)
    )
  ) &
  theme(
    plot.tag   = element_text(size = 9, face = "bold"),
    plot.caption = element_text(size = 6.5, color = "grey40", hjust = 0)
  )

# ── Salvar ────────────────────────────────────────────────────────────────────
# Dimensões: 174 mm × 100 mm (double-column journal standard)
W_in <- 174 / 25.4   # 6.85 pol
H_in <- 100 / 25.4   # 3.94 pol

cairo_pdf("kaks_figure.pdf", width = W_in, height = H_in, family = "Helvetica")
print(fig)
dev.off()
message("PDF publicação: kaks_figure.pdf (174 × 100 mm)")

png("kaks_figure.png", width = W_in, height = H_in,
    units = "in", res = 300, type = "cairo")
print(fig)
dev.off()
message("PNG 300 dpi: kaks_figure.png")

# ── Estatísticas resumo ───────────────────────────────────────────────────────
cat("\n=== Resumo Ka/Ks por tipo ===\n")
df %>%
  group_by(tipo) %>%
  summarise(n = n(), Ka_mean = mean(Ka), Ks_mean = mean(Ks),
            KaKs_mean = mean(KaKs), KaKs_min = min(KaKs), KaKs_max = max(KaKs)) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print()

cat("\n=== Pares (ordenados por Ka/Ks crescente) ===\n")
df %>%
  arrange(KaKs) %>%
  select(pair, tipo, Ka, Ks, KaKs, focal) %>%
  mutate(across(c(Ka,Ks,KaKs), ~ round(., 3))) %>%
  print(n = Inf)
