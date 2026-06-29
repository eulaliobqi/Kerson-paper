#!/usr/bin/env Rscript
# Arquitetura de domínios dos 49 LRR-RLPs (ou 7 selecionados)
# Uso: Rscript 02_plot_domain_architecture.R [hmmer_domains.tsv]
#
# Entrada: TSV com colunas:
#   gene_id, protein_len, domain, start, end, evalue, score, description
#
# Produz: domain_architecture.pdf (estilo TBtools Gene Structure View)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gggenes)       # install.packages("gggenes")
  library(RColorBrewer)
  library(scales)
})

# ── Metadados dos 7 genes focais (ordem no plot) ─────────────────────────────
gene_order <- c(
  "Solyc05g055190","Solyc03g112680","Solyc05g009990",
  "Solyc12g042760","Solyc02g072250","Solyc02g092040","Solyc10g007830"
)
gene_labels <- c(
  "Solyc05g055190"="SlRLP1 (CLV2)",  "Solyc03g112680"="SlRLP2",
  "Solyc05g009990"="SlRLP3 (RIC7)",  "Solyc12g042760"="SlRLP4 (TMM)",
  "Solyc02g072250"="SlRLP5 (SNC2/3)","Solyc02g092040"="SlRLP6",
  "Solyc10g007830"="SlRLP7"
)

# ── Paleta de domínios ────────────────────────────────────────────────────────
domain_colors <- c(
  "LRR_1"       = "#3498DB",
  "LRR_2"       = "#2980B9",
  "LRR_3"       = "#1A6BA0",
  "LRR_4"       = "#85C1E9",
  "LRR_8"       = "#AED6F1",
  "LRR_RI"      = "#5DADE2",
  "LRR_RII"     = "#76B7E8",
  "LRRNT"       = "#154360",
  "LRRCT"       = "#1B4F72",
  "Transmembrane"= "#E74C3C",
  "Signal_pep"  = "#F39C12",
  "PAN"         = "#27AE60",
  "EGF"         = "#8E44AD",
  "other"       = "#95A5A6"
)

# ── Carregar dados ────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "hmmer_domains.tsv"

if (!file.exists(input_file)) {
  message("[DEMO] Gerando dados sintéticos de domínios LRR-RLP...")
  set.seed(42)

  demo_rows <- list()
  lens <- c(890,930,480,690,810,770,950)
  names(lens) <- gene_order

  for (i in seq_along(gene_order)) {
    g <- gene_order[i]
    plen <- lens[g]
    # Signal peptide
    demo_rows[[length(demo_rows)+1]] <- tibble(
      gene_id=g, protein_len=plen, domain="Signal_pep", start=1, end=25,
      evalue=1e-5, description="Signal peptide"
    )
    # LRR domains (multiple copies)
    n_lrr <- c(NA,32,9,10,6,6,NA)[i]
    n_lrr[is.na(n_lrr)] <- sample(5:15,1)
    lrr_start <- 100
    lrr_size  <- 24
    lrr_gap   <- 3
    for (k in seq_len(n_lrr)) {
      s <- lrr_start + (k-1)*(lrr_size+lrr_gap)
      e <- s + lrr_size - 1
      if (e > plen - 50) break
      demo_rows[[length(demo_rows)+1]] <- tibble(
        gene_id=g, protein_len=plen,
        domain=sample(c("LRR_1","LRR_2","LRR_8","LRR_RI"),1),
        start=s, end=e, evalue=1e-8, description="Leucine-rich repeat"
      )
    }
    # Transmembrane
    tm_start <- plen - 80
    demo_rows[[length(demo_rows)+1]] <- tibble(
      gene_id=g, protein_len=plen, domain="Transmembrane",
      start=tm_start, end=tm_start+20, evalue=1e-4, description="Transmembrane domain"
    )
  }
  dom <- bind_rows(demo_rows)
} else {
  dom <- read_tsv(input_file, show_col_types = FALSE)
}

# ── Filtrar e preparar ────────────────────────────────────────────────────────
# Manter apenas os 7 genes focais (se mais genes presentes)
focal_only <- dom %>% filter(gene_id %in% gene_order)
if (nrow(focal_only) == 0) {
  message("Genes focais não encontrados — usando todos os genes do arquivo")
  focal_only <- dom
}

# Mapear nomes de domínios para paleta
focal_only <- focal_only %>%
  mutate(
    dom_group = case_when(
      str_detect(domain, "^LRR")      ~ domain,
      str_detect(domain, "ignal")     ~ "Signal_pep",
      str_detect(domain, "TM|Trans|Hydro") ~ "Transmembrane",
      TRUE                            ~ "other"
    ),
    gene_label = gene_labels[gene_id],
    gene_label = factor(gene_label, levels = rev(gene_labels[gene_order]))
  )

# Backbone das proteínas (linha base)
backbone <- focal_only %>%
  group_by(gene_id) %>%
  summarise(protein_len = first(protein_len), .groups="drop") %>%
  mutate(
    gene_label = gene_labels[gene_id],
    gene_label = factor(gene_label, levels = rev(gene_labels[gene_order]))
  )

used_domains <- unique(focal_only$dom_group)
domain_pal   <- c(domain_colors[intersect(names(domain_colors), used_domains)],
                  "other" = "#95A5A6")

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  # Backbone (linha cinza representando a proteína)
  geom_segment(data = backbone,
               aes(x = 0, xend = protein_len, y = gene_label, yend = gene_label),
               color = "grey70", linewidth = 1.5) +
  # Domínios (retângulos)
  geom_rect(data = focal_only,
            aes(xmin = start, xmax = end,
                ymin = as.numeric(gene_label) - 0.35,
                ymax = as.numeric(gene_label) + 0.35,
                fill = dom_group),
            color = "white", linewidth = 0.3) +
  scale_fill_manual(values = domain_pal, name = "Domínio") +
  scale_x_continuous(labels = comma, expand = c(0.01,0)) +
  labs(
    title  = "Arquitetura de domínios dos LRR-RLPs de S. lycopersicum",
    x      = "Posição na proteína (aa)",
    y      = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 10, face = "bold"),
    plot.title       = element_text(size = 12)
  )

ggsave("domain_architecture.pdf", p, width = 13, height = 5)
message("Figura salva: domain_architecture.pdf")

# Versão PNG para inserir no docx
ggsave("domain_architecture.png", p, width = 13, height = 5, dpi = 300)
message("PNG salvo: domain_architecture.png")

# ── Tabela de contagem de domínios ────────────────────────────────────────────
summary_tbl <- focal_only %>%
  group_by(gene_id, gene_label, dom_group) %>%
  summarise(count = n(), .groups="drop") %>%
  pivot_wider(names_from=dom_group, values_from=count, values_fill=0L)

write_tsv(summary_tbl, "domain_counts.tsv")
message("Tabela de domínios: domain_counts.tsv")
print(summary_tbl)
