#!/usr/bin/env Rscript
# Arquitetura de domínios dos 49 LRR-RLPs de S. lycopersicum
# Uso: Rscript 02_plot_domain_architecture.R [hmmer_domains.tsv]
#
# Entrada: TSV gerado por 01_run_hmmer.sh com colunas:
#   gene_id, isoform_id, protein_len, domain, pfam_acc,
#   ali_start, ali_end, env_start, env_end, ievalue, score, description
#
# Produz: domain_architecture.pdf  (todos os 49 genes, destacando os 7 focais)
#         domain_counts.tsv        (matriz gene × família de domínio)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
})

script_dir <- tryCatch(dirname(normalizePath(sys.frame(0)$ofile)), error = function(e) getwd())
setwd(script_dir)

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
  ),
  focal = c(
    rep(FALSE,15),
    FALSE,TRUE,TRUE,            # 16-18
    FALSE,FALSE,TRUE,           # 19-21
    FALSE,
    TRUE,FALSE,TRUE,            # 23-25
    rep(FALSE,12),
    TRUE,FALSE,                 # 38-39 (38 é focal? na lista original é Seca/Vírus)
    FALSE,
    rep(FALSE,4),TRUE,rep(FALSE,3)  # 45 é TMM/focal
  )
)

# Corrigir: 7 genes focais são Solyc05g055190, Solyc03g112680, Solyc05g009990,
#           Solyc12g042760, Solyc02g072250, Solyc02g092040, Solyc10g007830
focal_ids <- c("Solyc05g055190","Solyc03g112680","Solyc05g009990",
               "Solyc12g042760","Solyc02g072250","Solyc02g092040","Solyc10g007830")
gene_meta <- gene_meta %>% mutate(focal = gene_id %in% focal_ids)

# ── Paleta de domínios ────────────────────────────────────────────────────────
domain_colors <- c(
  "LRR_1"        = "#3498DB",
  "LRR_2"        = "#2980B9",
  "LRR_3"        = "#1A6BA0",
  "LRR_4"        = "#85C1E9",
  "LRR_8"        = "#AED6F1",
  "LRR_RI"       = "#5DADE2",
  "LRR_RII"      = "#76B7E8",
  "LRRNT"        = "#154360",
  "LRRCT"        = "#1B4F72",
  "Transmembrane" = "#E74C3C",
  "Signal_pep"   = "#F39C12",
  "PAN"          = "#27AE60",
  "EGF"          = "#8E44AD",
  "other"        = "#95A5A6"
)

# ── Carregar dados ────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "hmmer_out/hmmer_domains.tsv"

if (!file.exists(input_file)) {
  message("[DEMO] Gerando dados sintéticos para 49 LRR-RLPs...")
  set.seed(42)
  demo_rows <- list()
  for (i in seq_len(nrow(gene_meta))) {
    g    <- gene_meta$gene_id[i]
    plen <- sample(400:1100, 1)
    demo_rows[[length(demo_rows)+1]] <- tibble(
      gene_id=g, protein_len=plen, domain="Signal_pep",
      start=1, end=25, evalue=1e-5
    )
    n_lrr <- sample(3:20, 1)
    lrr_s <- 80
    for (k in seq_len(n_lrr)) {
      s <- lrr_s + (k-1)*27
      e <- s + 23
      if (e > plen - 50) break
      demo_rows[[length(demo_rows)+1]] <- tibble(
        gene_id=g, protein_len=plen,
        domain=sample(c("LRR_1","LRR_2","LRR_8","LRR_RI"), 1),
        start=s, end=e, evalue=1e-8
      )
    }
    demo_rows[[length(demo_rows)+1]] <- tibble(
      gene_id=g, protein_len=plen, domain="Transmembrane",
      start=plen-80, end=plen-58, evalue=1e-4
    )
  }
  dom <- bind_rows(demo_rows)
} else {
  dom <- read_tsv(input_file, show_col_types = FALSE)
  if ("ali_start" %in% colnames(dom)) dom <- dom %>% rename(start = ali_start, end = ali_end)
  if (!"protein_len" %in% colnames(dom) && "query_len" %in% colnames(dom))
    dom <- dom %>% rename(protein_len = query_len)
}

# ── Preparar dados para plot ──────────────────────────────────────────────────
# Usar todos os genes presentes no HMMER output, na ordem cromossômica
dom_prep <- dom %>%
  inner_join(gene_meta %>% select(gene_id, label, focal), by = "gene_id") %>%
  mutate(
    dom_group = case_when(
      str_detect(domain, "^LRR")           ~ domain,
      str_detect(domain, "ignal")          ~ "Signal_pep",
      str_detect(domain, "TM|Trans|Hydro") ~ "Transmembrane",
      TRUE                                 ~ "other"
    ),
    # Posição na ordem cromossômica (inverter para ggplot y de cima p/ baixo)
    gene_rank = match(gene_id, gene_meta$gene_id),
    gene_rank_inv = max(gene_rank) + 1 - gene_rank,
    label_f   = factor(label, levels = rev(gene_meta$label))
  )

# Genes sem nenhum domínio HMMER (presentes em gene_meta mas ausentes no dom)
missing_genes <- setdiff(gene_meta$gene_id, unique(dom_prep$gene_id))
if (length(missing_genes) > 0)
  message("Genes sem domínio HMMER: ", paste(missing_genes, collapse=", "))

backbone <- dom_prep %>%
  group_by(gene_id, label_f, focal) %>%
  summarise(protein_len = first(protein_len), .groups = "drop")

used_domains <- unique(dom_prep$dom_group)
domain_pal   <- c(domain_colors[intersect(names(domain_colors), used_domains)], "other"="#95A5A6")

n_genes  <- length(unique(dom_prep$gene_id))
fig_h    <- max(8, n_genes * 0.22 + 2)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  # Backbone — focal genes em cinza escuro, demais em cinza claro
  geom_segment(
    data = backbone,
    aes(x = 0, xend = protein_len,
        y = label_f, yend = label_f,
        color = focal),
    linewidth = 1.2
  ) +
  scale_color_manual(values = c("TRUE"="#555555","FALSE"="#CCCCCC"), guide = "none") +
  ggnewscale::new_scale_color() +
  # Domínios
  geom_rect(
    data = dom_prep,
    aes(xmin = start, xmax = end,
        ymin = as.numeric(label_f) - 0.38,
        ymax = as.numeric(label_f) + 0.38,
        fill = dom_group),
    color = "white", linewidth = 0.2
  ) +
  scale_fill_manual(values = domain_pal, name = "Domínio Pfam") +
  scale_x_continuous(labels = comma, expand = c(0.01, 0)) +
  labs(
    title   = "Arquitetura de domínios dos 49 LRR-RLPs de S. lycopersicum (ITAG4.0)",
    subtitle = "* = genes focais (SlRLP1–SlRLP7); ordenados por posição cromossômica",
    x       = "Posição na proteína (aa)",
    y       = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(
      size  = 7.5,
      face  = if_else(grepl("\\*", levels(dom_prep$label_f)[as.numeric(dom_prep$label_f)],
                             fixed=FALSE), "bold", "plain")[1]
    ),
    plot.title    = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

# Negrito nos genes focais via after_scale não é direto no theme;
# usar geom_text de um y-axis substituto para destacar os focais
focal_labels <- gene_meta %>% filter(focal) %>%
  mutate(label_f = factor(label, levels = rev(gene_meta$label)))

# ── Salvar ────────────────────────────────────────────────────────────────────
ggsave("domain_architecture.pdf", p, width = 14, height = fig_h)
message("Figura salva: domain_architecture.pdf")
ggsave("domain_architecture.png", p, width = 14, height = fig_h, dpi = 200)
message("PNG salvo: domain_architecture.png")

# ── Tabela de contagem — todos os 49 genes ────────────────────────────────────
summary_tbl <- dom_prep %>%
  group_by(gene_id, label, dom_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = dom_group, values_from = count, values_fill = 0L)

write_tsv(summary_tbl, "domain_counts.tsv")
message("Tabela de domínios (49 genes): domain_counts.tsv")
print(summary_tbl)
