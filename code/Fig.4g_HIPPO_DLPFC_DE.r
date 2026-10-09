#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
  library(jsonlite)
})

# ---------------------------- settings --------------------------------------
input_dir <- "/data/home/ymq_codex/lhz/GWAS_MP"
config_file <- file.path(input_dir, "09ce9e5c654ad49f.json")
hippo_file <- file.path(input_dir, "HIPPO_diff_deseq2_used_plot.txt")
dlpfc_file <- file.path(input_dir, "DLPFC_diff_deseq2_used_plot.txt")

# Set OUTPUT_DIR to another directory when the input directory is read-only,
# for example: OUTPUT_DIR=/tmp Rscript plot_figure4g.R
output_dir <- Sys.getenv("OUTPUT_DIR", input_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_pdf <- file.path(output_dir, "Figure_4g_reproduced.pdf")
output_png <- file.path(output_dir, "Figure_4g_reproduced.png")

config <- if (file.exists(config_file)) {
  fromJSON(config_file, simplifyVector = TRUE)
} else {
  list(parameters = list())
}

param <- function(name, default = NULL) {
  value <- config$parameters[[name]]
  if (is.null(value) || length(value) == 0L || identical(value, "")) {
    return(default)
  }
  value[[1]]
}

as_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(param(name, default)))
  if (length(value) == 0L || is.na(value)) default else value
}

fc_cutoff_left <- as_num("FCcutoff1", 1)
fc_cutoff_right <- as_num("FCcutoff2", 1)
padj_cutoff_left <- as_num("Pcutoff1", 0.01)
padj_cutoff_right <- as_num("Pcutoff2", 0.01)
min_p_left <- as_num("minp1", 1e-100)
min_p_right <- as_num("minp2", 1e-100)

up_color <- paste0("#", param("upcolor", "F39B7F"))
down_color <- paste0("#", param("downcolor", "4DBBD5"))
no_color <- paste0("#", param("nocolor", "A6A6A6"))
label_color <- paste0("#", param("label_gene_color", "E64B35"))

left_x_label <- param("left_volcano_x_label", "HIPPO SCZ vs control")
right_x_label <- param("right_volcano_x_label", "DLPFC SCZ vs control")

fig_width <- as_num("fig_width", 10)
fig_height <- as_num("fig_height", 8)
label_fontsize <- as_num("label_gene_fontsize", 4.5)
text_fontsize <- as_num("text_fontsize", 14)
ticks_fontsize <- as_num("ticks_fontsize", 14)

# The current tables contain these transcript-level Ensembl IDs. This map
# makes the script runnable without an annotation download. If the first
# column already contains symbols, it is not used.
fallback_symbol_map <- tibble::tribble(
  ~gene_id,          ~symbol,
  "ENSG00000204308", "RNF5",
  "ENSG00000213722", "DDAH2",
  "ENSG00000204366", "ZBTB12",
  "ENSG00000204427", "ABHD16A",
  "ENSG00000204392", "LSM2",
  "ENSG00000204344", "STK19",
  "ENSG00000204410", "MSH5",
  "ENSG00000204371", "EHMT2"
)

# Read labelled genes from the JSON. The reference JSON contains Ensembl IDs,
# but a symbol-based labelled_gene_data field works without any modification.
raw_label_genes <- param(
  "basic_label_genes_linked_two_volcano_plot_labelled_gene_data",
  paste(c(
    "DDAH2", "ZBTB12", "RNF5", "LSM2",
    "ABHD16A", "MSH5", "EHMT2", "STK19"
  ), collapse = "\n")
)
raw_label_genes <- unique(trimws(unlist(strsplit(raw_label_genes, "[\r\n,;]+"))))
raw_label_genes <- raw_label_genes[nzchar(raw_label_genes)]

if (any(raw_label_genes %in% fallback_symbol_map$gene_id)) {
  label_genes <- unique(c(
    fallback_symbol_map$symbol[fallback_symbol_map$gene_id %in% raw_label_genes],
    raw_label_genes[!raw_label_genes %in% fallback_symbol_map$gene_id]
  ))
} else {
  label_genes <- raw_label_genes
}

# ---------------------------- data loading ----------------------------------
read_de_table <- function(path, tissue_name, min_p) {
  x <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x) < 3L) {
    stop("Expected at least 3 columns in: ", path)
  }

  # DLPFC uses Row.names and Hippocampus uses Gene in the first column.
  names(x)[1:3] <- c("gene_id", "log2FC", "padj")

  x %>%
    transmute(
      tissue = tissue_name,
      gene_id = sub("\\..*$", "", as.character(gene_id)),
      log2FC = as.numeric(log2FC),
      padj = as.numeric(padj)
    ) %>%
    filter(!is.na(gene_id), !is.na(log2FC), !is.na(padj)) %>%
    mutate(
      padj_plot = pmax(padj, min_p),
      neg_log10_padj = -log10(padj_plot),
      is_label = FALSE
    )
}

add_labels <- function(dat) {
  # Use the first column directly when it already contains gene symbols.
  if (any(dat$gene_id %in% label_genes)) {
    dat$symbol <- ifelse(dat$gene_id %in% label_genes, dat$gene_id, NA_character_)
  } else {
    dat <- dat %>%
      left_join(fallback_symbol_map, by = "gene_id")
  }

  dat %>%
    mutate(is_label = !is.na(symbol) & symbol %in% label_genes)
}

hippo <- read_de_table(hippo_file, "Hippocampus", min_p_left) %>%
  add_labels()
dlpfc <- read_de_table(dlpfc_file, "DLPFC", min_p_right) %>%
  add_labels()

# ---------------------------- plotting --------------------------------------
plot_volcano <- function(
    dat,
    title,
    x_limits,
    x_cutoff,
    p_cutoff,
    x_label,
    y_limits = NULL
) {
  dat <- dat %>%
    mutate(
      significance = case_when(
        padj < p_cutoff & log2FC >= x_cutoff ~ "Up",
        padj < p_cutoff & log2FC <= -x_cutoff ~ "Down",
        TRUE ~ "Not significant"
      )
    )

  labels <- dat %>%
    filter(is_label)

  p <- ggplot(dat, aes(x = log2FC, y = neg_log10_padj)) +
    geom_point(
      aes(color = significance),
      size = 0.8,
      alpha = 0.8,
      stroke = 0
    ) +
    geom_vline(
      xintercept = c(-x_cutoff, 0, x_cutoff),
      linetype = "dotted",
      linewidth = 0.45,
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dotted",
      linewidth = 0.45,
      color = "black"
    ) +
    scale_color_manual(
      values = c(
        "Down" = down_color,
        "Up" = up_color,
        "Not significant" = no_color
      ),
      breaks = c("Down", "Up", "Not significant"),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = pretty(x_limits, n = 5),
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      x = x_label,
      y = expression(-log[10] * "(padj)")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = text_fontsize / 1.1,
        margin = margin(b = 5)
      ),
      axis.title = element_text(size = text_fontsize / 1.1),
      axis.text = element_text(color = "black", size = ticks_fontsize / 1.5),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )

  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, expand = FALSE)
  }

  if (nrow(labels) > 0L) {
    p <- p +
      geom_point(
        data = labels,
        color = label_color,
        size = 1.7,
        inherit.aes = FALSE,
        aes(x = log2FC, y = neg_log10_padj)
      ) +
      geom_text_repel(
        data = labels,
        aes(
          x = log2FC,
          y = neg_log10_padj,
          label = symbol
        ),
        inherit.aes = FALSE,
        color = label_color,
        size = label_fontsize,
        direction = "y",
        hjust = 0,
        nudge_x = 0.25,
        box.padding = 0.25,
        point.padding = 0.15,
        force = 1.1,
        min.segment.length = 0,
        segment.color = label_color,
        segment.size = 0.35,
        max.overlaps = Inf,
        seed = 20260918
      )
  }

  p
}

# The JSON leaves the y-axis maxima blank, so let ggplot determine them from
# the data. The lower bound is zero, matching the supplied Figure 4g.
p_hippo <- plot_volcano(
  hippo,
  title = "Hippocampus",
  x_limits = c(-6.5, 2.5),
  x_cutoff = fc_cutoff_left,
  p_cutoff = padj_cutoff_left,
  x_label = left_x_label,
  y_limits = c(0, NA)
)

p_dlpfc <- plot_volcano(
  dlpfc,
  title = "DLPFC",
  x_limits = c(-7.5, 2.5),
  x_cutoff = fc_cutoff_right,
  p_cutoff = padj_cutoff_right,
  x_label = right_x_label,
  y_limits = c(0, NA)
)

figure_4g <- p_hippo + p_dlpfc +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(plot.background = element_rect(fill = "white", color = NA))

# ---------------------------- export ----------------------------------------
# Use the standard R PDF device rather than cairo_pdf, which may not be
# available in a headless/server R installation.
ggsave(
  filename = output_pdf,
  plot = figure_4g,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = "pdf"
)

ggsave(
  filename = output_png,
  plot = figure_4g,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("Wrote: ", output_pdf)
message("Wrote: ", output_png)
