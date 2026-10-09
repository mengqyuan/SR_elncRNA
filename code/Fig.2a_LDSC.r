## Forest plot for TSS ±10kb S-LDSC enrichment (MVP SCZ EUR, baseline v1.2)


outdir <- "/mnt/results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---- Build data from our ±10kb results ----
d <- data.frame(
  lncRNA_set = c("WBG_lncRNA", "LBG_lncRNA", "Brain_elncRNA", "eQTL_lncRNA", "SR_elncRNA"),
  Enrichment = c(1.3944074, 11.5813393, 1.7080762, 63.8473081, 71.9465422),
  Enrichment_std_error = c(0.1237006, 0.7700634, 0.3107796, 11.0516527, 28.7914248),
  Enrichment_p = c(0.0017590, 8.12873e-30, 0.0231705, 6.14988e-08, 0.0144628),
  `#SNPs` = c(1202010, 134151, 236568, 5376, 1598),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

d$SNPs <- as.numeric(d[["#SNPs"]])

## Order: WBG, LBG, Brain, eQTL, SR
order_levels <- c("WBG_lncRNA", "LBG_lncRNA", "Brain_elncRNA", "eQTL_lncRNA", "SR_elncRNA")
d <- d[match(order_levels, d$lncRNA_set), ]

d$lower <- pmax(d$Enrichment - 1.96 * d$Enrichment_std_error, 0)
d$upper <- d$Enrichment + 1.96 * d$Enrichment_std_error

d$label <- c("WBG lncRNAs", "LBG lncRNAs", "Brain elncRNAs", "eQTL lncRNAs", "SR elncRNAs")

## Colors (same as user-provided code)
pal <- c(
  WBG_lncRNA     = "#619CFF",
  LBG_lncRNA     = "#B79F00",
  Brain_elncRNA  = "#00BFC4",
  eQTL_lncRNA    = "#00BA38",
  SR_elncRNA     = "#F564E3"
)

## Save source data
source_out <- file.path(outdir, "Figure2a_LDSC_enrichment_forest_tss10k_source_data.tsv")
write.table(d[, c("lncRNA_set", "label", "Enrichment", "Enrichment_std_error",
                   "lower", "upper", "Enrichment_p", "SNPs")],
            source_out, sep = "\t", quote = FALSE, row.names = FALSE)

## ---- Forest plot function (same style as user code) ----
plot_forest <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(family = "Arial", mar = c(3.8, 6.0, 1.0, 0.8), xpd = FALSE, lend = "butt")

  y <- rev(seq_len(nrow(d)))
  dplot <- d
  dplot$y <- y

  xlim <- c(0, ceiling(max(dplot$upper) / 20) * 20)
  ylim <- c(0.35, nrow(d) + 0.75)
  plot(NA, NA, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")

  ## Alternating row shading
  for (i in seq_len(nrow(dplot))) {
    if (i %% 2 == 1) {
      rect(xlim[1], dplot$y[i] - 0.39, xlim[2], dplot$y[i] + 0.39,
           col = "#F7F7F7", border = NA)
    }
  }

  ## Reference line at enrichment = 1
  abline(v = 1, lty = 2, lwd = 0.75, col = "#7F7F7F")

  ## X axis
  axis(1,
       at = seq(0, xlim[2], by = 20),
       labels = seq(0, xlim[2], by = 20),
       cex.axis = 0.76, lwd = 0.55, lwd.ticks = 0.55, padj = -0.25)
  mtext("LDSC enrichment score (95% CI)", side = 1, line = 2.3, cex = 0.86)

  ## Y axis (gene set labels)
  axis(2, at = dplot$y, labels = dplot$label, las = 1, tick = FALSE,
       cex.axis = 0.82, line = -0.2)

  ## CI bars and points
  for (i in seq_len(nrow(dplot))) {
    set <- dplot$lncRNA_set[i]
    col <- pal[[set]]
    segments(dplot$lower[i], dplot$y[i], dplot$upper[i], dplot$y[i],
             col = col, lwd = 2.6)
    segments(dplot$lower[i], dplot$y[i] - 0.12, dplot$lower[i], dplot$y[i] + 0.12,
             col = col, lwd = 1.1)
    segments(dplot$upper[i], dplot$y[i] - 0.12, dplot$upper[i], dplot$y[i] + 0.12,
             col = col, lwd = 1.1)
    points(dplot$Enrichment[i], dplot$y[i], pch = 21, bg = col, col = "black",
           lwd = 0.45, cex = ifelse(set == "SR_elncRNA", 1.35, 1.20))
  }

  ## Panel label
  text(xlim[1], nrow(d) + 0.48, "a", font = 2, cex = 1.0, adj = c(0, 0.5))
  text(xlim[1], nrow(d) + 0.16, "lncRNA set", font = 2, cex = 0.78, adj = c(0, 0.5))

  box(bty = "l", lwd = 0.55)
}

## ---- Output PDF (editable, Arial embedded) ----
base <- file.path(outdir, "Figure2a_LDSC_enrichment_forest_tss10k")

grDevices::cairo_pdf(paste0(base, ".pdf"), width = 3.55, height = 2.80,
                     family = "Arial", onefile = TRUE)
plot_forest()
dev.off()

cat("PDF=", paste0(base, ".pdf"), "\n", sep = "")
cat("SOURCE=", source_out, "\n", sep = "")
