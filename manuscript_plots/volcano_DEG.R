#!/usr/bin/Rscript
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library("DESeq2"))

# import DESeq2 results
load("DE_results/de_results.rData")
  
  # filter by expressionsasas
  de_results <- de_results[de_results$baseMean >= 100, ]
  
  # add significant
  de_results$threshold_DE <- de_results$padj < 0.05 & abs(de_results$log2FoldChange) >= 0.5849625
  
  DE_event <- cbind(de_results$log2FoldChange, de_results$padj, as.logical(de_results$threshold_DE))
  
  # data frame for all required information
  DE_event <- data.frame(
    ensemblid = rownames(de_results),
    log2FoldChange = de_results$log2FoldChange,
    padj = de_results$padj,
    threshold = de_results$threshold_DE
  )
  
  DE_event <- DE_event %>% mutate(pointcolor = case_when(
    log2FoldChange > 0 & threshold == 1 ~ "#F9B288",
    log2FoldChange < 0 & threshold == 1 ~ "#A2D9F9",
    threshold == 0 ~ "gray"
  ))
  DE_event <- DE_event %>% mutate(pointclass = case_when(
    log2FoldChange > 0 & threshold == 1 ~ "Increase at 31C",
    log2FoldChange < 0 & threshold == 1 ~ "Increase at 37C",
    threshold == 0 ~ "Not Significant"
  ))
  DE_event <- DE_event %>% mutate(pointsize = case_when(
    log2FoldChange > 0 & threshold == 1 ~ "0.5",
    log2FoldChange < 0 & threshold == 1 ~ "0.5",
    threshold == 0 ~ "0.2"
  ))
  
  
  g <- ggplot(DE_event) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), color = pointcolor, fill = pointclass), size = 0.1) +
    scale_colour_manual(
      name = "the colour",
      breaks = as.factor(DE_event$pointclass),
      values = c("#F9B288", "#A2D9F9", "gray"),
      labels = as.factor(DE_event$pointclass)
    ) +
    xlab(bquote(~ Log[2] ~ "fold change")) +
    ylab("-log10 BH adjusted p-value") +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "gray"
    ) +
    geom_vline(
      xintercept = c(-0.5849625, 0.5849625),
      linetype = "dashed",
      color = "light gray"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25))
    )
  
  g <- g + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7), legend.position = "bottom", legend.title = element_blank())
  
  ggsave(paste("DE_results/DEG_volcano.tiff"), plot = g, height = 6, width = 3, units = "in", dpi = 600)
  write(paste("DEG volcano_plot complete"), stdout())
