#!/usr/bin/Rscript

library(UpSetR)
library(xlsx)

# read in the DE gene ensembl ids
DE_file <- "DE_results/DESeq2_analysis.xlsx"
de_gene_list <- read.xlsx(DE_file, 1, header = TRUE, colIndex = 1, stringsAsFactors = FALSE)

AS_file <- "AS_results/rMats_analysis.xlsx"
se_gene_list <- read.xlsx(AS_file, 1, colIndex = 2, header = TRUE, stringsAsFactors = FALSE)
mxe_gene_list <- read.xlsx(AS_file, 2, colIndex = 2, header = TRUE, stringsAsFactors = FALSE)
ri_gene_list <- read.xlsx(AS_file, 3, colIndex = 2, header = TRUE, stringsAsFactors = FALSE)
a5ss_gene_list <- read.xlsx(AS_file, 4, colIndex = 2, header = TRUE, stringsAsFactors = FALSE)
a3ss_gene_list <- read.xlsx(AS_file, 5, colIndex = 2, header = TRUE, stringsAsFactors = FALSE)

listInput <- list(
  GENE = de_gene_list[, 1],
  SE = se_gene_list[, 1],
  MXE = mxe_gene_list[, 1],
  RI = ri_gene_list[, 1],
  A5SS = a5ss_gene_list[, 1],
  A3SS = a3ss_gene_list[, 1]
)

tiff("manuscript_plots/upset_plot/UpsetPlot.tiff",
     height = 20, width = 34, units = "cm",
     compression = "lzw", res = 300
)

upset(fromList(listInput),
      nsets = 6,
      order.by = "freq",
      empty.intersections = "on",
      sets.bar.color = "#227f74",
      main.bar.color = "#3b4154",
      text.scale = 2,
      mainbar.y.label = "Intersecting Genes",
      sets.x.label = "Total Genes",
      nintersects = 16
)
dev.off()


# upset(fromList(listInput), order.by = "freq")

save(listInput, file = "manuscript_plots/upset_plot/upset_data.rData")

se_symbol_list <- as.character(read.xlsx(AS_file, 2, header = TRUE)[, 7])
mxe_symbol_list <- as.character(read.xlsx(AS_file, 3, header = TRUE)[, 7])
ri_symbol_list <- as.character(read.xlsx(AS_file, 4, header = TRUE)[, 7])
a5ss_symbol_list <- as.character(read.xlsx(AS_file, 5, header = TRUE)[, 7])
a3ss_symbol_list <- as.character(read.xlsx(AS_file, 6, header = TRUE)[, 7])

as_symbols <- unique(c(se_symbol_list, mxe_symbol_list, ri_symbol_list, a5ss_symbol_list, a3ss_symbol_list))

de_gene_symbols <- read.table("differential_expression/DESeq2_results/TS_v_NTS.csv", sep = ",", header = T, stringsAsFactors = F)$Gene.Symbol
splicing_only_symbols <- setdiff(as_symbols, de_gene_symbols)

write.table(splicing_only_symbols, "manuscript_plots/splicing_only_symbols.txt")



