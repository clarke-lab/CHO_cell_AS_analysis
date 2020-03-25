#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: run_DEseq2.rep
##
## Purpose of script: To carry out count based differential expression
##
## Author: NIBRT
##
## Date Created: Dec-2020
##
## Email: colin.clarke@nibrt.ie
##
## -----------------------------------------------------------------------------
##
## Info: Called from the filter_rmats_results.R script
##
## -----------------------------------------------------------------------------

# load
suppressMessages(library("DESeq2"))
suppressMessages(library(biomaRt))


args <- commandArgs(TRUE)
count_dir <- args[1]
results_dir <- args[2]
#dir.create(results_dir, showWarnings = F)

count_file_names <- grep("counts", list.files(count_dir), value = T)

cell_type <- c("NTS", "NTS", "NTS", "NTS", "TS", "TS", "TS", "TS")

sample_information <- data.frame(
  sampleName = count_file_names,
  fileName = count_file_names,
  condition = cell_type
)

DESeq_data <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_information,
  directory = count_dir,
  design = ~condition
)

colData(DESeq_data)$condition <- factor(colData(DESeq_data)$condition,
  levels = c("TS", "NTS")
)

rld <- rlogTransformation(DESeq_data, blind = T)
pdf(file = paste(results_dir, "PCA_TS_v_NTS.pdf", sep = "/"))
plotPCA(rld, intgroup = "condition")
suppressMessages(dev.off())
print("PCA plot saved in results folder")

# set 37C as the comparator
DESeq_data$condition <- relevel(DESeq_data$condition, "NTS")

# calculate differential expression using the DESeq wrapper function
suppressMessages(DESeq_data <- DESeq(DESeq_data))

# set differential expression criteria
de_results <- results(DESeq_data, lfcThreshold = 0, independentFiltering = T)
# order results by padj value (most significant to least)
sig_de_results <- subset(
  de_results,
  abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >= 100
)

sig_de_results <- sig_de_results[
  order(sig_de_results$log2FoldChange, decreasing = T),
]

save(de_results, sig_de_results, file = paste(results_dir,
  "de_results.RData",
  sep = "/"
))

# write.csv(as.matrix(sig_de_results),paste(results_dir,"significant_transcript_changes.csv",sep="/"))

ensembl <- useMart("ensembl")
ensembl <- useDataset("cgcrigri_gene_ensembl", mart = ensembl)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

chok1_entrez <- readLines("reference_genome/chok1_ncbi_ids.txt")

biomart.out <- getBM(
  attributes = c(
    "ensembl_gene_id", "entrezgene_id",
    "external_gene_name", "description",
    "gene_biotype"
  ),
  filters = c("ensembl_gene_id"),
  values = rownames(sig_de_results),
  mart = ensembl,
  uniqueRows = TRUE
)

ncbi_annotated <- matrix(, length(biomart.out[, 2]), 5)
ncbi_annotated[, 1] <- biomart.out[, 1]
ncbi_annotated[, 2] <- biomart.out[, 2]
ncbi_annotated[, 3] <- biomart.out[, 3]
ncbi_annotated[, 4] <- biomart.out[, 4]
ncbi_annotated[, 5] <- biomart.out[, 5]

pb <- txtProgressBar(
  min = 0, max = length(biomart.out[, 2]), style = 3,
  title = "completing annotation"
)
for (i in 1:length(biomart.out[, 2])) {
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, i)
  if (biomart.out[i, 3] == "") {
    if (is.na(biomart.out[i, 2])) {
      ncbi_annotated[i, 2] <- "Novel"
      ncbi_annotated[i, 3] <- "Novel"
      ncbi_annotated[i, 4] <- "Novel"
    } else if (length(strsplit(chok1_entrez[grep(biomart.out[i, 2],
      chok1_entrez,
      fixed = T
    )], "\t")) == 0) {
      ncbi_annotated[i, 2] <- "Unannotated"
      ncbi_annotated[i, 3] <- "Unannotated"
      ncbi_annotated[i, 4] <- "Unannotated"
    } else {
      ncbi_annotated[i, 3] <- strsplit(chok1_entrez[grep(biomart.out[i, 2],
        chok1_entrez,
        fixed = T
      )], "\t")[[1]][3]
    }
  }
}
close(pb)

ncbi_annotated <- ncbi_annotated[!duplicated(ncbi_annotated[, 1]), ]
rownames(ncbi_annotated) <- ncbi_annotated[, 1]
ncbi_annotated <- ncbi_annotated[, -1]

sig_de_results <- sig_de_results[rownames(ncbi_annotated), ]
column_names <- c("NCBI.Gene.ID", "Gene.Symbol", "Gene.Name", "Ensembl.Biotype", colnames(sig_de_results))
sig_de_results <- cbind(ncbi_annotated, sig_de_results)
colnames(sig_de_results) <- column_names


sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing = T), ]

# select only protein coding genes
sig_de_results <- sig_de_results[sig_de_results$Ensembl.Biotype == "protein_coding", ]
print("Differentially Expressed Genes +/- 1.5 Fold and P < 0.05 and baseMean >= 100")
print(summary(sig_de_results))
print("10 most upregulated genes")
print(head(as.data.frame(sig_de_results), n = 10))
print("10 most downregulated genes")
print(tail(as.data.frame(sig_de_results), n = 10))


write.csv(sig_de_results, file = paste(results_dir,"/TS_v_NTS.csv", sep="/"), row.names = T)
save(list = ls(), file = paste(results_dir,"TS_v_NTS.RData", sep="/"))
save(sig_de_results, file = paste(results_dir,"sig_de_results.RData", sep="/"))
print("DE gene list saved in results folder")
print("PCA plot saved in results folder")
print("Data saved in results folder")
quit()


print("PCA plot saved in results folder")
print("Data saved in results folder")
quit()
######################   END OF SCRIPT ######################################
