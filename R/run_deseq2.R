#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: run_deseq2.r
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
## Info: Genes differentially expressed with +/- 1.5 fold difference, along with a
## BH adjusted p-value and baseMean >=100 are considered significant
##
## -----------------------------------------------------------------------------

# load
suppressMessages(library("DESeq2"))
suppressMessages(library("biomaRt"))
suppressMessages(library("dplyr"))
suppressMessages(library("xlsx"))

# input arguments
args <- commandArgs(TRUE)
count_dir <- args[1]
results_dir <- args[2]

print("The error - The query to the BioMart webservice returned an invalid result: biomaRt expected a character string of length 1 - can result from temporary unavailiability of BioMart. Rerun process_rmats.R or use a different host for the useMart function in rmats annotation.R")

ensembl <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "cgcrigri_gene_ensembl",
  host = "uswest.ensembl.org"
)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

count_file_names <- grep("counts", list.files(count_dir), value = T)

cell_type <- c("NTS", "NTS", "NTS", "NTS", 
               "TS", "TS", "TS", "TS")

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

# set 37C as the comparator condition
DESeq_data$condition <- relevel(DESeq_data$condition, "NTS")

# calculate differential expression using the DESeq wrapper function
suppressMessages(DESeq_data <- DESeq(DESeq_data))

# set differential expression criteria
de_results <- results(DESeq_data, 
                      lfcThreshold = 0, 
                      independentFiltering = T)

# retain significant results
sig_de_results <- subset(
  de_results,
  abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >= 100
)

sig_de_results <- sig_de_results[
  order(sig_de_results$log2FoldChange, decreasing = T),
]

biomart.out <- getBM(
  attributes = c(
    "ensembl_gene_id", "entrezgene_id",
    "external_gene_name", "description",
    "gene_biotype"
  ),
  filters = c(
    "ensembl_gene_id"
  ),
  values = rownames(sig_de_results),
  mart = ensembl,
  uniqueRows = TRUE
)

# finish annotation by determining missing symbols using NCBI
chok1_entrez <- readLines("reference_genome/chok1_ncbi_ids.txt")
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

# remove duplicates
ncbi_annotated <- ncbi_annotated[!duplicated(ncbi_annotated[, 1]), ]
rownames(ncbi_annotated) <- ncbi_annotated[, 1]
ncbi_annotated <- ncbi_annotated[, -1]

sig_de_results <- sig_de_results[rownames(ncbi_annotated), ]
column_names <- c("NCBI.Gene.ID", "Gene.Symbol", "Gene.Name", "Ensembl.Biotype", colnames(sig_de_results))
sig_de_results <- cbind(ncbi_annotated, sig_de_results)
colnames(sig_de_results) <- column_names

sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing = T), ]

# remove non-protein coding genes 
sig_de_results <-sig_de_results[sig_de_results$Ensembl.Biotype == "protein_coding",]

# save the de_results object
save(de_results, file = paste(results_dir, "/de_results.rData", sep=""))

fn <- paste(results_dir, "/DESeq2_analysis.xlsx",sep="")
suppressMessages(if (file.exists(fn)) {file.remove(fn)})

# save the results to Excel
write.xlsx(sig_de_results,
           file = fn,
           sheetName = "DE genes", append = TRUE
)

print(paste("Differential Expression results saved to", results_dir,"/DESeq2_analysis.xlsx"))
