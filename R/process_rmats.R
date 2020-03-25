#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: process_rmats.R
##
## Purpose of script: Filter the output from the rMats algorithm and annotate the results
##
## Author: NIBRT
##
## Date Created: Dec-2020
##
## Email: colin.clarke@nibrt.ie
##
## -----------------------------------------------------------------------------
##
## Info: each splicing event must have 10 reads spanning the junction, an absolute
## PSI >= 0.1 and an FDR < 0.05 to be considered signifciant
## Each rMats results is annotated using bioMart, to resolve transcripts with sense
## overlap the exon coordinates are searched
## ---------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(xlsx))
suppressMessages(library(biomaRt))
suppressMessages(library(stringr))

# input arguments
args <- commandArgs(TRUE)
rmats_dir <- args[1]
results_dir <- args[2]

source("R/rmats_annotation.R")

print("The error - The query to the BioMart webservice returned an invalid result: biomaRt expected a character string of length 1 - can result from temporary unavailiability of BioMart. Rerun process_rmats.R or use a different host for the useMart function in rmats annotation.R")

# filteres rMats AS events with mean junction coverage < 10
coverage_filter <- function(rmats_data_frame) {
  merged_counts <- paste(rmats_data_frame$IJC_SAMPLE_1,
    rmats_data_frame$SJC_SAMPLE_1,
    rmats_data_frame$IJC_SAMPLE_2,
    rmats_data_frame$SJC_SAMPLE_2,
    sep = ","
  )

  counts <- as.matrix(cbind(data.frame(do.call(
    "rbind",
    strsplit(as.character(merged_counts), ",", fixed = TRUE)
  ))))

  class(counts) <- "numeric"
  rownames(counts) <- rownames(rmats_data_frame)
  selected_IDs <- rownames(counts[rowMeans(counts) >= 10, ])
  rmats_data_frame <- rmats_data_frame[selected_IDs, ]
  return(rmats_data_frame)
}

# coverage, FDR and PSI
filter_rmats_results <- function(rmats_dir, event_type) {
  # import unfiltered results
  unfiltered <- read.table(
    file = paste(rmats_dir, "/", event_type, ".MATS.JC.txt", sep = ""),
    header = T,
    sep = "\t",
    row.names = 1,
    stringsAsFactors = FALSE
  )

  # retain events with 10 or more reads spanning the included and skipped junctions
  coverage_filtered <- coverage_filter(unfiltered)

  # keep all the events that have sufficent cover whether significant or not
  save(coverage_filtered, file=paste(rmats_dir,"/",event_type,".rData",sep=""))

  # retain events with an rMats FDR < 0.05 & a change in the PSI of 10%
  filtered <- coverage_filtered %>%
    filter(FDR < 0.05) %>%
    filter(abs(IncLevelDifference) >= 0.1)

  return(filtered)
}

# determine the significant splicing events and annotate
rmats_events <- c("SE", "MXE", "RI", "A3SS", "A5SS")

fn <- paste(results_dir, "/rMats_analysis.xlsx",sep="")
suppressMessages(if (file.exists(fn)) {file.remove(fn)})

for (i in 1:length(rmats_events)) {
  psi_filtered <- filter_rmats_results(rmats_dir, rmats_events[i])
  annotated_output <- rmats_annotation(psi_filtered, rmats_events[i])
  write.xlsx(annotated_output,
    file = paste(results_dir,"/rMats_analysis.xlsx",sep=""),
    sheetName = rmats_events[i], append = TRUE
  )
  print(paste("Processing:", rmats_events[i], ": complete", sep = ""))
}
print("Splicing events saved to AS_results/rMats_analysis.xlsx")
