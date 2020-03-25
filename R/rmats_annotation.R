#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: rmats_annotation.R
##
## Purpose of script: To annotate the rmats results
##
## Author: NIBRT
##
## Date Created: Dec-2020
##
## Email: colin.clarke@nibrt.ie
##
## -----------------------------------------------------------------------------
##
## Info: the scripts takes care of the additional coordinates outputs for MXE
## events. Where an Ensmebl annotation does not exsist NCBI is used
## host = "uswest.ensembl.org"
## -----------------------------------------------------------------------------

rmats_annotation <- function(psi_filtered, event_type) {
  ensembl <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "cgcrigri_gene_ensembl",
    host = "uswest.ensembl.org"
  )
  filters <- listFilters(ensembl)
  attributes <- listAttributes(ensembl)

  # remove the exsisting geneSymbol and ID.1 column
  psi_filtered <- psi_filtered[, !(colnames(psi_filtered) %in% c("geneSymbol", "ID.1"))]

  # deal with different rmats output for event type, MXE rMats coordinates are placed in different columns
  if (event_type == "RI") {
    coordinate_string1 <- "riExonStart_0base"
    coordinate_string2 <- "riExonEnd"
  } else if (event_type == "SE") {
    coordinate_string1 <- "exonStart_0base"
    coordinate_string2 <- "exonEnd"
  } else if (event_type == "A5SS" | event_type == "A3SS") {
    coordinate_string1 <- "longExonStart_0base"
    coordinate_string2 <- "longExonEnd"
  } else if (event_type == "MXE") {
    coordinate_string1 <- "X1stExonStart_0base"
    coordinate_string2 <- "X1stExonEnd"
  }

  # make the rMats results coordinates compatiable with Ensembl
  psi_filtered$chr <- gsub("chr", "", psi_filtered$chr, fixed = T)
  psi_filtered$strand <- gsub("+", "+1", psi_filtered$strand, fixed = T)
  psi_filtered$strand <- gsub("-", "-1", psi_filtered$strand, fixed = T)

  # split the rMats output depending on the ID
  # where an rMats has a single ENS gene id we can simply run bioMart R
  # where there is multiple ids, we will use the genomic coordinates to determine the gene id
  # any unannotated or non

  # row only has the MSTRG, no Ensembl annotation - removed from final results
  mstrg_only <- psi_filtered %>%
    filter(str_count(GeneID, "\\|") == 0 & str_count(GeneID, "MSTRG") == 1)

  ## Stage 1: single ensembl ID, query by Ensembl gene ID ##

  # 1) Ensembl ID only
  ensembl_without_mstrg <- psi_filtered %>%
    filter(str_count(GeneID, "\\|") == 0 & str_count(GeneID, "ENS") == 1) # row has the ensembl id

  ensembl_without_mstrg <- cbind.data.frame(ensembl_without_mstrg$GeneID, ensembl_without_mstrg, stringsAsFactors = F)

  colnames(ensembl_without_mstrg)[1] <- "ensembl_id"

  # 2) Ensembl ID with MSTRG
  single_ensembl_id <- psi_filtered %>%
    filter(str_count(GeneID, "\\|") == 1) # mstrg with only 1 ensembl id

  split <- strsplit(single_ensembl_id$GeneID, "|", fixed = TRUE)

  # stringtie_id <- sapply(split, "[", 1)
  ensembl_id <- sapply(split, "[", 2)

  single_ensembl_id <- cbind.data.frame(ensembl_id, single_ensembl_id, stringsAsFactors = F) # add the clean ensembl id

  # merge the rMats matrices
  rmats_single_ens_id <- rbind(single_ensembl_id, ensembl_without_mstrg)

  # combine IDs for bioMart
  biomart_query_1 <- c(rmats_single_ens_id$ensembl_id, ensembl_without_mstrg$ensembl_id)

  # query biomart
  biomart_output_1 <- getBM(
    attributes = c(
      "ensembl_gene_id", "entrezgene_id",
      "external_gene_name", "description",
      "gene_biotype"
    ),
    filters = c(
      "ensembl_gene_id"
    ),
    values = biomart_query_1,
    mart = ensembl,
    uniqueRows = TRUE
  )

  # some times we get duplicates due to 2 different entrez_gene_ids for the same gene
  biomart_output_1 <- biomart_output_1 %>% distinct(ensembl_gene_id, .keep_all = TRUE)

  # join the rmats results and biomart annotation
  biomart_output_1 <- merge(rmats_single_ens_id, biomart_output_1, by.x = "ensembl_id", by.y = "ensembl_gene_id")

  # where a gene symbol is not provided Ensembl, use NCBI to aquire the gene symol
  chok1_entrez <- readLines("reference_genome/chok1_ncbi_ids.txt")
  for (i in 1:dim(biomart_output_1)[1]) {
    if (biomart_output_1$external_gene_name[i] == "") {
      biomart_output_1$external_gene_name[i] <- strsplit(chok1_entrez[grep(biomart_output_1$entrezgene_id[i], chok1_entrez, fixed = T)], "\t")[[1]][3]
    }
  }

  # remove the non-protein coding genes
  annotated_single_ensid <- biomart_output_1 %>% filter(str_detect(gene_biotype, "protein_coding"))

  annotated_single_ensid <- annotated_single_ensid %>% dplyr::select(-gene_biotype)
  dim(annotated_single_ensid)

  ## Stage 2: multiple ensembl ID, query by exon_coordinates ##
  # now the second biomart query
  # this time because we have the appended Id with more than 1 ENSEMBL ID, the exon coordinates to are used to query biomart
  # a result is considered valid only if there is no sense overlap for the spliced exon

  # rMats results with multiple Ensembl IDs
  rmats_multiple_ens_id <- psi_filtered %>%
    filter(str_count(GeneID, "\\|") > 1)

  # make a frame with the coordinates, use the specifc string for the target exon
  coordinates <- cbind(
    as.character(rmats_multiple_ens_id$chr),
    rmats_multiple_ens_id[, coordinate_string1],
    rmats_multiple_ens_id[, coordinate_string2],
    as.character(rmats_multiple_ens_id$strand)
  )

  # data frame to hold results
  annotations <- data.frame(matrix(, dim(rmats_multiple_ens_id)[1], ncol = 7))
  rownames(annotations) <- rownames(rmats_multiple_ens_id)
  colnames(annotations) <- c(
    "GeneID",
    "ensembl_gene_id",
    "entrezgene_id",
    "external_gene_name",
    "description",
    "gene_biotype",
    "Sense Overlapping"
  )
  annotations$GeneID <- rmats_multiple_ens_id$GeneID
  if (dim(rmats_multiple_ens_id)[1] > 0) {
    #print(dim(rmats_multiple_ens_id)[1])
    for (i in 1:dim(rmats_multiple_ens_id)[1]) {
    
       biomart.out <- getBM(
        attributes = c(
          "ensembl_gene_id", "entrezgene_id",
          "external_gene_name", "description",
          "gene_biotype"
        ),
        filters = c(
          "chromosome_name", "start",
          "end", "strand"
        ),
        values = list(
          coordinates[i, 1], coordinates[i, 3],
          coordinates[i, 2], coordinates[i, 4]
        ),
        mart = ensembl,
        uniqueRows = TRUE
      )

      # deal with getting multiple hits per input
      if (dim(biomart.out)[1] == 1) {
        annotations[i, 2:7] <- cbind(as.matrix(biomart.out), "No")
      } else if (dim(biomart.out)[1] == 0) {
        annotations[i, 2:7] <- rep("Unannotated", 6)
      } else if (dim(biomart.out)[1] > 1) {
        if (biomart.out$ensembl_gene_id[1] == biomart.out$ensembl_gene_id[2]) {
          annotations[i, 2:7] <- cbind(as.matrix(biomart.out[1, ]), "No") # if there's two results returned with the same Ens gene ID keep the first
        } else if (is.na(biomart.out$entrezgene)[1] == 1) {
          annotations[i, 2:7] <- cbind(
            as.matrix(biomart.out[!is.na(biomart.out$entrezgene), ]),
            "No"
          )
        } else {
          annotations[i, 2:7] <- cbind(
            paste(biomart.out[, 1], collapse = "|"),
            paste(biomart.out[, 2], collapse = "|"),
            paste(biomart.out[, 3], collapse = "|"),
            paste(biomart.out[, 4], collapse = "|"),
            paste(biomart.out[, 5], collapse = "|"),
            "Yes"
          )
        }
      }
    }
    if (is.na(annotations[i, 4])) {
      annotations[i, 4] <- strsplit(chok1_entrez[grep(annotations[i, 3], chok1_entrez, fixed = T)], "\t")[[1]][3]
      if (is.na(annotations[i, 4])) {
        annotations[i, 2:6] <- "Unannotated"
      }
    }
  }

  # remove dupliates before merging
  annotations <- annotations %>%
    distinct(GeneID, .keep_all = TRUE)
  # merge
  biomart_output_2 <- merge(rmats_multiple_ens_id, annotations, by.x = "GeneID", by.y = "GeneID")

  # remove non-protien coding and sense overlapping
  annotated_multiple_ensid <- biomart_output_2 %>%
    filter(str_detect(gene_biotype, "protein_coding")) %>%
    filter(str_detect(`Sense Overlapping`, "No"))

  annotated_multiple_ensid <- annotated_multiple_ensid %>%
    dplyr::select(-`Sense Overlapping`) %>%
    dplyr::select(-gene_biotype) %>%
    dplyr::select(ensembl_gene_id, everything())
  colnames(annotated_multiple_ensid)[1] <- "ensembl_id"

  ## Stage 3: Merge annotated results

  # are columns and reformat column names
  annotated <- rbind(annotated_single_ensid, annotated_multiple_ensid)


  if (event_type == "MXE") {
    annotated <- annotated %>% dplyr::select(
      ensembl_id,
      entrezgene_id,
      external_gene_name,
      description,
      IncLevelDifference,
      PValue,
      FDR,
      IncLevel1,
      IncLevel2,
      IJC_SAMPLE_1,
      IJC_SAMPLE_2,
      SJC_SAMPLE_1,
      SJC_SAMPLE_2,
      IncFormLen,
      SkipFormLen,
      chr,
      strand,
      X1stExonStart_0base,
      X1stExonEnd,
      X2ndExonStart_0base,
      X2ndExonEnd,
      upstreamES,
      upstreamEE,
      downstreamES,
      downstreamEE
    )
  } else if (event_type == "SE") {
    annotated <- annotated %>% dplyr::select(
      ensembl_id,
      entrezgene_id,
      external_gene_name,
      description,
      IncLevelDifference,
      PValue,
      FDR,
      IncLevel1,
      IncLevel2,
      IJC_SAMPLE_1,
      IJC_SAMPLE_2,
      SJC_SAMPLE_1,
      SJC_SAMPLE_2,
      IncFormLen,
      SkipFormLen,
      chr,
      strand,
      exonStart_0base,
      exonEnd,
      upstreamES,
      upstreamEE,
      downstreamES,
      downstreamEE
    )
  } else if (event_type == "A5SS" | event_type == "A3SS") {
    annotated <- annotated %>% dplyr::select(
      ensembl_id,
      entrezgene_id,
      external_gene_name,
      description,
      IncLevelDifference,
      PValue,
      FDR,
      IncLevel1,
      IncLevel2,
      IJC_SAMPLE_1,
      IJC_SAMPLE_2,
      SJC_SAMPLE_1,
      SJC_SAMPLE_2,
      IncFormLen,
      SkipFormLen,
      chr,
      strand,
      longExonStart_0base,
      longExonEnd,
      shortES,
      shortEE,
      flankingES,
      flankingEE
    )
  } else {
    annotated <- annotated %>% dplyr::select(
      ensembl_id,
      entrezgene_id,
      external_gene_name,
      description,
      IncLevelDifference,
      PValue,
      FDR,
      IncLevel1,
      IncLevel2,
      IJC_SAMPLE_1,
      IJC_SAMPLE_2,
      SJC_SAMPLE_1,
      SJC_SAMPLE_2,
      IncFormLen,
      SkipFormLen,
      chr,
      strand,
      riExonStart_0base,
      riExonEnd,
      upstreamES,
      upstreamEE,
      downstreamES,
      downstreamEE
    )
  }

  # rename some columns
  annotated <- annotated %>% rename(
    Ensembl.Gene.ID = ensembl_id,
    Entrez.Gene.ID = entrezgene_id,
    Gene.Symbol = external_gene_name,
    NCBI.Gene.Name = description,
    Delta.PSI = IncLevelDifference,
    P.value = PValue,
    TS.norm.counts = IncLevel1,
    NTS.norm.counts = IncLevel2,
    Scaffold = chr,
    Strand = strand
  )
  # sort output by PSI
  annotated <- annotated %>% arrange(-Delta.PSI)

  return(annotated)
}
