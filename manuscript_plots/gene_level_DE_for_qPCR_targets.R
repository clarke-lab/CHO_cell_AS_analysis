#!/usr/bin/Rscript
suppressMessages(library(DESeq2))

load("differential_expression/DESeq2_results/de_results.RData")

goi<-data.frame(ensembl=c(
              "ENSCGRG00000017431",
              "ENSCGRG00000017026",
              "ENSCGRG00000017440",
              "ENSCGRG00000003601",
              "ENSCGRG00000016072",
              "ENSCGRG00000001449",
              "ENSCGRG00000016730",
              "ENSCGRG00000014828",
              "ENSCGRG00000011168",
              "ENSCGRG00000007344",
              "ENSCGRG00000013689",
              "ENSCGRG00000014636"),
              symbol=c(
              "Immp1l",
              "Mff",
              "Slirp",
              "Gpx8",
              "Dnm1l",
              "Dnajc15",
              "Hikeshi",
              "Rnf34",
              "Ubr2",
              "Fars2",
              "Cirbp",
              "Rbm3"
            ))

deseq2_output<-cbind(goi$symbol,data.frame(de_results[goi$ensembl, c(1,2,6)]))
