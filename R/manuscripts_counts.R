

SE_filtered<-read.csv("rmats/filtered/SE_filtered_annotated.csv")
annotated_protein_only_SE <-SE_filtered[SE_filtered$Ensembl.Biotype == "protein_coding" | SE_filtered$Ensembl.Biotype == "protein_coding|protein_coding",]
SE_up<-sum(annotated_protein_only_SE$Delta.PSI > 0)
SE_down<-sum(annotated_protein_only_SE$Delta.PSI < 0)

write(paste("Total annotated SE events",dim(annotated_protein_only_SE)[1],sep="="),stdout())
write(paste("SEs upregulated at 31C",SE_up,sep="="),stdout())
write(paste("SEs upregulated at 37C",SE_down,sep="="),stdout())


MXE_filtered<-read.csv("rmats/filtered/MXE_filtered_annotated.csv")
annotated_protein_only_MXE <-MXE_filtered[MXE_filtered$Ensembl.Biotype == "protein_coding" | MXE_filtered$Ensembl.Biotype == "protein_coding|protein_coding",]
MXE_up<-sum(annotated_protein_only_MXE$Delta.PSI > 0)
MXE_down<-sum(annotated_protein_only_MXE$Delta.PSI < 0)

write(paste("Total annotated MXE events",dim(annotated_protein_only_MXE)[1],sep="="),stdout())
write(paste("MXEs upregulated at 31C",MXE_up,sep="="),stdout())
write(paste("MXEs upregulated at 37C",MXE_down,sep="="),stdout())

RI_filtered<-read.csv("rmats/filtered/RI_filtered_annotated.csv")
annotated_protein_only_RI <-RI_filtered[RI_filtered$Ensembl.Biotype == "protein_coding" | RI_filtered$Ensembl.Biotype == "protein_coding|protein_coding",]
RI_up<-sum(annotated_protein_only_RI$Delta.PSI > 0)
RI_down<-sum(annotated_protein_only_RI$Delta.PSI < 0)

write(paste("Total annotated RI events",dim(annotated_protein_only_RI)[1],sep="="),stdout())
write(paste("RIs upregulated at 31C",RI_up,sep="="),stdout())
write(paste("RIs upregulated at 37C",RI_down,sep="="),stdout())

A5SS_filtered<-read.csv("rmats/filtered/A5SS_filtered_annotated.csv")
annotated_protein_only_A5SS <-A5SS_filtered[A5SS_filtered$Ensembl.Biotype == "protein_coding" | A5SS_filtered$Ensembl.Biotype == "protein_coding|protein_coding",]
A5SS_up<-sum(annotated_protein_only_A5SS$Delta.PSI > 0)
A5SS_down<-sum(annotated_protein_only_A5SS$Delta.PSI < 0)

write(paste("Total annotated A5SS events",dim(annotated_protein_only_A5SS)[1],sep="="),stdout())
write(paste("A5SS upregulated at 31C",A5SS_up,sep="="),stdout())
write(paste("A5SS upregulated at 37C",A5SS_down,sep="="),stdout())


A3SS_filtered<-read.csv("rmats/filtered/A3SS_filtered_annotated.csv")
annotated_protein_only_A3SS <-A3SS_filtered[A3SS_filtered$Ensembl.Biotype == "protein_coding" | A3SS_filtered$Ensembl.Biotype == "protein_coding|protein_coding",]
A3SS_up<-sum(annotated_protein_only_A3SS$Delta.PSI > 0)
A3SS_down<-sum(annotated_protein_only_A3SS$Delta.PSI < 0)

write(paste("Total annotated A3SS events",dim(annotated_protein_only_A3SS)[1],sep="="),stdout())
write(paste("A3SS upregulated at 31C",A3SS_up,sep="="),stdout())
write(paste("A3SS upregulated at 37C",A3SS_down,sep="="),stdout())

all_psi_values<-c(annotated_protein_only_SE$Delta.PSI,
      annotated_protein_only_MXE$Delta.PSI,
      annotated_protein_only_RI$Delta.PSI,
      annotated_protein_only_A3SS$Delta.PSI,
      annotated_protein_only_A5SS$Delta.PSI)

write(paste("maximum PSI", max(all_psi_values),sep="="),stdout())
write(paste("min PSI", min(all_psi_values),sep="="),stdout())

total_events<-sum(dim(annotated_protein_only_SE)[1],dim(annotated_protein_only_MXE)[1],dim(annotated_protein_only_RI)[1],dim(annotated_protein_only_A5SS)[1],dim(annotated_protein_only_A3SS)[1])
write(paste("total AS events detected", total_events,sep="="),stdout())



de_genes<-read.csv("differential_expression/DESeq2_results/TS_v_NTS.csv")
write(paste("Total DE genes",dim(de_genes)[1],sep="="),stdout())
