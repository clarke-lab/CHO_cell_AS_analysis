#!/usr/bin/Rscript

library(UpSetR)
library(xlsx)
de_gene_list<-read.table("manuscript_plots/de_gene.list", header =F, stringsAsFactors=F)$V1

AS_file="manuscript_plots/upset_plot/Table S5_final_cleaned.xlsx"
se_gene_list<-as.character(read.xlsx(AS_file, 2, header=TRUE)[,5])
mxe_gene_list<-as.character(read.xlsx(AS_file, 3, header=TRUE)[,5])
ri_gene_list<-as.character(read.xlsx(AS_file, 4, header=TRUE)[,5])
a5ss_gene_list<-as.character(read.xlsx(AS_file, 5, header=TRUE)[,5])
a3ss_gene_list<-as.character(read.xlsx(AS_file, 6, header=TRUE)[,5])

listInput <- list(
        GENE = de_gene_list,
        SE = se_gene_list,
        MXE = mxe_gene_list,
        RI = ri_gene_list,
        A5SS = a5ss_gene_list,
        A3SS =a3ss_gene_list)

#upset(fromList(listInput), order.by = "freq")



save(listInput,file="manuscript_plots/upset_plot/upset_data.rData")

se_symbol_list<-as.character(read.xlsx(AS_file, 2, header=TRUE)[,7])
mxe_symbol_list<-as.character(read.xlsx(AS_file, 3, header=TRUE)[,7])
ri_symbol_list<-as.character(read.xlsx(AS_file, 4, header=TRUE)[,7])
a5ss_symbol_list<-as.character(read.xlsx(AS_file, 5, header=TRUE)[,7])
a3ss_symbol_list<-as.character(read.xlsx(AS_file, 6, header=TRUE)[,7])

as_symbols<-unique(c(se_symbol_list,mxe_symbol_list,ri_symbol_list,a5ss_symbol_list,a3ss_symbol_list))

de_gene_symbols<-read.table("differential_expression/DESeq2_results/TS_v_NTS.csv",sep=",", header =T, stringsAsFactors=F)$Gene.Symbol
splicing_only_symbols<-setdiff(as_symbols, de_gene_symbols)

write.table(splicing_only_symbols,"manuscript_plots/splicing_only_symbols.txt")


tiff("manuscript_plots/upset_plot/UpsetPlot.tiff", height = 20, width = 34, units="cm",
     compression = "lzw", res = 300)

upset(fromList(listInput), nsets=6,
               order.by = "freq",
               empty.intersections = "on",
               sets.bar.color = "#227f74",
               main.bar.color ="#3b4154",
               text.scale=2,
               mainbar.y.label ="Intersecting Genes",
               sets.x.label = "Total Genes",
             nintersects=16)
dev.off()
