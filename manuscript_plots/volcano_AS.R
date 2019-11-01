#!/usr/bin/Rscript
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

load("rmats/filtered/r_data_files/SE_filtered_annotated.rData")
SE<-SE_coverage_filt

load("rmats/filtered/r_data_files/MXE_filtered_annotated.rData")
MXE<-MXE_coverage_filt

load("rmats/filtered/r_data_files/RI_filtered_annotated.rData")
RI<-RI_coverage_filt

load("rmats/filtered/r_data_files/A5SS_filtered_annotated.rData")
A5SS<-A5SS_coverage_filt

load("rmats/filtered/r_data_files/A3SS_filtered_annotated.rData")
A3SS<-A3SS_coverage_filt

threshold_SE <- SE$FDR < 0.05 & abs(SE$IncLevelDifference) >= 0.1
SE$threshold <- threshold_SE
SE$type<-"SE"

threshold_MXE <- MXE$FDR < 0.05 & abs(MXE$IncLevelDifference) >= 0.1
MXE$threshold <- threshold_MXE
MXE$type<-"MXE"

threshold_RI <- RI$FDR < 0.05 & abs(RI$IncLevelDifference) >= 0.1
RI$threshold <- threshold_RI
RI$type<-"RI"

threshold_A5SS <- A5SS$FDR < 0.05 & abs(A5SS$IncLevelDifference) >= 0.1
A5SS$threshold <- threshold_A5SS
A5SS$type<-"A5SS"

threshold_A3SS <- A3SS$FDR < 0.05 & abs(A3SS$IncLevelDifference) >= 0.1
A3SS$threshold <- threshold_A3SS
A3SS$type<-"A3SS"

AS_event<-rbind(
cbind(SE$type, as.numeric(SE$IncLevelDifference), as.numeric(SE$FDR), SE$threshold),
cbind(MXE$type, as.numeric(MXE$IncLevelDifference), as.numeric(MXE$FDR), MXE$threshold),
cbind(RI$type, as.numeric(RI$IncLevelDifference), as.numeric(RI$FDR), RI$threshold),
cbind(A5SS$type, as.numeric(A5SS$IncLevelDifference), as.numeric(A5SS$FDR), A5SS$threshold),
cbind(A3SS$type, as.numeric(A3SS$IncLevelDifference), as.numeric(A3SS$FDR), A3SS$threshold)
)

AS_event<-data.frame(type=c(SE$type,MXE$type,RI$type,A5SS$type, A3SS$type),
                            dPSI=as.numeric(c(SE$IncLevelDifference,
                            MXE$IncLevelDifference,
                            RI$IncLevelDifference,
                            A5SS$IncLevelDifference,
                            A3SS$IncLevelDifference)),
                            FDR=c(SE$FDR,
                            MXE$FDR,
                            RI$FDR,
                            A5SS$FDR,
                            A3SS$FDR),
                            threshold=c(
                            SE$threshold,
                            MXE$threshold,
                            RI$threshold,
                            A5SS$threshold,
                            A3SS$threshold))



AS_event=AS_event %>% mutate(pointcolor = case_when(dPSI>0 & threshold=="TRUE" ~ "#F9B288",
                                                 dPSI<0 & threshold=="TRUE" ~ "#A2D9F9",
                                                 threshold=="FALSE" ~ "gray"))
AS_event=AS_event %>% mutate(pointclass = case_when(dPSI>0 & threshold=="TRUE" ~ "Increase at 31C",
                                                 dPSI<0 & threshold=="TRUE" ~ "Increase at 37C ",
                                                 threshold=="FALSE" ~ "Not Significant"))
AS_event=AS_event %>% mutate(pointsize = case_when(dPSI>0 & threshold=="TRUE" ~ "0.5",
                                                 dPSI<0 & threshold=="TRUE" ~ "0.5",
                                                 threshold=="FALSE" ~ "0.2"))

g <- ggplot(AS_event) +
        geom_point(aes(x=dPSI,
                       y=-log10(FDR),
                       color=pointcolor,
                       fill=pointclass),
                       size=0.1)  +
        scale_colour_manual(name = 'the colour',
                            breaks = as.factor(AS_event$pointclass),
                            values = c("#F9B288","#A2D9F9","gray"),
                            labels = as.factor(AS_event$pointclass)) +
        xlab("dPSI") +
        ylab("-log10 FDR") +
        scale_x_continuous(limits = c(-0.8,0.8)) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "gray") +
        geom_vline(xintercept =c(-0.1,0.1), linetype="dashed", color = "light gray") +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))

AS_event$type_f = factor(AS_event$type, levels=c('SE','MXE','RI','A5SS','A3SS'))

g <- g + facet_grid(. ~ type) +theme(legend.position = "bottom", )
g <- g + theme_bw() + theme(panel.grid.major =element_blank(),
                            panel.grid.minor = element_blank(),
                            strip.background = element_blank(),
                            strip.text = element_text(face="bold", size=9),
                            legend.position="bottom",
                            legend.title=element_blank())

ggsave(paste("manuscript_plots/","AS_volcano_plot.tiff"),
              plot=g,
              height=3,
              width=6,
              units='in',
              dpi=600)

write(paste("AS volcano_plot complete"),stdout())
