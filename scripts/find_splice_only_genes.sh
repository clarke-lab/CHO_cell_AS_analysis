

awk -F, '{print $1}' differential_expression/DESeq2_results/TS_v_NTS.csv | sed -e 's/\"//g' | sed -i '/^$/d' > manuscript_plots/de_gene.list

#grep -v -wFf manuscript_plots/de_gene.list rmats/filtered/SE_filtered_annotated.csv | awk -F, '{print $5}'

rm manuscript_plots/upset_plot/*
grep -v -wFf manuscript_plots/de_gene.list rmats/filtered/SE_filtered_annotated.csv \
| egrep 'protein_coding' \
| grep -v 'lincRNA' \
| awk -F, '{print $5}' |  sort | uniq > manuscript_plots/upset_plot/SE_splice_only.list



awk -F, '{print $5}'  rmats/filtered/SE_filtered_annotated.csv | sort | uniq > manuscript_plots/upset_plot/SE_only.list

grep -v -wFf manuscript_plots/de_gene.list rmats/filtered/MXE_filtered_annotated.csv \
| egrep 'protein_coding' \
| grep -v 'lincRNA' \
| awk -F, '{print $5}' | sort | uniq > manuscript_plots/upset_plot/MXE_splice_only.list

grep -v -wFf manuscript_plots/de_gene.list rmats/filtered/RI_filtered_annotated.csv \
| egrep 'protein_coding' \
| grep -v 'lincRNA' \
| awk -F, '{print $5}' | sort | uniq > manuscript_plots/upset_plot/RI_only.list

grep -v -wFf manuscript_plots/de_gene.list rmats/filtered/A5SS_filtered_annotated.csv \
| egrep 'protein_coding' \
| grep -v 'lincRNA' \
| awk -F, '{print $5}' | sort | uniq > manuscript_plots/upset_plot/A5SS_only.list

grep -v -wFf manuscript_plots/de_gene.list rmats/filtered/A3SS_filtered_annotated.csv \
| egrep 'protein_coding' \
| grep -v 'lincRNA' \
| awk -F, '{print $5}' | sort | uniq > manuscript_plots/upset_plot/A3SS_only.list

cat manuscript_plots/upset_plot/MXE_splice_only.list \
manuscript_plots/upset_plot/SE_only.list \
manuscript_plots/upset_plot/RI_only.list \
manuscript_plots/upset_plot/A5SS_only.list \
manuscript_plots/upset_plot/A3SS_only.list  > manuscript_plots/splice_only_genes.list



awk -F, '{print $5}'  rmats/filtered/SE_filtered_annotated.csv | sort | uniq > manuscript_plots/upset_plot/SE_only.list
awk -F, '{print $5}'  rmats/filtered/MXE_filtered_annotated.csv | sort | uniq > manuscript_plots/upset_plot/MXE_only.list
awk -F, '{print $5}'  rmats/filtered/RI_filtered_annotated.csv | sort | uniq > manuscript_plots/upset_plot/RI_only.list
awk -F, '{print $5}'  rmats/filtered/A5SS_filtered_annotated.csv | sort | uniq > manuscript_plots/upset_plot/A5SS_only.list
awk -F, '{print $5}'  rmats/filtered/A3SS_filtered_annotated.csv | sort | uniq > manuscript_plots/upset_plot/A3SS_only.list
