# create the bed files for annotation track plotting

# first find the Dnm1l and Gpx8 in the GTF and parse

grep 'MSTRG.19199' stringtie/stringtie.gtf | convert2bed --input=gtf - > manuscript_plots/annotation_tracks/dnm1l.bed
grep 'MSTRG.18036' stringtie/stringtie.gtf | convert2bed --input=gtf - > manuscript_plots/annotation_tracks/gpx8.bed

grep 'MSTRG.19199' stringtie/stringtie.gtf | egrep -w 'exon' | sed 's/[";]//g;' | awk '{OFS="\t"; print $1, $4,$5,$12,$8,$7,$3}' > manuscript_plots/annotation_tracks/dnm1l.bed
grep 'MSTRG.18036' stringtie/stringtie.gtf | egrep -w 'exon' | sed 's/[";]//g;' | awk '{OFS="\t"; print $1, $4,$5,$12,$8,$7,$3}' > manuscript_plots/annotation_tracks/gpx8.bed
