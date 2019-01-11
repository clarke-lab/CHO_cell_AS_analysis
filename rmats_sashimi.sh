sample1_rep1=mapping/REP37_1Aligned.sortedByCoord.out.bam
sample1_rep2=mapping/REP37_2Aligned.sortedByCoord.out.bam
sample1_rep3=mapping/REP37_3Aligned.sortedByCoord.out.bam
sample1_rep4=mapping/REP37_4Aligned.sortedByCoord.out.bam
sample2_rep1=mapping/REP31_1Aligned.sortedByCoord.out.bam
sample2_rep2=mapping/REP31_2Aligned.sortedByCoord.out.bam
sample2_rep3=mapping/REP31_3Aligned.sortedByCoord.out.bam
sample2_rep4=mapping/REP31_4Aligned.sortedByCoord.out.bam

mkdir sashimi_plots

output_dir="sashimi_plots"

awk '$1~/^5617$/ {print $0}' rMats_results/unfiltered/MXE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Rnf34"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $23; $23 = $24; $24 = t; print; } ' OFS=$'\t' >  sashimi_plots/Rnf34.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Rnf34.txt -t MXE --l1 31 --l2 37 --exon_s 1 --intron_s 20 -o ./sashimi_plots/Rnf34 --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^1884$/ {print $0}' rMats_results/unfiltered/MXE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Fars2"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $23; $23 = $24; $24 = t; print; } ' OFS=$'\t' >  sashimi_plots/Fars2.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Fars2.txt -t MXE --l1 31 --l2 37 --exon_s 1 --intron_s 30 -o ./sashimi_plots/Fars2 --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^3173$/ {print $0}' rMats_results/unfiltered/MXE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Ubr2"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $23; $23 = $24; $24 = t; print; } ' OFS=$'\t' >  sashimi_plots/Ubr2.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Ubr2.txt -t MXE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Ubr2 --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^47046$/ {print $0}' rMats_results/unfiltered/SE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Dnm1l"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' >  sashimi_plots/Dnm1l.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Dnm1l.txt -t SE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Dnm1l --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^23432$/ {print $0}' rMats_results/unfiltered/SE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Gpx8"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' >  sashimi_plots/Gpx8.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Gpx8.txt -t SE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Gpx8 --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^6999$/ {print $0}' rMats_results/unfiltered/SE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Hikeshi"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' >  sashimi_plots/Hikeshi.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Hikeshi.txt -t SE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Hikeshi --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^21371$/ {print $0}' rMats_results/unfiltered/SE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Immp1l"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' >  sashimi_plots/Immp1l.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Immp1l.txt -t SE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Immp1l --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^11430$/ {print $0}' rMats_results/unfiltered/SE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Dnajc15"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' >  sashimi_plots/Dnajc15.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Dnajc15.txt -t SE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Dnajc15 --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^12207/ {print $0}' rMats_results/unfiltered/SE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Mff"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' >  sashimi_plots/Mff.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Mff.txt -t SE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Mff --group-info grp.txt --color "#ee6363","#436eee"

awk '$1~/^50250/ {print $0}' rMats_results/unfiltered/SE.MATS.JC.txt | sed 's/chr//g' | awk '$3="Slirp"' FS="\t" OFS="\t" | awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' >  sashimi_plots/Slirp.txt
rmats2sashimiplot --b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  --b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 -e sashimi_plots/Slirp.txt -t SE --l1 31 --l2 37 --exon_s 1 --intron_s 5 -o ./sashimi_plots/Slirp --group-info grp.txt --color "#ee6363","#436eee"

