## ====================
##  ggsashimi examples
## ====================

cp rmats_output_2/SAMPLE_1/REP_1/aligned.sorted.bam sashimi_plots/bam_files/REP31_1.bam&
cp rmats_output_2/SAMPLE_1/REP_2/aligned.sorted.bam sashimi_plots/bam_files/REP31_2.bam&
cp rmats_output_2/SAMPLE_1/REP_3/aligned.sorted.bam sashimi_plots/bam_files/REP31_3.bam&
cp rmats_output_2/SAMPLE_1/REP_4/aligned.sorted.bam sashimi_plots/bam_files/REP31_4.bam&
cp rmats_output_2/SAMPLE_2/REP_1/aligned.sorted.bam sashimi_plots/bam_files/REP37_1.bam&
cp rmats_output_2/SAMPLE_2/REP_2/aligned.sorted.bam sashimi_plots/bam_files/REP37_2.bam&
cp rmats_output_2/SAMPLE_2/REP_3/aligned.sorted.bam sashimi_plots/bam_files/REP37_3.bam&
cp rmats_output_2/SAMPLE_2/REP_4/aligned.sorted.bam sashimi_plots/bam_files/REP37_4.bam&


cp rmats_output_2/SAMPLE_1/REP_1/aligned.sorted.bam.bai sashimi_plots/bam_files/REP31_1.bam.bai&
cp rmats_output_2/SAMPLE_1/REP_2/aligned.sorted.bam.bai sashimi_plots/bam_files/REP31_2.bam.bai&
cp rmats_output_2/SAMPLE_1/REP_3/aligned.sorted.bam.bai sashimi_plots/bam_files/REP31_3.bam.bai&
cp rmats_output_2/SAMPLE_1/REP_4/aligned.sorted.bam.bai sashimi_plots/bam_files/REP31_4.bam.bai&
cp rmats_output_2/SAMPLE_2/REP_1/aligned.sorted.bam.bai sashimi_plots/bam_files/REP37_1.bam.bai&
cp rmats_output_2/SAMPLE_2/REP_2/aligned.sorted.bam.bai sashimi_plots/bam_files/REP37_2.bam.bai&
cp rmats_output_2/SAMPLE_2/REP_3/aligned.sorted.bam.bai sashimi_plots/bam_files/REP37_3.bam.bai&
cp rmats_output_2/SAMPLE_2/REP_4/aligned.sorted.bam.bai sashimi_plots/bam_files/REP37_4.bam.bai&


## Example #1. Overlay, intron shrinkage, gene annotation, PDF output, custom size and colors
python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/input_bams2.txt -c JH000890.1:567671-586157 -g stringtie_output/stringtie_merged.gtf -M 10 -C 3 -O 3 -A median --alpha 1 -F tiff -R 350 --base-size=16 --height=3 --width=18

sashimi_plots/bam_files/REP31_1.bam



python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/bam_files/REP31_2.bam -c JH000890.1:567671-586157 -g stringtie_output/stringtie_merged.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=2 --height=40 --width=18 -P sashimi_plots/palette.txt
## Example #2. Median coverage and number of reads supporting inclusion and exclusion, no gene annotation, TIFF output (350 PPI), custom size, default colors
../sashimi-plot.py -b input_bams.tsv -c chr10:27040584-27048100 -M 10 -C 3 -O 3 -A median --alpha 1 -F tiff -R 350 --base-size=16 --height=3 --width=18


python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/input_bams.tsv -c JH000890.1:567671-586157 -g stringtie_output/stringtie_merged.gtf -M 10 -C 3 -O 3 -A median --alpha 1 -F tiff -R 350 --base-size=16 --height=3 --width=18
 
python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/input_bams.tsv -c JH002291.1:7824-12636 -g stringtie_output/removed.overlapping.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.5 --base-size=20 --ann-height=4 --height=3 --width=18 -P sashimi_plots/palette.txt -F tiff -R 600 -o Dnm1l -A median_j
 
 

 
 52019
 
  python2 ~/bin/ggsashimi-master/sashimi-plot.py -b input_bams.tsv -c JH003517.1:58797-64053 -g stringtie_output/rmats_stringtie.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.5 --base-size=20 --ann-height=4 --height=3 --width=18 -P sashimi_plots/palette.txt -F tiff -R 600 -o sashimi_plots/Mff -A median_j
  
  
  python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/input_bams.tsv -c JH000890.1:567671-586157 -g stringtie_output/removed.overlapping.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.5 --base-size=20 --ann-height=4 --height=3 --width=18 -P sashimi_plots/palette.txt -F tiff -R 600 -o Hikeshi -A median_j
  
    python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/input_bams.tsv -c JH001768.1:241803-245625 -g stringtie_output/removed.overlapping.gtf -M 20 -C 3 -O 3 --shrink --alpha 0.5 --base-size=20 --ann-height=4 --height=3 --width=18 -P sashimi_plots/palette.txt -F tiff -R 600 -o Hikeshi -A median_j
  
  python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/input_bams.tsv -c JH000627.1:633914-640388 -g stringtie_output/stringtie_merged.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.5 --base-size=20 --ann-height=4 --height=3 --width=18 -P sashimi_plots/palette.txt -F tiff -R 600 -o Ubr2 -A median_j
  
  
  python2 ~/bin/ggsashimi-master/sashimi-plot.py -b sashimi_plots/input_bams.tsv -c JH000113.1:1892668-1909123 -g stringtie_output/stringtie_merged.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.5 --base-size=20 --ann-height=4 --height=3 --width=18 -P sashimi_plots/palette.txt -F tiff -R 600 -o Rnf34_a -A median_j