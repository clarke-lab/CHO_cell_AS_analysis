
export mapped_dir=mapped
export gtf_file=reference_genome/ensembl_chok1.gtf
export assembled_transcripts_dir=stringtie_output

readlink -f $assembled_transcripts_dir/*.gtf >> $assembled_transcripts_dir/mergelist.txt 

stringtie --merge $assembled_transcripts_dir/mergelist.txt -o $assembled_transcripts_dir/stringtie.gtf -f 0.1  -c 10 

gffcompare -o $assembled_transcripts_dir/gffcompare -r $gtf_file $assembled_transcripts_dir/stringtie.gtf -N -M