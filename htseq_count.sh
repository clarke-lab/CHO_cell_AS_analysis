SAMPLEID=$1
mkdir -p htseq_counts

mapped_dir="mapped"
outdir=htseq_counts
gtf_file=reference_genome/ensembl_chok1.gtf

htseq-count -r pos -f bam -i gene_id -s reverse \
            $mapped_dir/"$SAMPLEID"Aligned.sortedByCoord.out.bam \
            $gff_file > outdir/"$SAMPLEID".counts
            
# END 