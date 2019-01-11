INDIR=$1


mkdir -p rmats_output

# append ensembl gene ids to MSTRG GTF
perl mstrg_prep.pl stringtie_output/stringtie_merged.gtf > stringtie_output/merged.test.delete

# find instances where stringtie has asemebled transcripts from 2 or more overlaping loci and created a new "gene". The final field of the GTF file will contain and MSTRG ID
# these can be separated from new 
grep 'MSTRG.*|ENSCGRG.*|ENSC.*' stringtie_output/merged.test.delete| grep '\<transcript\>' | awk '$NF ~/MSTRG/ {print $NF}'  > deleted.transcripts

# remove superious transcripts
grep -v -F -f deleted.transcripts stringtie_output/merged.test.delete > stringtie_output/removed.overlapping.gtf

# run rmats 4.0.2
 python ~/bin/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --s1 reps31.txt --s2 reps37.txt \
 --gtf stringtie_output/removed.overlapping.gtf  --bi reference_genome/star_index \
 --od rmats_output_3 -t paired --readLength 125 --nthread 32 --libType fr-firststrand  
 
 
 python2 ~/bin/rMATS.3.2.5/RNASeq-MATS.py -s1 \
../data/rmats/trimmed/REP31_1_R1.fastq:../data/rmats/trimmed/REP31_1_R2.fastq,\
../data/rmats/trimmed/REP31_2_R1.fastq:../data/rmats/trimmed/REP31_2_R2.fastq,\
../data/rmats/trimmed/REP31_3_R1.fastq:../data/rmats/trimmed/REP31_3_R2.fastq,\
../data/rmats/trimmed/REP31_4_R1.fastq:../data/rmats/trimmed/REP31_4_R2.fastq \
-s2 \
../data/rmats/trimmed/REP37_1_R1.fastq:../data/rmats/trimmed/REP37_1_R2.fastq,\
../data/rmats/trimmed/REP37_2_R1.fastq:../data/rmats/trimmed/REP37_2_R2.fastq,\
../data/rmats/trimmed/REP37_3_R1.fastq:../data/rmats/trimmed/REP37_3_R2.fastq,\
../data/rmats/trimmed/REP37_4_R1.fastq:../data/rmats/trimmed/REP37_4_R2.fastq \
-gtf stringtie_output/stringtie_merged.gtf -o rmats_output_2 -t paired -len 125 -libType fr-firststrand -bi reference_genome/star_index

mkdir -p rmats_output_2/filtered_results

# skipped exons
awk 'function abs(x){return ((x < 0.0) ? -x : x)} { if (($20 < 0.05) && (abs($23) > 0.15 )) { print } }'    rmats_output_2/MATS_output/SE.MATS.JunctionCountOnly.txt > rmats_output_2/filtered_results/SE.MATS.JunctionCountOnly.txt

awk '{print $2}' rmats_output_2/filtered_results/SE.MATS.JunctionCountOnly.txt



# retained introns
awk 'function abs(x){return ((x < 0.0) ? -x : x)} { if (($20 < 0.05) && (abs($23) > 0.1 )) { print } }'   $input_dir/RI.MATS.JunctionCountOnly.txt > $output_dir/RI.MATS.JunctionCountOnly.txt

awk 'function abs(x){return ((x < 0.0) ? -x : x)} { if (($20 < 0.05) && (abs($23) > 0.1 )) { print } }'    $input_dir/A3SS.MATS.JunctionCountOnly.txt > $output_dir/A3SS.MATS.JunctionCountOnly.txt





#alternative 5' splice sites
awk 'function abs(x){return ((x < 0.0) ? -x : x)} { if (($20 < 0.05) && (abs($23) > 0.1 )) { print } }'    $input_dir/A5SS.MATS.JunctionCountOnly.txt > $output_dir/A5SS.MATS.JunctionCountOnly.txt

#MXE
awk 'function abs(x){return ((x < 0.0) ? -x : x)} { if (($22 < 0.05) && (abs($25) > 0.1 )) { print } }'  $input_dir/MXE.MATS.JunctionCountOnly.txt > $output_dir/MXE.MATS.JunctionCountOnly.txt




grep -e '\<MSTRG.849\>\|\<transcript\>\|ref_gene_id'  stringtie_output/stringtie_merged.gtf 

grep -E '\bMSTRG.849\b.*ref_gene_id'  stringtie_output/stringtie_merged.gtf | awk '{print $NF}' | uniq | sed 's/"//g; s/;//g'
while read p; do
 current_id=$(echo "$p" | awk '{print $2}' | sed 's/"//g')
 ensembl_id=$(grep -E "\<${current_id}\>"  stringtie_output/stringtie_merged.gtf | grep '\<ref_gene_id\>' | awk '{print $NF}' | uniq | sed 's/"//g; s/;//g')
 echo $ensembl_id  $current_id
done < rmats_output_2/MATS_output/SE.MATS.JunctionCountOnly.txt



rmats.se.df<-as.data.frame(temperature_top@SE_gr$exon_target)
search_list<-with(rmats.se.df[2,], list(gsub("chr", "", as.character(seqnames)), start, end,as.character(strand)))
if(search_list[[4]][1]=="-"){search_list[[4]][1]<-"-1"} else if(search_list[[4]][1]=="+"){search_list[[4]][1]<-"+1"}
biomart.out<-getBM(
    attributes=c('ensembl_gene_id','external_gene_name','entrezgene', 'chromosome_name', 'transcript_start','strand', 'transcript_end'), 
    filters = c('chromosome_name', 'start', 'end'), 
    values = with(rmats.se.df[3,], list(gsub("chr", "", as.character(seqnames)), start, end)), 
    mart = ensembl
)
biomart.out

cgr_symbols<-getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                                       filters = c("chromosomal_region","biotype"),
                                       values = list(chromosomal_region="JH000074.1:2053861:2053994",biotype="protein_coding"),
                                       uniqueRows = F,
                                       mart = ensembl)