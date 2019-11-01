SAMPLEID=$1
mkdir -p ../data/rmats/
mkdir -p ../data/rmats/trimmed
mkdir -p ../data/rmats/untrimmed

export rmats_fq_untrimmed_dir=../data/rmats/trimmed
export rmats_fq_trimmed_dir=../data/rmats/untrimmed
export preprocessed_dir=../data/preprocessed/paired

# FORWARD READ
cp $ preprocessed_dir/"$SAMPLEID"_1_sequence.txt.gz $rmats_fq_untrimmed_dir
gunzip $rmats_fq_untrimmed_dir/"$SAMPLEID"_1_sequence.txt.gz
python ~/bin/rMATS.3.2.5/bin/trimFastq.py  \
$rmats_fq_untrimmed_dir/"$SAMPLEID"_1_sequence.txt \
$rmats_fq_trimmed_dir/"$SAMPLEID"_1_sequence.txt 125
rm $rmats_fq_untrimmed_dir/"$SAMPLEID"_1_sequence.txt

# REVERSE READ
cp $preprocessed_data_directory/"$SAMPLEID"_2_sequence.txt.gz $rmats_fq_untrimmed_dir
gunzip $rmats_fq_untrimmed_dir/"$SAMPLEID"_2_sequence.txt.gz
python ~/bin/rMATS.3.2.5/bin/trimFastq.py  $rmats_fq_untrimmed_dir/"SAMPLEID"_2_sequence.txt $rmats_fq_trimmed_dir/"$SAMPLEID"_2_sequence.txt 125
rm $rmats_fq_untrimmed_dir/"$SAMPLEID"_2_sequence.txt
