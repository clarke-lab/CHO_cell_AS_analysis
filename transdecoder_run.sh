


~/bin/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl stringtie_output/rmats_stringtie.gtf reference_genome/ensembl_chok1_genome.fa > transdecoder_analysis/stringtie_transcripts.fasta
~/bin/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl stringtie_output/rmats_stringtie.gtf > transdecoder_analysis/transcripts.gff3

TransDecoder.LongOrfs -t transdecoder_analysis/stringtie_transcripts.fasta

mkdir pfam_db
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 
gunzip Pfam-A.hmm.gz 

hmmpress pfam_db/Pfam-A.hmm
hmmscan --cpu 28 --domtblout transdecoder_analysis/pfam.domtblout pfam_db/Pfam-A.hmm transdecoder_analysis/stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep -f transdecoder_analysis/pfam.std.out


mkdir transdecoder_analysis/uniref90
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz 

makeblastdb -in uniprot_sprot.fasta -dbtype prot

TransDecoder.Predict -t target_transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6