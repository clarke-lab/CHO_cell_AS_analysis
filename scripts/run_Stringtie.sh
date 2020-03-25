
#!/usr/bin/env bash
if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input bam directory"
        echo "-p = num processors"
        echo "-g = GTF file"
        echo "-o = path for output BAMs"
        exit 2
fi
while getopts s:i:p:g:o: option
  do
    case "${option}"
      in
      s) SAMPLEID=${OPTARG};;
      i) INDIR=${OPTARG};;
      p) THREADS=${OPTARG};;
      g) GTF=${OPTARG};;
      o) OUTDIR=${OPTARG};;
    esac
done

if [ ! -d $OUTDIR/individual_gtfs ]; then
mkdir -p $OUTDIR/individual_gtfs
fi

/mnt/HDD2/colin/bin/stringtie-2.0.3.Linux_x86_64/stringtie \
-p $THREADS \
-G $GTF \
--rf \
-j 5 \
-o $OUTDIR/individual_gtfs/"$SAMPLEID".gtf \
$INDIR/"$SAMPLEID"Aligned.sortedByCoord.out.bam
