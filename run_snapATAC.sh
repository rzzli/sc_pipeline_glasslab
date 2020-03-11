#!/bin/bash
#to run $./run_snapATAC.sh JB999
#if error ./configure : /bin/bash^M : bad interpreter [duplicate] occurs
#To fix, open your script with vi or vim and enter in vi command mode (key Esc), then type this:
#:set fileformat=unix
# :wq!
sample=$1
output="${sample}.snap"

samtools view possorted_bam.bam -H > possorted.header.sam
echo "done samtool view"
# create a bam file with the barcode embedded into the read name
cat <( cat possorted.header.sam ) \
<( samtools view possorted_bam.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
| samtools view -bS - > possorted_bam.snap.bam
echo "done embed barcode"

samtools sort -n -@ 10 -m 1G possorted_bam.snap.bam -o possorted_bam.snap.nsrt.bam
echo "sort done"

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
echo "finish downloading chrome size"
snaptools snap-pre  \
  --input-file=possorted_bam.snap.nsrt.bam  \
	--output-snap=${output}  \
	--genome-name=hg38  \
	--genome-size=hg38.chrom.sizes  \
	--min-mapq=30  \
	--min-flen=50  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=False  \
	--overwrite=True  \
	--max-num=20000  \
	--min-cov=500  \
	--verbose=True

echo "finish pre"

rm possorted_bam.bam
rm possorted.header.sam

echo "success"
