#!/bin/bash

# GET VARIABLES
helpFunction()
{
   echo ""
   echo "Usage:"
   echo ""
   echo "$0 -b Barcode -g gene reference; -r reads directory; -s SNPs"
   # echo - "\t-a Description of what is parameterA"
   exit 1 # Exit script after printing help
}

while getopts "b:f:g:p:r:s:a:i:q:?" opt
do
   case "$opt" in
      b ) currentbarcode="$OPTARG" ;;
      f ) generef="$OPTARG" ;;
      g ) gene="$OPTARG" ;;
      p ) primerfile="$OPTARG" ;;
      r ) readsdir="$OPTARG" ;;
      a ) maxlen="$OPTARG" ;;
      i ) minlen="$OPTARG" ;;
      q ) minqual="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$currentbarcode" ] || [ -z "$generef" ] || [ -z "$gene" ] || [ -z "$primerfile" ] || [ -z "$readsdir" ] || [ -z "$minlen" ] || [ -z "$maxlen" ] || [ -z "$minqual" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi


refid=$(head -n 1 $generef | perl -p -e 's/\>//g' | perl -p -e 's/[\r|\n]*//g')
echo "using contig $refid"






# INITIAL HOUSEKEEPING
rm -rf scriptfiles_$currentbarcode
mkdir scriptfiles_$currentbarcode

echo "target minlen maxlen
$gene $minlen $maxlen" > scriptfiles_$currentbarcode/cutofflength.txt


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )




# CONCATENATE AND FILTER READS
# cat $readsdir/barcode$currentbarcode/FAR*1*.fastq > scriptfiles_$currentbarcode/barcode$currentbarcode.fastq
 cat $readsdir/barcode$currentbarcode/*.fastq > scriptfiles_$currentbarcode/barcode$currentbarcode.fastq
# gzip -dc $readsdir/barcode$currentbarcode/*.fastq.gz > scriptfiles_$currentbarcode/barcode$currentbarcode.fastq


# APPLY HOMOPOLYMER COMPRESSION
# python3 $SCRIPT_DIR/hp_compr.py scriptfiles_$currentbarcode/barcode$currentbarcode.fastq scriptfiles_$currentbarcode/barcode$currentbarcode.compressed.fastq


# quality and length filtering
filtlong scriptfiles_$currentbarcode/barcode$currentbarcode.fastq --min_length $minlen --max_length $maxlen --min_mean_q $minqual --length_weight 10 > scriptfiles_$currentbarcode/barcode$currentbarcode.filtered.fastq

# filter for primers and more
SeekDeep extractor \
 --fastq scriptfiles_$currentbarcode/barcode$currentbarcode.filtered.fastq \
 --id $primerfile \
 --lenCutOffs scriptfiles_$currentbarcode/cutofflength.txt \
 --qualWindow 50,5,10 \
 --dout scriptfiles_$currentbarcode/sd_extractor \
 --overWriteDir \
 --barcodeErrors 2 \
 --primerWithinStart 100 \
 --checkRevComplementForPrimers \
 --checkRevComplementForMids \
 --primerNumOfMismatches 2 \
 --midWithinStart 100 \
 --primerCoverage 0.95




# MAP READS TO REFERENCE
minimap2 -a -t 8 -x map-ont --MD $generef scriptfiles_$currentbarcode/sd_extractor/$gene.fastq --sam-hit-only | samtools sort -T tmp -o scriptfiles_$currentbarcode/$gene.sorted.bam
samtools index scriptfiles_$currentbarcode/$gene.sorted.bam

# rm -rf scriptfiles_$currentbarcode/snp

# EXTRACT NUCLEOTIDES AT SNP POSITIONS
mkdir scriptfiles_$currentbarcode/snp

# loop through snp-file
# while read -r line; do
#   samtools tview scriptfiles_$currentbarcode/$gene.sorted.bam -d T -p $refid:$line | cut -c-1 | tail -n +4 > scriptfiles_$currentbarcode/snp/$line.txt
# done < $snpfile


# echo "please check if all reads contain the required SNPs"
# find scriptfiles_$currentbarcode/snp/*.txt | xargs wc -l


# paste scriptfiles_$currentbarcode/snp/*.txt | tr [:lower:] [:upper:] | grep -Po '^(?:[A|G|T|C][\t])*[A|G|T|C]$' | sort | uniq -c | sort -bgr | perl -p -e 's/(?<=\d) (?=\S)/\t/g' > ${gene}_${currentbarcode}.tsv
# SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# python3 $SCRIPT_DIR/variance_calc.py -c $refid scriptfiles_$currentbarcode/${gene}.sorted.bam -s $snpfile -q 10 | sort | uniq -c | sort -bgr | perl -p -e 's/(?<=\d) (?=\S)/\t/g' > python${gene}_${currentbarcode}10.tsv
# python3 $SCRIPT_DIR/variance_calc.py -c $refid scriptfiles_$currentbarcode/${gene}.sorted.bam -s $snpfile -q 9 | sort | uniq -c | sort -bgr | perl -p -e 's/(?<=\d) (?=\S)/\t/g' > python${gene}_${currentbarcode}9.tsv
# python3 $SCRIPT_DIR/variance_calc.py -c $refid scriptfiles_$currentbarcode/${gene}.sorted.bam -s $snpfile -q 8 | sort | uniq -c | sort -bgr | perl -p -e 's/(?<=\d) (?=\S)/\t/g' > python${gene}_${currentbarcode}8.tsv
# python3 $SCRIPT_DIR/variance_calc.py -c $refid scriptfiles_$currentbarcode/${gene}.sorted.bam -s $snpfile -q 7 | sort | uniq -c | sort -bgr | perl -p -e 's/(?<=\d) (?=\S)/\t/g' > python${gene}_${currentbarcode}7.tsv
# python3 $SCRIPT_DIR/variance_calc.py -c $refid scriptfiles_$currentbarcode/${gene}.sorted.bam -s $snpfile -q 6 | sort | uniq -c | sort -bgr | perl -p -e 's/(?<=\d) (?=\S)/\t/g' > python${gene}_${currentbarcode}6.tsv

python3 $SCRIPT_DIR/snp_selector.py scriptfiles_$currentbarcode/${gene}.sorted.bam -c $refid -q $minqual -n 0.05 -r $generef




# python3 $SCRIPT_DIR/variance_calc.py -c $refid scriptfiles_$currentbarcode/${gene}.sorted.bam -s $snpfile -q $minqual | sort | uniq -c | sort -bgr | perl -p -e 's/(?<=\d) (?=\S)/\t/g' > python${gene}_${currentbarcode}.tsv

