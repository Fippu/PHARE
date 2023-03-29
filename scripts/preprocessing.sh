#!/bin/bash

# GET VARIABLES
helpFunction()
{
   echo ""
   echo "Usage:"
   echo "See script source"
   echo "$0"
   # echo - "\t-a Description of what is parameterA"
   exit 1 # Exit script after printing help
}

while getopts "b:f:g:p:r:s:a:i:q:s:?" opt
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
      s ) seekdeep="$OPTARG" ;;
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


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )




# CONCATENATE AND FILTER READS
# cat $readsdir/barcode$currentbarcode/FAR*1*.fastq > scriptfiles_$currentbarcode/barcode$currentbarcode.fastq
 cat $readsdir/barcode$currentbarcode/*.fastq > scriptfiles_$currentbarcode/barcode$currentbarcode.fastq
# gzip -dc $readsdir/barcode$currentbarcode/*.fastq.gz > scriptfiles_$currentbarcode/barcode$currentbarcode.fastq


# quality and length filtering
filtlong scriptfiles_$currentbarcode/barcode$currentbarcode.fastq --min_length $minlen --max_length $maxlen --min_mean_q $minqual --length_weight 10 > scriptfiles_$currentbarcode/barcode$currentbarcode.filtered.fastq


if [ -z "$seekdeep" ]
then
   minimap2 -a -t 8 -x map-ont --MD $generef scriptfiles_$currentbarcode/barcode$currentbarcode.filtered.fastq --sam-hit-only | samtools sort -T tmp -o scriptfiles_$currentbarcode/$gene.sorted.bam
else
   echo "using SeekDeep with a target amplicon length of "$seekdeep
   # filter for primers and more
   SeekDeep extractor \
    --fastq scriptfiles_$currentbarcode/barcode$currentbarcode.filtered.fastq \
    --id $primerfile \
    --maxLen `expr $seekdeep \* 11 / 10` \
    --minLen `expr $seekdeep \* 9 / 10` \
    --qualWindow 50,5,$minqual \
    --dout scriptfiles_$currentbarcode/sd_extractor \
    --overWriteDir \
    --primerWithinStart $minlen \
    --checkRevComplementForPrimers \
    --primerNumOfMismatches 3

   minimap2 -a -t 8 -x map-ont --MD $generef scriptfiles_$currentbarcode/sd_extractor/$gene.fastq --sam-hit-only | samtools sort -T tmp -o scriptfiles_$currentbarcode/$gene.sorted.bam
fi




# MAP READS TO REFERENCE
samtools index scriptfiles_$currentbarcode/$gene.sorted.bam