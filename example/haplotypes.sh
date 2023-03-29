#!/bin/bash


## LOADING THE ENVIRONMENT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate phare



## FUNCTION TO RUN THE 4 SCRIPTS AND CREATE DIRECTORIES
run_phare () {
    # SETUP FOLDERS
    savedir=${gene}${savedirsup}

    cd $wdir
    mkdir $savedir
    mkdir $wdir/export
    cd ${wdir}/${savedir}


    # FIND SNPS
    for currentbarcode in ${barcodes[@]}
    do
        # bash ${pharedir}/haplotypes.sh -s $targetlen -b $currentbarcode -f ${pharedir}/templates/${gene}.fasta -g ${gene} -p ${pharedir}/templates/${gene}_primers.txt -r ${source} -i ${minlen} -a ${maxlen} -q ${minqual} &
        bash ${pharedir}/scripts/preprocessing.sh -b $currentbarcode -f ${pharedir}/templates/${gene}.fasta -g ${gene} -p ${pharedir}/templates/${gene}_primers.txt -r ${source} -i ${minlen} -a ${maxlen} -q ${minqual} &
        wait
        python3 ${pharedir}/scripts/snp_selector.py scriptfiles_$currentbarcode/${gene}.sorted.bam -c $contig -q $minqual -n 0.05 -r ${pharedir}/templates/${gene}.fasta -m 0.8 &
    done
    wait


    # create a snp list
    cat scriptfiles_*/*.sorted.bam.snps | sort -n | uniq > snp_list


    # FITLER SILENT MUTATIONS
    for currentbarcode in ${barcodes[@]}
    do
        python3 ${pharedir}/scripts/filter_silent.py scriptfiles_$currentbarcode/${gene}.sorted.bam -c $contig -g ${pharedir}/templates/${gene}.jsonl -s snp_list -q $minqual -o ${gene}_${currentbarcode}.tsv -r ${pharedir}/templates/${gene}.fasta &
    done
    wait

    # create a new snp list without the silent mutations
    cat scriptfiles_*/*.sorted.bam.non-silent.snps | sort -n | uniq > snp_list

    # calculate haplotypes
    for currentbarcode in ${barcodes[@]}
    do
        python3 ${pharedir}/scripts/variance_calc.py scriptfiles_$currentbarcode/${gene}.sorted.bam -c $contig -g ${pharedir}/templates/${gene}.jsonl -s snp_list -q $minqual -o ${gene}_${currentbarcode}.tsv &
    done
    wait

    mkdir $wdir/export/${gene}${savedirsup}
    cp $wdir/$savedir/$gene*.tsv $wdir/export/${gene}${savedirsup}/
}




#############################

# SETUP
minqual=15

# GENE SETUP
gene='dhps'
contig='NC_004329.3_dhps'
targetlen=2742
minlen=$((targetlen - 50))
maxlen=3100

# DIRECTORIES
wdir=$(pwd)
pharedir=$wdir/../
source=$wdir/input # sample directory

# SAMPLES
barcodes=(17 18 19 20 21 22 23 24)

run_phare

Rscript analysis.R