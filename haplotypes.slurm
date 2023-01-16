#!/bin/bash

#SBATCH --job-name=pf_haplotypes                  #This is the name of your job
#SBATCH --cpus-per-task=16                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=2G              #This is the memory reserved per core.

#SBATCH --time=00:30:00        #This is the time that your task will run
#SBATCH --qos=30min          #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/daubenbe/wagphi00/Documents/SNPs/slurm_out/haplotype_slurm.log    #This is the joined STDOUT and STDERR file

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=philipp.wagner@unibas.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################
module load Python/3.9.5-GCCcore-10.3.0

# module load Python/3.6.6-fosscuda-2018b
# module load fosscuda/2018b

#export your required environment variables below
#################################################
source ~/.bashrc
conda activate artic_new



# contig=NC_004318.2_dhfr
contig='NC_004329.3_dhps'
 # contig='XM_001350122.1_k13'
gene='dhps'
minqual='10'
minlen=2400
maxlen=3100
source='../dhps_refstrain_pass'

#add your command lines below
#############################

cd $HOME/Documents/SNPs/ref_exp


for currentbarcode in 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
do
	./../haplotypes.sh -b $currentbarcode -f ../${gene}.fasta -g ${gene} -p ../${gene}_primers.txt -r ${source} -i ${minlen} -a ${maxlen} -q ${minqual} &
done


wait




# concat snps
cat scriptfiles_*/*.sorted.bam.snps | sort -n | uniq > snp_list

# calculate haplotypes
# for currentbarcode in 13 14 15 16 17 18 19 20 21 22 23 24
for currentbarcode in 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
# for currentbarcode in 34 35 36 37 38 39
do
	python3 ../variance_calc.py scriptfiles_$currentbarcode/${gene}.sorted.bam -c $contig -g ../${gene}.jsonl -s snp_list -q $minqual -o ${gene}_${currentbarcode}.tsv  &
done

wait