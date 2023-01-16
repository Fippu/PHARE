##################################################
## A SCRIPT TO EXTRACT REGIONS OF INTEREST FROM A BAM FILE
##################################################
## FREE USE
##################################################
## Author: Philipp Wagner
## Copyright: Copyright 2022
## Version: 1.0.0
## Email: philipp.wagner@unibas.ch
## Status: dev
##
##
## This script finds regions of interest (ROI) in a bam file
## This works by counting the Nucleotides at each position
## and then setting a position as ROI if a one variant (base)
## reaches a frequency of at least the --noise_cutoff threshold
##
## outputs a file with all the ROI nucleotide positions, one per row
## 
##################################################


import argparse
import pysam
import sys, csv
import json
from Bio import SeqIO
from os import path
from numpy import loadtxt 



#----------------------------------------------------------------------------------------
# TERMINAL UI
#----------------------------------------------------------------------------------------

# SET UP ARGUMETNS TO BE PASSED IN TERMINAL
parser = argparse.ArgumentParser(description="haplotype extractor",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", help="input fastq file")
# parser.add_argument("output", help="output filename")
parser.add_argument("-c", "--contig", type=str, help="reference name/contig")
parser.add_argument("-q", "--min_quality", type=int, help="minimal quality score for base to qualify for haplotype")
parser.add_argument("-n", "--noise_cutoff", type=float, default=0.1, help="default 0.1 (10%) minumum proportion of reads at a given position to be considered a real nucleotide variant instead of sequencing noise (error)")
parser.add_argument("-r", "--reference", type=str, help="optional; if a reference is given, ROIs are called also if a difference to the reference is obsered")
args = parser.parse_args()
config = vars(args)


infile = args.input
contig = args.contig
noise_cutoff = args.noise_cutoff
reference = args.reference

minqual = args.min_quality




#----------------------------------------------------------------------------------------
# READ FILES
#----------------------------------------------------------------------------------------

# input bam file containing reads to be analyzed
samfile = pysam.AlignmentFile(infile, "rb")

contiglen = samfile.get_reference_length(contig)
# contig = 'NC_004318.2_dhfr'
# contig = 'NC_004329.3_dhps'

# contiglen = 3

#----------------------------------------------------------------------------------------
# LOOP THROUGH THE SEQUENCE
# and store the nucleotide of each read at the SNP site
#----------------------------------------------------------------------------------------



# allreads = samfile.fetch(contig)
# readcount = 0
# for read in allreads:
#     readcount += 1
res = [{'A': 0, 'T': 0, 'G': 0, 'C': 0} for i in range(contiglen)]
# for read in allreads:
#     res[read.query_name] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}


# the output of this loop is cached. read cache if not set to true
if not path.isfile(infile + '.json'):
    # create iterable object from fastq file using sample
    roi = samfile.pileup(contig, start=0, end=contiglen, truncate=True, irgnore_overlaps=False, min_base_quality=minqual)

    for col in roi:
        # print(col.pos)
        # print("coverage at base %s = %s (filtered). %s reads at this position" % (col.pos, col.get_num_aligned(), col.n))

        # iterate through all reads at this snp site
        for row in col.pileups:
            read = row.alignment
            if isinstance(row.query_position, int): # not empty

                if read.query_qualities[row.query_position] > minqual:  # apply quality filter

                    res[col.pos][read.query_sequence[row.query_position]] += 1


    with open(infile + '.json', 'w') as f:
        json.dump(res, f)



# read cached file
else:
    with open(infile + '.json') as f:
        res = json.load(f)
    # print(col.keys())






#----------------------------------------------------------------------------------------
# GET THE REFERENCE SEQUENCE
#----------------------------------------------------------------------------------------

refseq = False
if reference != '':
    for record in SeqIO.parse(reference, "fasta"):
        if record.id == contig:
            refseq = record.seq
            break
    if not refseq:
        raise Exception('The contig does not match any contig in the reference')


# get snp sites
snpsites = []
pos = 0
for col in res: # go through all nucleotides
    tot = col['A']+col['T']+col['G']+col['C']
    counter = 0

    # count a base if its frequency is above noise
    for base in col.items(): # go through the bases
        if (base[1] > tot * noise_cutoff): 
            counter += 1

    # count if the reference base is not present
    if (refseq and col[refseq[pos]] < noise_cutoff):
        print('gotcha')
        counter += 1

    if counter > 1: # we have more than one base not considered noise
        snpsites.append(pos + 1) # pisam uses a 0 index, but generally 1 is used for the first nucleotide
    pos += 1


#----------------------------------------------------------------------------------------
# OUTPUT THE RESULT AS CSV
#----------------------------------------------------------------------------------------

# print(dict(list(res.items())[0:50]))

with open(infile + '.snps', 'w') as f:
    for pos in snpsites:
        f.write(str(pos) + '\n')

# for pos in snpsites:
#     print(pos)