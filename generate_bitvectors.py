##########################################################
# Generate a bitvector for each alignment in the SAM file
# 
# Author: Kobie Kirven
##########################################################

# imports
import argparse
import subprocess
import os
import pysam
import re
import numpy as np
import matplotlib.pyplot as plt

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def parse_args():
    parser = argparse.ArgumentParser(description='Generate a bitvector for each alignment in the SAM file.')
    parser.add_argument('-s', '--sam', type=str, required=True, help='SAM file of alignments.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output bitvectors file.')
    return parser.parse_args()

def split_cigar(cigar):
    # Split the cigar string into a list of tuples
    # where the first element is the length and the
    # second element is the type of cigar operation
    # (M, I, D, N, S, H, P, =, X)
    cigar_list = []
    current = ''
    current_len = 0
    for i in range(len(cigar)):
        if cigar[i].isdigit():
            current_len = current_len * 10 + int(cigar[i])
        else:
            cigar_list.append((current_len, cigar[i]))
            current_len = 0
    return cigar_list

def complement_dna(dna):
    # Return the complement of a DNA sequence
    complement = ''
    for base in dna:
        if base == 'A':
            complement += 'T'
        elif base == 'T':
            complement += 'A'
        elif base == 'C':
            complement += 'G'
        elif base == 'G':
            complement += 'C'
    return complement

def get_query_seq(read):
    # Get the query sequence
    query_seq = str(read.query_sequence).upper()
    if read.is_reverse:
        query_seq = complement_dna(query_seq[::-1])
    return query_seq

def get_ref_seq(read):
    # Get the reference sequence
    ref_seq = str(read.get_reference_sequence()).upper()
    return ref_seq

def RepresentsInt(i):
    '''
    Helper function that returns True if `i` is an int.
    '''
    try:
        int(i)
        return True
    except ValueError:
        return False

def parseMD(mdstr):
    '''
    Returns a numpy array with positions for where a `fromMM` base is, in the
    reference genome. This can then be used to intersect with the query string
    to find all desired mutations.
    '''
    # Split MD String at every single [ACGT] or ^:
    mdSub = re.sub(r'([\\^]*[ACGT]+)[0]*', ' \\1 ', mdstr)
    mdSplit = re.split('[ ]+', mdSub)

    mutArr = np.array([]).astype(str)

    # Iterate over Array and replace all mutations from the MD string with the
    # letter of the corresponding reference.
    # eg: 2G1 will produce the numpy array: 'M M G M'
    # All ^[ACGT]* by "D" and the number with a corresponding stretch of "M"
    for i in range(len(mdSplit)):
        mdPos = mdSplit[i]
        if len(mdPos) > 0 and RepresentsInt(mdPos):
            mutArr = np.concatenate((mutArr, np.repeat("M", int(mdPos))))
        elif re.match('\\^', mdPos):
            mutArr = np.concatenate((mutArr, np.repeat('D', len(mdPos) - 1)))
        elif len(mdPos) == 1:
            mutArr = np.concatenate((mutArr, np.array([mdPos])))
        else:
            # I'm not yet quite sure, if this won't just break at some point.
            # In my BAM/SAM files I have seen rare cases with two consecutive
            # mismatches in the MD tag causing this series of ifs to report incorrect
            # positions if I don't catch this.
            mutArr = np.concatenate((mutArr, np.array(list(mdPos))))

    # Return mismatch positions
    return mutArr

def generate_bitvectors(sam, output):
    # Generate a bitvector for each alignment in the SAM file
    # To save space, the bitvectors will be represented in a 
    # similar way as a SAM cigar string. For example, the bitvector
    # 0001100100 would be represented as 3M2N2M2N where M represents
    # a match and N represents a mismatch.

    fn = open(output, 'w') # open output file

    # Parse the SAM file
    samfile = pysam.AlignmentFile(sam, check_sq=False)

    # Process each alignment in the SAM file
    for read in samfile.fetch():

        # Make sure that the read is mapped
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Get the MD tag for the alignment
        md = read.get_tag('MD')

        # get the positions of the mismatches
        mutArr = parseMD(md)

        # Write the mismatch string to a file
        fn.write(f'{read.qname}\t{"".join(list(mutArr))}\n')

    fn.close() # close output file
    samfile.close() # close SAM file

def main():
    args = parse_args()
    generate_bitvectors(args.sam, args.output)
    print(f"{bcolors.OKGREEN}Bitvectors generated.{bcolors.ENDC}")

if __name__ == '__main__':
    main()