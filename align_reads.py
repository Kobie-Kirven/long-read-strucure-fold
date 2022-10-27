###################################################################
# Align long reads to a reference transcriptome and output a 
# SAM file
#
# Author: Kobie Kirven
# Assmann Lab, Penn State University
###################################################################

# imports
import argparse
import subprocess

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
    parser = argparse.ArgumentParser(description='Align long reads to a reference transcriptome and output a SAM file.')
    parser.add_argument('-r', '--reads', type=str, required=True, help='Fasta file of long reads (may be in GZIP format)')
    parser.add_argument('-t', '--transcriptome', type=str, required=True, help='Fasta file of the transcriptome (may be in GZIP format)')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file name.')
    return parser.parse_args()

def align_reads(reads, transcriptome, output):
    # Align the reads to the transcriptome
    # Output a SAM file
    cmd = f"minimap2 -ax splice:hq -uf {transcriptome} {reads} > {output} 2> /dev/null"
    print(f"{bcolors.WARNING}Running minimap2...{bcolors.ENDC}")
    print(f"{bcolors.WARNING}Command --> {bcolors.ENDC}{bcolors.OKGREEN}{cmd}{bcolors.ENDC}")
    subprocess.run(cmd, shell=True)
    print(f"{bcolors.WARNING}Done!{bcolors.ENDC}")

def main():
    # Parse the arguments
    args = parse_args()

    # Align the reads to the transcriptome
    align_reads(args.reads, args.transcriptome, args.output)

if __name__ == "__main__":
    main()