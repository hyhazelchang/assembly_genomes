#!/usr/bin/python3

# vcf_to_fillgap.py

# Hsin-Ying Chang <dnx202138@gmail.com>
# v1 2023/07/16

# Usage: python3 /home/xinchang/pyscript_xin/vcf_to_fillgap.py --vcf_file=/scratch/xinchang/cyano13/cyano13.01/bwa_19/PE.vcf --len_file=/scratch/xinchang/cyano13/cyano13.01/bwa_19/cyano13.01.19.fasta.length --thread=100 --output=/scratch/xinchang/cyano13/cyano13.01/scaffolding/fill.17.gap


import argparse

def main():
    parser = argparse.ArgumentParser(
            description=("Create the fill gap by vcf calling"),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf_file",
                        type=str,
                        default=None,
                        help="input file vcf. Please give the absolute path.")
    parser.add_argument("--len_file",
                        type=str,
                        default=None,
                        help="input length file. Please give the absolute path.")
    parser.add_argument("--thread",
                        type=float,
                        default=None,
                        help="the cutoff")
    parser.add_argument("--output",
                        type=str,
                        default=None,
                        help="output file to fill gap. Please give the absolute path.")
    
    # Defining variables from input
    args = parser.parse_args()
    vcf_file = args.vcf_file
    len_file = args.len_file
    thread = args.thread
    output = args.output

    # parse input file
    lines = [line.rstrip("\n").split("\t") for line in open(vcf_file)]
    lengths = [line.rstrip("\n").split("\t") for line in open(len_file)]

    # get line
    for line in lines:
        if line[0] == '#CHROM':
            index = lines.index(line) + 1
            break
    
    out = open(output, "w")
    while index < len(lines):
        if float(lines[index][5]) > thread:
            chrom = lines[index][0]
            pos_from = lines[index][1]
            pos_to = str(int(lines[index][1]) + int(len(lines[index][3])) - 1)
            alt_len = str(len(lines[index][3]))
            alt = lines[index][4]
            for length in lengths:
                if length[0] == chrom:
                    out.write(chrom + "\t" + length[1] + "\t1\t" + pos_from + "\t" + pos_to + "\t" + alt_len + "\tY\tM\t" + alt + "\n")
                    break
        index += 1
    out.close()

if __name__ == "__main__":
    main()
