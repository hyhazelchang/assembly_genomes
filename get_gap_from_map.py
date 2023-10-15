#!/usr/bin/python3

# get_gap_from_map.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2023/02/22
# v2 2023/07/16

# Usage: python3 /home/xinchang/pyscript_xin/get_gap_from_map.py --input=/scratch/xinchang/cyano13/cyano13.01/bwa_14/uni.bed --prefix=cyano13.01 --seq_name=cyano13.01.01 --seq_len=2658240 > /scratch/xinchang/cyano13/cyano13.01/agp/cyano13.01.14.agp1

# if the system can't be installed the pandas module, get pip installer under user by this cmd: wget https://bootstrap.pypa.io/pip/3.6/get-pip.py (according to diiferent version of your python, this one is 3.6)
# then install pip under the user: python get-pip.py --user
# after installing pip, use pip to install pandas under the user. so go to the bin: cd .local/bin; ./pip install pandas --user

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
            description=("Get the gaps from alignment between velvet contigs and reference genome"),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input",
                        type=str,
                        default=None,
                        help="input file (bed file). Please give the absolute path.")
    parser.add_argument("--seq_name",
                        type=str,
                        default=None,
                        help="sequence name")
    parser.add_argument("--prefix",
                        type=str,
                        default=None,
                        help="prefix")
    parser.add_argument("--seq_len",
                        type=str,
                        default=None,
                        help="sequence length")           

    # Defining variables from input
    args = parser.parse_args()
    in_file = args.input
    name = args.seq_name
    prefix = args.prefix
    length = args.seq_len

    # Read the bed file into dataframe
    df = pd.read_csv(in_file, header = None, sep = "\t")

    # delete the repeat rows from top to bottom
    up, down = 0, 1
    while down < len(df.index):
        if df.iloc[up, 1] <= df.iloc[down, 1] and df.iloc[down, 1] < df.iloc[up, 2] and df.iloc[up, 1] < df.iloc[down, 2] and df.iloc[down, 2] < df.iloc[up, 2]:
            df.drop(df.index[down], inplace = True, errors = "raise")
            df = df.reset_index(drop = True)
        else:
            up = down
            down += 1

    # delete the repeat rows from bottom to top
    up, down = len(df.index)-2, len(df.index)-1
    while up >= 0:
        if df.iloc[down, 1] <= df.iloc[up, 1] and df.iloc[up, 1] < df.iloc[down, 2] and df.iloc[down, 1] < df.iloc[up, 2] and df.iloc[up, 2] < df.iloc[down, 2]:
            df.drop(df.index[up], inplace = True, errors = "raise")
            df = df.reset_index(drop = True)
        else:
            down = up
            up -= 1

    # create agp list for creating the Ns gap for further re-sequencing
    start, end = "1", 1
    contig = 1
    for nrow in range(1, len(df.index)):
        if df.iloc[nrow, 1] <= df.iloc[nrow-1, 2] and df.iloc[nrow, 2] > df.iloc[nrow-1, 2]:
            end = str(df.iloc[nrow-1, 2])
            if contig < 10:
                print(prefix + ".0" + str(contig) + "\t1\t" + length + "\t" + str(contig) + "\tW\t" + name + "\t" + start + "\t" + end + "\t+")
            else:
                print(prefix + "." + str(contig) + "\t1\t" + length + "\t" + str(contig) + "\tW\t" + name + "\t" + start + "\t" + end + "\t+")  
            start = str(df.iloc[nrow, 1])
            contig += 1
        if df.iloc[nrow, 2] == int(length):
            end = length
            if contig < 10:
                print(prefix + ".0" + str(contig) + "\t1\t" + length + "\t" + str(contig) + "\tW\t" + name + "\t" + start + "\t" + end + "\t+")
            else:
                print(prefix + "." + str(contig) + "\t1\t" + length + "\t" + str(contig) + "\tW\t" + name + "\t" + start + "\t" + end + "\t+")

if __name__ == '__main__':
    main()
