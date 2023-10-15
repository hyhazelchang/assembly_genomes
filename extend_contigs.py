#!/usr/bin/python3

# extend_contigs.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# notes: /home/xinchang/notes/python/pytest06_add_softclip_into_gaps.txt
# v1 2023/03/17

# Should install SE-MEI/extractSoftclipped first!
# Download to local
# Usage: python3 /home/xinchang/pyscript/pyscript_xin/extend_contigs.py --softclip_file="/scratch/xinchang/cyano22/cyano22.01/bwa_01/flashcut.fastq" --length_file="/scratch/xinchang/cyano22/cyano22.01/bwa_01/cyano22.01.fasta.length" --output="/scratch/xinchang/cyano22/cyano22.01/fill.01.gap"

import argparse

def main():
    parser = argparse.ArgumentParser(
        description=("Extend the velvet contigs by softclips"),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--softclip_file",
                        default=None,
                        type=str,
                        help="Softclip file. Please provide absolute path.")
    parser.add_argument("--length_file",
                        default=None,
                        type=str,
                        help="Length file. Please provide absolute path.")
    parser.add_argument("--output",
                        default="./fill.gap",
                        type=str,
                        help="Output files. Please provide absolute path.")

    args = parser.parse_args()

    seqfile = args.softclip_file
    lenfile = args.length_file
    fill = args.output
    nt = "ATCG"

    # parse out the needed lines
    merge = get_softclip(seqfile, nt)
    lens = get_length(lenfile)

    # create and write into a fill gaps file
    fillgap = open(fill, "w")
    for i in merge:
        for x in range(len(lens)):
            if i[0] == lens[x][0] and int(i[1]) == 0:
                fillgap.write(lens[x][0] + "\t" + lens[x][1] + "\t1\t" + "1\t1" + "\t" + str(len(i[2])) +"\tY\tM\t" + i[2] + "A\n")
            if i[0] == lens[x][0] and int(i[1]) == int(lens[x][1]):
                fillgap.write(lens[x][0] + "\t" + lens[x][1] + "\t1\t" + str(int(lens[x][1])) + "\t" + str(int(lens[x][1])) + "\t" + str(len(i[2])) +"\tY\tM\tA" + i[2] + "\n")
    fillgap.close()

def get_softclip(seqfile, nt):
    names, regions, clips = [], [], []
    for line in open(seqfile):
        if line[0] == "@":
            name = line.split(":")
            names.append(name[0])
            regions.append(name[1])
        if line[0] in nt:
            clips.append(line)

    # strip out useless characters
    names = [line.lstrip("@") for line in names]
    regions = [line.rstrip("\n") for line in regions]
    clips = [line.rstrip("\n") for line in clips]

    # merge regions and clips
    merge = [(names[i], regions[i], clips[i]) for i in range(len(clips))]
    merge.sort()

    # remove duplicates
    l = 0
    r = 1
    while r < len(merge):
        if merge[l][1] == merge[r][1]:
            if len(merge[l][2]) >= len(merge[r][2]):
                merge.pop(r)
            else:
                merge.pop(l)
        else:
            l = r
            r += 1
    return merge

def get_length(lenfile):
    names, length = [], []
    for line in open(lenfile):
        row = line.split("\t")
        names.append(row[0])
        length.append(row[1])

    length = [line.rstrip("\n") for line in length]

    # merge regions and clips
    lens = [(names[i], length[i]) for i in range(len(length))]
    return lens

if __name__ == '__main__':
    main()