#!/usr/bin/env python3


import argparse

"""
USAGE
    ./triple.py [-h] -f <input_file> -t <file_type> -o <output_file>
DESCRIPTION
    Script to triple sequences contain into a FASTA or FASTQ file.
PREREQUISITE
    - python3
"""


def parse_fasta(input_file):
    """
    Parse input FASTA file and keep it into a dict like : {seqId : {"seq" : }}
    """
    fastaDict = dict()
    seqId = None
    for line in input_file:
        if line.startswith(">"):
            seqId = line.strip().replace(">","").split()[0]
            fastaDict[seqId] = {}
            fastaDict[seqId]["seq"] = ""
        else:
            fastaDict[seqId]["seq"] += line.strip()
    return fastaDict


def parse_fastq(input_file):
    """
    Parse input FASTQ file and keep it into a dict like : {seqId : {"seq" : , "qual" : }}
    """
    fastqDict = dict()
    seqId = None
    for line in input_file:
        if line.startswith("@") and "~" not in line: # Rq : some quality lines start with @ this is why there is the second part of the if
            seqId = line.strip().replace("@","").split()[0]
            complementary_lines = []
            fastqDict[seqId] = {}
            fastqDict[seqId]["seq"] = ""
            fastqDict[seqId]["qual"] = ""
        else:
            complementary_lines.append(line.strip())
            if len(complementary_lines)==3:
                fastqDict[seqId]["seq"] += complementary_lines[0] 
                fastqDict[seqId]["qual"] += complementary_lines[2]
    return fastqDict


def write_output_files(fileDict, file_type, output_file):
    """
    Write the new FASTA or FASTQ file with the reads tripled
    """
    for seqId in fileDict.keys():
        if file_type=="fasta":
            print(">"+seqId, file=output_file)
            print(fileDict[seqId]["seq"]*3, file=output_file)
        elif file_type=="fastq":
            print("@"+seqId, file=output_file)
            print(fileDict[seqId]["seq"]*3, file=output_file)
            print("+", file=output_file)
            print(fileDict[seqId]["qual"]*3, file=output_file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to triple sequences")
    parser.add_argument('--input_file', '-f', required=True,help="Input file with all reads",type=argparse.FileType('r'))
    parser.add_argument('--file_type', '-t', required=True,help="Format of the input file : fasta or fastq",type=str)
    parser.add_argument('--output_file', '-o', required=True,help="Output file with all sequences tripled",type=argparse.FileType('w'))
    args = parser.parse_args()
    
    if args.file_type == "fasta" :
        file_dict = parse_fasta(args.input_file)
    elif args.file_type == "fastq" : 
        file_dict = parse_fastq(args.input_file)
    write_output_files(file_dict, args.file_type, args.output_file)
