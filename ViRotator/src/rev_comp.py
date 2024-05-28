#!/usr/bin/env python3


import argparse

"""
USAGE
    ./rev_comp.py [-h] -f <input_file> -t <file_type> -b <blast_file> -o <output_file> -l <output_log>
DESCRIPTION
    Script to reverse complement sequences on minus strand. fasta and fastq file are accepted.
PREREQUISITE
    - python3
    - A blast output file formated like 'qseqid sstrand'
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
        if line.startswith("@") and "~" not in line: # Rq : some quality lines start with @
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


def parse_blast_file(blast_file):
    """
    Parse input blast file and keep it into a dict like : {seqId : {"strand" : }}
    """
    blastDict = dict()
    for line in blast_file:
        line = line.strip().split("\t")
        seqId = line[0]
        strand = line[1]
        blastDict[seqId] = dict()
        blastDict[seqId]["strand"] = strand
    return blastDict


def reverse_complement(fileDict, blastDict, file_type):
    """
    Reverse complement sequences into fileDict if they're on the minus strand, info into blastDict. 
    Keep all id of sequences reverse complemented into a list idSeqRevComp.
    """
    idSeqRevComp = list()
    idNotBlast = list()
    for seqId in fileDict:
        if seqId in blastDict:
            if blastDict[seqId]["strand"] == "minus":
                idSeqRevComp.append(seqId)
                # Complement
                comp = fileDict[seqId]["seq"].maketrans('ATGC', 'TACG')
                # Reverse
                fileDict[seqId]["seq"] = fileDict[seqId]["seq"].translate(comp)[::-1]
                if file_type == "fastq":
                    fileDict[seqId]["qual"] = fileDict[seqId]["qual"][::-1]
        else:
            idNotBlast.append(seqId)
    for seqId in idNotBlast:
        fileDict.pop(seqId, None)
    return fileDict, idSeqRevComp, idNotBlast


def write_output_files(newFileDict, idSeqRevComp, idNotBlast, file_type, output_file, output_log):
    """
    Write the new FASTA or FASTQ file with all the reads in the plus strand 
    and the log file with all id of reads reverse complemented and where strand was not find by blast.
    """
    for seqId in newFileDict.keys():
        if file_type=="fasta":
            print(">"+seqId, file=output_file)
            print(newFileDict[seqId]["seq"], file=output_file)
        elif file_type=="fastq":
            print("@"+seqId, file=output_file)
            print(newFileDict[seqId]["seq"], file=output_file)
            print("+", file=output_file)
            print(newFileDict[seqId]["qual"], file=output_file)
    # Log file
    print("*---------- Sequences reverse complemented ("+str(len(idSeqRevComp))+") : ", file=output_log)
    print("\n".join(idSeqRevComp), file=output_log)
    print("\n*---------- Sequences where strand was not find ("+str(len(idNotBlast))+") : ", file=output_log)
    print("\n".join(idNotBlast), file=output_log)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to reverse complement sequences which are on the minus strand")
    parser.add_argument('--input_file', '-f', required=True,help="Input file with all reads",type=argparse.FileType('r'))
    parser.add_argument('--file_type', '-t', required=True,help="Format of the input file : fasta or fastq",type=str)
    parser.add_argument('--blast_file', '-b', required=True,help="Output file of blast formated with columns 'qseqid sstrand'",type=argparse.FileType('r'))
    parser.add_argument('--output_file', '-o', required=True,help="Output file with all sequences on plus strand",type=argparse.FileType('w'))
    parser.add_argument('--output_log', '-l', required=True,help="Output log",type=argparse.FileType('w'))
    args = parser.parse_args()
    
    if args.file_type == "fasta" :
        file_dict = parse_fasta(args.input_file)
    elif args.file_type == "fastq" : 
        file_dict = parse_fastq(args.input_file)
    blast_dict = parse_blast_file(args.blast_file)
    new_file_dict, id_seq_rev_comp, id_seq_not_blast = reverse_complement(file_dict, blast_dict, args.file_type)
    write_output_files(new_file_dict, id_seq_rev_comp, id_seq_not_blast, args.file_type, args.output_file, args.output_log)
   