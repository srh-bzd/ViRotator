#!/usr/bin/env python3


import os, argparse

"""
USAGE
    ./rotate.py [-h] -f <input_file> -t <file_type> -b <blast_file> -o <output_file> -l <output_log>
DESCRIPTION
    Script to rotate sequences.
PREREQUISITE
    - python3
    - A Blast file formated like 'qseqid qstart'
WARNING
    ! This script can be used only with tripled sequences !
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


def parse_blast_file(blast_file):
    """
    Parse input blast file and keep it into a dict like : {seqId : {"pos" : [] }}
    """
    blastDict = dict()
    for line in blast_file:
        line = line.strip().split("\t")
        seqId = line[0]
        pos = int(line[1])
        if seqId not in blastDict:
            blastDict[seqId] = dict()
            blastDict[seqId]["pos"] = [pos]
        else:
            blastDict[seqId]["pos"].append(pos)
    return blastDict


def rotate(fileDict, blastDict, file_type):
    """
    Rotate sequences according on the 2nd and 3rd starting position of the discovered gene in the ordered list.
    """
    idNotBlast = list()
    for seqId in fileDict:
        if seqId in blastDict:
            blastDict[seqId]["pos"].sort()
            start2 = blastDict[seqId]["pos"][1] - 1 # Minus 1 because Blast positions are on one base or Python in zero base
            start3 = blastDict[seqId]["pos"][2] - 1
            fileDict[seqId]["seq"] = fileDict[seqId]["seq"][start2:start3]
            if file_type == "fastq":
                fileDict[seqId]["qual"] = fileDict[seqId]["qual"][start2:start3]
        else:
            idNotBlast.append(seqId)
    for seqId in idNotBlast:
        fileDict.pop(seqId, None)
    return fileDict, idNotBlast


def write_output_files(newFileDict, idNotBlast, file_type, output_file, output_log):
    """
    Write the new FASTA or FASTQ file with all the sequences rotate 
    and the log file with all id of sequences where position of gene was not find.
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
    # Log files
    outputPath, outputName = os.path.split(output_log.name)
    outputName = os.path.splitext(outputName)[0]
    print("\n*---------- Sequences where position of the gene for rotation was not find ("+str(len(idNotBlast))+") : ", file=output_log)
    print("\n".join(idNotBlast), file=output_log)
    with open(outputPath+"/rejected.count.txt", "a") as output_count:
        print(outputName+"\t"+str(len(idNotBlast)), file=output_count)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to rotate sequences")
    parser.add_argument('--input_file', '-f', required=True,help="Input file with all reads tripled",type=argparse.FileType('r'))
    parser.add_argument('--file_type', '-t', required=True,help="Format of the input file : fasta or fastq",type=str)
    parser.add_argument('--blast_file', '-b', required=True,help="Output file of Blast formated with columns 'qseqid qstart'",type=argparse.FileType('r'))
    parser.add_argument('--output_file', '-o', required=True,help="Output file with all reads rotated",type=argparse.FileType('w'))
    parser.add_argument('--output_log', '-l', required=True,help="Output log",type=argparse.FileType('a'))
    args = parser.parse_args()
    
    if args.file_type == "fasta" :
        file_dict = parse_fasta(args.input_file)
    elif args.file_type == "fastq" : 
        file_dict = parse_fastq(args.input_file)
    blast_dict = parse_blast_file(args.blast_file)
    new_file_dict, id_seq_not_blast = rotate(file_dict, blast_dict, args.file_type)
    write_output_files(new_file_dict, id_seq_not_blast, args.file_type, args.output_file, args.output_log)
