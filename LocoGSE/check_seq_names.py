#!/usr/bin/env python3
import sys

def check_sequence_names(fastq_list):
    contains_dup = False
    
    for fastq in fastq_list:
        with open(fastq) as inf:
            seq_names = set()

            line = " "
            while line != "":
                # Get everything between the "@" and the first space
                line = inf.readline().strip("\n")
                header = line[1:].split(" ")[0]
            
                if header in seq_names:
                    contains_dup = True
                    break
                
                seq_names.add(header)

                # Skip sequence, comments and quality lines
                for i in range(0, 3):
                    inf.readline()

    if contains_dup:
        print("WARNING: Diamond reads sequence names up to the first space.")
        print("WARNING: Some of your sequences may have the same name.")
        print("WARNING: Please consider renaming your sequences for optimal results.")



def run():
    if "-h" in sys.argv:
        print("Usage: check_seq_names.py fastq1 fastq2 ...")
        exit()
    
    check_sequence_names(sys.argv[1:])

