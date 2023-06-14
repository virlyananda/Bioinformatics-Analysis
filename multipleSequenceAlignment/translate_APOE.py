#!/usr/bin/env python3

# Load modules
from Bio import SeqIO
from Bio.Seq import Seq

# Define fasta file
input_fasta = "APOE_refseq_transcript.fasta"

# Input 1 empty list to append the whole sequence(record)
seq = ''
list_seqrecords = []

# Use file handle to read the original fasta file
with open(input_fasta, "r") as apoe_nt:
    # Parse the fasta file after reading to get record. Seq record is needed for the file to be written on another fasta file
    for record in SeqIO.parse(apoe_nt, "fasta"):
        print(record.id)
        record.seq = record.seq.translate()
        print(record.seq)
        list_seqrecords.append(record)
# Use file handle to open output file and use SeqIO.write to put the translated file to the new fasta file        
with open("apoe_aa.fasta", 'w') as apoe_aa:
    SeqIO.write(list_seqrecords, apoe_aa, 'fasta')
