#!/usr/bin/env bash
# alignReads.sh
bwa mem -t 8 -R "@RG\tID:SRR6808334\tSM:bar" -p GRCh38_reference.fa SRR6808334_1.fastq SRR6808334_2.fastq \
    1>SRR6808334.sam 2>alignReads.err &

