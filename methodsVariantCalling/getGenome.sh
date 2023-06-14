#!/usr/bin/env bash
# getGenome.sh

# Get the GRCh38 reference genome
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz \
    -O GRCh38_reference.fa.gz \
    1>getGenome.log 2>getGenome.err &

