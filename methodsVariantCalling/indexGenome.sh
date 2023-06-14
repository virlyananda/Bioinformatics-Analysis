#!/usr/bin/env bash
# indexGenome.sh
bwa index -a bwtsw GRCh38_reference.fa \
    1>indexGenome.log 2>indexGenome.err &
