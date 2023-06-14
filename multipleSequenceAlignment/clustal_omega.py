#!/usr/bin/env python3

# Load module for Clustal Omega
from Bio.Align.Applications import ClustalOmegaCommandline

in_file = "apoe_aa.fasta"
out_file = "aligned_clustal.fasta"
clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
print(clustalomega_cline)
