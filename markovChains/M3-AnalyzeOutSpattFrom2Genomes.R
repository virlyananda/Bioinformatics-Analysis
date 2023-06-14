#!/usr/bin/env Rscript


options("scipen" = 30)
# Usage as a statement to print out input file
usage <- "\nUsage: AnalyzeOutSpattFrom2Genes.R <input file.txt>\n\n"
# Tell R to get arguments
args <- commandArgs(trailingOnly = TRUE)
# Test argument with if
if(length(args) == 0) {
    cat(prompt=usage)
    q(save="no")
}
# spattMetrics as data frame
spattMetrics <- read.table(args[1], sep="\t", header=FALSE)
# Names for each column on data frame. Put variable names.
names(spattMetrics) <- c("kmer1", "occurence1", "expected1", "pvalue1", "kmer2", "occurence2", "expected2", "pvalue2")
# Add new column to data frame by giving all 0
spattMetrics["probability_occurence"] <- 0
# Fill column
spattMetrics$probability_occurence <- (spattMetrics$expected1 / (spattMetrics$expected1 + spattMetrics$expected2))

spattMetrics["trials"] <- 0
spattMetrics$trials <- spattMetrics$occurence1 + spattMetrics$occurence2

spattMetrics["binomial_upper_tail"] <- 0
spattMetrics$binomial_upper_tail <- pbinom(spattMetrics$occurence1-1, spattMetrics$trials, spattMetrics$probability_occurence, lower.tail=FALSE)

# Write output to a file that can be opened on Excel
write.table(spattMetrics, file="out.txt", sep="\t", row.names = FALSE)
