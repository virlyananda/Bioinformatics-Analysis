#!/usr/bin/env Rscript

# Plot Data
m <- as.matrix(read.table("ibd_view.mibs"))
mds <- cmdscale(as.dist(1-m))
k <- c(rep("green", 45), rep("blue", 44))
plot(mds, pch=20, col=k)
