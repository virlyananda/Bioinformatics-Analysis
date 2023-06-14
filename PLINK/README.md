#### Author: Virly Y. Ananda
#### Affiliation: College of Science, Northeastern University

## Overview

In this project, we use PLINK to analyze WGAS and Population-Based
Linkage(Purcell et al. 2007). In this cluster plot, we analyze 2
different population from the green color representing Chinese
individuals and blue color representing Japanese individuals from
ibd\_view file generated.

## Results and Analysis

``` r
library(knitr)
```

``` r
lmiss_stat<-read.table("miss_stat.lmiss", header=T)
kable(head(lmiss_stat))
```

| CHR | SNP       | N\_MISS | N\_GENO | F\_MISS |
| --: | :-------- | ------: | ------: | ------: |
|   1 | rs6681049 |       0 |      89 |       0 |
|   1 | rs4074137 |       0 |      89 |       0 |
|   1 | rs7540009 |       0 |      89 |       0 |
|   1 | rs1891905 |       0 |      89 |       0 |
|   1 | rs9729550 |       0 |      89 |       0 |
|   1 | rs3813196 |       0 |      89 |       0 |

``` r
## Generate summary statistics of missing rates
## This table shows for every SNP, the number of missing individuals represented by (N_MISS) and the proportion of individuals missing (F_MISS) are shown.
```

``` r
imiss_stat<-read.table("miss_stat.imiss", header=T)
kable(head(imiss_stat))
```

| FID    | IID | MISS\_PHENO | N\_MISS | N\_GENO |  F\_MISS |
| :----- | --: | :---------- | ------: | ------: | -------: |
| HCB181 |   1 | N           |     671 |   83534 | 0.008033 |
| HCB182 |   1 | N           |    1156 |   83534 | 0.013840 |
| HCB183 |   1 | N           |     498 |   83534 | 0.005962 |
| HCB184 |   1 | N           |     412 |   83534 | 0.004932 |
| HCB185 |   1 | N           |     329 |   83534 | 0.003939 |
| HCB186 |   1 | N           |    1233 |   83534 | 0.014760 |

``` r
## The table above shows the actual genotyping rate for individual mentioned along with its high rate of genotyping rate.
```

``` r
freq_strats<-read.table("freq_stat.frq.strat", header=T, sep="\t")
kable(head(freq_strats))
```

| CHR………SNP…..CLST…A1…A2……MAF….MAC..NCHROBS |
| :---------------------------------------- |
| 1 rs6681049 1 1 2 0.2333 21 90            |
| 1 rs6681049 2 1 2 0.1932 17 88            |
| 1 rs4074137 1 1 2 0.1 9 90                |
| 1 rs4074137 2 1 2 0.05682 5 88            |
| 1 rs7540009 1 0 2 0 0 90                  |
| 1 rs7540009 2 0 2 0 0 88                  |

``` r
## Alleles Frequencies Summary Statistics
```

``` r
baa_snp<-read.table("as1.assoc", header=T, sep="\t")
kable(head(baa_snp))
```

| CHR………SNP………BP…A1……F\_A……F\_U…A2……..CHISQ…………P………..OR   |
| :------------------------------------------------------ |
| 1 rs6681049 1 1 0.1591 0.2667 2 3.067 0.07991 0.5203    |
| 1 rs4074137 2 1 0.07955 0.07778 2 0.001919 0.9651 1.025 |
| 1 rs7540009 3 0 0 0 2 NA NA NA                          |
| 1 rs1891905 4 1 0.4091 0.4 2 0.01527 0.9017 1.038       |
| 1 rs9729550 5 1 0.1705 0.08889 2 2.631 0.1048 2.106     |
| 1 rs3813196 6 1 0.03409 0.02222 2 0.2296 0.6318 1.553   |

``` r
## Basic Association Analysis: Disease trait for all single SNPs
```

``` r
adj_baa<-read.table("as2.assoc.adjusted", header=T, sep="\t")
kable(head(adj_baa))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY |
| :-------------------------------------------------------------------------- |
| 13 rs9585021 5.586e-06 4.994e-05 0.3839 0.3839 0.3188 0.3188 0.09719 1      |
| 2 rs2222162 5.918e-06 5.232e-05 0.4068 0.4067 0.3342 0.3342 0.09719 1       |
| 9 rs10810856 7.723e-06 6.483e-05 0.5308 0.5308 0.4118 0.4118 0.09719 1      |
| 2 rs4675607 8.05e-06 6.703e-05 0.5533 0.5533 0.4249 0.4249 0.09719 1        |
| 2 rs4673349 8.485e-06 6.994e-05 0.5832 0.5831 0.4419 0.4419 0.09719 1       |
| 2 rs1375352 8.485e-06 6.994e-05 0.5832 0.5831 0.4419 0.4419 0.09719 1       |

``` r
## Adjusted Basic Association Analysis
## The table above shows the adjusted version which contains the log file records the inflation factor calculated for the genomic control analysis.
## The mean chi-squared statistic (that should be 1 under the null)
```

``` r
## Genotypic and Other Association Model
## The table below shows the tests are being performed by Cochran-Armitage trend test.
## This shows different tests performed for each SNP.

gen1_model<-read.table("mod1.model", header=T, sep="\t")
kable(head(gen1_model))
```

| CHR………SNP…A1…A2…..TEST…………AFF……….UNAFF……..CHISQ…DF…………P |
| :------------------------------------------------------ |
| 2 rs2222162 1 2 GENO 3/19/22 17/22/6 NA NA NA           |
| 2 rs2222162 1 2 TREND 25/63 56/34 19.15 1 1.207e-05     |
| 2 rs2222162 1 2 ALLELIC 25/63 56/34 20.51 1 5.918e-06   |
| 2 rs2222162 1 2 DOM 22/22 39/6 NA NA NA                 |
| 2 rs2222162 1 2 REC 3/41 17/28 NA NA NA                 |

``` r
## Model test that produces values for genotypic test. Here --cell is being implemented as the default option.
gen2_model<-read.table("mod2.model", header=T, sep="\t")
kable(head(gen2_model))
```

| CHR………SNP…A1…A2…..TEST…………AFF……….UNAFF……..CHISQ…DF…………P |
| :------------------------------------------------------ |
| 2 rs2222162 1 2 GENO 3/19/22 17/22/6 19.15 2 6.932e-05  |
| 2 rs2222162 1 2 TREND 25/63 56/34 19.15 1 1.207e-05     |
| 2 rs2222162 1 2 ALLELIC 25/63 56/34 20.51 1 5.918e-06   |
| 2 rs2222162 1 2 DOM 22/22 39/6 13.87 1 0.0001958        |
| 2 rs2222162 1 2 REC 3/41 17/28 12.24 1 0.0004679        |

``` r
## Stratifiation Analysis by using whole genome data to cluster individuals into homogenous groups
clus_strat<-read.table("str1.cluster1", header=T, sep="\t")
kable(head(clus_strat))
```

| SOL.0 | HCB181\_1.JPT260\_1 |
| :---- | :------------------ |
| SOL-1 | HCB182\_1 HCB225\_1 |
| SOL-2 | HCB183\_1 HCB194\_1 |
| SOL-3 | HCB184\_1 HCB202\_1 |
| SOL-4 | HCB185\_1 HCB217\_1 |
| SOL-5 | HCB186\_1 HCB201\_1 |
| SOL-6 | HCB187\_1 HCB189\_1 |

``` r
## Association Analysis: Accounting for Clusters
## The table below shows association test conditional on the matching by using Cochran-Mantel-Haenszel (CMH) association statistic
## This test can help for SNP disease association conditional on the cluster

adj_clus<-read.table("aac1.cmh.adjusted", header=T, sep="\t")
kable(head(adj_clus))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY |
| :-------------------------------------------------------------------------- |
| 13 rs9585021 1.906e-06 4.418e-06 0.1274 0.1274 0.1196 0.1196 0.1274 1       |
| 21 rs3017432 2.209e-05 4.332e-05 1 1 0.7716 0.7716 0.7384 1                 |
| 2 rs2222162 4.468e-05 8.353e-05 1 1 0.9496 0.9495 0.8734 1                  |
| 17 rs3829612 7.177e-05 0.0001299 1 1 0.9918 0.9918 0.8734 1                 |
| 2 rs4673349 9.617e-05 0.0001707 1 1 0.9984 0.9984 0.8734 1                  |
| 2 rs1375352 9.617e-05 0.0001707 1 1 0.9984 0.9984 0.8734 1                  |

``` r
## Perform clustering with different constraints
## Here, we do not put the maximum size of the cluster. Instead, we want to make sure every cluster contains 1 case
## and 1 control

gen_cluster<-read.table("version2.cluster1", header=T, sep="\t")
kable(head(gen_cluster))
```

| SOL.0 | HCB181\_1.1..HCB182\_1.1..HCB186\_1.1..HCB201\_1.2..HCB217\_1.1..HCB205\_1.1..HCB194\_1.1..HCB199\_1.1..HCB221\_1.2..HCB225\_1.1..HCB197\_1.1..HCB214\_1.2..HCB212\_1.1..JPT227\_1.2..JPT229\_1.2..JPT243\_1.1..JPT266\_1.2..JPT240\_1.2..JPT246\_1.2..JPT230\_1.2..JPT262\_1.1..JPT267\_1.2..JPT245\_1.2..JPT254\_1.1..JPT264\_1.2. |
| :---- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| SOL-1 | HCB183\_1(2) HCB206\_1(1) HCB208\_1(1) HCB189\_1(1) HCB191\_1(2) HCB220\_1(1) HCB196\_1(1) HCB207\_1(2) HCB223\_1(1) JPT253\_1(1) HCB203\_1(1) HCB193\_1(1) HCB195\_1(1) HCB215\_1(2) HCB216\_1(1)                                                                                                                                   |
| SOL-2 | HCB184\_1(1) HCB219\_1(2) HCB218\_1(1) HCB200\_1(1) HCB190\_1(1) HCB202\_1(1) HCB224\_1(2) HCB185\_1(1) HCB198\_1(1) HCB210\_1(2) HCB222\_1(1) HCB188\_1(1) HCB192\_1(1) HCB187\_1(1) HCB209\_1(1) HCB211\_1(2) HCB213\_1(1)                                                                                                         |
| SOL-3 | HCB204\_1(1) JPT255\_1(2) JPT235\_1(2) JPT237\_1(2) JPT250\_1(1) JPT251\_1(2) JPT258\_1(2) JPT228\_1(1) JPT252\_1(2) JPT242\_1(2) JPT232\_1(2) JPT249\_1(2) JPT236\_1(2) JPT256\_1(1) JPT261\_1(2) JPT247\_1(2) JPT268\_1(2) JPT259\_1(2) JPT260\_1(2) JPT257\_1(2)                                                                  |
| SOL-4 | JPT226\_1(1) JPT244\_1(2) JPT238\_1(2) JPT233\_1(1) JPT248\_1(2) JPT241\_1(2) JPT234\_1(2) JPT265\_1(1) JPT269\_1(2) JPT231\_1(1) JPT239\_1(2) JPT263\_1(2)                                                                                                                                                                          |

``` r
## Changes in the adjusted association analysis: Disease SNP is now genome-wide significant
adj_clus2<-read.table("aac2.cmh.adjusted", header=T, sep="\t")
kable(head(adj_clus2))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY          |
| :----------------------------------------------------------------------------------- |
| 2 rs2222162 1.474e-08 2.276e-08 0.001013 0.001013 0.001013 0.001013 0.001013 0.01187 |
| 2 rs4675607 1.134e-05 1.479e-05 0.7794 0.7794 0.5413 0.5413 0.2779 1                 |
| 13 rs9585021 1.213e-05 1.58e-05 0.8338 0.8338 0.5656 0.5656 0.2779 1                 |
| 9 rs7046471 3.372e-05 4.278e-05 1 1 0.9015 0.9015 0.4296 1                           |
| 2 rs4673349 3.75e-05 4.746e-05 1 1 0.924 0.924 0.4296 1                              |
| 2 rs1375352 3.75e-05 4.746e-05 1 1 0.924 0.924 0.4296 1                              |

``` r
## Stratification Analysis: Number Specification on Cluster
## The table below shows a 2 class solution
two_c_sa<-read.table("version3.cluster1", header=T, sep="\t")
kable(head(two_c_sa))
```

| SOL.0 | HCB181\_1.HCB183\_1.HCB195\_1.HCB193\_1.HCB187\_1.HCB209\_1.HCB211\_1.HCB205\_1.HCB208\_1.HCB189\_1.HCB191\_1.HCB220\_1.HCB197\_1.HCB214\_1.HCB212\_1.HCB199\_1.HCB221\_1.HCB207\_1.HCB223\_1.HCB182\_1.HCB225\_1.HCB190\_1.HCB185\_1.HCB217\_1.HCB198\_1.HCB210\_1.HCB222\_1.HCB192\_1.HCB184\_1.HCB219\_1.HCB218\_1.HCB224\_1                                                                                                                                                                                                                                                           |
| :---- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| SOL-1 | HCB186\_1 HCB201\_1 HCB213\_1 HCB188\_1 HCB206\_1 HCB194\_1 HCB200\_1 HCB215\_1 HCB216\_1 HCB196\_1 JPT253\_1 HCB202\_1 HCB203\_1 HCB204\_1 JPT255\_1 JPT226\_1 JPT244\_1 JPT238\_1 JPT269\_1 JPT233\_1 JPT248\_1 JPT232\_1 JPT247\_1 JPT235\_1 JPT237\_1 JPT250\_1 JPT246\_1 JPT249\_1 JPT258\_1 JPT227\_1 JPT266\_1 JPT268\_1 JPT229\_1 JPT243\_1 JPT251\_1 JPT259\_1 JPT267\_1 JPT228\_1 JPT252\_1 JPT242\_1 JPT234\_1 JPT245\_1 JPT254\_1 JPT264\_1 JPT230\_1 JPT262\_1 JPT236\_1 JPT256\_1 JPT265\_1 JPT240\_1 JPT241\_1 JPT261\_1 JPT231\_1 JPT239\_1 JPT263\_1 JPT260\_1 JPT257\_1 |

``` r
## Obtaining a genome-wide significance
adj_clus3<-read.table("aac3.cmh.adjusted", header=T, sep="\t")
kable(head(adj_clus3))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY                 |
| :------------------------------------------------------------------------------------------ |
| 2 rs2222162 2.594e-10 2.594e-10 1.783e-05 1.783e-05 1.783e-05 1.783e-05 1.783e-05 0.0002089 |
| 2 rs4675607 4.03e-06 4.03e-06 0.277 0.277 0.2419 0.2419 0.1385 1                            |
| 2 rs4673349 1.204e-05 1.204e-05 0.8276 0.8276 0.5629 0.5629 0.1679 1                        |
| 2 rs1375352 1.204e-05 1.204e-05 0.8276 0.8276 0.5629 0.5629 0.1679 1                        |
| 13 rs9585021 1.222e-05 1.222e-05 0.8395 0.8395 0.5681 0.5681 0.1679 1                       |
| 2 rs2176427 4.309e-05 4.309e-05 1 1 0.9483 0.9483 0.4936 1                                  |

``` r
#!/usr/bin/env Rscript

# Plot Data
m <- as.matrix(read.table("ibd_view.mibs"))
mds <- cmdscale(as.dist(1-m))
k <- c(rep("green", 45), rep("blue", 44))
plot(mds, pch=20, col=k)
```

![](runPlink_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
## Quantitative Trait Association Analysis: Analyzing quantitative trait directly

quan_trait<-read.table("quant1.qassoc", header=T, sep="\t")
kable(head(quan_trait))
```

| CHR………SNP………BP….NMISS…….BETA………SE………R2……..T…………P        |
| :------------------------------------------------------ |
| 1 rs6681049 1 89 -0.2266 0.3626 0.004469 -0.6249 0.5336 |
| 1 rs4074137 2 89 -0.2949 0.6005 0.002765 -0.4911 0.6246 |
| 1 rs7540009 3 89 NA NA NA NA NA                         |
| 1 rs1891905 4 89 -0.1053 0.3165 0.001272 -0.3328 0.7401 |
| 1 rs9729550 5 89 0.5402 0.4616 0.0155 1.17 0.2451       |
| 1 rs3813196 6 89 0.8053 1.025 0.00705 0.7859 0.434      |

``` r
## Adjusted Quantitative Trait Association Analysis
adj_quan_trait<-read.table("quant1.qassoc.adjusted", header=T, sep="\t")
kable(head(adj_quan_trait))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY                |
| :----------------------------------------------------------------------------------------- |
| 2 rs2222162 5.273e-09 6.224e-08 0.0003624 0.0003624 0.0003623 0.0003623 0.0003624 0.004245 |
| 21 rs219746 1.095e-06 6.815e-06 0.07529 0.07529 0.07252 0.07252 0.03764 0.441              |
| 7 rs1922519 1.63e-05 7.167e-05 1 1 0.6739 0.6739 0.2527 1                                  |
| 2 rs2969348 2.886e-05 0.0001176 1 1 0.8624 0.8624 0.2527 1                                 |
| 3 rs6773558 3.581e-05 0.0001419 1 1 0.9146 0.9146 0.2527 1                                 |
| 10 rs3862003 3.716e-05 0.0001465 1 1 0.9222 0.9222 0.2527 1                                |

``` r
## PERMUTATION TESTING: A method that controls for any between-cluster association under permuted datasets
## The table shows how sorting adaptive permutation is being implemented to rank disease variant based on empirical significance value
adapt_perm<-read.table("quant2.qassoc.perm", header=T, sep="\t")
kable(head(adapt_perm))
```

| CHR………SNP………EMP1………..NP |
| :---------------------- |
| 1 rs6681049 0.5926 26   |
| 1 rs4074137 0.6154 25   |
| 1 rs7540009 NA NA       |
| 1 rs1891905 1 6         |
| 1 rs9729550 0.04461 806 |
| 1 rs3813196 0.3425 72   |

``` r
## The table below shows how continuous phenotype differentiate two populations
diff_pop<-read.table("quant3.qassoc.gxe", header=T, sep="\t")
kable(head(diff_pop))
```

| CHR………SNP…NMISS1……BETA1……..SE1…NMISS2……BETA2……..SE2….Z\_GXE……..P\_GXE |
| :-------------------------------------------------------------------- |
| 2 rs2222162 45 -2.271 0.2245 44 -1.997 0.1722 -0.9677 0.3332          |

``` r
d<-read.table("rec_snp1.raw", header=T)
summary(glm(PHENOTYPE-1 ~ rs2222162_1, data=d, family="binomial"))
```

    ## 
    ## Call:
    ## glm(formula = PHENOTYPE - 1 ~ rs2222162_1, family = "binomial", 
    ##     data = d)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.7690  -1.1042  -0.5848   0.6851   1.9238  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.3300     0.4107   3.238   0.0012 ** 
    ## rs2222162_1  -1.5047     0.3765  -3.997 6.42e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 123.37  on 88  degrees of freedom
    ## Residual deviance: 102.64  on 87  degrees of freedom
    ## AIC: 106.64
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
## Extracting an SNP of interst
## Here, we should be ablet to see the summary of the coefficient from rs2222162_1
```

## References

<div id="refs" class="references">

<div id="ref-PLINK">

Purcell, Shaun, Benjamin Neale, Kathe Todd-Brown, Lori Thomas, Manuel A.
R. Ferreira, David Bender, Julian Maller, et al. 2007. “PLINK: A Tool
Set for Whole-Genome Association and Population-Based Linkage Analyses.”
*American Journal of Human Genetics* 81 (3): 559–75.
<https://doi.org/10.1086/519795>.

</div>

</div>
