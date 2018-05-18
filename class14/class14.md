---
title: "class14"
author: "Amay Dankar"
date: "May 18, 2018"
output: 
  html_document: 
    keep_md: yes
---



loading packages


## Transcriptomics and analysis of RNA-seq data


```r
counts <- read.csv("data/airway_scaledcounts.csv", stringsAsFactors = FALSE)
metadata <-  read.csv("data/airway_metadata.csv", stringsAsFactors = FALSE)
```


```r
#View(counts)
#View(metadata)
```

Verifying metadata corresponds to correct column


```r
# excluding first column name "ensgene" which does not correspond to any metadata id
colnames(counts)[-1] == metadata$id
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```



## Exploratory differential gene expression analysis (done in an inadvisable way for demonstration purposes)


### Obtaining means of control samples


```r
control <- metadata[metadata[,"dex"]=="control",]
control.mean <- rowSums( counts[ ,control$id] )/nrow(control)
names(control.mean) <- counts$ensgene

head(control.mean)
```

```
## ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 
##          900.75            0.00          520.50          339.75 
## ENSG00000000460 ENSG00000000938 
##           97.25            0.75
```

### Obtaining means of treated samples


```r
treated <- metadata[metadata[,"dex"] == "treated",]
treated.mean <- rowSums(counts[ ,treated$id])/nrow(treated)
names(treated.mean) <- counts$ensgene

head(treated.mean)
```

```
## ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 
##          658.00            0.00          546.00          316.50 
## ENSG00000000460 ENSG00000000938 
##           78.75            0.00
```


Storing mean counts of both for convenience and summarizing number of samples in each (in an inadvisable fashion)


```r
meancounts <- data.frame(control.mean, treated.mean)
head(meancounts)
```

```
##                 control.mean treated.mean
## ENSG00000000003       900.75       658.00
## ENSG00000000005         0.00         0.00
## ENSG00000000419       520.50       546.00
## ENSG00000000457       339.75       316.50
## ENSG00000000460        97.25        78.75
## ENSG00000000938         0.75         0.00
```

```r
colSums(meancounts)
```

```
## control.mean treated.mean 
##     23005324     22196524
```

Plot mean counts for control and treated


```r
plot(control.mean, treated.mean, xlab = 'control', ylab = 'treated')
```

![](class14_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

Plotting with a log axis to better illustrate differences


```r
plot(control.mean, treated.mean, xlab = 'control', ylab = 'treated', log = 'xy')
```

```
## Warning in xy.coords(x, y, xlabel, ylabel, log): 15032 x values <= 0
## omitted from logarithmic plot
```

```
## Warning in xy.coords(x, y, xlabel, ylabel, log): 15281 y values <= 0
## omitted from logarithmic plot
```

![](class14_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Calculating log2foldchange to better highlight differences between data, and add it to our meancounts dataframe


```r
meancounts$log2fc <- log2(meancounts[,'treated.mean']/meancounts[,'control.mean'])
head(meancounts)
```

```
##                 control.mean treated.mean      log2fc
## ENSG00000000003       900.75       658.00 -0.45303916
## ENSG00000000005         0.00         0.00         NaN
## ENSG00000000419       520.50       546.00  0.06900279
## ENSG00000000457       339.75       316.50 -0.10226805
## ENSG00000000460        97.25        78.75 -0.30441833
## ENSG00000000938         0.75         0.00        -Inf
```

Small test of what "arr.ind" paramter does within **which()** function


```r
x <- matrix(1:10, ncol = 2, byrow = TRUE)
x
```

```
##      [,1] [,2]
## [1,]    1    2
## [2,]    3    4
## [3,]    5    6
## [4,]    7    8
## [5,]    9   10
```

```r
x[5,2] <- 0
x
```

```
##      [,1] [,2]
## [1,]    1    2
## [2,]    3    4
## [3,]    5    6
## [4,]    7    8
## [5,]    9    0
```

```r
which (x == 0)
```

```
## [1] 10
```

```r
which (x == 0, arr.ind = TRUE)
```

```
##      row col
## [1,]   5   2
```

Finding genes with 0 counts


```r
# labeling which values are 0
zero.vals <- which(meancounts[,1:2]==0, arr.ind=TRUE)
head(zero.vals)
```

```
##                 row col
## ENSG00000000005   2   1
## ENSG00000004848  65   1
## ENSG00000004948  70   1
## ENSG00000005001  73   1
## ENSG00000006059 121   1
## ENSG00000006071 123   1
```

```r
#assigning vector to remove unique 0 values
to.rm <- unique(zero.vals[,1])

#removed values equal to 0 to generate new dataframe

mycounts <- meancounts[-to.rm,]
head(mycounts)
```

```
##                 control.mean treated.mean      log2fc
## ENSG00000000003       900.75       658.00 -0.45303916
## ENSG00000000419       520.50       546.00  0.06900279
## ENSG00000000457       339.75       316.50 -0.10226805
## ENSG00000000460        97.25        78.75 -0.30441833
## ENSG00000000971      5219.00      6687.50  0.35769358
## ENSG00000001036      2327.00      1785.75 -0.38194109
```

Summing up non-zero values based on upregulation or downregulation


```r
up.ind <- mycounts$log2fc > 2
down.ind <- mycounts$log2fc < (-2)

paste("Up:", sum(up.ind))
```

```
## [1] "Up: 250"
```

```r
paste("Down:", sum(down.ind))
```

```
## [1] "Down: 367"
```

# Adding annotation data


```r
anno <- read.csv("data/annotables_grch38.csv")
head(anno)
```

```
##           ensgene entrez   symbol chr     start       end strand
## 1 ENSG00000000003   7105   TSPAN6   X 100627109 100639991     -1
## 2 ENSG00000000005  64102     TNMD   X 100584802 100599885      1
## 3 ENSG00000000419   8813     DPM1  20  50934867  50958555     -1
## 4 ENSG00000000457  57147    SCYL3   1 169849631 169894267     -1
## 5 ENSG00000000460  55732 C1orf112   1 169662007 169854080      1
## 6 ENSG00000000938   2268      FGR   1  27612064  27635277     -1
##          biotype
## 1 protein_coding
## 2 protein_coding
## 3 protein_coding
## 4 protein_coding
## 5 protein_coding
## 6 protein_coding
##                                                                                                  description
## 1                                                          tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
## 2                                                            tenomodulin [Source:HGNC Symbol;Acc:HGNC:17757]
## 3 dolichyl-phosphate mannosyltransferase polypeptide 1, catalytic subunit [Source:HGNC Symbol;Acc:HGNC:3005]
## 4                                               SCY1-like, kinase-like 3 [Source:HGNC Symbol;Acc:HGNC:19285]
## 5                                    chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
## 6                          FGR proto-oncogene, Src family tyrosine kinase [Source:HGNC Symbol;Acc:HGNC:3697]
```

Using **merge()** function to merge mycounts dataframe with annotation dataframe


```r
results <- merge(mycounts, anno, by.x = "row.names", by.y = "ensgene")
head(results)
```

```
##         Row.names control.mean treated.mean      log2fc entrez   symbol
## 1 ENSG00000000003       900.75       658.00 -0.45303916   7105   TSPAN6
## 2 ENSG00000000419       520.50       546.00  0.06900279   8813     DPM1
## 3 ENSG00000000457       339.75       316.50 -0.10226805  57147    SCYL3
## 4 ENSG00000000460        97.25        78.75 -0.30441833  55732 C1orf112
## 5 ENSG00000000971      5219.00      6687.50  0.35769358   3075      CFH
## 6 ENSG00000001036      2327.00      1785.75 -0.38194109   2519    FUCA2
##   chr     start       end strand        biotype
## 1   X 100627109 100639991     -1 protein_coding
## 2  20  50934867  50958555     -1 protein_coding
## 3   1 169849631 169894267     -1 protein_coding
## 4   1 169662007 169854080      1 protein_coding
## 5   1 196651878 196747504      1 protein_coding
## 6   6 143494811 143511690     -1 protein_coding
##                                                                                                  description
## 1                                                          tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
## 2 dolichyl-phosphate mannosyltransferase polypeptide 1, catalytic subunit [Source:HGNC Symbol;Acc:HGNC:3005]
## 3                                               SCY1-like, kinase-like 3 [Source:HGNC Symbol;Acc:HGNC:19285]
## 4                                    chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
## 5                                                     complement factor H [Source:HGNC Symbol;Acc:HGNC:4883]
## 6                                          fucosidase, alpha-L- 2, plasma [Source:HGNC Symbol;Acc:HGNC:4008]
```

Setting up DESeq and necessary database


```r
library(DESeq2)
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following object is masked from 'package:base':
## 
##     apply
```

```r
library(org.Hs.eg.db)
```

```
## Loading required package: AnnotationDbi
```

```
## 
```



```r
mycounts$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(mycounts),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```


# DESeq2 Analysis 



```r
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)
```

```
## converting counts to integer mode
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
## design formula are characters, converting to factors
```

```r
dds
```

```
## class: DESeqDataSet 
## dim: 38694 8 
## metadata(1): version
## assays(1): counts
## rownames(38694): ENSG00000000003 ENSG00000000005 ...
##   ENSG00000283120 ENSG00000283123
## rowData names(0):
## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
## colData names(4): id dex celltype geo_id
```


```r
dds <- DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

Getting results out of DESeq transformed dds 


```r
res <- results(dds)
res
```

```
## log2 fold change (MLE): dex treated vs control 
## Wald test p-value: dex treated vs control 
## DataFrame with 38694 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE       stat     pvalue
##                 <numeric>      <numeric> <numeric>  <numeric>  <numeric>
## ENSG00000000003 747.19420    -0.35070283 0.1682342 -2.0846111 0.03710462
## ENSG00000000005   0.00000             NA        NA         NA         NA
## ENSG00000000419 520.13416     0.20610652 0.1010134  2.0403876 0.04131173
## ENSG00000000457 322.66484     0.02452714 0.1451103  0.1690242 0.86577762
## ENSG00000000460  87.68263    -0.14714409 0.2569657 -0.5726216 0.56690095
## ...                   ...            ...       ...        ...        ...
## ENSG00000283115  0.000000             NA        NA         NA         NA
## ENSG00000283116  0.000000             NA        NA         NA         NA
## ENSG00000283119  0.000000             NA        NA         NA         NA
## ENSG00000283120  0.974916     -0.6682308  1.694063 -0.3944544  0.6932456
## ENSG00000283123  0.000000             NA        NA         NA         NA
##                      padj
##                 <numeric>
## ENSG00000000003 0.1630257
## ENSG00000000005        NA
## ENSG00000000419 0.1757326
## ENSG00000000457 0.9616577
## ENSG00000000460 0.8157061
## ...                   ...
## ENSG00000283115        NA
## ENSG00000283116        NA
## ENSG00000283119        NA
## ENSG00000283120        NA
## ENSG00000283123        NA
```

Ordering results by p value


```r
resOrdered <- res[order(res$pvalue),]
resOrdered
```

```
## log2 fold change (MLE): dex treated vs control 
## Wald test p-value: dex treated vs control 
## DataFrame with 38694 rows and 6 columns
##                  baseMean log2FoldChange      lfcSE      stat       pvalue
##                 <numeric>      <numeric>  <numeric> <numeric>    <numeric>
## ENSG00000152583  954.7709       4.368359 0.23713648  18.42129 8.867079e-76
## ENSG00000179094  743.2527       2.863888 0.17555825  16.31304 7.972621e-60
## ENSG00000116584 2277.9135      -1.034700 0.06505273 -15.90556 5.798513e-57
## ENSG00000189221 2383.7537       3.341544 0.21241508  15.73120 9.244206e-56
## ENSG00000120129 3440.7038       2.965211 0.20370277  14.55656 5.306416e-48
## ...                   ...            ...        ...       ...          ...
## ENSG00000283114         0             NA         NA        NA           NA
## ENSG00000283115         0             NA         NA        NA           NA
## ENSG00000283116         0             NA         NA        NA           NA
## ENSG00000283119         0             NA         NA        NA           NA
## ENSG00000283123         0             NA         NA        NA           NA
##                         padj
##                    <numeric>
## ENSG00000152583 1.342919e-71
## ENSG00000179094 6.037268e-56
## ENSG00000116584 2.927283e-53
## ENSG00000189221 3.500088e-52
## ENSG00000120129 1.607313e-44
## ...                      ...
## ENSG00000283114           NA
## ENSG00000283115           NA
## ENSG00000283116           NA
## ENSG00000283119           NA
## ENSG00000283123           NA
```

Adjusting p-value to 0.05


```r
res05 <- results(dds, alpha=0.05)
summary(res05)
```

```
## 
## out of 25258 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)     : 1237, 4.9% 
## LFC < 0 (down)   : 933, 3.7% 
## outliers [1]     : 142, 0.56% 
## low counts [2]   : 9033, 36% 
## (mean count < 6)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

How many are significant with an adjusted p-value < 0.05? How about 0.01? Save this last set of results as resSig01.

Accessing actual subset of data


```r
resSig05 <- subset(as.data.frame(res), padj < 0.05)
nrow(resSig05)
```

```
## [1] 2182
```

```r
resSig01 <- subset(as.data.frame(res), padj < 0.01)
nrow(resSig01)
```

```
## [1] 1437
```

Using either the previously generated anno object (annotations from the file annotables_grch38.csv file) or the mapIds() function (from the AnnotationDbi package) add annotation to your res01 results data.frame.

Using **mapIds()** function to add "SYMBOL" column form org.Hs.eg.db database and assign correct symbol based on "ENSEMBL" key. Multiple columns must be added separately, function only does one at a time. 


```r
resSig01$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(resSig01),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
#View(resSig01)
```

Arrange by adjusted p-value


```r
ord <- order( resSig01$padj )
#View(res01[ord,])
head(resSig01[ord,])
```

```
##                   baseMean log2FoldChange      lfcSE      stat
## ENSG00000152583   954.7709       4.368359 0.23713648  18.42129
## ENSG00000179094   743.2527       2.863888 0.17555825  16.31304
## ENSG00000116584  2277.9135      -1.034700 0.06505273 -15.90556
## ENSG00000189221  2383.7537       3.341544 0.21241508  15.73120
## ENSG00000120129  3440.7038       2.965211 0.20370277  14.55656
## ENSG00000148175 13493.9204       1.427168 0.10036663  14.21955
##                       pvalue         padj  symbol
## ENSG00000152583 8.867079e-76 1.342919e-71 SPARCL1
## ENSG00000179094 7.972621e-60 6.037268e-56    PER1
## ENSG00000116584 5.798513e-57 2.927283e-53 ARHGEF2
## ENSG00000189221 9.244206e-56 3.500088e-52    MAOA
## ENSG00000120129 5.306416e-48 1.607313e-44   DUSP1
## ENSG00000148175 6.929711e-46 1.749175e-42    STOM
```

Writing csv file of results


```r
write.csv(resSig01[ord,], "signif01_results.csv")
```


# Data Visualization

Getting gene of interest - CRISPLD2


```r
i <- grep("CRISPLD2", resSig01$symbol)

# printing value at index i 
resSig01[i,]
```

```
##                 baseMean log2FoldChange     lfcSE     stat       pvalue
## ENSG00000103196 3096.159       2.626034 0.2674705 9.818031 9.416441e-23
##                         padj   symbol
## ENSG00000103196 3.395524e-20 CRISPLD2
```

```r
# getting rowname (the ensemble gene id)
rownames(resSig01[i,])
```

```
## [1] "ENSG00000103196"
```

visualizing data with **plotCounts()** function


```r
# intgroup meaning interesting group 
# gene was pulled from previous function to locate CRISPLD2 in dataframe
plotCounts(dds, gene="ENSG00000103196", intgroup="dex")
```

![](class14_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

Returning data 


```r
d <- plotCounts(dds, gene="ENSG00000103196", intgroup="dex", returnData=TRUE)
head(d)
```

```
##                count     dex
## SRR1039508  774.5002 control
## SRR1039509 6258.7915 treated
## SRR1039512 1100.2741 control
## SRR1039513 6093.0324 treated
## SRR1039516  736.9483 control
## SRR1039517 2742.1908 treated
```

Building a box plot based on returned data from **plotCounts**


```r
boxplot(count ~ dex , data=d)
```

![](class14_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

Same result plotted with **ggplot()**


```r
library(ggplot2)
# color coded boxes, scaled y axis to be logarithmic, added a title
ggplot(d, aes(dex, count)) + geom_boxplot(aes(fill=dex)) + scale_y_log10() + ggtitle("CRISPLD2")
```

![](class14_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

## Which plot do you prefer? Ggplot offers a lot of aesthetic and functional parameters not available in base R. DataCamp course on ggplot can teach you

# MA and Volcano Plots

Adding a column of logicals to represent which samples  within res have a p-value less than 0.05


```r
res$sig <- res$padj<0.05

# How many of each?
table(res$sig)
```

```
## 
## FALSE  TRUE 
## 12963  2182
```



```r
sum(is.na(res$sig))
```

```
## [1] 23549
```

MA plot


```r
plotMA(res, ylim=c(-2,2))
```

![](class14_files/figure-html/unnamed-chunk-33-1.png)<!-- -->


Reducing noise (shrinking) and re-plotting


```r
resLFC <- lfcShrink(dds, coef=2)
resLFC
```

```
## log2 fold change (MAP): dex treated vs control 
## Wald test p-value: dex treated vs control 
## DataFrame with 38694 rows and 6 columns
##                  baseMean log2FoldChange      lfcSE       stat     pvalue
##                 <numeric>      <numeric>  <numeric>  <numeric>  <numeric>
## ENSG00000000003 747.19420    -0.31838595 0.15271739 -2.0846111 0.03710462
## ENSG00000000005   0.00000             NA         NA         NA         NA
## ENSG00000000419 520.13416     0.19883048 0.09744556  2.0403876 0.04131173
## ENSG00000000457 322.66484     0.02280238 0.13491699  0.1690242 0.86577762
## ENSG00000000460  87.68263    -0.11887370 0.20772938 -0.5726216 0.56690095
## ...                   ...            ...        ...        ...        ...
## ENSG00000283115  0.000000             NA         NA         NA         NA
## ENSG00000283116  0.000000             NA         NA         NA         NA
## ENSG00000283119  0.000000             NA         NA         NA         NA
## ENSG00000283120  0.974916    -0.05944174  0.1514839 -0.3944544  0.6932456
## ENSG00000283123  0.000000             NA         NA         NA         NA
##                      padj
##                 <numeric>
## ENSG00000000003 0.1630257
## ENSG00000000005        NA
## ENSG00000000419 0.1757326
## ENSG00000000457 0.9616577
## ENSG00000000460 0.8157061
## ...                   ...
## ENSG00000283115        NA
## ENSG00000283116        NA
## ENSG00000283119        NA
## ENSG00000283120        NA
## ENSG00000283123        NA
```


```r
plotMA(resLFC, ylim=c(-2,2))
```

![](class14_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

Make a volcano plot. Similarly, color-code by whether it’s significant or not.


```r
ggplot(as.data.frame(res), aes(log2FoldChange, -1*log10(pvalue), col=sig)) + 
    geom_point() + 
    ggtitle("Volcano plot")
```

```
## Warning: Removed 13578 rows containing missing values (geom_point).
```

![](class14_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

## Transforming data further for PCA analysis, clustering, etc. 

variance stabilizing transformation (VST) removes the dependence of the variance on the mean, particularly the high variance of the log counts when the mean is low.


```r
vsdata <- vst(dds, blind=FALSE)
```


Using DESeq provided function **plotPCA()**


```r
plotPCA(vsdata, intgroup="dex")
```

![](class14_files/figure-html/unnamed-chunk-38-1.png)<!-- -->

The sessionInfo() prints version information about R and any attached packages. It’s a good practice to always run this command at the end of your R session and record it for the sake of reproducibility in the future.


```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17134)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ggplot2_2.2.1              org.Hs.eg.db_3.5.0        
##  [3] AnnotationDbi_1.40.0       DESeq2_1.18.1             
##  [5] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
##  [7] matrixStats_0.53.1         Biobase_2.38.0            
##  [9] GenomicRanges_1.30.3       GenomeInfoDb_1.14.0       
## [11] IRanges_2.12.0             S4Vectors_0.16.0          
## [13] BiocGenerics_0.24.0       
## 
## loaded via a namespace (and not attached):
##  [1] locfit_1.5-9.1         Rcpp_0.12.16           lattice_0.20-35       
##  [4] rprojroot_1.3-2        digest_0.6.15          plyr_1.8.4            
##  [7] backports_1.1.2        acepack_1.4.1          RSQLite_2.1.1         
## [10] evaluate_0.10.1        pillar_1.2.1           zlibbioc_1.24.0       
## [13] rlang_0.2.0            lazyeval_0.2.1         rstudioapi_0.7        
## [16] data.table_1.11.2      annotate_1.56.2        blob_1.1.1            
## [19] rpart_4.1-13           Matrix_1.2-12          checkmate_1.8.5       
## [22] rmarkdown_1.9          labeling_0.3           splines_3.4.4         
## [25] BiocParallel_1.12.0    geneplotter_1.56.0     stringr_1.3.0         
## [28] foreign_0.8-69         htmlwidgets_1.2        RCurl_1.95-4.10       
## [31] bit_1.1-13             munsell_0.4.3          compiler_3.4.4        
## [34] pkgconfig_2.0.1        base64enc_0.1-3        htmltools_0.3.6       
## [37] nnet_7.3-12            tibble_1.4.2           gridExtra_2.3         
## [40] htmlTable_1.11.2       GenomeInfoDbData_1.0.0 Hmisc_4.1-1           
## [43] XML_3.98-1.11          bitops_1.0-6           grid_3.4.4            
## [46] xtable_1.8-2           gtable_0.2.0           DBI_1.0.0             
## [49] magrittr_1.5           scales_0.5.0           stringi_1.1.7         
## [52] XVector_0.18.0         genefilter_1.60.0      latticeExtra_0.6-28   
## [55] Formula_1.2-3          RColorBrewer_1.1-2     tools_3.4.4           
## [58] bit64_0.9-7            survival_2.41-3        yaml_2.1.18           
## [61] colorspace_1.3-2       cluster_2.0.6          memoise_1.1.0         
## [64] knitr_1.20
```


































