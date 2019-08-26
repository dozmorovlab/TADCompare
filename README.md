`TADCompare` is an R package providing tools for differential Topologically Associated Domain (TAD) detection and TAD boundary calling across multiple replicates. It has three main functions, `TADCompare` for differential TAD analysis, `TimeCompare` for time course analysis and `ConsensusTADs` for consensus boundary identification. 

```
install.packages(c('dplyr', 'PRIMME', 'cluster',
    'Matrix',
    'magrittr',
    'HiCcompare'))

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
```
The latest version of `TADCompare` can be directly installed from Github:

```
devtools::install_github('cresswellkg/TADCompare', build_vignettes = TRUE)
library(TADCompare)
```

Alternatively, the package can be installed from Bioconductor (to be submitted):

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TADCompare")
library(TADCompare)
```

## Input

There are three types of input accepted:

1. n x n contact matrices
2. n x (n+3) contact matrices
3. 3-column sparse contact matrices

It is required that the same format be used for each of the inputs to a given function or an error will occur. These formats are explained in depth in the [vignette](vignettes/TADCompare.Rmd).

# Usage

## TADcompare

Differential TAD detection (`TADCompare`) involves identifying differing TAD boundaries between two contact matrices. Accordingly, the input is the two contact matrices that we would like to find these boundaries in and their corresponding resolution. The output is a dataset containing boundary scores for all regions and a dataset containing all differential and non-differential TAD boundaries. 

```
#Load example contact matrices
data("rao_chr22_prim")
data("rao_chr22_rep")
#Find differential TADs
tads = TADCompare(rao_chr22_prim, rao_chr22_rep, resolution = 50000)
head(tads$TAD_Frame)
```
We can then print the set of regions with at least one TAD boundary:

```
  Boundary Gap_Score TAD_Score1 TAD_Score2 Differential Enriched_In            Type
16800000 -8.616720   1.447318   4.923084 Differential    Matrix 2    Differential
17200000  3.006105   3.037070   2.012388 Differential    Matrix 1 Strength Change
17250000  2.259350   6.080978   5.512345      Shifted    Matrix 1         Shifted
18700000 15.333645  10.385900   4.888309 Differential    Matrix 1 Strength Change
18750000 -6.981617   3.360815   6.293102      Shifted    Matrix 2         Shifted
18800000 -8.491078   1.021460   4.425071      Shifted    Matrix 2         Shifted
```

Additionally, we can plot the proportion of different types of boundaries:

```

```

