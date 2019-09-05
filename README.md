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

# Input

There are three types of input accepted:

1. n x n contact matrices
2. n x (n+3) contact matrices
3. 3-column sparse contact matrices

It is required that the same format be used for each of the inputs to a given function or an error will occur. These formats are explained in depth in the [vignette](vignettes/Data_Input.Rmd).

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
## TimeCompare

Another function of the TADCompare package is time course analysis of TAD boundary change. 

```
data("time_mats")
Time_Mats = TimeCompare(time_mats, resolution = 50000)
head(Time_Mats$Time_Bounds)
```
The resulting output is:

```
 Coordinate   Sample 1   Sample 2   Sample 3   Sample 4 Consensus_Score
1   16900000 -0.6733709 -0.7751516 -0.7653696 15.1272253     -0.71937026
2   17350000  3.6406563  2.3436229  3.0253018  0.7840556      2.68446232
3   18850000  0.6372268  6.3662245 -0.7876844  6.9255446      3.50172563
4   20700000  1.5667926  3.0968633  2.9130479  2.8300136      2.87153075
5   22000000 -1.0079676 -0.7982571  0.6007264  3.1909178     -0.09876534
6   22050000 -1.0405532 -0.9892242 -0.2675822  4.2737511     -0.62840320
                Category
1     Late Appearing TAD
2 Early Disappearing TAD
3    Early Appearing TAD
4            Dynamic TAD
5     Late Appearing TAD
6     Late Appearing TAD
```

## ConsensusTADs

ConsensusTADs uses a novel metric called the consensus boundary score to identify TAD boundaries summarized across multiple contact matrices. It can operate on an unlimited number of replicates, time points or conditions.

```
data("time_mats")
con_tads = ConsensusTADs(time_mats, resolution = 50000)
head(con_tads)
```

```
 Coordinate  Sample 1 Sample 2   Sample 3 Sample 4 Consensus_Score
1   18850000 0.6372268 6.366224 -0.7876844 6.925545        3.501726
2   28450000 3.1883107 3.313883  3.1711743 3.913620        3.251097
3   30000000 2.0253285 3.477652  2.9314321 3.926354        3.204542
4   32350000 2.6978488 2.455860  3.5131909 3.550942        3.105520
5   36900000 3.0731406 3.153978  3.1861296 4.489285        3.170054
```

# Contributions and Support

Suggestions for new features and bug reports are welcome. Please, create a new issue for any of these or contact the author directly: @cresswellkg (cresswellkg[at]vcu[dot]edu)

# Contributors

Authors: @cresswellkg (cresswellkg[at]vcu[dot]edu) & @mdozmorov (mikhail.dozmorov[at]vcuhealth[dot]org)



