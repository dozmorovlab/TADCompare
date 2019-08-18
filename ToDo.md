# ToDo, delete after completion

- Add a function that returns the violin plot of TAD width. If one sample is provided, there will be one plot. If two, or multiple, then several. E.g., for two samples the width distributions may look like https://seaborn.pydata.org/_images/seaborn-violinplot-5.png. They should be comparable, t-test or Wilcoxon. Goal - to visualize and answer the question, is the TAD width changes between the conditions?
    - Besides plotting, the function should output basic TAD statistics (for each sample). These parameters may change, so should be tested for differences:
        - Total number of TADs, gaps
        - Width distribution (min, mean, median, max) of TADs, gaps
        - Boundary strength statistics (min, mean, median, max) of TADs
            - Variance of boundary strength (also may change)
 
 - For differential detection and time course analysis, add a function that returns the stacked bar plot, like https://www.smartsheet.com/sites/default/files/ic-excel-stacked-bar-charts-part-to-hole.png. Goal - to see the proportions of types of TAD boundary changes, both for two-group and time course analyses. 
 
- Add `README.md`, formatted similar to `SpectralTAD` package. See "Add to README.md" section at https://github.com/mdozmorov/Programming_notes#r-packages

- Installation error: `remotes::install_github("cresswellkg/TADCompare", build_vignettes = TRUE)`
```
E  creating vignettes (18.6s)
   --- re-building ‘TADCompare.Rmd’ using rmarkdown
   Warning in data("rao_chr22_prim") : data set 'rao_chr22_prim' not found
   Quitting from lines 60-63 (TADCompare.Rmd) 
   Error: processing vignette 'TADCompare.Rmd' failed with diagnostics:
   object 'rao_chr22_prim' not found
   --- failed re-building ‘TADCompare.Rmd’
   
   SUMMARY: processing the following file failed:
     ‘TADCompare.Rmd’
   
   Error: Vignette re-building failed.
   Execution halted
```
