# ToDo, delete after completion

- Add a function that returns the violin plot of TAD width. If one sample is provided, there will be one plot. If two, or multiple, then several. E.g., for two samples the width distributions may look like https://seaborn.pydata.org/_images/seaborn-violinplot-5.png. They should be comparable, t-test or Wilcoxon. Goal - to visualize and answer the question, is the TAD width changes between the conditions?
    - Besides plotting, the function should output basic TAD statistics (for each sample). These parameters may change, so should be tested for differences:
        - Total number of TADs, gaps
        - Width distribution (min, mean, median, max) of TADs, gaps
        - Boundary strength statistics (min, mean, median, max) of TADs
            - Variance of boundary strength (also may change)
 
 - For differential detection and time course analysis, add a function that returns the stacked bar plot, like https://www.smartsheet.com/sites/default/files/ic-excel-stacked-bar-charts-part-to-hole.png. Goal - to see the proportions of types of TAD boundary changes, both for two-group and time course analyses. 
 
- Add `README.md`, formatted similar to `SpectralTAD` package. See "Add to README.md" section at https://github.com/mdozmorov/Programming_notes#r-packages

## Vignette

- Address current ???
- Split the current vignette into four, titles quoted:
    1) "Input data". Includes file format description (beginning of the current vignette) and benchmarking of file format conversions (end)
    2) "TAD comparison between two conditions"
    3) "TAD comparison across time course"
    4) "Consensus TAD calling"
    - Write each vignette to be completely self-standing. Assume zero prior knowkedge, be verbose. Assume a user never read the paper, clarify everything. Make your vignettes separate papers, that is, reasonably follow IMRaD format.

- Add gene enrichment analysis into vignette using rGREAT
