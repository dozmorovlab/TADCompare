# ToDo, delete after completion

- Add `README.md`, formatted similar to `SpectralTAD` package. See "Add to README.md" section at https://github.com/mdozmorov/Programming_notes#r-packages

- Upon installation, the following message appears - doesn't look good. Possible to fix?
```
Registered S3 method overwritten by 'R.oo':
  method        from       
  throw.default R.methodsS3
```



## Vignette

- Address current ???
- Split the current vignette into four, titles quoted:
    1) "Input data". Includes file format description (beginning of the current vignette) and benchmarking of file format conversions (end)
    2) "TAD comparison between two conditions"
    3) "TAD comparison across time course"
    4) "Consensus TAD calling"
    - Write each vignette to be completely self-standing. Assume zero prior knowkedge, be verbose. Assume a user never read the paper, clarify everything. Make your vignettes separate papers, that is, reasonably follow IMRaD format.

- Add gene enrichment analysis into vignette using rGREAT
