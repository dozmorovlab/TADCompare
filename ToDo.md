# ToDo, delete after completion

- Address current ???

- Make `README.md` formatted similar to `SpectralTAD` package. 
    - Add this to both `SpectralTAD` and `TADcompare`: The developmental version is available at https://github.com/jstansfield0/multiHiCcompare, the stable version is available at https://github.com/dozmorovlab/multiHiCcompare.
    - Add Travis CI and other bells and whistles to README.md, like in https://github.com/dozmorovlab/SpectralTAD
    - See "Add to README.md" section at https://github.com/mdozmorov/Programming_notes#r-packages
    - Take it seriously, README is the first thing people see when google your package

- Upon installation, the following message appears - doesn't look good. Possible to fix?
```
Registered S3 method overwritten by 'R.oo':
  method        from       
  throw.default R.methodsS3
```

- Add gene enrichment analysis into vignette using rGREAT

- Rename vignettes as: `1_Input_Data`, `2_TADCompare` etc. Adjust links to them in the TADcompare vignette.

- Remove `*.bak` and `*.sav` files from `vignettes` folder. And, don't commit at the first place 
