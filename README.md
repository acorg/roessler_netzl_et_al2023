# roessler_netzl_et_al2023
This repository contains the code for the manuscript "Characterizing SARS-CoV-2 neutralization profiles after bivalent boosting using antigenic cartography" by Annika Rössler, Antonia Netzl, et. al., 2023. Please cite the original publication in Nature Communications if any data or code from this repository is used. 

The repository's DOI
[![DOI](https://zenodo.org/badge/494721761.svg)](https://zenodo.org/badge/latestdoi/494721761)
was created with Zenodo (https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)

Raw data can be found in the `data` directory. The code for the analyses shown in the main manuscript can be found in the `code` directory. To obtain a titer table for antigenic map construction, execute the `excel_to_titertable.R` script. To construct maps, first execute the `make_map.R` followed by the `reactivity_adjusted_maps.R` scripts. Antibody landscapes are fit in the `ablandscapes_fit_multi.R` and  `ablandscapes_fit_idvl.R` scripts. The `3D_landscapes.Rmd` creates an html script with 3D illustrations of the landscapes. 

All SOM analyses can be found in the `som` directory. 

The `function` directory contains utility functions, including two bash scripts to screenshot and crop html landscapes. For cropping ImageMagick was used: 
ImageMagick Studio LLC. (2023). ImageMagick. Retrieved from https://imagemagick.org

All analyses were performed in R version 4.2.2 (2022-10-31).
R Core Team (2022). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna,
  Austria. URL https://www.R-project.org/.
  
Antigenic maps were constructed using the Racmacs package, Version 1.1.35:
Wilks S (2022). _Racmacs: R Antigenic Cartography Macros_. https://acorg.github.io/Racmacs,
  https://github.com/acorg/Racmacs.
  
Antibody landscapes were constructed using the ablandscapes package, Version 1.1.0: 
Wilks S (2021). _ablandscapes: Making Antibody landscapes Using R_. R package
  version 1.1.0, <https://github.com/acorg/ablandscapes>.

Geometric mean titers and fold changes were calculated using the titertools package, Version 0.0.0.9001:
Wilks SH (2022). _titertools: A statistical toolkit for the annalysis of censored
  titration data_. R package version 0.0.0.9001.
