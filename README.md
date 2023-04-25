# roessler_netzl_et_al2023
This repository contains the code for the manuscript "Characterizing SARS-CoV-2 neutralization profiles after bivalent boosting using antigenic cartography" by Annika Rössler, Antonia Netzl, et. al., 2023

The raw data will be added at the stage of acceptance of the manuscript in the `data` directory. The code for the analyses shown in the main manuscript can be found in the `code` directory. 
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
