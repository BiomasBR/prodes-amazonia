Prodes-Amazônia experiments and pipelines developed using the SITS package
================

<img src="./inst/extdata/sticker/biomasbr_logo.jpeg" alt="RESTORE+ icon" align="right" height="150" width="150"/>

This repository brings together reproducible experiments and processing pipelines
from the Prodes-Amazônia project, developed using the SITS package. Its purpose
is to clearly and systematically document the adopted workflows,
providing references for experimentation, validation, and methodological
improvements, while supporting reproducibility and the continuous
evolution of project's analyses.

# Getting started

To use the scripts in this repository, clone the project to
your local machine using the command below:

``` sh
git clone https://github.com/BiomasBR/prodes-amazonia.git
```

After cloning, open the sits-prodes directory in RStudio and install the
package using the command:

``` r
devtools::install(".")
```

# Documentation
Visit our website to explore the sampling protocol and learn more about the data used:

``` sh
https://migualex.github.io/sits-prodes/
```

# Repository structure

- `data/`: Datasets used and generated throughout the analyses
- `inst/`: Supplementary resources required for the package to run
- `R/`: Package functions
- `scripts/`: Processing and experimentation routines

# License

The data and results available in this repository are licensed under the
terms of the Creative Commons license:

<img style="display: inline-block; vertical-align: middle; margin-right: 5px;" src="./inst/extdata/licenses/Cc-by-nc-sa_icon.png" alt="CC BY Icon" width="70">[Link](https://creativecommons.org/licenses/by-nc-sa/4.0/).

## Support

For questions, suggestions, or issues, please use the **Issues** section or
contact the repository maintainers.
