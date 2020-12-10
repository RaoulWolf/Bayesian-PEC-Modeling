Bayesian PEC modeling
================

This repository contains the data and the R scripts for Wolf & Tollefsen
2021, Environmental Science & Technology. *This README will be updated
with the final details of the publication as soon as they are
available.*

## Requirements

To analyze the data, the statistical software [R 4.0.1 or
higher](https://cloud.r-project.org/) needs to be installed, as well as
[Rtools 40](https://cloud.r-project.org/bin/windows/Rtools/) for users
operating from Windows.

The package {brms} needs to be installed as well. Within R, execute the
following command to install:

``` r
install.packages("brms")
```

## Information

All scripts are inside the `scripts/` folder of this repository. For
each of the three campaigns, a separate R file exists: `sorfjord.R` for
the Sørfjord campaign, `kaldvellfjord.R` for the Kaldvellfjord campaign,
and `oslofjord.R` for the Oslofjord campaign.

Additionally, the file `empirical_prior.R` contains a custom function to
calculate empirical priors. Detailed information on the modelling
procedure and the calculations of the empirical priors are given in the
Supporting Information of the publication.
