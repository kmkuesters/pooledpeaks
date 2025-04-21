<!-- badges: start -->
  [![R-CMD-check](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Package Overview

pooledpeaks is designed for analyzing genetic data obtained from Fragment
Analysis output files (`.fsa`) of pooled biological samples. It provides functions
for a comprehensive analysis pipeline from processing `.fsa` files, to cleaning
the peak data, and conducting population genetic analyses. Some features are
listed below and a usage example of the entire pipeline is included as a vignette.
Please check out the [Contributing Guidelines](inst/CONTRIBUTING.md) for
information on how to add to this package.

## Installation Instructions

You can install the package directly from GitHub using the following instructions:

Open R and copy the following code into your console

### Install devtools and pooledpeaks from GitHub

``` r
install.packages("devtools")
devtools::install_github("kmkuesters/pooledpeaks")
```



### Install pooledpeaks directly from CRAN

```r
install.packages("pooledpeaks")
```

## Features

For a detailed example of how to apply the functions contained in this package
please see the [Introduction to Using the pooledpeaks
Workflow](https://CRAN.R-project.org/package=pooledpeaks/vignettes/workflow.html).
Example data can be found on GitHub under the inst/extdata folder including `.fsa`
files and a formatted "Multiplex_frequencies.txt" file for the Genetic Analysis
portion.

  * **Peak Scoring:** Process `.fsa` files and score peaks contained therein.

```r
check_fsa_v_batch()
fsa_metadata()
fsa_batch_imp()
associate_dyes()
score_markers_rev3()
```

  * **Data Manipulation:** Clean and prepare peak data for downstream analyses.

```r
clean_scores()
lf_to_tdf()
data_manipulation()
Rep_check()
PCDM()
LoadData()
```

  * **Population Genetics Analysis:**

    * Calculate Gene Identity Matrix and Genetic Distance Matrix

    * Calculate diversity indices

    * Calculate differentiation indices

    * Perform cluster analysis

```r
TypedLoci()
GeneIdentityMatrix()
GeneticDistanceMatrix()
GST()
JostD()
cluster()
```

  * **Visualization:** Visualize the peak scoring and genetic analysis results.

```r
MDSplot()
```

## Sample Data


The sample `.fsa` files included in this package are provided for demonstration
purposes and originate from two sources:

- **Schistosoma haematobium laboratory isolates**, used for preliminary testing
of the `pooledpeaks` workflow. These data contain no identifiable or human
subject information.
- **De-identified Schistosoma mansoni samples** from a three studies conducted in
Brazil, extracted from discarded human waste. These files were originally used
for genetic analysis and are shared here in anonymized form to illustrate
compatibility with additional species and data sources.These studies are
described in detail by Long et al. (2022), available at
https://www.nature.com/articles/s41598-022-04776-0:

These files are intended solely to demonstrate the functionality of the
`pooledpeaks` package and are not for diagnostic or clinical use.To access the
example `.fsa` files included with the package, use the following path within R:

```r
system.file("extdata", package = "pooledpeaks")
```

The pooledpeaks package was developed by the Blanton Lab as part of Kathleen
Kuesters' dissertation.


## References:

* Covarrubias-Pazaran et al. (2016) [doi:10.1186/s12863-016-0365-6](https://doi.org/10.1186/s12863-016-0365-6)
* Long et al. (2022) [doi:10.1038/s41598-022-04776-0](https://doi.org/10.1038/s41598-022-04776-0)
* Jost (2008) [doi:10.1111/j.1365-294x.2008.03887.x](https://doi.org/10.1111/j.1365-294x.2008.03887.x)
* Nei (1973) [doi:10.1073/pnas.70.12.3321](https://doi.org/10.1073/pnas.70.12.3321)
* Foulley et al. (2006) [doi:10.1016/j.livprodsci.2005.10.021](https://doi.org/10.1016/j.livprodsci.2005.10.021)
* Chao et al. (2008) [doi:10.1111/j.1541-0420.2008.01010.x](https://doi.org/10.1111/j.1541-0420.2008.01010.x)
