# PCA with Partial SVD Implementation

[![R-CMD-check](https://github.com/dereklei12/pca_R/workflows/R-CMD-check/badge.svg)](https://github.com/dereklei12/pca_R/actions)

## Overview

Efficient PCA implementations for large datasets, providing drop-in replacements for `mixOmics::pca` with improved performance using partial SVD methods.

- **pca_pSVD**: Partial SVD for general use (fastest when `k << min(n,p)`)

## Key Features

- High performance on large datasets (tested on 30,000+ cell scRNA-seq data)
- Native sparse matrix support (`dgCMatrix`, `dgRMatrix`)
- Full compatibility with `mixOmics::pca` output format
- Cross-platform tested (Ubuntu, Windows, macOS)
- 50+ comprehensive test cases

## Installation

### Core dependencies
```r
install.packages(c("testthat", "ggplot2", "dplyr", "tidyr", "knitr"))
install.packages(c("RSpectra", "rARPACK", "Matrix"))
BiocManager::install("mixOmics")
```

## Quick Start

```r
# Load function
source("R/pca_pSVD.R")

# Basic usage
pca_result <- pca_pSVD(
    X = your_data,
    ncomp = 10,
    center = TRUE,
    scale = TRUE
)

# Access results
pca_result$x                    # Scores
pca_result$rotation             # Loadings
pca_result$sdev                 # Standard deviations
pca_result$prop_expl_var$X      # Variance explained
```

## Running Tests

```r
library(testthat)
test_dir("test/testthat")
```

## Project Structure

```
.
├── .github/workflows/
│   └── R-CMD-check.yml         # CI/CD configuration
├── R/
│   ├── pca_pSVD.R              # Partial SVD (main)
│   ├── pca_asym.R              # Asymmetric optimization (Experimental)
│   └── pca_em_woodbury.R       # EM-PPCA(Experimental Undone)
├── test/testthat/
│   ├── test-pca_pSVD.R         # Main tests
│   └── test-pca_woodbury.R     # EM-PPCA tests
├── Data/
│   └── HaffinaCovidPBMC_30000cells_dense.rds
├── pca_pSVD.Rmd                # Benchmark report
└── README.md
```

## Test Coverage

- ✅ Basic functionality and output validation
- ✅ Comparison with `prcomp` and `mixOmics::pca`
- ✅ Sparse matrix support
- ✅ Edge cases (wide/tall matrices, zero-variance)
- ✅ Multilevel/repeated measures
- ✅ Cross-platform numerical stability

See `pca_pSVD.Rmd` for detailed benchmarks.
