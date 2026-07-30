# Developer Guide for modulr

Welcome to the developer guide for **modulr** (v0.3.1), an R package in
the `foundr` ecosystem that provides Weighted Gene Co-expression Network
Analysis (WGCNA) module creation, soft-thresholding topology analysis,
topological overlap matrix (TOM) calculations, and eigentrait extraction
for multiparent genetic cross data (specifically Collaborative Cross
mouse studies).

This document details package architecture, data models, developer
environment setup, function references, organizing WGCNA objects for
downstream `foundr` routines, testing procedures, and release
guidelines.

------------------------------------------------------------------------

## 1. Ecosystem & Package Architecture

`modulr` provides the network module construction layer in the `foundr`
ecosystem. It transforms wide expression/trait data frames or harmonized
long-format datasets from `foundrHarmony` into structured WGCNA
co-expression trait modules (`wgcnaModules` and `listof_wgcnaModules` S3
objects).

### Role in the Ecosystem

    Raw Data (Excel/CSV)
           │
           ▼
    foundrHarmony  [normalize & standardize]
           │
           ├─────────────────────┐
           ▼                     ▼
    foundr (Core)  ◄────────  modulr  [WGCNA trait modules & eigentraits]
      [partition & stats]        │
           │                     │
           └──────────┬──────────┘
                      ▼
                 foundrShiny   [visualize & interact]

| Package / Repository | Role / Description |
|----|----|
| [`foundrHarmony`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/foundrHarmony/README.md) | Data normalization, standardizing raw matrices, specialized dataset parsers, and multi-study consolidation |
| [modulr](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/README.md) | WGCNA co-expression module creation, soft-thresholding topology estimation, TOM dissimilarity matrices, and eigentrait extraction |
| `foundr` | Core statistical models, orthogonal variance partitioning (`partition()`), summary stats (`strainstats()`), and visualization generics |
| `foundrShiny` | Interactive Shiny web application for exploratory data visual analytics |

------------------------------------------------------------------------

## 2. Developer Environment & Setup

### System Requirements

- **R Version**: `>= 3.5.0`
- **Core Dependencies**: Listed in
  [`DESCRIPTION`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/DESCRIPTION)
  (`dplyr`, `tidyr`, `purrr`, `tibble`, `WGCNA`, `fastcluster`,
  `dynamicTreeCut`, `stringr`).

### Local Development Workflow

To work on `modulr` locally:

``` r

# 1. Open the project in RStudio / Posit Workbench or VS Code

# 2. Install devtools if needed
install.packages("devtools")

# 3. Load package into active session during development
devtools::load_all()

# 4. Generate updated documentation & NAMESPACE
devtools::document()

# 5. Run local package check
devtools::check(cran = FALSE, vignettes = FALSE)
```

------------------------------------------------------------------------

## 3. WGCNA Object Structures & Technical Organization for `foundr`

In `modulr`, “modules” refer to WGCNA clusters of highly correlated
phenotypic traits or genes. `modulr` organizes network results into
structured S3 objects designed for seamless pass-through into `foundr`
routines (such as
[`foundr::partition()`](https://rdrr.io/pkg/foundr/man/partition.html)
or
[`foundr::strainstats()`](https://rdrr.io/pkg/foundr/man/strainstats.html)).

### `wgcnaModules` S3 Class

A `wgcnaModules` object is a structured S3 list containing:

``` r

wgcnaModules_object <- list(
  ID      = ID_df,    # Data frame containing strain, sex, condition, and animal IDs
  dendro  = dendro,   # hclust object from hierarchical clustering of TOM dissimilarity
  eigen   = eigen_df, # Data frame of module eigentraits (ME) with sample IDs as rownames
  modules = kME_df    # Long-format data frame of traits, assigned module colors, and kME correlations
)
class(wgcnaModules_object) <- c("wgcnaModules", class(wgcnaModules_object))
attr(wgcnaModules_object, "params") <- params
```

#### Fields Breakdown

1.  **`ID`** (`data.frame` or `tbl_df`): Standardized metadata
    containing columns `ID` (`strain_sex_condition`), `animal`,
    `strain`, `sex`, `condition`/`diet`. Generated via
    [`wgcna_ID()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgnca_ID.R).
2.  **`dendro`** (`hclust`): Hierarchical clustering tree produced by
    [`fastcluster::hclust()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R#L33-L35)
    using average linkage on the TOM dissimilarity matrix.
3.  **`eigen`** (`data.frame`): Matrix/data frame of module eigentraits
    (first principal component of module trait expressions), with
    leading `"ME"` prefix removed. Rownames match sample `ID`
    identifiers.
4.  **`modules`** (`data.frame`): Data frame generated by
    [`module_factors()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/module_factors.R)
    containing:
    - `trait`: Measurement or gene symbol.
    - `module`: Assigned module color factor ordered by module size
      (number of traits).
    - `kME`: Module membership correlation coefficient
      ($`r \in [-1, 1]`$) with the eigentrait.

### `listof_wgcnaModules` S3 Class

When network analysis is executed across subsets (e.g. split by `sex`,
`strain`, or variance partition components `cellmean`, `signal`, `rest`,
`noise`), `modulr` packages the results into a `listof_wgcnaModules` S3
object:

``` r

listof_object <- purrr::map(split(traitContr, traitContr$sex), wgcnaModules, params)
class(listof_object) <- c("listof_wgcnaModules", class(listof_object))
attr(listof_object, "params") <- attr(listof_object[[1]], "params")
```

------------------------------------------------------------------------

## 4. Package Function Catalog

`modulr` consists of 10 primary R source files in
[`R/`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R):

| Function | Source File | Description |
|----|----|----|
| [`load_wgcna()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/load_wgcna.R) | [`load_wgcna.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/load_wgcna.R) | Loads legacy WGCNA `.Rdata` objects and converts them into standard `wgcnaModules` or `listof_wgcnaModules` objects |
| [`wgcna_harmonize()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_harmonize.R) | [`wgcna_harmonize.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_harmonize.R) | Harmonizes `.Rdata` files into cached `.rds` serialized objects with parameter encoding |
| [`wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R) | [`wgcna_modules.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R) | Core pipeline runner: pivots data, calculates TOM similarity, performs dynamic tree cut, merges close modules, computes kME correlations |
| [`listof_wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R#L101-L108) | [`wgcna_modules.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R) | Splits trait dataset by factor (e.g. sex) and maps [`wgcnaModules()`](https://byandell-sysgen.github.io/modulr/reference/wgcnaModules.md) across subsets |
| [`plot_wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R#L142-L157) | [`wgcna_modules.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R) | Visualizes gene/trait dendrogram with dynamic tree cut module colors using [`WGCNA::plotDendroAndColors`](https://rdrr.io/pkg/WGCNA/man/plotDendroAndColors.html) |
| [`wgcna_pivot()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_pivot.R) | [`wgcna_pivot.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_pivot.R) | Pivots long-format harmonized data into wide numeric matrix format (traits as columns, sample ID as rownames) |
| [`wgcna_ID()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgnca_ID.R) | [`wgnca_ID.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgnca_ID.R) | Standardizes sample identifiers (`strain`, `sex`, `condition`, `animal`), supporting annotation mappings and fixing non-standard strain tags (`A129`/`X129` -\> `129`) |
| [`wgcna_dist()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_dist.R) | [`wgcna_dist.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_dist.R) | Calculates Topological Overlap Matrix (TOM) dissimilarity: `1 - TOMsimilarity(adjacency(...))` |
| [`wgcna_params()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_params.R) | [`wgcna_params.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_params.R) | Validates and fills default WGCNA parameters (`power = 6`, `signType = "unsigned"`, `minSize = 4`, `cutHeight = 0.995`, `thresholdMEDiss = 0.25`) |
| [`wgcna_topology()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R#L14-L32) | [`wgcna_topology.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R) | Multi-threaded soft-thresholding power estimation using [`WGCNA::pickSoftThreshold`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html) |
| [`ggplot_wgcna_topology()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R#L45-L65) | [`wgcna_topology.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R) | Two-panel ggplot2 diagnostic plot for scale-free topology fit $`R^2`$ and mean connectivity |
| [`module_factors()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/module_factors.R) | [`module_factors.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/module_factors.R) | Extracts trait-to-module assignments and kME correlation values from WGCNA results |
| [`plot_null()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/plot_null.R) | [`plot_null.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/plot_null.R) | Utility function displaying fallback void ggplot message when dataset is empty |

------------------------------------------------------------------------

## 5. Testing & Verification

To verify `modulr` package integrity:

``` r

# 1. Regenerate docs
devtools::document()

# 2. Test execution of core workflow on mock matrix
params <- modulr::wgcna_params(list(power = 6, minSize = 4))
```

------------------------------------------------------------------------

## 6. Release & Versioning Guidelines

1.  **Version Numbering**: Follow Semantic Versioning
    (`MAJOR.MINOR.PATCH`). Update
    [`DESCRIPTION`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/DESCRIPTION).
2.  **Roxygen Documentation**: Update docstrings in `R/*.R` source files
    before running `devtools::document()`.
3.  **No Automatic Commits**: Per repository rules, perform all local
    testing and verification commands, leaving `git commit` and
    `git push` for manual developer execution.
