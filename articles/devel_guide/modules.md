# modulr Function Index & WGCNA Trait Module Breakdown

## Function Index & WGCNA Trait Module Breakdown

This document provides a comprehensive technical reference for all
functions in **modulr** and details how WGCNA trait modules and
`wgcnaModules` objects are organized for downstream `foundr` analysis
routines.

> **NOTE:** In `modulr`, **“modules”** refer specifically to WGCNA
> co-expression clusters of highly correlated phenotypic traits or gene
> expressions summarized by eigentraits.

------------------------------------------------------------------------

### 1. Technical Organization of WGCNA Objects for `foundr` Routines

The primary objective of `modulr` is to ingest trait data, execute
network clustering, and structure WGCNA objects so that extracted
eigentraits and module membership factors can be passed seamlessly into
core `foundr` functions (such as variance partitioning via
[`foundr::partition()`](https://rdrr.io/pkg/foundr/man/partition.html)
or strain summary statistics via
[`foundr::strainstats()`](https://rdrr.io/pkg/foundr/man/strainstats.html)).

#### Integration Flow with `foundr`

``` mermaid
flowchart LR
    subgraph modulrPkg ["modulr Package"]
        wgcnaMod["wgcnaModules() / load_wgcna()"]
        eigenExtract["out$eigen (Eigentraits Data Frame)"]
        kMEExtract["out$modules (Trait-Module kME Data Frame)"]
        wgcnaMod --> eigenExtract
        wgcnaMod --> kMEExtract
    end

    subgraph foundrPkg ["foundr Package"]
        partition["foundr::partition() Orthogonal Variance Decomposition"]
        stats["foundr::strainstats() Strain Linear Models & F-Stats"]
    end

    subgraph shinyPkg ["foundrShiny Package"]
        dashboard["Interactive Module & Eigentrait Visualization"]
    end

    eigenExtract --> partition
    eigenExtract --> stats
    kMEExtract --> dashboard
    partition --> dashboard
```

#### Component Object Specifications

A harmonized `wgcnaModules` S3 object comprises four principal
components:

1.  **`out$ID`**: Sample metadata table produced by
    [`wgcna_ID()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgnca_ID.R)
    maintaining alignment across `strain`, `sex`, `condition`, and
    `animal`.
2.  **`out$dendro`**: The `hclust` hierarchical tree object generated
    from TOM dissimilarity, allowing visual inspection of cluster
    branches via
    [`plot_wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R#L142-L157).
3.  **`out$eigen`**: Module eigentrait matrix where columns correspond
    to module colors (e.g. `blue`, `turquoise`, `brown`) and rownames
    correspond to standardized sample `ID`s. Eigentraits represent the
    1st principal component of each module.
4.  **`out$modules`**: Tidy data frame output from
    [`module_factors()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/module_factors.R)
    containing mapping between individual traits, assigned `module`
    color factors (ordered by module trait count), and `kME`
    correlations ($`r \in [-1, 1]`$).

------------------------------------------------------------------------

### 2. Exhaustive Function Breakdown by Functional Area

#### Category 1: Ingestion & Legacy Object Loading

##### [`load_wgcna()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/load_wgcna.R)

- **Source File**:
  [`R/load_wgcna.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/load_wgcna.R)
- **Signature**:
  `load_wgcna(moddir, modobj = "WGCNA_objects_ms10.Rdata", params = list(signType = "unsigned", power = 12, minSize = 4), listof = TRUE, annot = NULL)`
- **Description**: Loads an `.Rdata` object containing pre-calculated
  WGCNA structures (`merge`, `kMEs`), strips the `"ME"` prefix from
  eigentraits, computes standardized sample IDs via
  [`wgcna_ID()`](https://byandell-sysgen.github.io/modulr/reference/wgcna_ID.md),
  updates parameters, and assigns the `wgcnaModules` or
  `listof_wgcnaModules` S3 class.

##### [`wgcna_harmonize()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_harmonize.R)

- **Source File**:
  [`R/wgcna_harmonize.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_harmonize.R)
- **Signature**:
  `wgcna_harmonize(mod = "WGCNAmodule", moddir = ".", modRdata = NULL, params = NULL, annot = NULL, harmonizeddir = ".", force = FALSE)`
- **Description**: Automation wrapper that checks for existing RDS
  representations, invokes
  [`load_wgcna()`](https://byandell-sysgen.github.io/modulr/reference/load_wgcna.md)
  if needed, and caches the resulting `wgcnaModules` object as an `.rds`
  file with parameter-encoded filenames.

------------------------------------------------------------------------

#### Category 2: Data Pivoting & Sample ID Standardization

##### [`wgcna_pivot()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_pivot.R)

- **Source File**:
  [`R/wgcna_pivot.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_pivot.R)
- **Signature**: `wgcna_pivot(object)`
- **Description**: Transforms long-format trait data frames into wide
  matrices with traits as columns and unified `ID` strings as rownames.
  Returns a list with `$matrix` and `$ID`.

##### [`wgcna_ID()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgnca_ID.R)

- **Source File**:
  [`R/wgnca_ID.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgnca_ID.R)
- **Signature**:
  `wgcna_ID(object, condition_under = TRUE, annot = NULL)`
- **Description**: Parses sample identifiers from matrix rownames or
  data frames. Unifies strain nomenclature (e.g. converting `A129` or
  `X129` to `129`), separates animal IDs, and integrates external
  annotation mappings when provided.

------------------------------------------------------------------------

#### Category 3: Topology & Soft-Thresholding Analysis

##### [`wgcna_topology()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R#L14-L32)

- **Source File**:
  [`R/wgcna_topology.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R)
- **Signature**:
  `wgcna_topology(object, powers = c(c(1:10), seq(from = 12, to = 20, by = 2)), cores = parallel::detectCores(), verbose = 0, ...)`
- **Description**: Enables multi-threaded processing via
  [`WGCNA::enableWGCNAThreads()`](https://rdrr.io/pkg/WGCNA/man/allowWGCNAThreads.html)
  and calls
  [`WGCNA::pickSoftThreshold()`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html)
  to evaluate scale-free topology model fit ($`R^2`$) and mean
  connectivity across candidate soft-thresholding powers $`\beta`$.

##### [`ggplot_wgcna_topology()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R#L45-L65)

- **Source File**:
  [`R/wgcna_topology.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_topology.R)
- **Signature**: `ggplot_wgcna_topology(object, cutoff = 0.9, ...)`
- **Description**: Creates a side-by-side two-panel ggplot2 figure
  showing (1) signed $`R^2`$ vs soft threshold power with a cutoff line,
  and (2) mean connectivity vs soft threshold power.

------------------------------------------------------------------------

#### Category 4: TOM Distance & Module Construction

##### [`wgcna_dist()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_dist.R)

- **Source File**:
  [`R/wgcna_dist.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_dist.R)
- **Signature**: `wgcna_dist(object, params)`
- **Description**: Computes the Topological Overlap Matrix (TOM)
  dissimilarity matrix:
  `1 - WGCNA::TOMsimilarity(WGCNA::adjacency(object, type, power))`.

##### [`wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R)

- **Source File**:
  [`R/wgcna_modules.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R)
- **Signature**: `wgcnaModules(object, params = NULL, ...)`
- **Description**: The core function orchestrating network construction.
  Pivots data, calculates TOM dissimilarity, performs hierarchical
  clustering and dynamic tree cutting, merges close modules, calculates
  signed kME correlations, and returns a `wgcnaModules` object.

##### [`listof_wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R#L101-L108)

- **Source File**:
  [`R/wgcna_modules.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R)
- **Signature**: `listof_wgcnaModules(traitContr, params = NULL)`
- **Description**: Maps
  [`wgcnaModules()`](https://byandell-sysgen.github.io/modulr/reference/wgcnaModules.md)
  across subsets of data (e.g., split by `sex` or variance partition
  components) using
  [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html).

##### [`module_factors()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/module_factors.R)

- **Source File**:
  [`R/module_factors.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/module_factors.R)
- **Signature**: `module_factors(kME, colors)`
- **Description**: Merges trait names, assigned module colors, and kME
  correlations into a long-format data frame, ordering module color
  factors by size.

------------------------------------------------------------------------

#### Category 5: Configuration & Utilities

##### [`wgcna_params()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_params.R)

- **Source File**:
  [`R/wgcna_params.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_params.R)
- **Signature**: `wgcna_params(params = NULL)`
- **Description**: Sets and validates parameter defaults (`signType`,
  `power`, `minSize`, `method`, `cutHeight`, `split`, `thresholdKME`,
  `thresholdMEDiss`, `verbose`).

##### [`plot_wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R#L142-L157)

- **Source File**:
  [`R/wgcna_modules.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R)
- **Signature**:
  `plot_wgcnaModules(x, main = "Gene dendrogram and module colors", ...)`
- **Description**: Plots the hierarchical clustering dendrogram with
  module color bars using
  [`WGCNA::plotDendroAndColors`](https://rdrr.io/pkg/WGCNA/man/plotDendroAndColors.html).

##### [`plot_null()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/plot_null.R)

- **Source File**:
  [`R/plot_null.R`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/plot_null.R)
- **Signature**: `plot_null(msg = "no data")`
- **Description**: Fallback void ggplot display showing centered message
  text when input data is empty.
