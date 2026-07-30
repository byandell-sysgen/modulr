# modulr Data Pipeline & WGCNA Harmonization Methodology

## Data Pipeline & WGCNA Harmonization Methodology

This document details the underlying numerical data models, network
transformation algorithms, TOM distance formulas, dynamic tree cut
module merging, and S3 class structures in **modulr**.

------------------------------------------------------------------------

### 1. Long-to-Wide Matrix Pivot & Sample ID Standardization

Input datasets in the `foundr` ecosystem typically arrive in long-format
data frames containing sample metadata (`strain`, `sex`, `condition`,
`animal`), trait identifiers (`trait`), and numeric values (`value`).

WGCNA algorithms require wide numeric matrices where rows represent
individual sample observations and columns represent individual
phenotypic traits or genes.

#### Pivot Mechanism ([`wgcna_pivot()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_pivot.R))

[`wgcna_pivot()`](https://byandell-sysgen.github.io/modulr/reference/wgcna_pivot.md)
performs wide transformation:

``` r

# Extract metadata columns
IDcols <- c("dataset", "strain", "sex", "condition", "animal")

# Pivot traits into columns
wide_df <- tidyr::pivot_wider(
  tidyr::unite(object, ID, tidyr::any_of(IDcols)),
  names_from = "trait",
  values_from = "value"
)
```

#### Sample ID Standardization ([`wgcna_ID()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgnca_ID.R))

To ensure sample alignment,
[`wgcna_ID()`](https://byandell-sysgen.github.io/modulr/reference/wgcna_ID.md)
processes matrix rownames or input tables to construct unified `ID` keys
of the form `strain_sex_condition` alongside distinct `animal`
identifiers:

- **Strain Normalization**: Converts non-standard strain designations
  (e.g. `A129` or `X129`) to standard founder strain name `129`.
- **Delimiters**: Converts underscores to colons or hyphens as required
  during string parsing and separates `animal` numbers from composite
  strain strings.
- **External Annotations**: When an annotation table (`annot`) is
  provided, performs a `left_join` against `mouse_id` to map `strain`,
  `sex`, `diet`, and animal `number`.

------------------------------------------------------------------------

### 2. Topological Overlap Matrix (TOM) Dissimilarity & Soft Thresholding

#### Mathematical Formulation

Given a matrix $`X \in \mathbb{R}^{N \times P}`$ with $`N`$ samples and
$`P`$ traits:

1.  **Pearson Correlation Matrix**: Compute pairwise Pearson
    correlations $`S_{ij} = \text{cor}(x_i, x_j)`$.
2.  **Adjacency Matrix**: Apply soft-thresholding power $`\beta`$:
    ``` math
    A_{ij} = |S_{ij}|^\beta \quad (\text{for unsigned networks})
    ```
3.  **Topological Overlap Matrix (TOM)**:
    ``` math
    \text{TOM}_{ij} = \frac{l_{ij} + A_{ij}}{\min(k_i, k_j) + 1 - A_{ij}}
    ```
    where $`l_{ij} = \sum_u A_{iu} A_{uj}`$ measures shared
    connectivity, and $`k_i = \sum_u A_{iu}`$ is node degree.
4.  **TOM Dissimilarity Matrix**:
    ``` math
    \text{dissTOM}_{ij} = 1 - \text{TOM}_{ij}
    ```

#### Implementation ([`wgcna_dist()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_dist.R))

``` r

wgcna_dist <- function(object, params) {
  if(!is.matrix(object))
    object <- wgcna_pivot(object)$matrix

  params <- wgcna_params(params)

  1 - WGCNA::TOMsimilarity(
    WGCNA::adjacency(
      object,
      type = params$signType,
      power = params$power),
    TOMType = params$signType,
    verbose = params$verbose)
}
```

------------------------------------------------------------------------

### 3. Dynamic Tree Cutting, Module Merging, & kME Calculations

The core pipeline execution in
[`wgcnaModules()`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/R/wgcna_modules.R)
proceeds through four distinct stages:

``` mermaid
flowchart TD
    step1["1. TOM Dissimilarity Calculation (wgcna_dist)"]
    step2["2. Average Linkage Hierarchical Clustering (fastcluster::hclust)"]
    step3["3. Dynamic Branch Cut (dynamicTreeCut::cutreeDynamic)"]
    step4["4. Merge Similar Modules (WGCNA::mergeCloseModules)"]
    step5["5. Module Membership Correlation (WGCNA::signedKME)"]
    step6["6. Format Output wgcnaModules S3 Object"]

    step1 --> step2
    step2 --> step3
    step3 --> step4
    step4 --> step5
    step5 --> step6

    classDef proc fill:#1f77b4,stroke:#333,stroke-width:2px,color:#fff
    classDef merge fill:#ff7f0e,stroke:#333,stroke-width:2px,color:#fff
    classDef out fill:#2ca02c,stroke:#333,stroke-width:2px,color:#fff

    class step1,step2,step3 proc
    class step4,step5 merge
    class step6 out
```

#### Module Identification & Merging Rules

1.  **Initial Clustering**: Fast average-linkage hierarchical clustering
    is performed on $`\text{dissTOM}`$.
2.  **Adaptive Branch Cut**: Branches are dynamically cut based on shape
    parameters (`deepSplit = 2`, `minClusterSize = 4`).
3.  **Eigentrait Merging**: Module eigentraits (first principal
    component) are calculated for initial modules. Modules whose
    eigentraits exhibit correlation distance less than `thresholdMEDiss`
    (default `0.25`, corresponding to $`r \ge 0.75`$) are merged into a
    single module using
    [`WGCNA::mergeCloseModules()`](https://rdrr.io/pkg/WGCNA/man/mergeCloseModules.html).
4.  **kME Correlation**: Pairwise correlations between every trait
    $`x_j`$ and every module eigentrait $`E^{(m)}`$ are computed via
    [`WGCNA::signedKME()`](https://rdrr.io/pkg/WGCNA/man/signedKME.html).

------------------------------------------------------------------------

### 4. S3 Class Specifications & Downstream Integration

#### `wgcnaModules` Structure

``` r

out <- list(
  ID      = ID,      # Standardized sample metadata
  dendro  = dendro,  # Hierarchical clustering tree
  eigen   = eigen,   # Eigentraits data frame (rownames = sample IDs)
  modules = kME      # Long data frame of trait, module color factor, and kME
)
class(out) <- c("wgcnaModules", class(out))
attr(out, "params") <- params
```

#### `listof_wgcnaModules` Structure

``` r

out <- purrr::map(split(traitContr, traitContr$sex), wgcnaModules, params)
class(out) <- c("listof_wgcnaModules", class(out))
attr(out, "params") <- attr(out[[1]], "params")
```

The resulting eigentraits (`out$eigen`) can be directly consumed by
[`foundr::partition()`](https://rdrr.io/pkg/foundr/man/partition.html)
or
[`foundr::strainstats()`](https://rdrr.io/pkg/foundr/man/strainstats.html)
for downstream variance decomposition across founder strains, sex, and
diet conditions.
