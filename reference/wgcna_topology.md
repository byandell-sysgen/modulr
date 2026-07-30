# WGCNA topology

WGCNA topology

GGplot of WGCNA Topology

## Usage

``` r
wgcna_topology(
  object,
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2)),
  cores = parallel::detectCores(),
  verbose = 0,
  ...
)

ggplot_wgcna_topology(object, cutoff = 0.9, ...)
```

## Arguments

- object:

  object of class \`wgcna_topology\`

- powers:

  vector of beta values

- cores:

  number of cores to enable

- verbose:

  level of verbose messages

- ...:

  additional parameters (ignored)

- cutoff:

  explained variation cutoff

## Value

data frame of class \`wgcna_topology\`

ggplot object
