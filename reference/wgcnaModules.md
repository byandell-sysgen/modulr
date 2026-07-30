# Create WGCNA Modules

Data in \`object\` are assumed to be in long format with columns a
subset of "dataset", "strain", "sex", "animal", "condition"

## Usage

``` r
wgcnaModules(object, params = NULL, ...)

listof_wgcnaModules(traitContr, params = NULL)

plot_wgcnaModules(x, main = "Gene dendrogram and module colors", ...)

# S3 method for class 'wgcnaModules'
plot(x, ...)
```

## Arguments

- object:

  harmonized data frame from routine \`foundr\`

- params:

  list of parameters for WGCNA routines

- ...:

  additional parameters

- traitContr:

  data frame of trait contributions

- x:

  object of class \`wgcnaModules\`

- main:

  title for dendrogram plot

## Value

object of class wgcnaModules

object of class \`listof_wgcnaModules\`
