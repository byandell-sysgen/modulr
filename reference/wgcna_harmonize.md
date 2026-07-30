# Harmonize WGCNA module

This loads \`geneTree\`, \`kMEs\`, \`merge\`, \`params\` and
\`dynamicColors_2\` from an \`Rdata\` object and creates and \`RDS\`
object for later reuse.

## Usage

``` r
wgcna_harmonize(
  mod = "WGCNAmodule",
  moddir = ".",
  modRdata = NULL,
  params = NULL,
  annot = NULL,
  harmonizeddir = ".",
  force = FALSE
)
```

## Arguments

- mod:

  name for module object

- moddir:

  director

- modRdata:

  file name ending with \`.Rdata\`

- params:

  list of parameters

- annot:

  annotation file (ignored if \`NULL\`)

- harmonizeddir:

  name of directory to save \`RDS\` object in its \`moddir\`

- force:

  force creation if \`TRUE\`

## Value

invisible file name for created object
