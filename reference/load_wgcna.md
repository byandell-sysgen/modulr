# Load WGCNA object

Load WGCNA object

## Usage

``` r
load_wgcna(
  moddir,
  modobj = "WGCNA_objects_ms10.Rdata",
  params = list(signType = "unsigned", power = 12, minSize = 4),
  listof = TRUE,
  annot = NULL
)
```

## Arguments

- moddir:

  directory name containing module object

- modobj:

  name of module object, ending in \`.Rdata\` or \`.RData\`

- params:

  non-default parameters for WGCNA (see \`wgcna_params\`)

- listof:

  embed object as \`listof_wgcnaModules\` if \`TRUE\`

- annot:

  annotation file (ignored if \`NULL\`)

## Value

list object that has all WGCNA components
