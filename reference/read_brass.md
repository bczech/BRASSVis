# Read BRASS BEDPE output

Parse a BEDPE file produced by the BRASS structural variant caller into
a data.frame suitable for visualization.

## Usage

``` r
read_brass(path)
```

## Arguments

- path:

  Path to a `.bedpe` file.

## Value

A data.frame with columns: gene1, gene2, strand1, strand2, direction1,
direction2, breakpoint1, breakpoint2, type, transcript_id1,
transcript_id2, fusion_flag, and others.

## Examples

``` r
bedpe <- system.file("extdata", "example.bedpe", package = "BRASSVis")
if (nzchar(bedpe)) fusions <- read_brass(bedpe)
```
