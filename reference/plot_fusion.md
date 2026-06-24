# Plot a single fusion event

Draw a detailed visualization of one structural variant / gene fusion,
including exon structure, ideograms, breakpoint labels, circos plot, and
protein domain retention.

## Usage

``` r
plot_fusion(
  fusion,
  fusions,
  fusion_idx,
  exons,
  cytobands,
  protein_domains = NULL,
  color1 = "#e5a5a5",
  color2 = "#a7c4e5",
  font_size = 1,
  transcript_selection = "provided"
)
```

## Arguments

- fusion:

  A single-row data.frame from
  [`read_brass()`](https://bczech.github.io/BRASSVis/reference/read_brass.md)
  (preprocessed with breakpoints parsed and directions set).

- fusions:

  The full fusions data.frame (needed for circos context).

- fusion_idx:

  Row index of this fusion in `fusions`.

- exons:

  Parsed exon annotation from
  [`load_annotations()`](https://bczech.github.io/BRASSVis/reference/load_annotations.md).

- cytobands:

  Parsed cytobands from
  [`load_annotations()`](https://bczech.github.io/BRASSVis/reference/load_annotations.md).

- protein_domains:

  Parsed protein domains (or NULL).

- color1, color2:

  Colors for the two genes.

- font_size:

  Font size multiplier.

- transcript_selection:

  One of `"provided"` or `"canonical"`.

## Value

Invisible NULL. Draws to the current graphics device.
