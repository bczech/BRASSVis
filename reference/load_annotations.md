# Load annotation files for fusion visualization

Parse GTF, cytobands, and protein domain files needed by
[`plot_fusion()`](https://bczech.github.io/BRASSVis/reference/plot_fusion.md)
and
[`brass_report()`](https://bczech.github.io/BRASSVis/reference/brass_report.md).

## Usage

``` r
load_annotations(
  gtf_path,
  cytobands_path = NULL,
  protein_domains_path = NULL,
  print_exon_labels = TRUE
)
```

## Arguments

- gtf_path:

  Path to a GTF annotation file (`.gtf` or `.gtf.gz`).

- cytobands_path:

  Path to a cytobands TSV file. If NULL, uses the built-in GRCh38
  cytobands from `inst/ref/`.

- protein_domains_path:

  Optional path to a protein domains GFF3 file. If NULL, protein domains
  are not drawn.

- print_exon_labels:

  If TRUE, parse and include exon numbers.

## Value

A list with elements: `exons`, `cytobands`, `protein_domains`.

## Examples

``` r
if (FALSE) { # \dontrun{
ann <- load_annotations("Homo_sapiens.GRCh38.gtf.gz")
} # }
```
