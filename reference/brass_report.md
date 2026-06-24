# Generate a full BRASS fusion report

Read a BEDPE file, load annotations, and produce a multi-page PDF with
one page per fusion event showing exon structure, ideograms, circos
plots, and protein domain retention.

## Usage

``` r
brass_report(
  bedpe_path,
  gtf_path,
  output_pdf = "brass_report.pdf",
  cytobands_path = NULL,
  protein_domains_path = NULL,
  fusion_flag_threshold = 720,
  transcript_selection = "provided",
  color1 = "#e5a5a5",
  color2 = "#a7c4e5",
  font_size = 1,
  pdf_width = 11.692,
  pdf_height = 8.267,
  print_exon_labels = TRUE
)
```

## Arguments

- bedpe_path:

  Path to a BRASS `.bedpe` file.

- gtf_path:

  Path to a GTF annotation file.

- output_pdf:

  Path for the output PDF file.

- cytobands_path:

  Path to cytobands TSV. If NULL, uses built-in GRCh38 cytobands.

- protein_domains_path:

  Path to protein domains GFF3. If NULL, protein domains are not drawn.

- fusion_flag_threshold:

  Minimum fusion flag score to include a fusion (default 720).

- transcript_selection:

  One of `"provided"` or `"canonical"`.

- color1, color2:

  Colors for the two genes in each fusion.

- font_size:

  Font size multiplier.

- pdf_width, pdf_height:

  PDF page dimensions in inches.

- print_exon_labels:

  Include exon numbers on the plot.

## Value

Invisible path to the output PDF.

## Examples

``` r
if (FALSE) { # \dontrun{
brass_report(
  bedpe_path = "results.bedpe",
  gtf_path = "Homo_sapiens.GRCh38.gtf.gz",
  output_pdf = "fusions_report.pdf"
)
} # }
```
