# BRASSVis

<!-- badges: start -->
[![R-CMD-check](https://github.com/bczech/BRASSVis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bczech/BRASSVis/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Visualize gene fusions and structural variants detected by
[BRASS](https://github.com/cancerit/BRASS) (BReakpoint AnalySiS) from
whole-genome sequencing data. Based on the
[Arriba](https://github.com/suhrig/arriba) visualization approach
adapted for WGS structural variant BEDPE output.

Each fusion is rendered as a multi-panel page showing:
- Chromosome ideograms with breakpoint positions
- Exon structure with intron squishing
- Circos plot for genomic context
- Retained protein domains in the fusion transcript
- Supporting read counts and assembly scores

## Installation

```r
# From GitHub
pak::pak("bczech/BRASSVis")

# Or with remotes
remotes::install_github("bczech/BRASSVis")
```

### Dependencies

BRASSVis requires Bioconductor packages `GenomicRanges`, `IRanges`, and
`S4Vectors`. These are installed automatically by `pak`. If installing
manually:

```r
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors"))
```

## Quick start

### From R

```r
library(BRASSVis)

# Generate a full PDF report
brass_report(
  bedpe_path = "results.bedpe",
  gtf_path = "Homo_sapiens.GRCh38.gtf.gz",
  output_pdf = "fusions_report.pdf"
)
```

### With protein domains and custom cytobands

```r
brass_report(
  bedpe_path = "results.bedpe",
  gtf_path = "Homo_sapiens.GRCh37.75.gtf.gz",
  output_pdf = "report.pdf",
  cytobands_path = system.file("ref", "cytobands_hg19_hs37d5_GRCh37_v2.1.0.tsv",
                               package = "BRASSVis"),
  protein_domains_path = system.file("ref", "protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3",
                                     package = "BRASSVis")
)
```

### Step by step

```r
# Read BEDPE
fusions <- read_brass("results.bedpe")

# Load annotations
ann <- load_annotations(
  gtf_path = "Homo_sapiens.GRCh38.gtf.gz",
  protein_domains_path = system.file("ref",
    "protein_domains_hg38_GRCh38_v2.1.0.gff3", package = "BRASSVis")
)

# Plot a single fusion
pdf("single_fusion.pdf", width = 11.7, height = 8.3)
plot_fusion(fusions[1, ], fusions, 1, ann$exons, ann$cytobands,
            ann$protein_domains)
dev.off()
```

### From the command line (legacy CLI)

```bash
Rscript inst/scripts/BRASSVis.R \
  -i results.bedpe \
  -e Homo_sapiens.GRCh38.gtf.gz \
  -o report.pdf \
  --proteinDomainsFile ref/protein_domains_hg38.gff3 \
  --cytobandsFile ref/cytobands_hg38.tsv
```

### With Docker

```bash
docker build -t brassvis .
docker run brassvis Rscript -e 'BRASSVis::brass_report(
  "input.bedpe", "annotation.gtf", "output.pdf")'
```

## Built-in reference data

The package includes cytobands and protein domain annotations for:

| Genome | Cytobands | Protein domains |
|---|---|---|
| Human GRCh38 (hg38) | Yes | Yes |
| Human GRCh37 (hg19) | Yes | Yes |
| Mouse GRCm38 (mm10) | Yes | Yes |

Access with `system.file("ref", "<filename>", package = "BRASSVis")`.

## Input format

BRASSVis reads annotated BEDPE files produced by BRASS. The file must
contain columns: `chr1`, `start1`, `chr2`, `start2`, `strand1`,
`strand2`, `svclass`, `gene1`, `gene2`, `fusion_flag`, and others.

## Example output

An example report generated from test data is available in
[`tst/report.pdf`](tst/report.pdf).

## Exported functions

| Function | Description |
|---|---|
| `read_brass()` | Parse BRASS BEDPE file |
| `load_annotations()` | Load GTF, cytobands, and protein domains |
| `plot_fusion()` | Visualize a single fusion event |
| `brass_report()` | Generate a complete multi-page PDF report |

## Authors

- **Bartosz Czech**
- **Pawel Sztromwasser**
- **Marzena Wojtaszewska**

## License

GPL-3
