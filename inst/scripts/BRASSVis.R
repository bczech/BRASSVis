#!/usr/bin/env Rscript
suppressMessages(library("docopt"))

"usage: BRASSVis.R [options]

options:
-i --inputFile=<file>       BRASS BEDPE file path.
-e --exonsFile=<file>       GTF annotation file.
-o --outputFile=<file>      Output PDF file name [default: ./output.pdf].
--fusionFlagThreshold=<n>   Fusion flag threshold [default: 720].
--transcriptSelection=<s>   Transcript selection [default: provided].
--pdfWidth=<n>              PDF width [default: 11.692].
--pdfHeight=<n>             PDF height [default: 8.267].
--cytobandsFile=<file>      Cytobands file.
--proteinDomainsFile=<file> Protein domains file.
--color1=<s>                Color 1 [default: #e5a5a5].
--color2=<s>                Color 2 [default: #a7c4e5].
--printExonLabels=<l>       Print exon labels [default: TRUE].
--fontSize=<n>              Font size [default: 1]." -> doc

opt <- docopt(doc)

library(BRASSVis)

brass_report(
  bedpe_path = opt$inputFile,
  gtf_path = opt$exonsFile,
  output_pdf = opt$outputFile,
  cytobands_path = opt$cytobandsFile,
  protein_domains_path = opt$proteinDomainsFile,
  fusion_flag_threshold = as.numeric(opt$fusionFlagThreshold),
  transcript_selection = opt$transcriptSelection,
  color1 = opt$color1,
  color2 = opt$color2,
  font_size = as.numeric(opt$fontSize),
  pdf_width = as.numeric(opt$pdfWidth),
  pdf_height = as.numeric(opt$pdfHeight),
  print_exon_labels = as.logical(opt$printExonLabels)
)
