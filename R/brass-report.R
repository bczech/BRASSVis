#' Generate a full BRASS fusion report
#'
#' Read a BEDPE file, load annotations, and produce a multi-page PDF
#' with one page per fusion event showing exon structure, ideograms,
#' circos plots, and protein domain retention.
#'
#' @param bedpe_path Path to a BRASS `.bedpe` file.
#' @param gtf_path Path to a GTF annotation file.
#' @param output_pdf Path for the output PDF file.
#' @param cytobands_path Path to cytobands TSV. If NULL, uses built-in
#'   GRCh38 cytobands.
#' @param protein_domains_path Path to protein domains GFF3. If NULL,
#'   protein domains are not drawn.
#' @param fusion_flag_threshold Minimum fusion flag score to include
#'   a fusion (default 720).
#' @param transcript_selection One of `"provided"` or `"canonical"`.
#' @param color1,color2 Colors for the two genes in each fusion.
#' @param font_size Font size multiplier.
#' @param pdf_width,pdf_height PDF page dimensions in inches.
#' @param print_exon_labels Include exon numbers on the plot.
#' @return Invisible path to the output PDF.
#' @export
#' @examples
#' \dontrun{
#' brass_report(
#'   bedpe_path = "results.bedpe",
#'   gtf_path = "Homo_sapiens.GRCh38.gtf.gz",
#'   output_pdf = "fusions_report.pdf"
#' )
#' }
brass_report <- function(bedpe_path,
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
                         print_exon_labels = TRUE) {

  checkmate::assert_string(transcript_selection,
                           pattern = "provided|canonical")

  fusions <- read_brass(bedpe_path)
  ann <- load_annotations(gtf_path, cytobands_path, protein_domains_path,
                          print_exon_labels)

  fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1,
                                 perl = TRUE)
  fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2,
                                 perl = TRUE)
  fusions$contig1 <- gsub(":.*", "", fusions$breakpoint1)
  fusions$contig2 <- gsub(":.*", "", fusions$breakpoint2)
  fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1,
                                        perl = TRUE))
  fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2,
                                        perl = TRUE))
  fusions$site1 <- rep("exon", nrow(fusions))
  fusions$site2 <- rep("exon", nrow(fusions))
  fusions$confidence <- rep("high", nrow(fusions))

  pdf(output_pdf, onefile = TRUE, width = pdf_width, height = pdf_height,
      title = "BRASSVis Report")
  on.exit(dev.off(), add = TRUE)

  if (nrow(fusions) == 0) {
    plot(0, 0, type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(0, 0, "Error: empty input file\n")
    message("No fusions found in input file.")
    return(invisible(output_pdf))
  }

  exons <- ann$exons

  if (any(fusions$site1 == "intergenic" | fusions$site2 == "intergenic")) {
    inter <- rbind(
      setNames(fusions[fusions$site1 == "intergenic",
        c("gene1", "strand1", "contig1", "breakpoint1")],
        c("gene", "strand", "contig", "breakpoint")),
      setNames(fusions[fusions$site2 == "intergenic",
        c("gene2", "strand2", "contig2", "breakpoint2")],
        c("gene", "strand", "contig", "breakpoint"))
    )
    exons <- rbind(exons, data.frame(
      contig = inter$contig, type = "intergenic",
      start = inter$breakpoint - 1000, end = inter$breakpoint + 1000,
      strand = ".", attributes = "", geneName = inter$gene,
      geneID = inter$gene, transcript = inter$gene,
      exonNumber = "intergenic", stringsAsFactors = FALSE
    ))
  }

  fusions$tandem <- ifelse(fusions$type == "tandem-duplication", 2, 1)
  fusions <- fusions[rep(seq_len(nrow(fusions)), fusions$tandem), ]
  dup_idx <- which(duplicated(fusions))
  fusions[dup_idx, "direction1"] <- "downstream"
  fusions[dup_idx, "direction2"] <- "upstream"

  for (i in seq_len(nrow(fusions))) {
    message(paste0("Drawing fusion #", i, ": ",
                   fusions[i, "gene1"], ":", fusions[i, "gene2"]))

    if (fusions[i, "fusion_flag"] < fusion_flag_threshold) {
      message(sprintf("Fusion flag lower than %s. Skipping...",
                      fusion_flag_threshold))
      next
    }

    plot_fusion(fusions[i, ], fusions, i, exons, ann$cytobands,
      ann$protein_domains, color1, color2, font_size,
      transcript_selection)
  }

  message("Done")
  invisible(output_pdf)
}
