#' Load annotation files for fusion visualization
#'
#' Parse GTF, cytobands, and protein domain files needed by
#' [plot_fusion()] and [brass_report()].
#'
#' @param gtf_path Path to a GTF annotation file (`.gtf` or `.gtf.gz`).
#' @param cytobands_path Path to a cytobands TSV file. If NULL, uses
#'   the built-in GRCh38 cytobands from `inst/ref/`.
#' @param protein_domains_path Optional path to a protein domains GFF3
#'   file. If NULL, protein domains are not drawn.
#' @param print_exon_labels If TRUE, parse and include exon numbers.
#' @return A list with elements: `exons`, `cytobands`, `protein_domains`.
#' @export
#' @examples
#' \dontrun{
#' ann <- load_annotations("Homo_sapiens.GRCh38.gtf.gz")
#' }
load_annotations <- function(gtf_path,
                             cytobands_path = NULL,
                             protein_domains_path = NULL,
                             print_exon_labels = TRUE) {

  if (!file.exists(gtf_path)) stop("GTF file not found: ", gtf_path)

  message("Loading ideograms")
  if (is.null(cytobands_path)) {
    cytobands_path <- system.file("ref",
      "cytobands_hg38_GRCh38_v2.1.0.tsv", package = "BRASSVis")
  }
  cytobands <- read.table(cytobands_path, header = TRUE, sep = "\t")
  colnames(cytobands)[1] <- "contig"
  names(cytobands)[5] <- "giemsa"
  cytobands <- cytobands[order(cytobands$contig, cytobands$start,
                               cytobands$end), ]

  message("Loading annotation")
  exons <- read.table(gtf_path, header = FALSE, sep = "\t",
    comment.char = "#", quote = "", stringsAsFactors = FALSE)[,
    c(1, 3, 4, 5, 7, 9)]
  colnames(exons) <- c("contig", "type", "start", "end", "strand",
                        "attributes")
  exons <- exons[exons$type %in% c("exon", "CDS"), ]
  exons$contig <- .remove_chr(exons$contig)

  exons$geneID <- .parse_gtf_attribute("gene_id", exons)
  exons$geneName <- .parse_gtf_attribute("gene_name", exons)
  exons$geneName <- ifelse(exons$geneName == "", exons$geneID,
                           exons$geneName)
  exons$transcript <- .parse_gtf_attribute("transcript_id", exons)
  exons$exonNumber <- if (print_exon_labels) {
    .parse_gtf_attribute("exon_number", exons)
  } else {
    rep("", nrow(exons))
  }

  protein_domains <- NULL
  if (!is.null(protein_domains_path)) {
    message("Loading protein domains")
    protein_domains <- read.table(protein_domains_path, header = FALSE,
      sep = "\t", comment.char = "", quote = "",
      stringsAsFactors = FALSE)[, c(1, 4, 5, 7, 9)]
    colnames(protein_domains) <- c("contig", "start", "end", "strand",
                                   "attributes")
    protein_domains$color <- sub(";.*", "",
      sub(".*color=", "", protein_domains$attributes, perl = TRUE),
      perl = TRUE)
    protein_domains$proteinDomainName <- sapply(
      sub(";.*", "",
        sub(".*Name=", "", protein_domains$attributes, perl = TRUE),
        perl = TRUE),
      URLdecode)
    protein_domains$proteinDomainID <- sub(";.*", "",
      sub(".*protein_domain_id=", "", protein_domains$attributes,
          perl = TRUE),
      perl = TRUE)
    protein_domains <- protein_domains[,
      colnames(protein_domains) != "attributes"]
  }

  list(exons = exons, cytobands = cytobands,
       protein_domains = protein_domains)
}

.parse_gtf_attribute <- function(attribute, exons) {
  parsed <- gsub(paste0(".*", attribute, " \"?([^;\"]+)\"?;.*"),
                 "\\1", exons$attributes)
  failed <- parsed == exons$attributes
  if (any(failed)) {
    warning(paste0("Failed to parse '", attribute, "' attribute of ",
                   sum(failed), " GTF record(s)."))
    parsed <- ifelse(failed, "", parsed)
  }
  parsed
}
