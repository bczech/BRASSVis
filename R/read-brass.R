#' Read BRASS BEDPE output
#'
#' Parse a BEDPE file produced by the BRASS structural variant caller
#' into a data.frame suitable for visualization.
#'
#' @param path Path to a `.bedpe` file.
#' @return A data.frame with columns: gene1, gene2, strand1, strand2,
#'   direction1, direction2, breakpoint1, breakpoint2, type,
#'   transcript_id1, transcript_id2, fusion_flag, and others.
#' @export
#' @examples
#' bedpe <- system.file("extdata", "example.bedpe", package = "BRASSVis")
#' if (nzchar(bedpe)) fusions <- read_brass(bedpe)
read_brass <- function(path) {
  checkmate::assert_string(path)
  if (!file.exists(path)) stop("File not found: ", path)

  bras_input <- data.table::fread(path, sep = "\t", fill = TRUE,
                                  skip = "# chr1")

  data.frame(
    gene1 = bras_input$gene1,
    gene2 = bras_input$gene2,
    strand1 = paste0(bras_input$strand1, "/", bras_input$strand1),
    strand2 = paste0(bras_input$strand2, "/", bras_input$strand2),
    direction1 = ifelse(bras_input$strand1 == "+", "downstream", "upstream"),
    direction2 = ifelse(bras_input$strand2 == "-", "downstream", "upstream"),
    breakpoint1 = paste0(bras_input$`# chr1`, ":", bras_input$start1),
    breakpoint2 = paste0(bras_input$chr2, ":", bras_input$start2),
    site2 = "splice-site",
    type = bras_input$svclass,
    transcript_id1 = bras_input$transcript_id1,
    transcript_id2 = bras_input$transcript_id2,
    fusion_flag = bras_input$fusion_flag,
    total_region_count1 = bras_input$total_region_count1,
    total_region_count2 = bras_input$total_region_count2,
    assembly_score = bras_input$assembly_score,
    exonsSelected1 = bras_input$region_number1,
    exonsSelected2 = bras_input$region_number2,
    stringsAsFactors = FALSE
  )
}
