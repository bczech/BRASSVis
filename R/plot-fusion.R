#' Plot a single fusion event
#'
#' Draw a detailed visualization of one structural variant / gene fusion,
#' including exon structure, ideograms, breakpoint labels, circos plot,
#' and protein domain retention.
#'
#' @param fusion A single-row data.frame from [read_brass()] (preprocessed
#'   with breakpoints parsed and directions set).
#' @param fusions The full fusions data.frame (needed for circos context).
#' @param fusion_idx Row index of this fusion in `fusions`.
#' @param exons Parsed exon annotation from [load_annotations()].
#' @param cytobands Parsed cytobands from [load_annotations()].
#' @param protein_domains Parsed protein domains (or NULL).
#' @param color1,color2 Colors for the two genes.
#' @param font_size Font size multiplier.
#' @param transcript_selection One of `"provided"` or `"canonical"`.
#' @return Invisible NULL. Draws to the current graphics device.
#' @export
plot_fusion <- function(fusion, fusions, fusion_idx, exons, cytobands,
                        protein_domains = NULL,
                        color1 = "#e5a5a5", color2 = "#a7c4e5",
                        font_size = 1,
                        transcript_selection = "provided") {

  dark1 <- .get_dark_color(color1)
  dark2 <- .get_dark_color(color2)

  exons1 <- .find_exons(exons, fusion$contig1, fusion$gene1,
    fusion$direction1, fusion$breakpoint1, NULL,
    fusion$transcript_id1, transcript_selection)
  exons2 <- .find_exons(exons, fusion$contig2, fusion$gene2,
    fusion$direction2, fusion$breakpoint2, NULL,
    fusion$transcript_id2, transcript_selection)
  if (nrow(exons1) == 0 || nrow(exons2) == 0) return(invisible(NULL))
  if (length(unique(exons1$geneID)) > 0 &&
      length(unique(exons2$geneID)) > 0 &&
      unique(exons1$geneID)[1] == unique(exons2$geneID)[1])
    return(invisible(NULL))

  exons1 <- exons1[order(exons1$start, -rank(exons1$type)), ]
  exons2 <- exons2[order(exons2$start, -rank(exons2$type)), ]

  bp1 <- fusion$breakpoint1
  bp2 <- fusion$breakpoint2

  if (bp1 < min(exons1$start)) {
    exons1 <- rbind(
      data.frame(contig = exons1[1, "contig"], type = "dummy",
        start = bp1 - 1000, end = bp1 - 1000,
        strand = exons1[1, "strand"], attributes = "",
        geneID = "dummy", geneName = exons1[1, "geneName"],
        transcript = exons1[1, "transcript"], exonNumber = "",
        stringsAsFactors = FALSE),
      exons1)
  } else if (bp1 > max(exons1$end)) {
    exons1 <- rbind(exons1,
      data.frame(contig = exons1[1, "contig"], type = "dummy",
        start = bp1 + 1000, end = bp1 + 1000,
        strand = exons1[1, "strand"], attributes = "",
        geneID = "dummy", geneName = exons1[1, "geneName"],
        transcript = exons1[1, "transcript"], exonNumber = "",
        stringsAsFactors = FALSE))
  }
  if (bp2 < min(exons2$start)) {
    exons2 <- rbind(
      data.frame(contig = exons2[1, "contig"], type = "dummy",
        start = bp2 - 1000, end = bp2 - 1000,
        strand = exons2[1, "strand"], attributes = "",
        geneID = "dummy", geneName = exons2[1, "geneName"],
        transcript = exons2[1, "transcript"], exonNumber = "",
        stringsAsFactors = FALSE),
      exons2)
  } else if (bp2 > max(exons2$end)) {
    exons2 <- rbind(exons2,
      data.frame(contig = exons2[1, "contig"], type = "dummy",
        start = bp2 + 1000, end = bp2 + 1000,
        strand = exons2[1, "strand"], attributes = "",
        geneID = "dummy", geneName = exons2[1, "geneName"],
        transcript = exons2[1, "transcript"], exonNumber = "",
        stringsAsFactors = FALSE))
  }

  exons1$start <- as.integer(exons1$start)
  exons1$end <- as.integer(exons1$end)
  exons2$start <- as.integer(exons2$start)
  exons2$end <- as.integer(exons2$end)
  exons1$left <- exons1$start
  exons1$right <- exons1$end
  exons2$left <- exons2$start
  exons2$right <- exons2$end

  squished <- 200
  cum <- 0; prev <- -squished
  for (ex in seq_len(nrow(exons1))) {
    if (bp1 > prev + 1 && bp1 < exons1[ex, "left"])
      bp1 <- (bp1 - prev) / (exons1[ex, "left"] - prev) * squished +
        prev - cum
    if (exons1[ex, "left"] > prev) {
      cum <- cum + exons1[ex, "left"] - prev - squished
      prev <- exons1[ex, "right"]
    }
    if (bp1 >= exons1[ex, "left"] && bp1 <= exons1[ex, "right"] + 1)
      bp1 <- bp1 - cum
    exons1[ex, "left"] <- exons1[ex, "left"] - cum
    exons1[ex, "right"] <- exons1[ex, "right"] - cum
  }

  cum <- 0; prev <- -squished
  for (ex in seq_len(nrow(exons2))) {
    if (bp2 > prev + 1 && bp2 < exons2[ex, "left"])
      bp2 <- (bp2 - prev) / (exons2[ex, "left"] - prev) * squished +
        prev - cum
    if (exons2[ex, "left"] > prev) {
      cum <- cum + exons2[ex, "left"] - prev - squished
      prev <- exons2[ex, "right"]
    }
    if (bp2 >= exons2[ex, "left"] && bp2 <= exons2[ex, "right"] + 1)
      bp2 <- bp2 - cum
    exons2[ex, "left"] <- exons2[ex, "left"] - cum
    exons2[ex, "right"] <- exons2[ex, "right"] - cum
  }

  sf <- max(exons1$right) + max(exons2$right)
  exons1$left <- exons1$left / sf
  exons1$right <- exons1$right / sf
  exons2$left <- exons2$left / sf
  exons2$right <- exons2$right / sf
  bp1 <- bp1 / sf
  bp2 <- bp2 / sf

  g2_off <- max(exons1$right) + 0.05
  f_off1 <- (max(exons1$right) + g2_off) / 2 -
    ifelse(fusion$direction1 == "downstream", bp1,
           max(exons1$right) - bp1)
  f_off2 <- f_off1 + ifelse(fusion$direction1 == "downstream", bp1,
                             max(exons1$right) - bp1)

  layout(matrix(c(1, 1, 1, 2, 3, 4), 2, 3, byrow = TRUE),
         widths = c(0.9, 1.2, 0.9))
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, type = "l", xlim = c(-0.12, 1.12), ylim = c(0.4, 1.1),
       bty = "n", xaxt = "n", yaxt = "n")

  y_ideo <- 0.84; y_bp_lab <- 0.76; y_exons <- 0.67
  y_gene <- 0.58; y_fus <- 0.5; y_scale <- 0.407
  y_traj_bp <- y_bp_lab - 0.035
  y_traj_ex_top <- y_exons + 0.03; y_traj_ex_bot <- y_exons - 0.055
  y_traj_fus <- y_fus + 0.03

  if (!is.null(cytobands)) {
    .draw_ideogram("left", min(exons1$left), max(exons1$right), y_ideo,
      cytobands, fusion$contig1, fusion$breakpoint1, font_size)
    .draw_ideogram("right", g2_off, g2_off + max(exons2$right), y_ideo,
      cytobands, fusion$contig2, fusion$breakpoint2, font_size)
  }

  text(max(exons1$right) / 2, y_gene, fusion$gene1, font = 2,
       cex = font_size, adj = c(0.5, 0))
  if (fusion$site1 != "intergenic")
    text(max(exons1$right) / 2, y_gene - 0.01, head(exons1$transcript, 1),
         cex = 0.9 * font_size, adj = c(0.5, 1))
  text(g2_off + max(exons2$right) / 2, y_gene, fusion$gene2, font = 2,
       cex = font_size, adj = c(0.5, 0))
  if (fusion$site2 != "intergenic")
    text(g2_off + max(exons2$right) / 2, y_gene - 0.01,
         head(exons2$transcript, 1), cex = 0.9 * font_size,
         adj = c(0.5, 1))

  text(bp1 + 0.01, y_bp_lab - 0.03,
    paste0("breakpoint\n", fusion$display_contig1, ":",
           fusion$breakpoint1),
    adj = c(1, 0), cex = font_size)
  text(g2_off + bp2 - 0.01, y_bp_lab - 0.03,
    paste0("breakpoint\n", fusion$display_contig2, ":",
           fusion$breakpoint2),
    adj = c(0, 0), cex = font_size)

  lines(c(min(exons1$left), max(exons1$right)), c(y_exons, y_exons),
        col = dark1)
  for (g in unique(exons1$geneName))
    .draw_strand(min(exons1[exons1$geneName == g, "left"]),
      max(exons1[exons1$geneName == g, "right"]), y_exons, dark1,
      head(exons1[exons1$geneName == g, "strand"], 1))
  for (ex in seq_len(nrow(exons1)))
    .draw_exon(exons1[ex, "left"], exons1[ex, "right"], y_exons, color1,
      exons1[ex, "exonNumber"], exons1[ex, "type"], font_size)

  lines(c(g2_off, g2_off + max(exons2$right)), c(y_exons, y_exons),
        col = dark2)
  for (g in unique(exons2$geneName))
    .draw_strand(g2_off + min(exons2[exons2$geneName == g, "left"]),
      g2_off + max(exons2[exons2$geneName == g, "right"]), y_exons, dark2,
      head(exons2[exons2$geneName == g, "strand"], 1))
  for (ex in seq_len(nrow(exons2)))
    .draw_exon(g2_off + exons2[ex, "left"], g2_off + exons2[ex, "right"],
      y_exons, color2, exons2[ex, "exonNumber"], exons2[ex, "type"],
      font_size)

  real_scale <- max(exons1$end - exons1$start, exons2$end - exons2$start)
  map_scale <- max(exons1$right - exons1$left,
                   exons2$right - exons2$left)
  desired <- 0.2
  real_scale <- desired / map_scale * real_scale
  map_scale <- desired
  real_opt <- signif(real_scale, 1)
  map_opt <- real_opt / real_scale * map_scale
  lines(c(0, map_opt), c(y_scale, y_scale))
  lines(c(0, 0), c(y_scale - 0.007, y_scale + 0.007))
  lines(c(map_opt, map_opt), c(y_scale - 0.007, y_scale + 0.007))
  thousands <- max(0, min(3, floor(log10(real_opt) / 3)))
  units <- c("bp", "kbp", "Mbp", "Gbp")
  text(map_opt / 2, y_scale + 0.005,
    paste(real_opt / max(1, 1000^thousands), units[thousands + 1]),
    adj = c(0.5, 0), cex = font_size * 0.9)
  text(map_opt, y_scale, "  introns not to scale",
       adj = c(0, 0.5), cex = font_size * 0.9, font = 3)

  par(mar = c(0, 0, 0, 0))
  .draw_circos(fusion_idx, fusions, cytobands, font_size)
  par(mar = c(0, 0, 0, 0))

  plot(0, 0, type = "l", xlim = c(-0.1, 1.1), ylim = c(0, 1),
       bty = "n", xaxt = "n", yaxt = "n")
  par(xpd = NA)
  if (!is.null(protein_domains))
    .draw_protein_domains(fusion, exons1, exons2, protein_domains,
                          color1, color2, font_size)
  par(xpd = FALSE)

  plot(0, 0, type = "l", xlim = c(0, 1), ylim = c(0, 1),
       bty = "n", xaxt = "n", yaxt = "n")
  text(0, 0.575, "SUPPORTING READ COUNT", font = 2, adj = c(0, 0),
       cex = font_size)
  text(0, 0.525, paste0(
    "Total region count in ",
    unique(exons1$geneName[exons1$type != "dummy"]),
    " = ", fusion$total_region_count1, "\n",
    "Total region count in ",
    unique(exons2$geneName[exons2$type != "dummy"]),
    " = ", fusion$total_region_count2, "\n",
    "Assembly score = ", fusion$assembly_score, "\n",
    "Fusion flag = ", fusion$fusion_flag),
    adj = c(0, 1), cex = font_size)

  invisible(NULL)
}
