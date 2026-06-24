.draw_vertical_gradient <- function(left, right, y, color, selection = NULL) {
  if (!is.null(selection)) {
    y <- y[selection]
    left <- left[selection]
    right <- right[selection]
  }
  for (i in seq_len(length(y))) {
    polygon(
      c(left[seq_len(i)], right[seq_len(i)]),
      c(y[seq_len(i)], y[seq_len(i)]),
      border = NA,
      col = rgb(col2rgb(color)["red", ], col2rgb(color)["green", ],
                col2rgb(color)["blue", ],
                col2rgb(color, alpha = TRUE)["alpha", ] * (1 / length(y)),
                max = 255)
    )
  }
}

.draw_curly_brace <- function(left, right, top, bottom, tip) {
  smoothness <- 20
  x <- cumsum(exp(-seq(-2.5, 2.5, len = smoothness)^2))
  x <- x / max(x)
  y <- seq(top, bottom, len = smoothness)
  lines(left + (tip - left) + x * (left - tip), y)
  lines(tip + x * (right - tip), y)
}

.draw_ideogram <- function(adjust, left, right, y, cytobands, contig,
                           breakpoint, font_size = 1) {
  band_colors <- setNames(
    rgb(100:0, 100:0, 100:0, maxColorValue = 100),
    paste0("gpos", 0:100)
  )
  band_colors <- c(band_colors, gneg = "#ffffff", acen = "#ec4f4f",
                   stalk = "#0000ff")
  arc_steps <- 30
  curly_brace_height <- 0.03
  ideogram_height <- 0.04
  ideogram_width <- 0.4

  bands <- cytobands[cytobands$contig == contig, ]
  if (nrow(bands) == 0) {
    warning(paste("Ideogram of contig", contig,
                  "cannot be drawn, no Giemsa staining info available."))
    return(invisible(NULL))
  }

  bands$color <- band_colors[bands$giemsa]
  bands$left <- bands$start / max(cytobands$end) * ideogram_width
  bands$right <- bands$end / max(cytobands$end) * ideogram_width
  offset <- ifelse(adjust == "left", left,
                   right - max(bands$right))
  bands$left <- bands$left + offset
  bands$right <- bands$right + offset

  tip <- min(bands$left) + (max(bands$right) - min(bands$left)) /
    (max(bands$end) - min(bands$start)) * breakpoint
  .draw_curly_brace(left, right, y - 0.05 + curly_brace_height,
                    y - 0.05, tip)

  text((max(bands$right) + min(bands$left)) / 2, y + 0.07,
       paste("chromosome", contig), font = 2, cex = font_size,
       adj = c(0.5, 0))

  band_name <- bands[which(.between(breakpoint, bands$start, bands$end)),
                     "name"]
  tryCatch(
    text(tip, y + 0.03, band_name, cex = font_size, adj = c(0.5, 0)),
    error = function(e) NULL
  )

  left_arc_x <- bands[1, "left"] +
    (1 + cos(seq(pi / 2, 1.5 * pi, len = arc_steps))) *
    (bands[1, "right"] - bands[1, "left"])
  left_arc_y <- y + sin(seq(pi / 2, 1.5 * pi, len = arc_steps)) *
    (ideogram_height / 2)
  polygon(left_arc_x, left_arc_y, col = bands[1, "color"])

  centromere_start <- NULL
  centromere_end <- NULL
  for (band in 2:(nrow(bands) - 1)) {
    if (bands[band, "giemsa"] != "acen") {
      rect(bands[band, "left"], y - ideogram_height / 2,
           bands[band, "right"], y + ideogram_height / 2,
           col = bands[band, "color"])
    } else {
      if (is.null(centromere_start)) {
        polygon(c(bands[band, "left"], bands[band, "right"],
                  bands[band, "left"]),
                c(y - ideogram_height / 2, y, y + ideogram_height / 2),
                col = bands[band, "color"])
        centromere_start <- bands[band, "left"]
      } else {
        polygon(c(bands[band, "right"], bands[band, "left"],
                  bands[band, "right"]),
                c(y - ideogram_height / 2, y, y + ideogram_height / 2),
                col = bands[band, "color"])
        centromere_end <- bands[band, "right"]
      }
    }
  }

  band <- nrow(bands)
  right_arc_x <- bands[band, "right"] -
    (1 + cos(seq(1.5 * pi, pi / 2, len = arc_steps))) *
    (bands[band, "right"] - bands[band, "left"])
  right_arc_y <- y + sin(seq(pi / 2, 1.5 * pi, len = arc_steps)) *
    ideogram_height / 2
  polygon(right_arc_x, right_arc_y, col = bands[band, "color"])

  if (is.null(centromere_start) || is.null(centromere_end)) {
    centromere_start <- bands[1, "right"]
    centromere_end <- bands[1, "right"]
  }

  .draw_vertical_gradient(left_arc_x, rep(centromere_start, arc_steps),
    left_arc_y, rgb(0, 0, 0, 0.8), seq_len(round(arc_steps * 0.4)))
  .draw_vertical_gradient(left_arc_x, rep(centromere_start, arc_steps),
    left_arc_y, rgb(1, 1, 1, 0.7),
    round(arc_steps * 0.4):round(arc_steps * 0.1))
  .draw_vertical_gradient(left_arc_x, rep(centromere_start, arc_steps),
    left_arc_y, rgb(1, 1, 1, 0.7),
    round(arc_steps * 0.4):round(arc_steps * 0.6))
  .draw_vertical_gradient(left_arc_x, rep(centromere_start, arc_steps),
    left_arc_y, rgb(0, 0, 0, 0.9),
    arc_steps:round(arc_steps * 0.5))
  .draw_vertical_gradient(right_arc_x, rep(centromere_end, arc_steps),
    right_arc_y, rgb(0, 0, 0, 0.8), seq_len(round(arc_steps * 0.4)))
  .draw_vertical_gradient(right_arc_x, rep(centromere_end, arc_steps),
    right_arc_y, rgb(1, 1, 1, 0.7),
    round(arc_steps * 0.4):round(arc_steps * 0.1))
  .draw_vertical_gradient(right_arc_x, rep(centromere_end, arc_steps),
    right_arc_y, rgb(1, 1, 1, 0.7),
    round(arc_steps * 0.4):round(arc_steps * 0.6))
  .draw_vertical_gradient(right_arc_x, rep(centromere_end, arc_steps),
    right_arc_y, rgb(0, 0, 0, 0.9),
    arc_steps:round(arc_steps * 0.5))
}

.draw_strand <- function(left, right, y, color, strand) {
  if (strand %in% c("+", "-")) {
    lines(c(left + 0.001, right - 0.001), c(y, y), col = color, lwd = 2)
    lines(c(left + 0.001, right - 0.001), c(y, y),
          col = rgb(1, 1, 1, 0.1), lwd = 1)
    if (right - left > 0.01) {
      for (i in seq(left + 0.005, right - 0.005,
                    by = sign(right - left - 2 * 0.005) * 0.01)) {
        arrows(i, y, i + 0.001 * ifelse(strand == "+", 1, -1), y,
               col = color, length = 0.05, lwd = 2, angle = 60)
        arrows(i, y, i + 0.001 * ifelse(strand == "+", 1, -1), y,
               col = rgb(1, 1, 1, 0.1), length = 0.05, lwd = 1, angle = 60)
      }
    }
  }
}

.draw_exon <- function(left, right, y, color, title, type, font_size = 1) {
  gradient_steps <- 10
  exon_height <- 0.03
  dark_color <- .get_dark_color(color)

  if (type == "CDS") {
    rect(left, y + exon_height, right, y + exon_height / 2 - 0.001,
         col = color, border = NA)
    rect(left, y - exon_height, right, y - exon_height / 2 + 0.001,
         col = color, border = NA)
    lines(c(left, left, right, right),
          c(y + exon_height / 2, y + exon_height, y + exon_height,
            y + exon_height / 2), col = dark_color, lend = 2)
    lines(c(left, left, right, right),
          c(y - exon_height / 2, y - exon_height, y - exon_height,
            y - exon_height / 2), col = dark_color, lend = 2)
    .draw_vertical_gradient(rep(left, gradient_steps),
      rep(right, gradient_steps),
      seq(y + 0.03, y + 0.015, len = gradient_steps),
      rgb(0, 0, 0, 0.2))
    .draw_vertical_gradient(rep(left, gradient_steps),
      rep(right, gradient_steps),
      seq(y - 0.03, y - 0.015, len = gradient_steps),
      rgb(0, 0, 0, 0.3))
  } else if (type == "exon") {
    rect(left, y + exon_height / 2, right, y - exon_height / 2,
         col = color, border = dark_color)
    .draw_vertical_gradient(rep(left, gradient_steps),
      rep(right, gradient_steps),
      seq(y, y + exon_height / 2, len = gradient_steps),
      rgb(1, 1, 1, 0.6))
    .draw_vertical_gradient(rep(left, gradient_steps),
      rep(right, gradient_steps),
      seq(y, y - exon_height / 2, len = gradient_steps),
      rgb(1, 1, 1, 0.6))
    text((left + right) / 2, y, title, cex = 0.9 * font_size)
  }
}

.draw_circos <- function(fusion, fusions, cytobands, font_size = 1) {
  for (contig in unlist(fusions[fusion, c("contig1", "contig2")])) {
    if (!any(cytobands$contig == contig)) {
      warning(paste0("Circos plot cannot be drawn, no Giemsa staining ",
                     "info for contig ", contig, "."))
      return(invisible(NULL))
    }
  }

  circlize::circos.clear()
  circlize::circos.initializeWithIdeogram(cytoband = cytobands,
    labels.cex = 1.1, axis.labels.cex = 0.6)

  gene_labels <- data.frame(
    contig = c(fusions[fusion, "contig1"], fusions[fusion, "contig2"]),
    start = c(fusions[fusion, "breakpoint1"],
              fusions[fusion, "breakpoint2"])
  )
  gene_labels$end <- gene_labels$start + 1
  gene_labels$gene <- c(fusions[fusion, "gene1"],
                        fusions[fusion, "gene2"])
  gene_labels$gene <- ifelse(
    c(fusions[fusion, "site1"], fusions[fusion, "site2"]) == "intergenic",
    paste0(c(fusions[fusion, "display_contig1"],
             fusions[fusion, "display_contig2"]), ":",
           gene_labels$start),
    gene_labels$gene
  )

  circlize::circos.genomicLabels(gene_labels, labels.column = 4,
    side = "inside", cex = font_size - 0.1)

  for (contig in unique(cytobands$contig)) {
    circlize::set.current.cell(track.index = 2, sector.index = contig)
    circlize::circos.text(circlize::CELL_META$xcenter,
      circlize::CELL_META$ycenter, contig, cex = 0.25)
  }
}
