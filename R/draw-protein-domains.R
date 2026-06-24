.draw_protein_domains <- function(fusion, exons1, exons2,
                                  protein_domains, color1, color2,
                                  font_size = 1) {
  exon_height <- 0.2
  exons_y <- 0.5
  gene_names_y <- exons_y - exon_height / 2 - 0.05

  coding1 <- exons1[exons1$type == "CDS" & fusion$site1 != "intergenic", ]
  coding2 <- exons2[exons2$type == "CDS" & fusion$site2 != "intergenic", ]

  if (fusion$direction1 == "upstream") {
    coding1 <- coding1[coding1$end >= fusion$breakpoint1, ]
    coding1$start <- ifelse(coding1$start < fusion$breakpoint1,
                            fusion$breakpoint1, coding1$start)
  } else {
    coding1 <- coding1[coding1$start <= fusion$breakpoint1, ]
    coding1$end <- ifelse(coding1$end > fusion$breakpoint1,
                          fusion$breakpoint1, coding1$end)
  }
  if (fusion$direction2 == "upstream") {
    coding2 <- coding2[coding2$end >= fusion$breakpoint2, ]
    coding2$start <- ifelse(coding2$start < fusion$breakpoint2,
                            fusion$breakpoint2, coding2$start)
  } else {
    coding2 <- coding2[coding2$start <= fusion$breakpoint2, ]
    coding2$end <- ifelse(coding2$end > fusion$breakpoint2,
                          fusion$breakpoint2, coding2$end)
  }

  gr1 <- GenomicRanges::GRanges(coding1$contig,
    IRanges::IRanges(coding1$start, coding1$end), strand = coding1$strand)
  gr2 <- GenomicRanges::GRanges(coding2$contig,
    IRanges::IRanges(coding2$start, coding2$end), strand = coding2$strand)
  dom_gr <- GenomicRanges::GRanges(protein_domains$contig,
    IRanges::IRanges(protein_domains$start, protein_domains$end),
    strand = protein_domains$strand)
  dom_gr$proteinDomainName <- protein_domains$proteinDomainName
  dom_gr$proteinDomainID <- protein_domains$proteinDomainID
  dom_gr$color <- protein_domains$color
  hits <- GenomicRanges::findOverlaps(dom_gr,
    GenomicRanges::union(gr1, gr2))
  dom_gr <- dom_gr[unique(S4Vectors::queryHits(hits))]

  dom_list <- GenomicRanges::GRangesList(lapply(
    unique(dom_gr$proteinDomainID),
    function(x) dom_gr[dom_gr$proteinDomainID == x]))

  trim_fn <- function(dl, egr) {
    do.call("rbind", lapply(dl, function(x) {
      inter <- as.data.frame(GenomicRanges::reduce(
        suppressWarnings(GenomicRanges::intersect(x, egr))))
      if (nrow(inter) > 0) {
        inter$proteinDomainName <- utils::head(x$proteinDomainName, 1)
        inter$proteinDomainID <- utils::head(x$proteinDomainID, 1)
        inter$color <- utils::head(x$color, 1)
      } else {
        inter$proteinDomainName <- character()
        inter$proteinDomainID <- character()
        inter$color <- character()
      }
      inter
    }))
  }
  retained1 <- trim_fn(dom_list, gr1)
  retained2 <- trim_fn(dom_list, gr2)

  coding1$length <- coding1$end - coding1$start + 1
  coding2$length <- coding2$end - coding2$start + 1

  if (sum(exons1$type == "CDS") + sum(exons2$type == "CDS") == 0) {
    text(0.5, 0.5, "Genes are not protein-coding.")
    return(invisible(NULL))
  }
  len1 <- sum(coding1$length)
  len2 <- sum(coding2$length)
  if (len1 + len2 == 0) {
    text(0.5, 0.5, "No coding regions retained in fusion transcript.")
    return(invisible(NULL))
  }
  anti1 <- sub("/.*", "", fusion$strand1) != sub(".*/", "", fusion$strand1)
  anti2 <- sub("/.*", "", fusion$strand2) != sub(".*/", "", fusion$strand2)
  if ((len1 == 0 || anti1) && (len2 == 0 || anti2)) {
    text(0.5, 0.5, "No coding regions due to antisense transcription.")
    return(invisible(NULL))
  }

  remove_introns <- function(coding, retained) {
    if (nrow(coding) == 0) return(NULL)
    cum_intron <- 0
    prev_end <- 0
    for (ex in seq_len(nrow(coding))) {
      if (coding[ex, "start"] > prev_end)
        cum_intron <- cum_intron + coding[ex, "start"] - prev_end
      idx <- which(.between(retained$start,
        coding[ex, "start"], coding[ex, "end"]))
      retained[idx, "start"] <- retained[idx, "start"] - cum_intron
      idx2 <- which(.between(retained$end,
        coding[ex, "start"], coding[ex, "end"]))
      retained[idx2, "end"] <- retained[idx2, "end"] - cum_intron
      prev_end <- coding[ex, "end"]
    }
    do.call("rbind", lapply(unique(retained$proteinDomainID), function(x) {
      d <- retained[retained$proteinDomainID == x, ]
      m <- GenomicRanges::reduce(GenomicRanges::GRanges(
        d$seqnames, IRanges::IRanges(d$start, d$end), strand = d$strand))
      m$proteinDomainName <- utils::head(d$proteinDomainName, 1)
      m$proteinDomainID <- utils::head(d$proteinDomainID, 1)
      m$color <- utils::head(d$color, 1)
      as.data.frame(m)
    }))
  }
  retained1 <- remove_introns(coding1, retained1)
  retained2 <- remove_introns(coding2, retained2)

  if (is.null(retained1) && is.null(retained2)) {
    text(0.5, 0.5, "No protein domains retained in fusion.")
    return(invisible(NULL))
  }

  merge_similar <- function(domains) {
    if (is.null(domains)) return(domains)
    merged <- domains[FALSE, ]
    domains <- domains[order(domains$end - domains$start,
                             decreasing = FALSE), ]
    for (i in seq_len(nrow(domains))) {
      overlap <- (abs(merged$start - domains[i, "start"]) +
        abs(merged$end - domains[i, "end"])) /
        (domains[i, "end"] - domains[i, "start"])
      if (!any(overlap <= 0.1))
        merged <- rbind(merged, domains[i, ])
    }
    merged
  }
  retained1 <- merge_similar(retained1)
  retained2 <- merge_similar(retained2)

  coding1$length <- coding1$length / (len1 + len2)
  coding2$length <- coding2$length / (len1 + len2)
  if (!is.null(retained1)) {
    retained1$start <- retained1$start / (len1 + len2)
    retained1$end <- retained1$end / (len1 + len2)
  }
  if (!is.null(retained2)) {
    retained2$start <- retained2$start / (len1 + len2)
    retained2$end <- retained2$end / (len1 + len2)
  }

  rect(0, exons_y - exon_height / 2, sum(coding1$length),
       exons_y + exon_height / 2, col = color1, border = NA)
  rect(sum(coding1$length), exons_y - exon_height / 2,
       sum(coding1$length) + sum(coding2$length),
       exons_y + exon_height / 2, col = color2, border = NA)

  boundaries <- cumsum(c(coding1$length, coding2$length))
  if (length(boundaries) > 1) {
    boundaries <- boundaries[-length(boundaries)]
    for (b in boundaries)
      lines(c(b, b), c(exons_y - exon_height, exons_y + exon_height),
            col = "white", lty = 3)
  }

  nest_domains <- function(domains) {
    if (length(unlist(domains)) == 0) return(domains)
    domains <- domains[order(domains$end - domains$start,
                             decreasing = TRUE), ]
    rownames(domains) <- seq_len(nrow(domains))
    domains$parent <- 0
    for (d in rownames(domains))
      domains[domains$start >= domains[d, "start"] &
              domains$end <= domains[d, "end"] &
              rownames(domains) != d, "parent"] <- d
    max_overlap <- max(IRanges::coverage(
      IRanges::IRanges(domains$start * 10e6, domains$end * 10e6)))
    padding <- 1 / max_overlap * 0.4
    domains$y <- 0
    domains$height <- 0
    adjust_fn <- function(parent_d, y_pos, h, pad, e) {
      for (d in which(e$domains$parent == parent_d)) {
        overlapping <- which(
          (.between(e$domains$start, e$domains[d, "start"],
                    e$domains[d, "end"]) |
           .between(e$domains$end, e$domains[d, "start"],
                    e$domains[d, "end"])) &
          e$domains$parent == parent_d)
        e$domains[d, "height"] <- h / length(overlapping) -
          pad * (length(overlapping) - 1) / length(overlapping)
        e$domains[d, "y"] <- y_pos +
          (which(d == overlapping) - 1) *
          (e$domains[d, "height"] + pad)
        adjust_fn(d, e$domains[d, "y"] + pad,
                  e$domains[d, "height"] - 2 * pad, pad, e)
      }
    }
    adjust_fn(0, 0, 1, padding, environment())
    domains[order(domains$height, decreasing = TRUE), ]
  }
  retained1 <- nest_domains(retained1)
  retained2 <- nest_domains(retained2)

  if (!is.null(retained1) && nrow(retained1) > 0) {
    retained1$y <- exons_y - exon_height / 2 + 0.025 +
      (exon_height - 2 * 0.025) * retained1$y
    retained1$height <- retained1$height * (exon_height - 2 * 0.025)
  }
  if (!is.null(retained2) && nrow(retained2) > 0) {
    retained2$y <- exons_y - exon_height / 2 + 0.025 +
      (exon_height - 2 * 0.025) * retained2$y
    retained2$height <- retained2$height * (exon_height - 2 * 0.025)
  }

  draw_rect <- function(l, b, r, t, color) {
    rect(l, b, r, t, col = color, border = .get_dark_color(color))
    gs <- 20
    .draw_vertical_gradient(rep(l, gs), rep(r, gs),
      seq(t, b, len = gs), rgb(1, 1, 1, 0.7))
    .draw_vertical_gradient(rep(l, gs), rep(r, gs),
      seq(b, b + (t - b) * 0.4, len = gs), rgb(0, 0, 0, 0.1))
  }

  if (!is.null(retained1) && nrow(retained1) > 0)
    for (d in seq_len(nrow(retained1)))
      draw_rect(retained1[d, "start"], retained1[d, "y"],
                retained1[d, "end"],
                retained1[d, "y"] + retained1[d, "height"],
                retained1[d, "color"])
  if (!is.null(retained2) && nrow(retained2) > 0)
    for (d in seq_len(nrow(retained2)))
      draw_rect(sum(coding1$length) + retained2[d, "start"],
                retained2[d, "y"],
                sum(coding1$length) + retained2[d, "end"],
                retained2[d, "y"] + retained2[d, "height"],
                retained2[d, "color"])

  if (len1 > 0)
    text(sum(coding1$length) / 2, gene_names_y, fusion$gene1,
         font = 2, cex = font_size)
  if (len2 > 0)
    text(sum(coding1$length) + sum(coding2$length) / 2, gene_names_y,
         fusion$gene2, font = 2, cex = font_size)

  count_unique <- function(domains) {
    if (length(unlist(domains)) == 0) return(0)
    n <- 1
    if (nrow(domains) > 1) {
      prev <- domains[1, "proteinDomainID"]
      for (d in 2:nrow(domains)) {
        if (prev != domains[d, "proteinDomainID"]) n <- n + 1
        prev <- domains[d, "proteinDomainID"]
      }
    }
    n
  }

  if (!is.null(retained1) && nrow(retained1) > 0)
    retained1 <- retained1[order(retained1$start), ]
  unique1 <- count_unique(retained1)
  if (!is.null(retained2) && nrow(retained2) > 0)
    retained2 <- retained2[order(retained2$end, decreasing = TRUE), ]
  unique2 <- count_unique(retained2)

  title_y <- exons_y + exon_height / 2 + (unique1 + 2) * 0.05
  text(0.5, title_y + 0.01, "RETAINED PROTEIN DOMAINS",
       adj = c(0.5, 0), font = 2, cex = font_size)

  .draw_domain_labels_gene1(retained1, coding1, exons_y, exon_height,
                            unique1, font_size)
  .draw_domain_labels_gene2(retained2, coding1, exons_y, exon_height,
                            unique2, font_size)
}

.draw_domain_labels_gene1 <- function(retained, coding, exons_y,
                                      exon_height, unique_n, font_size) {
  if (is.null(retained) || nrow(retained) == 0) return(invisible(NULL))
  prev_cx <- -1
  prev_lx <- -1
  label_y <- exons_y + exon_height / 2 + unique_n * 0.05
  for (d in seq_len(nrow(retained))) {
    cx <- min(retained[d, "start"] + 0.01,
              (retained[d, "start"] + retained[d, "end"]) / 2)
    if (cx - prev_cx < 0.01 && retained[d, "end"] > prev_cx + 0.01)
      cx <- prev_cx + 0.01
    lx <- max(cx, prev_lx) + 0.02
    adj_same <- d + 1 <= nrow(retained) &&
      retained[d + 1, "proteinDomainID"] == retained[d, "proteinDomainID"]
    if (adj_same) {
      lx <- retained[d + 1, "start"] + 0.015
    } else {
      text(lx, label_y, retained[d, "proteinDomainName"],
           adj = c(0, 0.5),
           col = .get_dark_color(retained[d, "color"]), cex = font_size)
    }
    lines(c(lx - 0.005, cx, cx),
          c(label_y, label_y, retained[d, "y"] + retained[d, "height"]),
          col = .get_dark_color(retained[d, "color"]))
    if (!adj_same) label_y <- label_y - 0.05
    prev_cx <- cx
    prev_lx <- lx
  }
}

.draw_domain_labels_gene2 <- function(retained, coding1, exons_y,
                                      exon_height, unique_n, font_size) {
  if (is.null(retained) || nrow(retained) == 0) return(invisible(NULL))
  code_len <- sum(coding1$length)
  prev_cx <- 100
  prev_lx <- 100
  label_y <- exons_y - exon_height / 2 - (unique_n + 1) * 0.05
  for (d in seq_len(nrow(retained))) {
    cx <- code_len + max(retained[d, "end"] - 0.01,
      (retained[d, "start"] + retained[d, "end"]) / 2)
    if (prev_cx - cx < 0.01 &&
        code_len + retained[d, "start"] < prev_cx - 0.01)
      cx <- prev_cx - 0.01
    lx <- min(cx, prev_lx) - 0.02
    adj_same <- d + 1 <= nrow(retained) &&
      retained[d + 1, "proteinDomainID"] == retained[d, "proteinDomainID"]
    if (adj_same) {
      lx <- code_len + retained[d + 1, "end"] - 0.015
    } else {
      text(lx, label_y, retained[d, "proteinDomainName"],
           adj = c(1, 0.5),
           col = .get_dark_color(retained[d, "color"]), cex = font_size)
    }
    lines(c(lx + 0.005, cx, cx),
          c(label_y, label_y, retained[d, "y"]),
          col = .get_dark_color(retained[d, "color"]))
    if (!adj_same) label_y <- label_y + 0.05
    prev_cx <- cx
    prev_lx <- lx
  }
}
