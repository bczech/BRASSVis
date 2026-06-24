.find_exons <- function(exons, contig, gene, direction, breakpoint,
                       coverage = NULL, transcript_id = ".",
                       transcript_selection = "provided") {
  if (transcript_selection == "provided" && transcript_id != "." &&
      transcript_id != "") {
    candidate <- exons[exons$transcript == transcript_id, ]
    if (nrow(candidate) > 0) return(candidate)
    warning(paste0("Unknown transcript (", transcript_id,
                   "), selecting a different one"))
  }

  if (transcript_selection == "canonical") {
    candidate <- exons[exons$geneName == gene & exons$contig == contig, ]
  } else {
    transcripts <- exons[
      exons$geneName == gene & exons$contig == contig &
      exons$type == "exon" &
      (direction == "downstream" & abs(exons$end - breakpoint) <= 2 |
       direction == "upstream" & abs(exons$start - breakpoint) <= 2),
      "transcript"
    ]
    candidate <- exons[exons$transcript %in% transcripts, ]
    if (nrow(candidate) == 0) {
      candidate <- exons[exons$geneName == gene &
                         exons$contig == contig, ]
      if (length(unique(candidate$geneID)) > 0) {
        dist <- aggregate(seq_len(nrow(candidate)),
          by = list(candidate$geneID),
          function(x) min(abs(candidate[x, "start"] - breakpoint),
                          abs(candidate[x, "end"] - breakpoint)))
        closest <- head(dist[dist[, 2] == min(dist[, 2]), 1], 1)
        candidate <- candidate[candidate$geneID == closest, ]
      }
    }

    if (!is.null(coverage) && length(unique(candidate$transcript)) > 1) {
      best_cov <- -1
      best_transcript <- NULL
      best_len <- 0
      for (tr in unique(candidate$transcript)) {
        ex <- candidate[candidate$transcript == tr, ]
        ex$start <- sapply(ex$start, max, min(IRanges::start(coverage)))
        ex$end <- sapply(ex$end, min, max(IRanges::end(coverage)))
        tr_len <- sum(ex$end - ex$start + 1)
        cov_sum <- sum(as.numeric(
          coverage[IRanges::IRanges(ex$start, ex$end)]))
        diff <- (1 - min(tr_len, best_len) /
                 max(tr_len, best_len)) / 10
        if ((tr_len > best_len && cov_sum * (1 - diff) > best_cov) ||
            (tr_len < best_len && cov_sum > best_cov * (1 - diff))) {
          best_cov <- cov_sum
          best_transcript <- tr
          best_len <- tr_len
        }
      }
      candidate <- candidate[candidate$transcript == best_transcript, ]
    }

    if (length(unique(candidate$transcript)) > 1) {
      tr_start <- aggregate(candidate$start,
        by = list(candidate$transcript), min)
      rownames(tr_start) <- tr_start[, 1]
      tr_end <- aggregate(candidate$end,
        by = list(candidate$transcript), max)
      rownames(tr_end) <- tr_end[, 1]
      candidate <- candidate[
        .between(breakpoint, tr_start[candidate$transcript, 2],
                 tr_end[candidate$transcript, 2]), ]
    }
  }

  if (length(unique(candidate$transcript)) > 1) {
    score <- ifelse(grepl("appris_principal_1", candidate$attributes), 12,
      ifelse(grepl("appris_principal_2", candidate$attributes), 11,
      ifelse(grepl("appris_principal_3", candidate$attributes), 10,
      ifelse(grepl("appris_principal_4", candidate$attributes), 9,
      ifelse(grepl("appris_principal_5", candidate$attributes), 8,
      ifelse(grepl("appris_principal", candidate$attributes), 7,
      ifelse(grepl("appris_candidate_longest", candidate$attributes), 6,
      ifelse(grepl("appris_candidate", candidate$attributes), 5,
      ifelse(grepl("appris_alternative_1", candidate$attributes), 4,
      ifelse(grepl("appris_alternative_2", candidate$attributes), 3,
      ifelse(grepl("appris_alternative", candidate$attributes), 2,
      ifelse(grepl("CCDS", candidate$attributes), 1, 0))))))))))))
    candidate <- candidate[score == max(score), ]
  }

  if (length(unique(candidate$transcript)) > 1) {
    cds_len <- ifelse(candidate$type == "CDS",
                      candidate$end - candidate$start, 0)
    total <- aggregate(cds_len, by = list(candidate$transcript), sum)
    rownames(total) <- total[, 1]
    candidate <- candidate[total[candidate$transcript, 2] ==
                           max(total[, 2]), ]
  }

  if (length(unique(candidate$transcript)) > 1) {
    ex_len <- candidate$end - candidate$start
    total <- aggregate(ex_len, by = list(candidate$transcript), sum)
    rownames(total) <- total[, 1]
    candidate <- candidate[total[candidate$transcript, 2] ==
                           max(total[, 2]), ]
  }

  unique(candidate[candidate$transcript ==
                   head(unique(candidate$transcript), 1), ])
}
