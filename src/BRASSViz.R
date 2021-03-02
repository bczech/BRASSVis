#!/usr/bin/env Rscript

### <PACKAGES> ###
suppressMessages(library("circlize"))
suppressMessages(library("ggplot2"))
suppressMessages(library("grDevices"))
suppressMessages(library("RSQLite"))
suppressMessages(library("D3GB"))
suppressMessages(library("glue"))
suppressMessages(library("docopt"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("checkmate"))

### </PACKAGES> ###


### <MAN> ###
"usage: BRASSViz.R [options]

options:
-i --inputFile=<file> BRASS file path.\n
-e --exonsFile=<file> annnotation file.\n
-o --outputFile=<file> output file name.\n
--confidence=<character> confidence of fusions [default: high].\n
--fusionFlagThreshold=<numeric> fusion flag threshold [default: 720].\n
--transcriptSelection=<character> transcript selection [default: provided].\n
--pdfWidth=<numeric> pdf width [default: 11.692].\n
--pdfHeight=<numeric> pdf height [default: 8.267].\n
--cytobandsFile=<file> cytobans file.\n
--proteinDomainsFile=<character> protein domains file.\n
--color1=<character> color 1 used in the visualization [default: #e5a5a5].\n
--color2=<character> color 2 used in the visualization [default: #a7c4e5].\n
--printExonLabels=<logical> print exon labels [default: TRUE].\n
--minConfidenceForCircosPlot=<character> minimal confidence for circos plot [default: medium].\n
--mergeDomainsOverlappingBy=<numeric> merge domains overlapping by [default: 0.9].\n
--optimizeDomainColors=<logical> optimize domain colors [default: FALSE]
--fontSize=<numeric> font size [default: 1]" -> doc
opt <- docopt(doc)

opt$fusionFlagThreshold <- as.numeric(opt$fusionFlagThreshold)
opt$pdfWidth <- as.numeric(opt$pdfWidth)
opt$pdfHeight <- as.numeric(opt$pdfHeight)
opt$printExonLabels <- as.logical(opt$printExonLabels)
opt$mergeDomainsOverlappingBy <- as.numeric(opt$mergeDomainsOverlappingBy)
opt$optimizeDomainColors <- as.logical(opt$optimizeDomainColors)
opt$fontSize <- as.numeric(opt$fontSize)

if(is.null(opt$inputFile)) stop("Input BEDPE file neeed to be defined!")
if(is.null(opt$outputFile)) stop("Output file neeed to be defined!")
if(is.null(opt$exonsFile)) stop("Annotation file neeed to be defined!")

if(is.null(opt$cytobandsFile)) message("Cytobands file not provided. Using default GRCh38 ideograms...")
if(is.null(opt$proteinDomainsFile)) message("Protein domains file not provided.")

### </MAN> ###


### <FUNCTIONS> ###
# read and transform BRASS bedpe
readBrass <- function(path) {
  checkmate::assert_string(path,
                           pattern = "bedpe")

  brasInput <- data.table::fread(path,
                                 sep = "\t",
                                 fill = TRUE,
                                 skip = "# chr1")

  data.frame(`gene1` = brasInput$gene1,
             `gene2` = brasInput$gene2,
             `strand1` = paste0(brasInput$strand1, "/", brasInput$strand1),
             `strand2` = paste0(brasInput$strand2, "/", brasInput$strand2),
             `direction1` = ifelse(brasInput$strand1 == "+", "downstream", "upstream"),
             `direction2` = ifelse(brasInput$strand2 == "-", "downstream", "upstream"),
             `breakpoint1` = paste0(brasInput$`# chr1`, ":", brasInput$start1),
             `breakpoint2` = paste0(brasInput$chr2, ":", brasInput$start2),
             `site2` = 'splice-site',
             `type` = brasInput$svclass,
             `transcript_id1` = brasInput$transcript_id1,
             `transcript_id2` = brasInput$transcript_id2,
             `fusion_flag` = brasInput$fusion_flag,
             `total_region_count1` = brasInput$total_region_count1,
             `total_region_count2` = brasInput$total_region_count2,
             `assembly_score` = brasInput$assembly_score,
             `exonsSelected1` = brasInput$region_number1,
             `exonsSelected2` = brasInput$region_number2,
             `stringsAsFactors` = FALSE)

}

# remove 'chr' from contig name
removeChr <- function(contig) {
  sub("^chr", "", sub("^chrM", "MT", contig), perl = TRUE)
}

# convenience function to check if a value is between two others
between <- function(value, start, end) {
  value >= start & value <= end
}

# define colors
changeColorBrightness <- function(color, delta) {
  rgb(
    min(255, max(0,col2rgb(color)["red",] + delta)),
    min(255, max(0,col2rgb(color)["green",] + delta)),
    min(255, max(0,col2rgb(color)["blue",] + delta)),
    maxColorValue = 255
  )
}

# get dark color of a given RGB
getDarkColor <- function(color) {
  changeColorBrightness(color, - 100)
}

# draw vertical gradient
drawVerticalGradient <- function(left, right, y, color, selection=NULL) {
  # check if gradient should only be drawn in part of the region
  if (!is.null(selection)) {
    y <- y[selection]
    left <- left[selection]
    right <- right[selection]
  }
  # draw gradient
  for (i in seq_len(length(y))) {
    polygon(
      c(left[seq_len(i)], right[seq_len(i)]),
      c(y[seq_len(i)], y[seq_len(i)]),
      border = NA,
      col = rgb(col2rgb(color)["red",], col2rgb(color)["green",], col2rgb(color)["blue",], col2rgb(color, alpha = T)["alpha",] * (1 / length(y)), max = 255)
    )
  }
}

# parse gtf file
parseGtfAttribute <- function(attribute, exons) {
  parsed <- gsub(paste0(".*", attribute, " \"?([^;\"]+)\"?;.*"), "\\1", exons$attributes)
  failedToParse <- parsed == exons$attributes
  if (any(failedToParse)) {
    warning(paste0("Failed to parse '", attribute, "' attribute of ", sum(failedToParse), " GTF record(s)."))
    parsed <- ifelse(failedToParse, "", parsed)
  }
  return(parsed)
}

# draw curly brace on plot
drawCurlyBrace <- function(left, right, top, bottom, tip) {
  smoothness <- 20
  x <- cumsum(exp(-seq(-2.5, 2.5, len = smoothness) ^ 2))
  x <- x / max(x)
  y <- seq(top, bottom, len = smoothness)
  lines(left + (tip-left) + x * (left - tip), y)
  lines(tip + x * (right - tip), y)
}

# draw ideogram
drawIdeogram <- function(adjust, left, right, y, cytobands, contig, breakpoint) {
  # define design of ideogram
  bandColors <- setNames(rgb(100:0, 100:0, 100:0, maxColorValue = 100), paste0("gpos", 0:100))
  bandColors <- c(bandColors, gneg="#ffffff", acen="#ec4f4f", stalk="#0000ff")
  cytobands$color <- bandColors[cytobands$giemsa]
  arcSteps <- 30 # defines roundness of arc
  curlyBraceHeight <- 0.03
  ideogramHeight <- 0.04
  ideogramWidth <- 0.4
  # extract bands of given contig
  bands <- cytobands[cytobands$contig == contig,]
  if (nrow(bands) == 0) {
    warning(paste("Ideogram of contig", contig, "cannot be drawn, because no Giemsa staining information is available."))
    return(NULL)
  }
  # scale width of ideogram to fit inside given region
  bands$left <- bands$start / max(cytobands$end) * ideogramWidth
  bands$right <- bands$end / max(cytobands$end) * ideogramWidth
  # left/right-align cytobands
  offset <- ifelse(adjust=="left", left, right - max(bands$right))
  bands$left <- bands$left + offset
  bands$right <- bands$right + offset
  # draw curly braces
  tip <- min(bands$left) + (max(bands$right)-min(bands$left)) / (max(bands$end) - min(bands$start)) * breakpoint
  drawCurlyBrace(left, right, y-0.05+curlyBraceHeight, y - 0.05, tip)
  # draw title of chromosome
  text((max(bands$right)+min(bands$left)) / 2, y + 0.07, paste("chromosome", contig), font = 2, cex = as.numeric(opt$fontSize), adj = c(0.5,0))
  # draw name of band
  bandName <- bands[which(between(breakpoint, bands$start, bands$end)), "name"]
  tryCatch({
    text(tip, y + 0.03, bandName, cex = as.numeric(opt$fontSize), adj = c(0.5,0))
  }, error = function(e) {
    print("Lack of name of band")
  })

  # draw start of chromosome
  leftArcX <- bands[1,"left"] + (1+cos(seq(pi / 2,1.5 * pi,len=arcSteps))) * (bands[1,"right"] - bands[1,"left"])
  leftArcY <- y + sin(seq(pi / 2,1.5 * pi, len = arcSteps)) * (ideogramHeight / 2)
  polygon(leftArcX, leftArcY, col = bands[1,"color"])
  # draw bands
  centromereStart <- NULL
  centromereEnd <- NULL
  for (band in 2:(nrow(bands)-1)) {
    if (bands[band,"giemsa"] != "acen") {
      rect(bands[band,"left"], y - ideogramHeight / 2, bands[band,"right"], y + ideogramHeight / 2, col = bands[band,"color"])
    } else { # draw centromere
      if (is.null(centromereStart)) {
        polygon(c(bands[band,"left"], bands[band,"right"], bands[band,"left"]), c(y - ideogramHeight / 2, y, y + ideogramHeight / 2), col = bands[band,"color"])
        centromereStart <- bands[band,"left"]
      } else {
        polygon(c(bands[band,"right"], bands[band,"left"], bands[band,"right"]), c(y - ideogramHeight / 2, y, y + ideogramHeight / 2), col = bands[band,"color"])
        centromereEnd <- bands[band,"right"]
      }
    }
  }
  # draw end of chromosome
  band <- nrow(bands)
  rightArcX <- bands[band,"right"] - (1 + cos(seq(1.5 * pi,pi / 2, len = arcSteps))) * (bands[band,"right"] - bands[band,"left"])
  rightArcY <- y + sin(seq(pi / 2,1.5 * pi, len = arcSteps)) * ideogramHeight/2
  polygon(rightArcX, rightArcY, col = bands[band,"color"])
  # if there is no centromere, make an artificial one with length zero
  if (is.null(centromereStart) || is.null(centromereEnd)) {
    centromereStart <- bands[1,"right"]
    centromereEnd <- bands[1,"right"]
  }
  # draw gradients for 3D effect
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.8), seq_len(round(arcSteps * 0.4))) # black from top on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps * 0.4):round(arcSteps * 0.1)) # white to top on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps * 0.4):round(arcSteps * 0.6)) # white to bottom on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps * 0.5)) # black from bottom on p-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.8), seq_len(round(arcSteps * 0.4))) # black from top on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps * 0.4):round(arcSteps * 0.1)) # white to top on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps * 0.4):round(arcSteps * 0.6)) # white to bottom on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps * 0.5)) # black from bottom on q-arm
}

# draw strand
drawStrand <- function(left, right, y, color, strand) {
  if (strand %in% c("+", "-")) {
    # draw strand
    lines(c(left + 0.001, right - 0.001), c(y, y), col = color, lwd = 2)
    lines(c(left + 0.001, right - 0.001), c(y, y), col = rgb(1,1,1,0.1), lwd = 1)
    # indicate orientation
    if (right - left > 0.01)
      for (i in seq(left + 0.005, right - 0.005, by = sign(right - left - 2 * 0.005) * 0.01)) {
        arrows(i, y, i + 0.001 * ifelse(strand == "+", 1, -1), y, col = color, length = 0.05, lwd = 2, angle = 60)
        arrows(i, y, i + 0.001 * ifelse(strand == "+", 1, -1), y, col = rgb(1,1,1,0.1), length = 0.05, lwd = 1, angle = 60)
      }
  }
}

# draw exons
drawExon <- function(left, right, y, color, title, type) {
  gradientSteps <- 10 # defines smoothness of gradient
  exonHeight <- 0.03
  if (type == "CDS") {
    # draw coding regions as thicker bars
    rect(left, y + exonHeight, right, y + exonHeight / 2 - 0.001, col = color, border = NA)
    rect(left, y - exonHeight, right, y - exonHeight / 2 + 0.001, col = color, border = NA)
    # draw border
    lines(c(left, left, right, right), c(y + exonHeight / 2, y + exonHeight, y + exonHeight, y + exonHeight / 2), col = getDarkColor(color), lend = 2)
    lines(c(left, left, right, right), c(y - exonHeight / 2, y - exonHeight, y - exonHeight, y - exonHeight / 2), col = getDarkColor(color), lend = 2)
    # draw gradients for 3D effect
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y + 0.03, y + 0.015, len = gradientSteps), rgb(0,0,0,0.2))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y - 0.03, y - 0.015, len = gradientSteps), rgb(0,0,0,0.3))
  } else if (type == "exon") {
    rect(left, y + exonHeight / 2, right, y - exonHeight / 2, col = color, border = getDarkColor(color))
    # draw gradients for 3D effect
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y + exonHeight / 2, len = gradientSteps), rgb(1,1,1,0.6))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y - exonHeight / 2, len = gradientSteps), rgb(1,1,1,0.6))
    # add exon label
    text((left + right) / 2, y, title, cex = 0.9 * as.numeric(opt$fontSize))
  }
}

# draw circos
drawCircos <- function(fusion, fusions, cytobands, minConfidenceForCircosPlot) {
  # check if Giemsa staining information is available
  for (contig in unlist(fusions[fusion,c("contig1", "contig2")])) {
    if (!any(cytobands$contig == contig)) {
      warning(paste0("Circos plot cannot be drawn, because no Giemsa staining information is available for contig ", contig, "."))
      return(NULL)
    }
  }
  # initialize with empty circos plot
  circos.clear()
  circos.initializeWithIdeogram(cytoband=cytobands, labels.cex = 1.1, axis.labels.cex = 0.6)
  # use gene names as labels or <contig>:<position> for intergenic breakpoints
  geneLabels <- data.frame(
    contig = c(fusions[fusion,"contig1"], fusions[fusion,"contig2"]),
    start = c(fusions[fusion,"breakpoint1"], fusions[fusion,"breakpoint2"])
  )
  geneLabels$end <- geneLabels$start + 1
  geneLabels$gene <- c(fusions[fusion,"gene1"], fusions[fusion,"gene2"])
  geneLabels$gene <- ifelse(c(fusions[fusion,"site1"], fusions[fusion,"site2"]) == "intergenic", paste0(c(fusions[fusion,"display_contig1"], fusions[fusion,"display_contig2"]), ":", geneLabels$start), geneLabels$gene)
  # draw gene labels
  circos.genomicLabels(geneLabels, labels.column = 4, side ="inside", cex = as.numeric(opt$fontSize) - 0.1)
  # draw chromosome labels in connector plot
  for (contig in unique(cytobands$contig)) {
    set.current.cell(track.index = 2, sector.index = contig) # draw in gene label connector track (track.index=2)
    circos.text(CELL_META$xcenter, CELL_META$ycenter, contig, cex = 0.25)
  }
}

drawProteinDomains <- function(fusion, exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors) {

  exonHeight <- 0.2
  exonsY <- 0.5
  geneNamesY <- exonsY - exonHeight / 2 - 0.05

  # find coding exons
  codingExons1 <- exons1[exons1$type == "CDS" & fusion$site1 != "intergenic",]
  codingExons2 <- exons2[exons2$type == "CDS" & fusion$site2 != "intergenic",]

  # cut off coding regions beyond breakpoint
  if (fusion$direction1 == "upstream") {
    codingExons1 <- codingExons1[codingExons1$end >= fusion$breakpoint1, ]
    codingExons1$start <- ifelse(codingExons1$start < fusion$breakpoint1, fusion$breakpoint1, codingExons1$start)
  } else {
    codingExons1 <- codingExons1[codingExons1$start <= fusion$breakpoint1, ]
    codingExons1$end <- ifelse(codingExons1$end > fusion$breakpoint1, fusion$breakpoint1, codingExons1$end)
  }
  if (fusion$direction2 == "upstream") {
    codingExons2 <- codingExons2[codingExons2$end >= fusion$breakpoint2, ]
    codingExons2$start <- ifelse(codingExons2$start < fusion$breakpoint2, fusion$breakpoint2, codingExons2$start)
  } else {
    codingExons2 <- codingExons2[codingExons2$start <= fusion$breakpoint2, ]
    codingExons2$end <- ifelse(codingExons2$end > fusion$breakpoint2, fusion$breakpoint2, codingExons2$end)
  }

  # find overlapping domains
  exonsGRanges1 <- GRanges(codingExons1$contig, IRanges(codingExons1$start, codingExons1$end), strand = codingExons1$strand)
  exonsGRanges2 <- GRanges(codingExons2$contig, IRanges(codingExons2$start, codingExons2$end), strand = codingExons2$strand)
  domainsGRanges <- GRanges(proteinDomains$contig, IRanges(proteinDomains$start, proteinDomains$end), strand = proteinDomains$strand)
  domainsGRanges$proteinDomainName <- proteinDomains$proteinDomainName
  domainsGRanges$proteinDomainID <- proteinDomains$proteinDomainID
  domainsGRanges$color <- proteinDomains$color
  domainsGRanges <- domainsGRanges[suppressWarnings(unique(queryHits(findOverlaps(domainsGRanges, GenomicRanges::union(exonsGRanges1, exonsGRanges2)))))]

  # group overlapping domains by domain ID
  domainsGRangesList <- GRangesList(lapply(unique(domainsGRanges$proteinDomainID), function(x) { domainsGRanges[domainsGRanges$proteinDomainID == x] }))

  # trim protein domains to exon boundaries
  trimDomains <- function(domainsGRangesList, exonsGRanges) {
    do.call(
      "rbind",
      lapply(
        domainsGRangesList,
        function(x) {
          intersected <- as.data.frame(reduce(suppressWarnings(GenomicRanges::intersect(x, exonsGRanges))))
          if (nrow(intersected) > 0) {
            intersected$proteinDomainName <- head(x$proteinDomainName, 1)
            intersected$proteinDomainID <- head(x$proteinDomainID, 1)
            intersected$color <- head(x$color, 1)
          } else {
            intersected$proteinDomainName <- character()
            intersected$proteinDomainID <- character()
            intersected$color <- character()
          }
          return(intersected)
        }
      )
    )
  }
  retainedDomains1 <- trimDomains(domainsGRangesList, exonsGRanges1)
  retainedDomains2 <- trimDomains(domainsGRangesList, exonsGRanges2)

  # calculate length of coding exons
  codingExons1$length <- codingExons1$end - codingExons1$start + 1
  codingExons2$length <- codingExons2$end - codingExons2$start + 1

  # abort, if there are no coding regions
  if (sum(exons1$type == "CDS") + sum(exons2$type == "CDS") == 0) {
    text(0.5, 0.5, "Genes are not protein-coding.")
    return(NULL)
  }
  codingLength1 <- sum(codingExons1$length)
  codingLength2 <- sum(codingExons2$length)
  if (codingLength1 + codingLength2 == 0) {
    text(0.5, 0.5, "No coding regions retained in fusion transcript.")
    return(NULL)
  }
  antisenseTranscription1 <- sub("/.*", "", fusion$strand1) != sub(".*/", "", fusion$strand1)
  antisenseTranscription2 <- sub("/.*", "", fusion$strand2) != sub(".*/", "", fusion$strand2)
  if ((codingLength1 == 0 || antisenseTranscription1) && (codingLength2 == 0 || antisenseTranscription2)) {
    text(0.5, 0.5, "No coding regions due to antisense transcription.")
    return(NULL)
  }

  # remove introns from protein domains
  removeIntronsFromProteinDomains <- function(codingExons, retainedDomains) {
    if (nrow(codingExons) == 0) return(NULL)
    cumulativeIntronLength <- 0
    previousExonEnd <- 0
    for (exon in seq_len(nrow(codingExons))) {
      if (codingExons[exon, "start"] > previousExonEnd)
        cumulativeIntronLength <- cumulativeIntronLength + codingExons[exon,"start"] - previousExonEnd
      domainsInExon <- which(between(retainedDomains$start, codingExons[exon,"start"], codingExons[exon, "end"]))
      retainedDomains[domainsInExon, "start"] <- retainedDomains[domainsInExon,"start"] - cumulativeIntronLength
      domainsInExon <- which(between(retainedDomains$end, codingExons[exon,"start"], codingExons[exon, "end"]))
      retainedDomains[domainsInExon, "end"] <- retainedDomains[domainsInExon,"end"] - cumulativeIntronLength
      previousExonEnd <- codingExons[exon, "end"]
    }
    # merge adjacent domains
    retainedDomains <- do.call(
      "rbind",
      lapply(
        unique(retainedDomains$proteinDomainID),
        function(x) {
          domain <- retainedDomains[retainedDomains$proteinDomainID == x,]
          merged <- reduce(GRanges(domain$seqnames, IRanges(domain$start, domain$end), strand=domain$strand))
          merged$proteinDomainName <- head(domain$proteinDomainName, 1)
          merged$proteinDomainID <- head(domain$proteinDomainID, 1)
          merged$color <- head(domain$color, 1)
          return(as.data.frame(merged))
        }
      )
    )
    return(retainedDomains)
  }
  retainedDomains1 <- removeIntronsFromProteinDomains(codingExons1, retainedDomains1)
  retainedDomains2 <- removeIntronsFromProteinDomains(codingExons2, retainedDomains2)

  # abort, if no domains are retained
  if (is.null(retainedDomains1) && is.null(retainedDomains2)) {
    text(0.5, 0.5, "No protein domains retained in fusion.")
    return(NULL)
  }

  # merge domains with similar coordinates
  mergeSimilarDomains <- function(domains, mergeDomainsOverlappingBy) {
    if (is.null(domains)) return(domains)
    merged <- domains[FALSE, ] # create empty data frame
    domains <- domains[order(domains$end - domains$start, decreasing = FALSE), ] # start with bigger domains => bigger domains are retained
    for (domain in rownames(domains)) {
      if (!any((abs(merged$start - domains[domain, "start"]) + abs(merged$end - domains[domain,"end"])) / (domains[domain,"end"] - domains[domain,"start"]) <= 1 - mergeDomainsOverlappingBy))
        merged <- rbind(merged, domains[domain,])
    }
    return(merged)
  }
  retainedDomains1 <- mergeSimilarDomains(retainedDomains1, as.logical(opt$mergeDomainsOverlappingBy))
  retainedDomains2 <- mergeSimilarDomains(retainedDomains2, as.logical(opt$mergeDomainsOverlappingBy))

  # if desired, reassign colors to protein domains to maximize contrast
  if (as.logical(opt$optimizeDomainColors)) {
    uniqueDomains <- unique(c(retainedDomains1$proteinDomainID, retainedDomains2$proteinDomainID))
    # make rainbow of pretty pastell colors
    colors <- rainbow(length(uniqueDomains))
    colors <- apply(col2rgb(colors), 2, function(x) { 0.3 + y / 255 * 0.7 }) # make pastell colors
    colors <- apply(colors, 2, function(x) {rgb(x["red"], x["green"], x["blue"])}) # convert back to rgb
    # reassign colors
    names(colors) <- uniqueDomains
    retainedDomains1$color <- colors[retainedDomains1$proteinDomainID]
    retainedDomains2$color <- colors[retainedDomains2$proteinDomainID]
  }

  # normalize length to 1
  codingExons1$length <- codingExons1$length / (codingLength1 + codingLength2)
  codingExons2$length <- codingExons2$length / (codingLength1 + codingLength2)
  retainedDomains1$start <- retainedDomains1$start / (codingLength1 + codingLength2)
  retainedDomains1$end <- retainedDomains1$end / (codingLength1 + codingLength2)
  retainedDomains2$start <- retainedDomains2$start / (codingLength1 + codingLength2)
  retainedDomains2$end <- retainedDomains2$end / (codingLength1 + codingLength2)

  # draw coding regions
  rect(0, exonsY - exonHeight / 2, sum(codingExons1$length), exonsY + exonHeight / 2, col = opt$color1, border = NA)
  rect(sum(codingExons1$length), exonsY - exonHeight / 2, sum(codingExons1$length) + sum(codingExons2$length), exonsY + exonHeight / 2, col = opt$color2, border = NA)

  # indicate exon boundaries as dotted lines
  exonBoundaries <- cumsum(c(codingExons1$length, codingExons2$length))
  if (length(exonBoundaries) > 1) {
    exonBoundaries <- exonBoundaries[1:(length(exonBoundaries)-1)]
    for (exonBoundary in exonBoundaries)
      lines(c(exonBoundary, exonBoundary), c(exonsY - exonHeight, exonsY + exonHeight), col = "white", lty = 3)
  }

  # find overlapping domains
  # nest if one is contained in another
  # stack if they overlap partially
  nestDomains <- function(domains) {
    if (length(unlist(domains)) == 0) return(domains)
    domains <- domains[order(domains$end - domains$start, decreasing = TRUE), ]
    rownames(domains) <- seq_len(nrow(domains))
    # find nested domains and make tree structure
    domains$parent <- 0
    for (domain in rownames(domains))
      domains[domains$start >= domains[domain, "start"] & domains$end <= domains[domain, "end"] & rownames(domains) != domain, "parent"] <- domain
    # find partially overlapping domains
    maxOverlappingDomains <- max(coverage(IRanges(domains$start*10e6, domains$end*10e6)))
    padding <- 1 / maxOverlappingDomains * 0.4
    domains$y <- 0
    domains$height <- 0
    adjustPositionAndHeight <- function(parentDomain, y, height, padding, e) {
      for (domain in which(e$domains$parent == parentDomain)) {
        overlappingDomains <- which((between(e$domains$start, e$domains[domain, "start"], e$domains[domain, "end"]) |
                                       between(e$domains$end, e$domains[domain, "start"], e$domains[domain, "end"])) &
                                      e$domains$parent == parentDomain)
        e$domains[domain, "height"] <- height / length(overlappingDomains) - padding * (length(overlappingDomains) - 1) / length(overlappingDomains)
        e$domains[domain, "y"] <- y + (which(domain == overlappingDomains) - 1) * (e$domains[domain,"height"] + padding)
        adjustPositionAndHeight(domain, e$domains[domain, "y"] + padding, e$domains[domain, "height"] - 2 * padding, padding, e)
      }
    }
    adjustPositionAndHeight(0, 0, 1, padding, environment())
    domains <- domains[order(domains$height, decreasing=TRUE), ] # draw nested domains last
    return(domains)
  }
  retainedDomains1 <- nestDomains(retainedDomains1)
  retainedDomains2 <- nestDomains(retainedDomains2)
  retainedDomains1$y <- exonsY - exonHeight / 2 + 0.025 + (exonHeight - 2 * 0.025) * retainedDomains1$y
  retainedDomains2$y <- exonsY - exonHeight / 2 + 0.025 + (exonHeight - 2 * 0.025) * retainedDomains2$y
  retainedDomains1$height <- retainedDomains1$height * (exonHeight - 2 * 0.025)
  retainedDomains2$height <- retainedDomains2$height * (exonHeight - 2 * 0.025)

  # draw domains
  drawProteinDomainRect <- function(left, bottom, right, top, color) {
    rect(left, bottom, right, top, col = color, border = getDarkColor(color))
    # draw gradients for 3D effect
    gradientSteps <- 20
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(top, bottom, len=gradientSteps), rgb(1,1,1,0.7))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(bottom, bottom+(top-bottom) * 0.4, len = gradientSteps), rgb(0,0,0,0.1))
  }
  if (length(unlist(retainedDomains1)) > 0)
    for (domain in seq_len(nrow(retainedDomains1)))
      drawProteinDomainRect(retainedDomains1[domain, "start"], retainedDomains1[domain, "y"], retainedDomains1[domain, "end"], retainedDomains1[domain, "y"] + retainedDomains1[domain, "height"], retainedDomains1[domain, "color"])
  if (length(unlist(retainedDomains2)) > 0)
    for (domain in seq_len(nrow(retainedDomains2)))
      drawProteinDomainRect(sum(codingExons1$length)+retainedDomains2[domain,"start"], retainedDomains2[domain,"y"], sum(codingExons1$length)+retainedDomains2[domain,"end"], retainedDomains2[domain,"y"]+retainedDomains2[domain,"height"], retainedDomains2[domain,"color"])

  # draw gene names, if there are coding exons
  if (codingLength1 > 0)
    text(sum(codingExons1$length) / 2, geneNamesY, fusion$gene1, font = 2, cex = opt$fontSize)
  if (codingLength2 > 0)
    text(sum(codingExons1$length) + sum(codingExons2$length) / 2, geneNamesY, fusion$gene2, font = 2, cex = opt$fontSize)

  # calculate how many non-adjacent unique domains there are
  # we need this info to know where to place labels vertically
  countUniqueDomains <- function(domains) {
    uniqueDomains <- 0
    if (length(unlist(domains)) > 0) {
      uniqueDomains <- 1
      if (nrow(domains) > 1) {
        previousDomain <- domains[1, "proteinDomainID"]
        for (domain in 2:nrow(domains)) {
          if (previousDomain != domains[domain, "proteinDomainID"])
            uniqueDomains <- uniqueDomains + 1
          previousDomain <- domains[domain, "proteinDomainID"]
        }
      }
    }
    return(uniqueDomains)
  }
  if (length(unlist(retainedDomains1)) > 0)
    retainedDomains1 <- retainedDomains1[order(retainedDomains1$start), ]
  uniqueDomains1 <- countUniqueDomains(retainedDomains1)
  if (length(unlist(retainedDomains2)) > 0)
    retainedDomains2 <- retainedDomains2[order(retainedDomains2$end, decreasing = TRUE), ]
  uniqueDomains2 <- countUniqueDomains(retainedDomains2)

  # draw title of plot
  titleY <- exonsY + exonHeight / 2 + (uniqueDomains1 + 2) * 0.05
  text(0.5, titleY + 0.01, "RETAINED PROTEIN DOMAINS", adj = c(0.5, 0), font = 2, cex = opt$fontSize)
  # draw domain labels for gene1
  if (length(unlist(retainedDomains1)) > 0) {
    previousConnectorX <- -1
    previousLabelX <- -1
    labelY <- exonsY + exonHeight / 2 + uniqueDomains1 * 0.05
    for (domain in seq_len(nrow(retainedDomains1))) {
      # if possible avoid overlapping lines of labels
      connectorX <- min(retainedDomains1[domain, "start"] + 0.01, (retainedDomains1[domain, "start"] + retainedDomains1[domain, "end"]) / 2)
      if (connectorX - previousConnectorX < 0.01 && retainedDomains1[domain, "end"] > previousConnectorX + 0.01)
        connectorX <- previousConnectorX + 0.01
      labelX <- max(connectorX, previousLabelX) + 0.02
      # use a signle label for adjacent domains of same type
      adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains1) && retainedDomains1[domain + 1, "proteinDomainID"] == retainedDomains1[domain, "proteinDomainID"]
      if (adjacentDomainsOfSameType) {
        labelX <- retainedDomains1[domain + 1, "start"] + 0.015
      } else {
        text(labelX, labelY, retainedDomains1[domain, "proteinDomainName"], adj = c(0,0.5), col = getDarkColor(retainedDomains1[domain, "color"]), cex = opt$fontSize)
      }
      lines(c(labelX - 0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains1[domain, "y"] + retainedDomains1[domain, "height"]), col = getDarkColor(retainedDomains1[domain, "color"]))
      if (!adjacentDomainsOfSameType)
        labelY <- labelY - 0.05
      previousConnectorX <- connectorX
      previousLabelX <- labelX
    }
  }

  # draw domain labels for gene2
  if (length(unlist(retainedDomains2)) > 0) {
    previousConnectorX <- 100
    previousLabelX <- 100
    labelY <- exonsY - exonHeight / 2 - (uniqueDomains2+1) * 0.05
    for (domain in seq_len(nrow(retainedDomains2))) {
      # if possible avoid overlapping connector lines of labels
      connectorX <- sum(codingExons1$length) + max(retainedDomains2[domain, "end"] - 0.01, (retainedDomains2[domain, "start"] + retainedDomains2[domain, "end"]) / 2)
      if (previousConnectorX - connectorX < 0.01 && sum(codingExons1$length) + retainedDomains2[domain, "start"] < previousConnectorX - 0.01)
        connectorX <- previousConnectorX - 0.01
      labelX <- min(connectorX, previousLabelX) - 0.02
      # use a signle label for adjacent domains of same type
      adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains2) && retainedDomains2[domain + 1, "proteinDomainID"] == retainedDomains2[domain, "proteinDomainID"]
      if (adjacentDomainsOfSameType) {
        labelX <- sum(codingExons1$length) + retainedDomains2[domain + 1, "end"] - 0.015
      } else {
        text(labelX, labelY, retainedDomains2[domain, "proteinDomainName"], adj = c(1,0.5), col = getDarkColor(retainedDomains2[domain, "color"]), cex = opt$fontSize)
      }
      lines(c(labelX + 0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains2[domain, "y"]), col = getDarkColor(retainedDomains2[domain, "color"]))
      if (!adjacentDomainsOfSameType)
        labelY <- labelY + 0.05
      previousConnectorX <- connectorX
      previousLabelX <- labelX
    }
  }

}

findExons <- function(exons,  contig, gene, direction, breakpoint, coverage, transcriptId, transcriptSelection) {
  # use the provided transcript if desired
  if (transcriptSelection == "provided" && transcriptId != "." && transcriptId != "") {
    candidateExons <- exons[exons$transcript == transcriptId,]
    if (nrow(candidateExons) == 0) {
      warning(paste0("Unknown transcript given in fusions file (", transcriptId, "), selecting a different one"))
    } else {
      return(candidateExons)
    }
  }

  if (transcriptSelection == "canonical") {
    candidateExons <- exons[exons$geneName == gene & exons$contig == contig,]
  } else {
    # look for exon with breakpoint as splice site
    transcripts <- exons[exons$geneName == gene & exons$contig == contig & exons$type == "exon" & (direction == "downstream" & abs(exons$end - breakpoint) <= 2 | direction == "upstream" & abs(exons$start - breakpoint) <= 2), "transcript"]
    candidateExons <- exons[exons$transcript %in% transcripts, ]
    # if none was found, use all exons of the gene closest to the breakpoint
    if (nrow(candidateExons) == 0) {
      candidateExons <- exons[exons$geneName == gene & exons$contig == contig,]
      if (length(unique(candidateExons$geneID)) > 0) { # more than one gene found with the given name => use the closest one
        distanceToBreakpoint <- aggregate(seq_len(nrow(candidateExons)), by = list(candidateExons$geneID), function(x) { min(abs(candidateExons[x, "start"] - breakpoint), abs(candidateExons[x, "end"] - breakpoint)) })
        closestGene <- head(distanceToBreakpoint[distanceToBreakpoint[, 2] == min(distanceToBreakpoint[, 2]), 1], 1)
        candidateExons <- candidateExons[candidateExons$geneID == closestGene,]
      }
    }
    # if we have coverage information, use the transcript with the highest coverage if there are multiple hits
    if (!is.null(coverage)) {
      highestCoverage <- -1
      transcriptWithHighestCoverage <- NULL
      lengthOfTranscriptWithHighestCoverage <- 0
      for (transcript in unique(candidateExons$transcript)) {
        exonsOfTranscript <- candidateExons[candidateExons$transcript == transcript,]
        exonsOfTranscript$start <- sapply(exonsOfTranscript$start, max, min(start(coverage)))
        exonsOfTranscript$end <- sapply(exonsOfTranscript$end, min, max(end(coverage)))
        lengthOfTranscript <- sum(exonsOfTranscript$end - exonsOfTranscript$start + 1)
        coverageSum <- sum(as.numeric(coverage[IRanges(exonsOfTranscript$start, exonsOfTranscript$end)]))
        # we prefer shorter transcripts over longer ones, because otherwise there is a bias towards transcripts with long UTRs
        # => a longer transcript must have substantially higher coverage to replace a shorter one
        substantialDifference <- (1 - min(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage) / max(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage)) / 10
        if (lengthOfTranscript > lengthOfTranscriptWithHighestCoverage && coverageSum * (1 - substantialDifference) > highestCoverage ||
            lengthOfTranscript < lengthOfTranscriptWithHighestCoverage && coverageSum > highestCoverage * (1 - substantialDifference)) {
          highestCoverage <- coverageSum
          transcriptWithHighestCoverage <- transcript
          lengthOfTranscriptWithHighestCoverage <- lengthOfTranscript
        }
      }
      candidateExons <- candidateExons[candidateExons$transcript==transcriptWithHighestCoverage,]
    }
    # if the gene has multiple transcripts, search for transcripts which encompass the breakpoint
    if (length(unique(candidateExons$transcript)) > 1) {
      transcriptStart <- aggregate(candidateExons$start, by = list(candidateExons$transcript), min)
      rownames(transcriptStart) <- transcriptStart[, 1]
      transcriptEnd <- aggregate(candidateExons$end, by = list(candidateExons$transcript), max)
      rownames(transcriptEnd) <- transcriptEnd[, 1]
      candidateExons <- candidateExons[between(breakpoint, transcriptStart[candidateExons$transcript, 2], transcriptEnd[candidateExons$transcript, 2]), ]
    }
  }

  # find the consensus transcript, if there are multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    consensusTranscript <-
      ifelse(grepl("appris_principal_1", candidateExons$attributes), 12,
             ifelse(grepl("appris_principal_2", candidateExons$attributes), 11,
                    ifelse(grepl("appris_principal_3", candidateExons$attributes), 10,
                           ifelse(grepl("appris_principal_4", candidateExons$attributes), 9,
                                  ifelse(grepl("appris_principal_5", candidateExons$attributes), 8,
                                         ifelse(grepl("appris_principal", candidateExons$attributes), 7,
                                                ifelse(grepl("appris_candidate_longest", candidateExons$attributes), 6,
                                                       ifelse(grepl("appris_candidate", candidateExons$attributes), 5,
                                                              ifelse(grepl("appris_alternative_1", candidateExons$attributes), 4,
                                                                     ifelse(grepl("appris_alternative_2", candidateExons$attributes), 3,
                                                                            ifelse(grepl("appris_alternative", candidateExons$attributes), 2,
                                                                                   ifelse(grepl("CCDS", candidateExons$attributes), 1,
                                                                                          0
                                                                                   ))))))))))))
    candidateExons <- candidateExons[consensusTranscript == max(consensusTranscript), ]
  }
  # use the transcript with the longest coding sequence, if there are still multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    codingSequenceLength <- ifelse(candidateExons$type == "CDS", candidateExons$end - candidateExons$start, 0)
    totalCodingSequenceLength <- aggregate(codingSequenceLength, by=list(candidateExons$transcript), sum)
    rownames(totalCodingSequenceLength) <- totalCodingSequenceLength[, 1]
    candidateExons <- candidateExons[totalCodingSequenceLength[candidateExons$transcript,2] == max(totalCodingSequenceLength[, 2]), ]
  }
  # use the transcript with the longest overall sequence, if there are still multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    exonLength <- candidateExons$end - candidateExons$start
    totalExonLength <- aggregate(exonLength, by=list(candidateExons$transcript), sum)
    rownames(totalExonLength) <- totalExonLength[, 1]
    candidateExons <- candidateExons[totalExonLength[candidateExons$transcript, 2] == max(totalExonLength[, 2]), ]
  }
  # if there are still multiple hits, select the first one
  candidateExons <- unique(candidateExons[candidateExons$transcript == head(unique(candidateExons$transcript), 1), ])
  return(candidateExons)
}
### </FUNCTIONS> ###

### <MAIN> ###
darkColor1 <- getDarkColor(opt$color1)
darkColor2 <- getDarkColor(opt$color2)


fusions <- readBrass(opt$inputFile)
fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1, perl = TRUE)
fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2, perl = TRUE)
fusions$contig1 <- gsub(":.*", "", fusions$breakpoint1)
fusions$contig2 <- gsub(":.*", "", fusions$breakpoint2)
fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl = TRUE))
fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl = TRUE))
fusions$site1 <- rep("exon", nrow(fusions))
fusions$site2 <- rep("exon", nrow(fusions))
fusions$confidence <- rep(opt$confidence, nrow(fusions))

pdf(opt$outputFile, onefile = TRUE, width = as.numeric(opt$pdfWidth), height = as.numeric(opt$pdfHeight), title = "Title")
if (nrow(fusions) == 0) {
  plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
  text(0, 0, "Error: empty input file\n")
  dev.off()
  quit("no")
}
message("Loading ideograms")
if(is.null(opt$cytobandsFile)){
  cytobands <- GRCh38.bands
} else {
  cytobands <- read.table(opt$cytobandsFile, header = TRUE, sep = "\t")
}
colnames(cytobands)[1] <- "contig"
names(cytobands)[5] <- "giemsa"
cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end),]


message("Loading annotation")
exons <- read.table(opt$exonsFile, header = FALSE, sep = "\t", comment.char = "#", quote = "", stringsAsFactors = FALSE)[,c(1, 3, 4, 5, 7, 9)]
colnames(exons) <- c("contig", "type", "start", "end", "strand", "attributes")
exons <- exons[exons$type %in% c("exon", "CDS"),]
exons$contig <- removeChr(exons$contig)


exons$geneID <- parseGtfAttribute("gene_id", exons)
exons$geneName <- parseGtfAttribute("gene_name", exons)
exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
exons$transcript <- parseGtfAttribute("transcript_id", exons)
exons$exonNumber <- ifelse(rep(as.logical(opt$printExonLabels), nrow(exons)), parseGtfAttribute("exon_number", exons), "")

proteinDomains <- NULL
if (!is.null(opt$proteinDomainsFile)) {
  message("Loading protein domains")
  proteinDomains <- read.table(opt$proteinDomainsFile, header = FALSE, sep = "\t", comment.char = "", quote = "", stringsAsFactors = FALSE)[,c(1,4,5,7,9)]
  colnames(proteinDomains) <- c("contig", "start", "end", "strand", "attributes")
  proteinDomains$color <- sub(";.*", "", sub(".*color=", "", proteinDomains$attributes, perl = TRUE), perl = TRUE)
  proteinDomains$proteinDomainName <- sapply(sub(";.*", "", sub(".*Name=", "", proteinDomains$attributes, perl = TRUE), perl = TRUE), URLdecode)
  proteinDomains$proteinDomainID <- sub(";.*", "", sub(".*protein_domain_id=", "", proteinDomains$attributes, perl = TRUE), perl = TRUE)
  proteinDomains <- proteinDomains[ , colnames(proteinDomains) != "attributes"]
}


# insert dummy annotations for intergenic breakpoints
if (any(fusions$site1 == "intergenic" | fusions$site2 == "intergenic")) {
  intergenicBreakpoints <- rbind(
    setNames(fusions[fusions$site1 == "intergenic", c("gene1", "strand1", "contig1", "breakpoint1")], c("gene", "strand", "contig", "breakpoint")),
    setNames(fusions[fusions$site2 == "intergenic", c("gene2", "strand2", "contig2", "breakpoint2")], c("gene", "strand", "contig", "breakpoint"))
  )
  exons <- rbind(exons, data.frame(
    contig = intergenicBreakpoints$contig,
    type = "intergenic",
    start = intergenicBreakpoints$breakpoint-1000,
    end = intergenicBreakpoints$breakpoint+1000,
    strand = ".",
    attributes = "",
    geneName = intergenicBreakpoints$gene,
    geneID = intergenicBreakpoints$gene,
    transcript = intergenicBreakpoints$gene,
    exonNumber = "intergenic"
  ))
}

fusions$tandem <- ifelse(fusions$type == "tandem-duplication", 2, 1)
fusions <- fusions[rep(seq_len(nrow(fusions)), fusions$tandem), ]
duplicatedFusion <- which(duplicated(fusions))
fusions[duplicatedFusion, "direction1"] <- "downstream"
fusions[duplicatedFusion, "direction2"] <- "upstream"

for (fusion in seq_len(nrow(fusions))) {
  message(paste0("Drawing fusion #", fusion, ": ", fusions[fusion, "gene1"], ":", fusions[fusion, "gene2"]))

  if(fusions[fusion, "fusion_flag"] < opt$fusionFlagThreshold) {
    message(glue("Fusion flag lower than {opt$fusionFlagThreshold}. Omiting..."))
    next
  }
  # compute coverage from alignments file
  coverage1 <- NULL
  coverage2 <- NULL
  # find all exons belonging to the fused genes
  exons1 <- findExons(exons, fusions[fusion, "contig1"], fusions[fusion, "gene1"], fusions[fusion, "direction1"], fusions[fusion, "breakpoint1"], coverage1, fusions[fusion, "transcript_id1"], opt$transcriptSelection)
  exons2 <- findExons(exons, fusions[fusion, "contig2"], fusions[fusion, "gene2"], fusions[fusion, "direction2"], fusions[fusion, "breakpoint2"], coverage2, fusions[fusion, "transcript_id2"], opt$transcriptSelection)
  if (nrow(exons1) == 0 || nrow(exons2) == 0) next

  if(unique(exons1$geneID) == unique(exons2$geneID)) next
  # sort coding exons last, such that they are drawn over the border of non-coding exons
  exons1 <- exons1[order(exons1$start, -rank(exons1$type)), ]
  exons2 <- exons2[order(exons2$start, -rank(exons2$type)), ]

  # insert dummy exons, if breakpoints are outside the gene (e.g., in UTRs)
  # this avoids plotting artifacts
  breakpoint1 <- fusions[fusion, "breakpoint1"]
  breakpoint2 <- fusions[fusion, "breakpoint2"]
  if (breakpoint1 < min(exons1$start)) {
    exons1 <- rbind(c(exons1[1, "contig"], "dummy", breakpoint1 - 1000, breakpoint1 - 1000, exons1[1, "strand"], "", "dummy", exons1[1, "geneID"], exons1[1, "transcript"], ""), exons1)
  } else if (breakpoint1 > max(exons1$end)) {
    exons1 <- rbind(exons1, c(exons1[1, "contig"], "dummy", breakpoint1 + 1000, breakpoint1 + 1000, exons1[1, "strand"], "", "dummy", exons1[1, "geneID"], exons1[1, "transcript"], ""))
  }
  if (breakpoint2 < min(exons2$start)) {
    exons2 <- rbind(c(exons2[1, "contig"], "dummy", breakpoint2 - 1000, breakpoint2 - 1000, exons2[1, "strand"], "", "dummy", exons2[1, "geneID"], exons2[1, "transcript"], ""), exons2)
  } else if (breakpoint2 > max(exons2$end)) {
    exons2 <- rbind(exons2, c(exons2[1, "contig"], "dummy", breakpoint2 + 1000, breakpoint2 + 1000, exons2[1, "strand"], "", "dummy", exons2[1, "geneID"], exons2[1, "transcript"], ""))
  }

  exons1$start <- as.integer(exons1$start)
  exons1$end <- as.integer(exons1$end)
  exons2$start <- as.integer(exons2$start)
  exons2$end <- as.integer(exons2$end)

  exons1$left <- exons1$start
  exons1$right <- exons1$end
  exons2$left <- exons2$start
  exons2$right <- exons2$end

  squishedIntronSize <- 200
  # hide introns in gene1
  cumulativeIntronLength <- 0
  previousExonEnd <- -squishedIntronSize
  for (exon in seq_len(nrow(exons1))) {
    if (breakpoint1 > previousExonEnd + 1 && breakpoint1 < exons1[exon, "left"])
        breakpoint1 <- (breakpoint1 - previousExonEnd) / (exons1[exon, "left"] - previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
    if (exons1[exon, "left"] > previousExonEnd) {
        cumulativeIntronLength <- cumulativeIntronLength + exons1[exon, "left"] - previousExonEnd - squishedIntronSize
        previousExonEnd <- exons1[exon, "right"]
      }
      if (breakpoint1 >= exons1[exon, "left"] && breakpoint1 <= exons1[exon, "right"] + 1)
        breakpoint1 <- breakpoint1 - cumulativeIntronLength
      exons1[exon, "left"] <- exons1[exon, "left"] - cumulativeIntronLength
      exons1[exon, "right"] <- exons1[exon, "right"] - cumulativeIntronLength
    }

    # hide introns in gene2
    cumulativeIntronLength <- 0
    previousExonEnd <- -squishedIntronSize
    for (exon in seq_len(nrow(exons2))) {
      if (breakpoint2 > previousExonEnd + 1 && breakpoint2 < exons2[exon, "left"])
        breakpoint2 <- (breakpoint2 - previousExonEnd) / (exons2[exon, "left"] - previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
      if (exons2[exon, "left"] > previousExonEnd) {
        cumulativeIntronLength <- cumulativeIntronLength + exons2[exon, "left"] - previousExonEnd - squishedIntronSize
        previousExonEnd <- exons2[exon, "right"]
      }
      if (breakpoint2 >= exons2[exon, "left"] && breakpoint2 <= exons2[exon, "right"] + 1)
        breakpoint2 <- breakpoint2 - cumulativeIntronLength
      exons2[exon, "left"] <- exons2[exon, "left"] - cumulativeIntronLength
      exons2[exon, "right"] <- exons2[exon, "right"] - cumulativeIntronLength
    }
  # scale exon sizes to 1
  scalingFactor <- max(exons1$right) + max(exons2$right)
  exons1$left <- exons1$left / scalingFactor
  exons1$right <- exons1$right / scalingFactor
  exons2$left <- exons2$left / scalingFactor
  exons2$right <- exons2$right / scalingFactor
  breakpoint1 <- breakpoint1 / scalingFactor
  breakpoint2 <- breakpoint2 / scalingFactor

  # shift gene2 to the right of gene1 with a little bit of padding
  gene2Offset <- max(exons1$right) + 0.05

  # center fusion horizontally
  fusionOffset1 <- (max(exons1$right) + gene2Offset) / 2 - ifelse(fusions[fusion, "direction1"] == "downstream", breakpoint1, max(exons1$right) - breakpoint1)
  fusionOffset2 <- fusionOffset1 + ifelse(fusions[fusion, "direction1"] == "downstream", breakpoint1, max(exons1$right) - breakpoint1)

  # layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
  layout(matrix(c(1, 1, 1, 2, 3, 4), 2, 3, byrow = TRUE), widths = c(0.9, 1.2, 0.9))
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, type = "l", xlim = c(-0.12, 1.12), ylim = c(0.4, 1.1), bty = "n", xaxt = "n", yaxt = "n")

  # vertical coordinates of layers
  yIdeograms <- 0.84
  yBreakpointLabels <- 0.76
  yCoverage <- 0.72
  yExons <- 0.67
  yGeneNames <- 0.58
  yFusion <- 0.5
  yTranscript <- 0.45
  yScale <- 0.407

  yTrajectoryBreakpointLabels <- yBreakpointLabels - 0.035
  yTrajectoryExonTop <- yExons + 0.03
  yTrajectoryExonBottom <- yExons - 0.055
  yTrajectoryFusion <- yFusion + 0.03

  # draw ideograms
  if (!is.null(cytobands)) {
    drawIdeogram("left", min(exons1$left), max(exons1$right), yIdeograms, cytobands, fusions[fusion, "contig1"], fusions[fusion, "breakpoint1"])
    drawIdeogram("right", gene2Offset, gene2Offset + max(exons2$right), yIdeograms, cytobands, fusions[fusion, "contig2"], fusions[fusion, "breakpoint2"])
  }

  # draw gene & transcript names
  text(max(exons1$right) / 2, yGeneNames, fusions[fusion, "gene1"], font = 2, cex = opt$fontSize, adj = c(0.5, 0))
  if (fusions[fusion,"site1"] != "intergenic")
    text(max(exons1$right)/2, yGeneNames-0.01, head(exons1$transcript,1), cex = 0.9 * opt$fontSize, adj = c(0.5, 1))
  text(gene2Offset+max(exons2$right) / 2, yGeneNames, fusions[fusion,"gene2"], font = 2, cex = opt$fontSize, adj = c(0.5, 0))
  if (fusions[fusion,"site2"] != "intergenic")
    text(gene2Offset+max(exons2$right) / 2, yGeneNames - 0.01, head(exons2$transcript,1), cex = 0.9 * opt$fontSize, adj = c(0.5, 1))

  # if multiple genes in the vicinity are shown, label them
  if (fusions[fusion,"site1"] == "intergenic")
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene & exons1$type != "dummy",]
      if (any(exonsOfGene$type == "exon"))
        text(mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yExons - 0.04, gene, cex = 0.9 * opt$fontSize, adj = c(0.5, 1))
    }
  if (fusions[fusion,"site2"] == "intergenic")
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene & exons2$type != "dummy",]
      if (any(exonsOfGene$type == "exon"))
        text(gene2Offset + mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yExons - 0.04, gene, cex = 0.9 * opt$fontSize, adj = c(0.5, 1))
    }

  # label breakpoints
  text(breakpoint1 + 0.01, yBreakpointLabels - 0.03, paste0("breakpoint\n", fusions[fusion, "display_contig1"], ":", fusions[fusion, "breakpoint1"]), adj = c(1, 0), cex = opt$fontSize)
  text(gene2Offset + breakpoint2 - 0.01, yBreakpointLabels - 0.03, paste0("breakpoint\n", fusions[fusion, "display_contig2"], ":", fusions[fusion, "breakpoint2"]), adj = c(0, 0), cex = opt$fontSize)

  # plot gene 1
  lines(c(min(exons1$left), max(exons1$right)), c(yExons, yExons), col = darkColor1)
  for (gene in unique(exons1$geneName))
    drawStrand(min(exons1[exons1$geneName == gene,"left"]), max(exons1[exons1$geneName == gene, "right"]), yExons, darkColor1, head(exons1[exons1$geneName == gene, "strand"], 1))
  for (exon in seq_len(nrow(exons1)))
    drawExon(exons1[exon, "left"], exons1[exon, "right"], yExons, opt$color1, exons1[exon, "exonNumber"], exons1[exon, "type"])

  # plot gene 2
  lines(c(gene2Offset, gene2Offset+max(exons2$right)), c(yExons, yExons), col = darkColor2)
  for (gene in unique(exons2$geneName))
    drawStrand(gene2Offset + min(exons2[exons2$geneName == gene,"left"]), gene2Offset + max(exons2[exons2$geneName == gene,"right"]), yExons, darkColor2, head(exons2[exons2$geneName == gene,"strand"], 1))
  for (exon in seq_len(nrow(exons2)))
    drawExon(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yExons, opt$color2, exons2[exon,"exonNumber"], exons2[exon,"type"])

  # plot gene1 of fusion
  if (fusions[fusion,"direction1"] == "downstream") {
    # plot strands
    lines(c(fusionOffset1, fusionOffset1 + breakpoint1), c(yFusion, yFusion), col = darkColor1)
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene,]
      if (min(exonsOfGene$start) <= fusions[fusion, "breakpoint1"])
        drawStrand(fusionOffset1 + min(exonsOfGene$left), fusionOffset1 + min(breakpoint1, max(exonsOfGene$right)), yFusion, col = darkColor1, exonsOfGene$strand[1])
    }
    # plot exons
    for (exon in seq_len(nrow(exons1)))
      if (exons1[exon,"start"] <= fusions[fusion,"breakpoint1"])
        drawExon(fusionOffset1 + exons1[exon, "left"], fusionOffset1 + min(breakpoint1, exons1[exon,"right"]), yFusion, opt$color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
    # plot trajectories
    lines(c(0, 0, fusionOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col = "red", lty = 2)
    lines(c(breakpoint1, breakpoint1, fusionOffset1 + breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col = "red", lty = 2)
  } else if (fusions[fusion, "direction1"] == "upstream")
    # plot strands
    lines(c(fusionOffset1, fusionOffset2), c(yFusion, yFusion), col = darkColor1)
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene, ]
      if (max(exonsOfGene$end + 1) >= fusions[fusion, "breakpoint1"])
        drawStrand(fusionOffset2 - max(exonsOfGene$right) + breakpoint1, min(fusionOffset2, fusionOffset2 - min(exonsOfGene$left) + breakpoint1), yFusion, col = darkColor1, chartr("+-", "-+", exonsOfGene$strand[1]))
    }
    # plot exons
    for (exon in seq_len(nrow(exons1))) {
      if (exons1[exon, "end"] + 1 >= fusions[fusion, "breakpoint1"])
        drawExon(fusionOffset1 + max(exons1$right) - exons1[exon, "right"], min(fusionOffset2, fusionOffset1 + max(exons1$right) - exons1[exon, "left"]), yFusion, opt$color1, exons1[exon, "exonNumber"], exons1[exon, "type"])
    # plot trajectories
    lines(c(max(exons1$right), max(exons1$right), fusionOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col = "red", lty = 2)
    lines(c(breakpoint1, breakpoint1, fusionOffset1 + max(exons1$right) - breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col = "red", lty = 2)
  }

   # plot gene2 of fusion
  if (fusions[fusion, "direction2"] == "downstream") {
    # plot strands
    lines(c(fusionOffset2, fusionOffset2 + breakpoint2), c(yFusion, yFusion), col = darkColor2)
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene,]
      if (min(exonsOfGene$start) <= fusions[fusion, "breakpoint2"])
        drawStrand(max(fusionOffset2, fusionOffset2 + breakpoint2 - max(exonsOfGene$right)), fusionOffset2 + breakpoint2 - min(exonsOfGene$left), yFusion, col = darkColor2, chartr("+-", "-+", exonsOfGene$strand[1]))
    }
    # plot exons
    for (exon in seq_len(nrow(exons2)))
      if (exons2[exon,"start"] <= fusions[fusion, "breakpoint2"])
        drawExon(max(fusionOffset2, fusionOffset2+breakpoint2 - exons2[exon, "right"]), fusionOffset2 + breakpoint2 - exons2[exon, "left"], yFusion, opt$color2, exons2[exon, "exonNumber"], exons2[exon, "type"])
    # plot trajectories
    lines(c(gene2Offset, gene2Offset, fusionOffset2+breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col = "red", lty = 2)
    lines(c(gene2Offset + breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col = "red", lty = 2)
  } else if (fusions[fusion, "direction2"] == "upstream") {
    # plot strands
    lines(c(fusionOffset2, fusionOffset2 + max(exons2$right) - breakpoint2), c(yFusion, yFusion), col = darkColor2)
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene,]
      if (max(exonsOfGene$end+1) >= fusions[fusion, "breakpoint2"])
        drawStrand(max(fusionOffset2, fusionOffset2 + min(exonsOfGene$left) - breakpoint2), fusionOffset2 + max(exonsOfGene$right) - breakpoint2, yFusion, col = darkColor2, exonsOfGene$strand[1])
    }
  }
    # plot exons
    for (exon in seq_len(nrow(exons2))) {
      if (exons2[exon,"end"] + 1 >= fusions[fusion, "breakpoint2"])
        drawExon(max(fusionOffset2, fusionOffset2 + exons2[exon, "left"] - breakpoint2), fusionOffset2 + exons2[exon, "right"] - breakpoint2, yFusion, opt$color2, exons2[exon, "exonNumber"], exons2[exon, "type"])
    # plot trajectories
    lines(c(gene2Offset + max(exons2$right), gene2Offset + max(exons2$right), fusionOffset2 + max(exons2$right) - breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col = "red", lty = 2)
    lines(c(gene2Offset + breakpoint2, gene2Offset + breakpoint2, fusionOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col= "red", lty = 2)
  }

  # draw scale
  realScale <- max(exons1$end - exons1$start, exons2$end - exons2$start)
  mapScale <- max(exons1$right - exons1$left, exons2$right - exons2$left)
  # choose scale which is closest to desired scale length
  desiredScaleSize <- 0.2
  realScale <- desiredScaleSize / mapScale * realScale
  mapScale <- desiredScaleSize
  realScaleOptimalFit <- signif(realScale, 1) # round to most significant digit
  mapScaleOptimalFit <- realScaleOptimalFit / realScale * mapScale
  # draw scale line
  lines(c(0, mapScaleOptimalFit), c(yScale, yScale)) # scale line
  lines(c(0, 0), c(yScale - 0.007, yScale + 0.007)) # left whisker
  lines(c(mapScaleOptimalFit, mapScaleOptimalFit), c(yScale - 0.007, yScale + 0.007)) # right whisker
  # draw units above scale line
  realScaleThousands <- max(0, min(3, floor(log10(realScaleOptimalFit) / 3)))
  scaleUnits <- c("bp", "kbp", "Mbp", "Gbp")
  scaleLabel <- paste(realScaleOptimalFit / max(1, 1000 ^ realScaleThousands), scaleUnits[realScaleThousands + 1])
  text(mapScaleOptimalFit / 2, yScale + 0.005, scaleLabel, adj = c(0.5, 0), cex = opt$fontSize * 0.9)
  text(mapScaleOptimalFit, yScale, "  introns not to scale", adj = c(0, 0.5), cex = opt$fontSize * 0.9, font = 3)

  # draw circos plot
  par(mar = c(0, 0, 0, 0))
  drawCircos(fusion, fusions, cytobands, opt$minConfidenceForCircosPlot)
  par(mar = c(0, 0, 0, 0))

  # draw protein domains
  plot(0, 0, type = "l", xlim = c(-0.1, 1.1), ylim = c(0, 1), bty = "n", xaxt = "n", yaxt = "n")
  par(xpd = NA)
  if (!is.null(proteinDomains))
    drawProteinDomains(fusions[fusion, ], exons1, exons2, proteinDomains, opt$color1, opt$color2, opt$mergeDomainsOverlappingBy, opt$optimizeDomainColors)
  par(xpd = FALSE)

  # print statistics about supporting alignments
  plot(0, 0, type = "l", xlim = c(0, 1), ylim = c(0, 1), bty = "n", xaxt = "n", yaxt = "n")
  text(0, 0.575, "SUPPORTING READ COUNT", font = 2, adj = c(0, 0), cex = opt$fontSize)
  text(0, 0.525, paste0("Total region count in ", unique(exons1$geneName[exons1$type != "dummy"]), " = ", fusions[fusion, "total_region_count1"], "\n",
                        "Total region count in ", unique(exons2$geneName[exons2$type != "dummy"]), " = ", fusions[fusion, "total_region_count2"], "\n",
                        "Assembly score", " = ", fusions[fusion, "assembly_score"], "\n",
                        "Fusion flag = ", fusions[fusion, "fusion_flag"]), adj = c(0, 1), cex = opt$fontSize)

  }
devNull <- dev.off()
message("Done")
### </MAIN> ###
