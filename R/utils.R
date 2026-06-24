.remove_chr <- function(contig) {
  sub("^chr", "", sub("^chrM", "MT", contig), perl = TRUE)
}

.between <- function(value, start, end) {
  value >= start & value <= end
}

.change_color_brightness <- function(color, delta) {
  rgb(
    min(255, max(0, col2rgb(color)["red", ] + delta)),
    min(255, max(0, col2rgb(color)["green", ] + delta)),
    min(255, max(0, col2rgb(color)["blue", ] + delta)),
    maxColorValue = 255
  )
}

.get_dark_color <- function(color) {
  .change_color_brightness(color, -100)
}
