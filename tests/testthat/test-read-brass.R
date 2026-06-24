test_that("read_brass parses example BEDPE", {
  f <- system.file("extdata", "example.bedpe", package = "BRASSVis")
  skip_if(f == "")
  result <- read_brass(f)
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true(all(c("gene1", "gene2", "breakpoint1", "breakpoint2",
                     "fusion_flag") %in% names(result)))
})

test_that("read_brass errors on missing file", {
  expect_error(read_brass("nonexistent.bedpe"))
})
