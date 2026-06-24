test_that(".remove_chr strips chr prefix", {
  expect_equal(BRASSVis:::.remove_chr("chr1"), "1")
  expect_equal(BRASSVis:::.remove_chr("chrX"), "X")
  expect_equal(BRASSVis:::.remove_chr("chrM"), "MT")
  expect_equal(BRASSVis:::.remove_chr("1"), "1")
})

test_that(".between works as expected", {
  expect_true(BRASSVis:::.between(5, 1, 10))
  expect_false(BRASSVis:::.between(15, 1, 10))
})

test_that(".get_dark_color returns valid color", {
  col <- BRASSVis:::.get_dark_color("#e5a5a5")
  expect_true(grepl("^#", col))
})
