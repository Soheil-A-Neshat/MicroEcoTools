test_that("Hill_diversity works with count data (one column)", {
  data("NMA_data")
  result <- Hill_diversity(counts = NMA_data[1507], q = 2)
  # Check that the output is a data.frame and has expected columns.
  expect_true(is.numeric(result))
})

test_that("Hill_diversity works with count data (multiple column)", {
  data("NMA_data")
  result <- Hill_diversity(counts = NMA_data[4:1507], q = 2)
  # Check that the output is a data.frame and has expected columns.
  expect_true(is.numeric(result))
})