test_that("Null model analysis works", {
  nma_result <- NMA(DaTa = NMA_data, NSim = 10)
  expect_true(is.list(nma_result))
  expect_true(is.data.frame(nma_result$Standard_effect_size_data))
  expect_true(is.data.frame(nma_result$Hill_table))
  expect_true(is.list(NMA_plot(nma_result, plot_type = "SES")))
  expect_true(is.list(NMA_plot(nma_result, plot_type = "SES_species_number")))
})

test_that("NMA returns results with expected components", {
  data("NMA_data")
  nma_results <- NMA(DaTa = NMA_data, NSim = 10, InDeX = "invsimpson", Deb = TRUE)
  # Check that essential components are in the result
  expect_true("Input_data" %in% names(nma_results))
  expect_true("Standard_effect_size_data" %in% names(nma_results))
  expect_true("Hill_table" %in% names(nma_results))
})

test_that("NMA returns results with expected components", {
  data("NMA_data")
  nma_results <- NMA(DaTa = NMA_data, NSim = 10, ObSsIm = TRUE)
  # Check that essential components are in the result
  expect_true("Input_data" %in% names(nma_results))
  expect_true("Standard_effect_size_data" %in% names(nma_results))
  expect_true("Hill_table" %in% names(nma_results))
})

test_that("NMA returns results with expected components", {
  nma_results <- NMA(NSim = 10)
  # Check that essential components are in the result
  expect_true("Input_data" %in% names(nma_results))
  expect_true("Standard_effect_size_data" %in% names(nma_results))
  expect_true("Hill_table" %in% names(nma_results))
})

test_that("NMA errors on invalid input", {
  expect_warning(NMA(DaTa = c(1,2,3)))
})