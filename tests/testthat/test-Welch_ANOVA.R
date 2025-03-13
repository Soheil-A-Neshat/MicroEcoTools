test_that("Welch_ANOVA works with valid input (parallel)", {
  data("CSR_IP2G_data")
  result <- Welch_ANOVA(dAtA = CSR_IP2G_data, var.name = "Traits", p_adj = "BH", Parallel = TRUE)
  # Check that the output is a data.frame and has expected columns.
  expect_true(is.data.frame(result))
  expect_true(all(c("Traits", "Welch-ANOVA p-value", "Method", "Remarks_WA") %in% colnames(result)))
})

test_that("Welch_ANOVA works with valid input (sequential)", {
  data("CSR_IP2G_data")
  result <- Welch_ANOVA(dAtA = CSR_IP2G_data, var.name = "Traits", p_adj = "BH", Parallel = FALSE)
  expect_true(is.data.frame(result))
})

test_that("Welch_ANOVA errors with missing data", {
  expect_warning (Welch_ANOVA())
})
