test_that("CSR_Simulation runs and returns expected output", {
  data("CSR_IP2G_data")
  # Use a low number of simulations for testing speed.
  sim_results <- CSR_Simulation(DaTa = CSR_IP2G_data, NSim = 3, p_adj = "BH", var.name = "Trait", NuLl.test = FALSE, Keep_data = FALSE, Parallel = TRUE)
  # Check that the final verdict has the proper class
  expect_true(inherits(sim_results$Summary, "Final_verdict"))
  # Check that key elements are present in the returned list
  expect_true("data" %in% names(sim_results))
  expect_true("CSR" %in% names(sim_results))
})

test_that("CSR_Simulation errors with missing data", {
  expect_warning(CSR_Simulation())
})
