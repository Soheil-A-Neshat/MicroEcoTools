test_that("CSR_assign returns CSR_assignments with proper class", {
  data("CSR_IP2G_data")
  csr_results <- CSR_assign(dAtA = CSR_IP2G_data, var.name = "Trait", p_adj = "BH", Parallel = TRUE)
  # Check that the returned list contains the CSR_assignments element and that it has the right class.
  expect_true(is.list(csr_results))
})

test_that("CSR_assign returns CSR_assignments with proper class", {
  data("CSR_IP2G_data")
  csr_results <- CSR_assign(dAtA = CSR_IP2G_data, var.name = "Trait", p_adj = "BH", Parallel = TRUE, v.equal = TRUE)
  # Check that the returned list contains the CSR_assignments element and that it has the right class.
  expect_true(is.list(csr_results))
})

test_that("CSR_assign returns CSR_assignments with proper class", {
  data("CSR_IP2G_data")
  csr_results <- CSR_assign(dAtA = CSR_IP2G_data, var.name = "Trait", p_adj = "BH", Parallel = FALSE)
  # Check that the returned list contains the CSR_assignments element and that it has the right class.
  expect_true(is.list(csr_results))
})

test_that("CSR_assign handles missing data gracefully", {
  expect_warning(CSR_assign())
})

# Optionally, test specific branches by setting Verbose = "all"
test_that("CSR_assign functions running under verbose mode", {
  data("CSR_IP2G_data")
  csr_results <- CSR_assign(dAtA = CSR_IP2G_data, var.name = "Trait", p_adj = "BH", Parallel = TRUE, Verbose = TRUE)
  expect_true (is.list(csr_results))
})
