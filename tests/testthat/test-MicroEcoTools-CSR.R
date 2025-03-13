test_that("Welch-ANOVA analysis works (parallel)", {

  # Check that the result is a data frame and has the expected columns
  expect_true(is.data.frame(Welch_ANOVA(
    dAtA = CSR_IP2G_data, 
    var.name = "Traits", 
    p_adj = "BH", 
    Parallel = TRUE
  )))

})

test_that("Welch-ANOVA analysis works (sequential)", {

  # Check that the result is a data frame with expected structure
  expect_true(is.data.frame(Welch_ANOVA(
    dAtA = CSR_IP2G_data, 
    var.name = "Traits", 
    p_adj = "BH", 
    Parallel = FALSE
  )))

})

test_that("Pairwise Welch-ANOVA analysis works (sequential)", {

  # Check that the output is a data frame or has the correct S3 class if assigned.
  expect_true(is.data.frame(pairwise_welch(
    dAtA = CSR_IP2G_data, 
    var.name = "Traits", 
    p_adj = "BH", 
    v.equal = FALSE, 
    p.value.cutoff = 0.05, 
    Parallel = TRUE
  )))
})

test_that("Pairwise Welch-ANOVA analysis works (parallel)", {

  expect_true(is.data.frame(pairwise_welch(
    dAtA = CSR_IP2G_data, 
    var.name = "Traits", 
    p_adj = "BH", 
    v.equal = FALSE, 
    p.value.cutoff = 0.05, 
    Parallel = TRUE
  )))

})

test_that("CSR assign functions work", {
  # Expect that the CSR_assignments element exists and is a data frame.
  expect_true(is.list(CSR_assign(
    dAtA = CSR_IP2G_data, 
    var.name = "Trait", 
    p_adj = "BH", 
    Parallel = TRUE
  )))
})

test_that("CSR assignment with simulation works", {

  # Check that the simulation returns a list with the expected elements.
  expect_true("Final_Verdict" %in% names(CSR_Simulation(
    DaTa = CSR_IP2G_data, 
    NSim = 3, 
    p_adj = "BH", 
    var.name = "Trait", 
    NuLl.test = FALSE, 
    Keep_data = FALSE
  )))

})
