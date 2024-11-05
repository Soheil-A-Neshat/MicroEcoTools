test_that("Welch-ANOVA analysis works", {

  Welch_ANOVA(dAtA = CSR_IP2G_data, var.name = "taxa", p_adj = "BH", Parallel = TRUE)
  
})

test_that("Welch-ANOVA analysis works", {
  
  Welch_ANOVA(dAtA = CSR_IP2G_data, var.name = "taxa", p_adj = "BH", Parallel = FALSE)
  
})

test_that("Pairwise Welch-ANOVA analysis works", {
  pairwise_welch(dAtA = CSR_IP2G_data, Parallel = TRUE)
  
})

test_that("Pairwise Welch-ANOVA analysis works", {
  pairwise_welch(dAtA = CSR_TAXA_data, var.name = "taxa", p_adj = "BH", v.equal = FALSE, p.value.cutoff = 0.05, Parallel = FALSE)
  
})

test_that("CSR assign function works", {
  CSR_assign(dAtA = CSR_IP2G_data)
  
})

test_that("CSR assignment with simulation works", {
  CSR_Simulation(DaTa = CSR_IP2G_data, NSim = 3)
  
})