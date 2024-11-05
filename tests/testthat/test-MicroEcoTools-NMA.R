test_that("Null model analysis works", {
  data("NMA_data")
  data("3CA_microcosm_Metagenomics_IP2G_traits_abundances")
  
  NMA(DaTa = NMA_data, NSim = 100)
  
})

test_that("Hill_diversity analysis works", {
  data("NMA_data")
  data("3CA_microcosm_Metagenomics_IP2G_traits_abundances")
  
  Hill_diversity(dAtA = CSR_TAXA_data[3:1509])
  
})
