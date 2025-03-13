test_that("Testing convert to MicroEcoTools df function", {
  # Generate a toy OTU matrix.
  # In phyloseq, rows can be taxa and columns samples (if taxa_are_rows is TRUE).
  otu_mat <- matrix(sample(1:50, 30, replace = TRUE), nrow = 5, ncol = 6)
  rownames(otu_mat) <- paste0("Taxa", 1:5)    # Taxa names
  colnames(otu_mat) <- paste0("Sample", 1:6)    # Sample names
  
  # Create the OTU table.
  # Set taxa_are_rows = TRUE if rows correspond to taxa.
  OTU <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  
  # Create sample metadata as a data.frame.
  sample_metadata_ps <- data.frame(
    Exp.Grp = rep(c("Control", "Treatment"), each = 3),  # Two groups with 3 samples each
    Replicate = rep(1:3, times = 2)                       # Replicate numbers for each group
  )
  rownames(sample_metadata_ps) <- paste0("Sample", 1:6)    # Make sure rownames match column names in OTU
  
  # Create the sample_data object.
  SD <- phyloseq::sample_data(sample_metadata_ps)
  
  # Finally, create the phyloseq object using the OTU table and sample metadata.
  ps_toy <- phyloseq::phyloseq(OTU, SD)
  expect_true(class(ps_toy) == "phyloseq")
  df_ps <- convert_to_MicroEcoTools_df(ps_toy,
                                       group_col = "Exp.Grp",
                                       replicate_col = "Replicate")
  print("Converted phyloseq object:")
  print(head(df_ps))
  
  expect_true(is.data.frame(df_ps))
  
  nma_result <- NMA(DaTa = df_ps, NSim = 10)
  expect_true(is.list(nma_result))
  })
