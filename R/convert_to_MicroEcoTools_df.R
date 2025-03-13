#' Convert SummarizedExperiment or phyloseq Object to a MicroEcoTools-Compatible Data Frame
#'
#' This function takes in either a [SummarizedExperiment::SummarizedExperiment()] or a
#' [phyloseq::phyloseq()] object and converts it to a data frame formatted for MicroEcoTools analyses.
#' The returned data frame will have two metadata columns (by default, `Exp.Grp` and `Replicate`)
#' followed by the taxa/feature count data.
#'
#' @param x A SummarizedExperiment or phyloseq object.
#' @param group_col A character string indicating the column name in the sample metadata that
#'   represents the experimental group. Default is `"Exp.Grp"`.
#' @param replicate_col A character string indicating the column name in the sample metadata that
#'   represents replicates. Default is `"Replicate"`.
#' @param assay_name For SummarizedExperiment objects, the name of the assay to use. Default is `"counts"`.
#'
#' @return A data frame with the first two columns corresponding to the experimental group and replicate,
#'   and the remaining columns containing the taxa/feature counts.
#'
#' @examples
#' \dontrun{
#'   ## For a SummarizedExperiment object:
#'   library(SummarizedExperiment)
#'   # Suppose 'se_obj' is a SummarizedExperiment with colData containing columns "Exp.Grp" and "Replicate"
#'   df <- convert_to_microecotools_df(se_obj, group_col = "Exp.Grp", replicate_col = "Replicate")
#'
#'   ## For a phyloseq object:
#'   library(phyloseq)
#'   # Suppose 'ps_obj' is a phyloseq object with sample_data containing columns "Exp.Grp" and "Replicate"
#'   df <- convert_to_microecotools_df(ps_obj, group_col = "Exp.Grp", replicate_col = "Replicate")
#'   
#'     # Generate a toy OTU matrix.
#'   In phyloseq, rows can be taxa and columns samples (if taxa_are_rows is TRUE).
#'   otu_mat <- matrix(sample(1:50, 30, replace = TRUE), nrow = 5, ncol = 6)
#'   rownames(otu_mat) <- paste0("Taxa", 1:5)    # Taxa names
#'   colnames(otu_mat) <- paste0("Sample", 1:6)    # Sample names
#' 
#'   # Create the OTU table.
#'   # Set taxa_are_rows = TRUE if rows correspond to taxa.
#'   OTU <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
#' 
#'   # Create sample metadata as a data.frame.
#'   sample_metadata_ps <- data.frame(
#'   Exp.Grp = rep(c("Control", "Treatment"), each = 3),  # Two groups with 3 samples each
#'   Replicate = rep(1:3, times = 2)                       # Replicate numbers for each group
#'   )
#'   rownames(sample_metadata_ps) <- paste0("Sample", 1:6)    # Make sure rownames match column names in OTU
#' 
#'   # Create the sample_data object.
#'   SD <- phyloseq::sample_data(sample_metadata_ps)
#' 
#'   # Finally, create the phyloseq object using the OTU table and sample metadata.
#'   ps_toy <- phyloseq::phyloseq(OTU, SD)
#'   df_ps <- convert_to_MicroEcoTools_df(ps_toy,
#'                                      group_col = "Exp.Grp",
#'                                      replicate_col = "Replicate")
#'   print("Converted phyloseq object:")
#'   print(head(df_ps))
# 
#'  # Test the toy dataset with NMA function
#'  nma_result <- NMA(DaTa = df_ps, NSim = 10)
#'  }
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom phyloseq sample_data otu_table taxa_are_rows
#' @export
convert_to_MicroEcoTools_df <- function(x,
                                        group_col = "Exp.Grp",
                                        replicate_col = "Replicate",
                                        assay_name = "counts") {
  if (inherits(x, "SummarizedExperiment")) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      # Extract sample metadata
      sample_df <- as.data.frame(SummarizedExperiment::colData(x))
      if (!(group_col %in% colnames(sample_df))) {
        stop(
          sprintf(
            "group_col '%s' not found in colData of the SummarizedExperiment object.",
            group_col
          )
        )
      }
      if (!(replicate_col %in% colnames(sample_df))) {
        stop(
          sprintf(
            "replicate_col '%s' not found in colData of the SummarizedExperiment object.",
            replicate_col
          )
        )
      }
      # Extract the assay (e.g., counts)
      assay_mat <- SummarizedExperiment::assay(x, assay_name)
      # In SummarizedExperiment, columns are samples and rows are features.
      # We need a data frame where each row is a sample.
      counts_df <- as.data.frame(t(assay_mat))
      # Combine the selected metadata columns and the counts data
      final_df <- cbind(sample_df[, c(group_col, replicate_col), drop = FALSE], counts_df)
      rownames(final_df) <- NULL
      return(final_df)
    } else {
      warning(
        "Package 'SummarizedExperiment' is not installed; cannot convert input to MicroEcoTools data frame format."
      )
      return(invisible(NULL))
    }
  } else if (inherits(x, "phyloseq")) {
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      # Extract sample metadata from the phyloseq object
      sample_df <- as.data.frame(phyloseq::sample_data(x))
      if (!(group_col %in% colnames(sample_df))) {
        stop(
          sprintf(
            "group_col '%s' not found in sample_data of the phyloseq object.",
            group_col
          )
        )
      }
      if (!(replicate_col %in% colnames(sample_df))) {
        stop(
          sprintf(
            "replicate_col '%s' not found in sample_data of the phyloseq object.",
            replicate_col
          )
        )
      }
      # Extract the OTU (or feature) table
      otu <- phyloseq::otu_table(x)
      otu_mat <- as(otu, "matrix")
      # If taxa are stored as rows, transpose so that rows correspond to samples
      if (phyloseq::taxa_are_rows(x)) {
        otu_df <- as.data.frame(t(otu_mat))
      } else {
        otu_df <- as.data.frame(otu_mat)
      }
      # Reorder rows (samples) if needed so that the row names match between sample_df and otu_df
      common_samples <- intersect(rownames(sample_df), rownames(otu_df))
      sample_df <- sample_df[common_samples, , drop = FALSE]
      otu_df <- otu_df[common_samples, , drop = FALSE]
      # Combine the selected metadata columns and the OTU count data
      final_df <- cbind(sample_df[, c(group_col, replicate_col), drop = FALSE], otu_df)
      rownames(final_df) <- NULL
      return(final_df)
    } else {
      warning(
        "Package 'phyloseq' is not installed; cannot convert input to MicroEcoTools data frame format."
      )
      return(invisible(NULL))
    }
  } else {
    warning("Input object must be either a SummarizedExperiment or a phyloseq object.")
    return(invisible(NULL))
  }
}
