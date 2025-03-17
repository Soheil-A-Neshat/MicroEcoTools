#' @title Welch-ANOVA comparison (MicroEcoTools)
#'
#' @description Statistical comparison for traits/taxa between experimental groups using Welch-ANOVA (or ANOVA when variance assumptions fail).
#'              This function uses a oneway test with the assumption of non-equal variances and automatically falls back to the equal-variance test if needed.
#'
#' @param dAtA Data frame containing the experimental group, replicates, and measurement columns.
#' @param var.name Character string specifying the variable name to be used in the output table. Defaults to \code{"Variable"}.
#' @param p_adj Character string specifying the p-value adjustment method (as per \code{\link[stats]{p.adjust}}). Defaults to \code{"BH"}.
#' @param Parallel Logical value. If TRUE, the function uses parallel processing (80% of available CPU cores) to speed up calculations.
#'                 Defaults to \code{TRUE}.
#'
#' @return A data frame (table) containing the following columns:
#'         - The name of each trait/taxa,
#'         - Its Welch-ANOVA p-value,
#'         - The test method used,
#'         - A remarks column indicating if the fallback equal variance test was applied,
#'         - An adjusted p-value column.
#'
#' @details
#' The input data sample is expected to be structured as follows:
#'
#' | Exp.Grp | Replicate | TAXA1 | TAXA2 | TAXA3 |
#' |---------|:---------:|------:|------:|------:|
#' | X0      | 1         | 10    | 91    | 68    |
#' | X0      | 2         | 11    | 86    | 70    |
#' | X0      | 3         | 12    | 84    | 70    |
#' | X1      | 1         | 15    | 30    | 3452  |
#' | X1      | 2         | 19    | 45    | 3274  |
#' | X1      | 3         | 13    | 25    | 3601  |
#'
#' @examples
#' \dontrun{
#'   # Using the CSR_IP2G_data dataset included in the package to generate a comparison table.
#'   wa_results <- Welch_ANOVA(dAtA = CSR_IP2G_data, var.name = "Traits", p_adj = "BH", Parallel = TRUE)
#' }
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom progress progress_bar
#' @importFrom stats oneway.test p.adjust
#' @export
Welch_ANOVA <- function(dAtA, var.name = "Variable", p_adj = "BH", Parallel = TRUE) {
  
  # Check if the input data is missing
  if (missing(dAtA)) {
    warning("No data input!")
    return(invisible(NULL))
  } else {
    # -------------------------------------------------------------------------
    # 1. Convert Input Data if Required
    # -------------------------------------------------------------------------
    # If the input is a SummarizedExperiment or phyloseq object, convert it to a compatible data frame.
    if (inherits(dAtA, "SummarizedExperiment") || inherits(dAtA, "phyloseq")) {
      message("Detected SummarizedExperiment/phyloseq input. Converting to MicroEcoTools data frame format...")
      dAtA <- convert_to_microecotools_df(
        dAtA,
        group_col = "Exp.Grp",    # Change if your metadata column name differs
        replicate_col = "Replicate",
        assay_name = "counts"     # Adjust if necessary
      )
    }
    
    # -------------------------------------------------------------------------
    # 2. Set Up Parameters for Analysis
    # -------------------------------------------------------------------------
    # Determine the number of iterations based on the number of measurement columns.
    total_iterations <- length(dAtA[1,]) - 2
    
    # -------------------------------------------------------------------------
    # 3. Define a Helper Function for Parallel Processing
    # -------------------------------------------------------------------------
    # This function performs Welch's ANOVA for each measurement column using parallel processing.
    run_parallel_welch_anova <- function(dAtA, total_iterations) {
      # Initialize an empty data frame to store results.
      a <- data.frame()
      
      # Function to determine the number of CPU cores to use.
      get_cores <- function() {
        if (as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE"))) {
          return(2)
        } else {
          ncores <- max(1, floor(0.8 * parallel::detectCores()))
          return(ncores)
        }
      }
      ncores <- get_cores()
      cl <- parallel::makeCluster(ncores)
      
      # Register the cluster for parallel processing.
      doParallel::registerDoParallel(cl)
      
      # Ensure the cluster is stopped and reset to sequential processing on exit.
      on.exit({
        suppressWarnings(parallel::stopCluster(cl))
        foreach::registerDoSEQ()
        #print("Closing the parallel tasks ...")
      }, add = TRUE)
      
      #message(paste0("Parallel Welch-ANOVA comparisons using ", ncores, " CPU cores. Please be patient..."))
      
      # Perform the oneway tests in parallel.
      a_result <- foreach::foreach(i = 1:total_iterations, .combine = rbind) %dopar% {
        # Create a placeholder for the result.
        result <- data.frame(matrix(nrow = 1, ncol = 4))
        # Dynamically construct the formula: "Taxon ~ Exp.Grp"
        f <- as.formula(noquote(paste0("`", names(dAtA)[i + 2], "` ~ ", paste(names(dAtA)[1], collapse = " + "))))
        
        # Save the taxon name (or trait name) as the first column.
        result[1, 1] <- names(dAtA)[i + 2]
        
        # Perform the test with non-equal variance assumption.
        test_result <- oneway.test(f, data = dAtA, var.equal = FALSE)
        if (is.na(test_result$p.value)) {
          # If p-value is NA, perform the test assuming equal variance.
          test_result <- oneway.test(f, data = dAtA, var.equal = TRUE)
          result[1, 2] <- test_result$p.value
          result[1, 3] <- test_result$method
          result[1, 4] <- "*"  # Mark with an asterisk if fallback test is used.
        } else {
          result[1, 2] <- test_result$p.value
          result[1, 3] <- test_result$method
          result[1, 4] <- "-"  # No fallback needed.
        }
        
        # Simulate a slight delay (can be removed if not needed).
        Sys.sleep(1 / total_iterations)
        
        # Return the result row.
        return(result)
      }
      
      return(a_result)
    }
    
    # -------------------------------------------------------------------------
    # 4. Run Welch's ANOVA (Parallel or Sequential)
    # -------------------------------------------------------------------------
    if (Parallel) {
      a <- run_parallel_welch_anova(dAtA = dAtA, total_iterations = total_iterations)
    } else {
      # Sequential processing with a progress bar.
      a <- data.frame()
      pb <- progress::progress_bar$new(format = "Welch ANOVA [:bar] :percent in :elapsed", clear = FALSE, total = total_iterations)
      for (i in 1:total_iterations) {
        f <- as.formula(noquote(paste0("`", names(dAtA)[i + 2], "` ~ ", paste(names(dAtA)[1], collapse = " + "))))
        a[i, 1] <- names(dAtA)[i + 2]
        test_result <- oneway.test(f, data = dAtA, var.equal = FALSE)
        if (is.na(test_result$p.value)) {
          test_result <- oneway.test(f, data = dAtA, var.equal = TRUE)
          a[i, 2] <- test_result$p.value
          a[i, 3] <- test_result$method
          a[i, 4] <- "*"
        } else {
          a[i, 2] <- test_result$p.value
          a[i, 3] <- test_result$method
          a[i, 4] <- "-"
        }
        Sys.sleep(1 / total_iterations)
        pb$tick()
      }
      pb$terminate()
    }
    
    # -------------------------------------------------------------------------
    # 5. Adjust p-values and Prepare Output Table
    # -------------------------------------------------------------------------
    colnames(a) <- c(var.name, "Welch-ANOVA p-value", "Method", "Remarks_WA")
    a$`adjusted Welch-ANOVA p-value` <- p.adjust(a$`Welch-ANOVA p-value`, method = p_adj)
    
    # Return the final results table.
    return(a)
  }
}


#' Summary method for Welch_Anova objects
#' @param object An object of class "Welch_anova".
#'
#' @return Invisibly returns a table summarizing the Welch Anova results.
#' @examples
#' \dontrun{
#'   WA_result <- Welch_ANOVA(dAtA = CSR_IP2G_data, var.name = "Traits", p_adj = "BH", Parallel = TRUE)
#'   summary(WA_result)
#' }
#'
#' @export
summary.Welch_ANOVA <- function(object) {
  cat("Summary of Welch-ANOVA Results:\n")
  # Prints the summary of Welch-ANOVA analysis
  print(object[,c(-3,-4)])
  invisible(object)
}

#' @title Pairwise Welch-ANOVA (MicroEcoTools)
#'
#' @description Pairwise comparison between experimental groups using Welch-ANOVA or ANOVA. This function uses oneway test assuming non-equal variance to perform pairwise comparison.
#' In cases where the assumptions are not valid it chooses equal variance by automatically setting var.equal parameter to TRUE it can be forced to use equal/not equal variance by setting v.equal to TRUE or FALSE.
#'
#' @param dAtA Dataframe that contains experimental group, replicates, and parameters.
#' @note
#' The input data can be relative abundance or count data. The user must ensure the total count across groups is similar.
#' @param var.name Variable name used for generating output table, for example, taxa, trait, function. Defaults to \code{"Variable"}.
#' @param p_adj P-value adjustment method. All methods mentioned in \code{\link[stats]{p.adjust}} function can be used ("BH", "BY", "holm", "hommel" or "none"). Defaults to \code{"BH"}.
#' @param v.equal Assumption of equal variance can be forced using this parameter. By default it is set to FALSE where it cannot be set to FALSE it will switch to TRUE automatically; a note will be shown in the remarks column in those cases.Defaults to \code{FALSE}.
#' @param p.value.cutoff A cut-off value for calling the P-value significant can be set using this parameter. Defaults to \code{0.05}.
#' @param Parallel Using this parameter you can make use of more CPU cores to decrease run time for the calculation. By default it is set to TRUE and uses 80% of the available CPU cores to perform the calculations.Defaults to \code{TRUE}
#' @return This function returns a table containing the pairwise statistical comparison results. The output table can be fed into the CSR_assign function to assign CSR categories.
#' @details
#'The input data sample:
#'
#' | Exp.Grp | Replicate | TAXA1 | TAXA2 | TAXA3 |
#'|-----------|:-----------:|-----------:|:-----------:|-----------:|
#'  | X0 | 1 | 10 | 91 | 68 |
#'  | X0 | 2 | 11 | 86 | 70 |
#'  | X0 | 3 | 12 | 84 | 70 |
#'  | X1 | 1 | 15 | 30 | 3452 |
#'  | X1 | 2 | 19 | 45 | 3274 |
#'  | X1 | 3 | 13 | 25 | 3601 |
#' @md
#'
#' @examples
#' # Generate a pairwise comparison table for traits across experimental groups.
#' # This example uses abundance data from a perturbation experiment described in Santillan et al. (2019).
#' # The dataset \code{\link[MicroEcoTools]{CSR_IP2G_data}} is included in the MicroEcoTools package.
#' # Running the command below produces a pairwise comparison table for IP2G functions.
#' # Setting the Parallel to TRUE forces the function to use all available cpu cores to speed up the process.
#' \dontrun{
#' pwa_reults <- pairwise_welch(dAtA = CSR_IP2G_data, var.name = "Traits", p_adj = "BH", v.equal = FALSE, p.value.cutoff = 0.05, Parallel = TRUE)
#' }
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom progress progress_bar
#' @importFrom effectsize hedges_g interpret_hedges_g
#' @importFrom stats oneway.test as.formula sd p.adjust
#' @importFrom utils combn
#' @export
pairwise_welch <- function(dAtA,
                           var.name = "Variable",
                           p_adj = "BH",
                           v.equal = FALSE,
                           p.value.cutoff = 0.05,
                           Parallel = TRUE) {
# ---------------------------------------------------------------------------
# 1. Input Validation and Conversion
# ---------------------------------------------------------------------------
  
  if (missing(dAtA))
    print("No data input!")
  else {
    # Check if the input object is a SummarizedExperiment or phyloseq object.
    # If so, convert it to a MicroEcoTools-compatible data frame.
    if (inherits(dAtA, "SummarizedExperiment") ||
        inherits(dAtA, "phyloseq")) {
      message(
        "Detected SummarizedExperiment/phyloseq input. Converting to MicroEcoTools data frame format..."
      )
      dAtA <- convert_to_microecotools_df(
        dAtA,
        group_col = "Exp.Grp",
        # change if your metadata column name is different
        replicate_col = "Replicate",
        assay_name = "counts"
      )    # adjust if needed for SummarizedExperiment objects
    }
    
# -------------------------------------------------------------------------
# 2. Data Preprocessing
# -------------------------------------------------------------------------
    # If data appears to be in relative abundance format (small counts), convert to count data.
    
    if (rowSums(dAtA[1, c(-1, -2)] <= 1) |
        rowSums(dAtA[1, c(-1, -2)]) <= 100) {
      message(
        "  Relative abundance data detected. Converting to count data by multiplying the relative abundances by 100,000 ..."
      )
    # Remove the second column (assumed to be 'Replicate') as it is not needed in this analysis.
      dAtA[, c(-1, -2)] <- dAtA[, c(-1, -2)] * 100000
    }
    dAtA <- dAtA[, -2]
    
    # Determine the number of unique pairwise comparisons (n) among the experimental groups.
    n <- (factorial(length(unique(dAtA[, 1])))) / (2 * factorial(length(unique(dAtA[, 1])) - 2))
    
    # Generate all pairwise combinations of experimental group names.
    c <- t(combn(unique(dAtA[, 1]), 2))
    
    # Create an initial results data frame 'a' containing:
    # - The variable (trait/taxon) names repeated for each pairwise comparison.
    a <- data.frame()
    a <- rep(names(dAtA)[2:length(dAtA[1, ])], each = n)
    a <- as.data.frame(a)
    # Bind the first group name from each pair.
    a <- cbind(a, unlist(rep(as.character(
      c[, 1], length(names(dAtA)[2:length(dAtA[1, ])])
    ))))
    # Bind the second group name from each pair.
    a <- cbind(a, unlist(rep(as.character(
      c[, 2], length(names(dAtA)[2:length(dAtA[1, ])])
    ))))
# -------------------------------------------------------------------------
# 3. Determine Number of Cores for Parallel Processing
# -------------------------------------------------------------------------
    
    get_cores <- function() {
      if (as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE"))) {
        return(2)
      } else {
        ncores <- max(1, floor(0.8 * parallel::detectCores()))
        return(ncores)
      }
    }
# -------------------------------------------------------------------------
# 4. Parallel Processing Block
# -------------------------------------------------------------------------
    
    if (Parallel) {
      ncores <- get_cores()
      message(
        paste(
          "  Pairwise comparison of variables with variable name",
          var.name,
          ", P-value cutoff of" ,
          p.value.cutoff,
          "in parallel mode using ",
          ncores,
          "cores."
        )
      )
      
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      #on.exit({parallel::stopCluster(cl) foreach::registerDoSEQ()}, add = TRUE)
      
      
      # Create a progress bar
      
      message("  Parallel processing of pairwise comparisons step 1/4. Please be patient...")
      
      # Run the pairwise comparisons in parallel.
      a_result <- foreach::foreach (i = 1:length(a[, 1]), .combine = rbind) %dopar% {
        if (!v.equal) {
          if (sd(dAtA[dAtA[1] == a[i, 2], a[i, 1]]) == 0 |
              sd(dAtA[dAtA[1] == a[i, 3], a[i, 1]]) == 0)
            v.equal <- TRUE
        }
        # Compute the difference in means between the two groups.
        a[i, 4] <- mean(dAtA[dAtA[1] == a[i, 2], a[i, 1]]) - mean(dAtA[dAtA[1] == a[i, 3], a[i, 1]])
        
        # Build the formula for the oneway test.
        f <- as.formula(noquote(paste0(
          "`", a[i, 1], "`", "~", paste(names(dAtA)[1], collapse = " + ")
        )))
        
        a[i, 5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = v.equal)$statistic
        a[i, 6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = v.equal)$parameter[1]
        a[i, 7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = v.equal)$parameter[2]
        a[i, 8] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = v.equal)$p.value
        a[i, 9] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = v.equal)$method
        
        # If the p-value is not estimable, rerun using equal variance.
        if (is.na(a[i, 5])) {
          a[i, 5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = TRUE)$statistic
          a[i, 6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = TRUE)$parameter[1]
          a[i, 7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = TRUE)$parameter[2]
          a[i, 8] <- NA
          a[i, 9] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA)[1], a[i, 1])], var.equal = TRUE)$method
        }
        
        # Set column 10 to the p-value for later p-value adjustment.
        a[i, 10] <- a[i, 8]
        
        return(a[i, , drop = FALSE])
      }
      
      # Stop the parallel cluster and revert to sequential processing.
      suppressWarnings(parallel::stopCluster(cl))
      foreach::registerDoSEQ()
      
      # Save the parallel result.
      a <- a_result
# -----------------------------------------------------------------------
# 5. P-value Adjustment for Parallel Processing
# -----------------------------------------------------------------------
      
      message("  Pairwise Welch ANOVA P-value adjustment step 2/4")
      for (i in unlist(unique(a[1]))) {
        a[a[1] == i , 10] <- p.adjust(a[a[1] == i , 8], method = p_adj)
      }
      
# -----------------------------------------------------------------------
# 6. Calculate Effect Sizes and Significance Markers (Parallel)
# -----------------------------------------------------------------------
      
      total_iterations <- length(a[, 1])
      pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA effect size calculation step 3/4 [:bar] :percent in :elapsed",
                                       clear = FALSE,
                                       total = total_iterations)
      for (i in 1:length(a[, 1])) {
        if (!is.na(a[i, 8])) {
          if (a[i, 10] < p.value.cutoff)
            a[i, 11] <- "*"
          else
            a[i, 11] <- "ns"
        } else
          a[i, 11] <- "ns"
        pb$tick()
      }
      pb$terminate()
      

# -----------------------------------------------------------------------
# 7. Interpret Effect Sizes using Hedges' g (Parallel)
# -----------------------------------------------------------------------
      
      message("  Pairwise Welch ANOVA effect size interpretation step 4/4. Please be patient...")
      
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      a <- foreach::foreach (i = 1:length(a[, 1]), .combine = rbind) %dopar% {
        tryCatch({
          a[i, 12] <- effectsize::interpret_hedges_g(effectsize::hedges_g(dAtA[dAtA[1] == a[i, 2], a[i, 1]], dAtA[dAtA[1] == a[i, 3], a[i, 1]], pooled_sd = FALSE)$Hedges_g)
        }, error = function(e) {
          a[i, 12] <- "ne"
        })
        a[i, 12][is.na(a[i, 12])] <- "ne"
        
        return (a[i, , drop = FALSE])
      }
      
      suppressWarnings(parallel::stopCluster(cl))
      foreach::registerDoSEQ()
      
    } else {
      message(
        paste(
          "Pairwise comparison of variables with variable name",
          var.name,
          ", P-value cutoff of" ,
          p.value.cutoff,
          "using a single core. Please consider using parallel mode by setting Parallel parameter to TRUE: Parallel=TRUE"
        )
      )
      
# ---------------------------------------------------------------------
# 8. Sequential Processing Block (if Parallel = FALSE)
# ---------------------------------------------------------------------
      
      total_iterations <- length(a[, 1])
      pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA 1/4 [:bar] :percent in :elapsed",
                                       clear = FALSE,
                                       total = total_iterations)
      
      for (i in 1:length(a[, 1])) {
        if (!v.equal) {
          if (sd(dAtA[dAtA[1] == a[i, 2], a[i, 1]]) == 0 |
              sd(dAtA[dAtA[1] == a[i, 3], a[i, 1]]) == 0)
            v.equal <- TRUE
        }
        
        # Calculate difference in means.
        a[i, 4] <- mean(dAtA[dAtA[1] == a[i, 2], a[i, 1]]) - mean(dAtA[dAtA[1] == a[i, 3], a[i, 1]])
        f <- as.formula(noquote(paste0(
          "`", a[i, 1], "`", "~", paste(names(dAtA)[1], collapse = " + ")
        )))
        
        # Run oneway test.
        a[i, 5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = v.equal)$statistic
        a[i, 6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = v.equal)$parameter[1]
        a[i, 7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = v.equal)$parameter[2]
        a[i, 8] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = v.equal)$p.value
        a[i, 9] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = v.equal)$method
        # Fallback to equal variance test.
        if (is.na(a[i, 5])) {
          a[i, 5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = TRUE)$statistic
          a[i, 6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = TRUE)$parameter[1]
          a[i, 7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = TRUE)$parameter[2]
          a[i, 8] <- NA
          a[i, 9] <- oneway.test(f, data = dAtA[dAtA[1] == a[i, 2] |
                                                  dAtA[1] == a[i, 3], c(names(dAtA[1]), a[i, 1])], var.equal = TRUE)$method
        }
        
        pb$tick()
      }
      pb$terminate()
      
      # Adjust p-values for sequential mode. 
      a[, 10] <- a[, 8]
      
      
      message("  Pairwise Welch ANOVA P-value adjustment 2/4")
      for (i in unlist(unique(a[1]))) {
        a[a[1] == i , 10] <- p.adjust(a[a[1] == i , 8], method = p_adj)
      }
      
      # Compute significance markers based on the p-value cutoff.
      total_iterations <- length(a[, 1])
      pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA effect size calculation 3/4 [:bar] :percent in :elapsed",
                                       clear = FALSE,
                                       total = total_iterations)
      for (i in 1:length(a[, 1])) {
        if (!is.na(a[i, 8])) {
          if (a[i, 10] < 0.05)
            a[i, 11] <- "*"
          else
            a[i, 11] <- "ns"
        } else
          a[i, 11] <- "ns"
        pb$tick()
      }
      pb$terminate()
      
      # Interpret effect sizes using Hedges' g in sequential mode.
      total_iterations <- length(a[, 1])
      pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA effect size interpretation 4/4 [:bar] :percent in :elapsed",
                                       clear = FALSE,
                                       total = total_iterations)
      for (i in 1:length(a[, 1])) {
        tryCatch({
          a[i, 12] <- effectsize::interpret_hedges_g(effectsize::hedges_g(dAtA[dAtA[1] == a[i, 2], a[i, 1]], dAtA[dAtA[1] == a[i, 3], a[i, 1]], pooled_sd = FALSE)$Hedges_g)
        }, error = function(e) {
          a[i, 12] <- "ne"
        })
        a[i, 12][is.na(a[i, 12])] <- "ne"
        pb$tick()
      }
      pb$terminate()
      
# -------------------------------------------------------------------------
# 9. Finalize Output Table
# -------------------------------------------------------------------------
      
      colnames(a) <- c(
        var.name,
        "P1",
        "P2",
        "Distance",
        "Statistic",
        "dF1",
        "dF2",
        "P-value",
        "Method",
        "Adj.P-value",
        "Significance",
        "Effect_size"
      )
      #print(Sys.time())
      return(a)
    }
    
    colnames(a) <- c(
      var.name,
      "P1",
      "P2",
      "Distance",
      "Statistic",
      "dF1",
      "dF2",
      "P-value",
      "Method",
      "Adj.P-value",
      "Significance",
      "Effect_size"
    )
    
    # Return the final pairwise comparison results table.
    return(a)
  }
}

#' Summary method for Pairwise_Welch objects
#' @param object An object of class "Pariwise_Welch_ANOVA".
#'
#' @return Invisibly returns a table summarizing the Pairwise Welch Anova results.
#' @examples
#' \dontrun{
#' pw_results <- pairwise_welch(dAtA = CSR_IP2G_data, var.name = "Traits", p_adj = "BH", 
#'                       v.equal = FALSE, p.value.cutoff = 0.05, Parallel = TRUE)
#' summary(pw_results)
#' }
#'
#' @export
summary.Pairwise_welch <- function(object) {
  cat(
    "Summary of Pairwise Welch-ANOVA Results. Comparing a trait actross two groups shown under P1 and P2 columns:\n"
  )
  # Prints the summary of Welch-ANOVA analysis
  print(object[, c(-4, -5, -6, -7)])
  invisible(object)
}



#' @title Competitor - Stress tolerant - Ruderal assignment (MicroEcoTools)
#' @description Assigning CSR categories from Grime's framework to the input variable.
#' This function assigns CSR categories to the input variables using the output table from the \code{\link[MicroEcoTools]{pairwise_welch}} function.
#' To do so, it looks at the pairwise comparisons and assigns C to the ones abundant in no disturbance group, S to the ones abundant in high disturbance group, and R to the ones abundant in intermediate disturbance groups.
#' It also assigns intermediate categories, CR, CS, and SR to groups that are abundant in no-disturbance and intermediate-disturbance, no-disturbance and high-disturbance, and intermediate-disturbance and high-disturbance, respectively.
#' Traits that are not differentially abundant in any groups will be assigned to the CSR category.
#' Please note that if the variability is lacking the group will be categorised as NA.
#' @param dAtA Dataframe that contains experimental group, replicates, and parameters.
#' @note
#' The input data can be relative abundance or count data. The user must ensure the total count across groups is similar.
#' @param var.name Variable name used for generating output table, for example, taxa, trait, function. Defaults to \code{"Trait"}.
#' @param p_adj P-value adjustment method. All methods mentioned in \code{\link[stats]{p.adjust}} function can be used ("BH", "BY", "holm", "hommel" or "none"). Defaults to \code{"BH"}.
#' @param v.equal Assumption of equal variance can be forced using this parameter. By default it is set to FALSE where it cannot be set to FALSE it will switch to TRUE automatically; a note will be shown in the remarks column in those cases. Defaults to \code{FALSE}.
#' @param p.value.cutoff A cut-off value for calling the P-value significant can be set using this parameter.  Defaults to \code{0.05}.
#' @param Parallel Using this parameter you can make use of more CPU cores to decrease run time for the calculation. By default it is set to TRUE and uses half of the available CPU cores to perform the calculations. Defaults to \code{TRUE}.
#' @param Verbose Run CSR_assign function in verbose mode to explain the assignment. Defaults to \code{FALSE}.
#' @return This function returns a table containing the pairwise statistical comparison results. The output table can be fed into the CSR_assign function to assign CSR categories.
#' 
#'@details
#'The input data sample:
#'
#' | Exp.Grp | Replicate | TAXA1 | TAXA2 | TAXA3 |
#'|-----------|:-----------:|-----------:|:-----------:|-----------:|
#'  | X0 | 1 | 10 | 91 | 68 |
#'  | X0 | 2 | 11 | 86 | 70 |
#'  | X0 | 3 | 12 | 84 | 70 |
#'  | X1 | 1 | 15 | 30 | 3452 |
#'  | X1 | 2 | 19 | 45 | 3274 |
#'  | X1 | 3 | 13 | 25 | 3601 |
#'@md
#'  
#' @examples
#' # To assign the CSR categories to the traits in the IP2G abundance dataset.
#' # This example uses abundance data from a perturbation experiment described in Santillan et al. (2019).
#' # The dataset \code{\link[MicroEcoTools]{CSR_IP2G_data}} is included in the MicroEcoTools package.
#' # Running the command below runs the welch_anova and pairwise_welch_anova functions.
#' # This command  automatically assigns C, S, and R traits to traits abundant in undisturbed, highly diturbed and intermediately disturbed experimental groups.
#' # In addition, it assigns intermediate groups CS, CR, and SR to groups prevalent in two categories, e.g., abundant in undisturbed and intermediately disturbed in compare to highly distrubed --> CR. 
#' \dontrun{
#' csr_assignment <- CSR_assign(dAtA = CSR_IP2G_data, var.name = "Trait", p_adj = "BH", Parallel = TRUE)
#' }
#' @seealso \code{vignette("MicroEcoTools_vignette")}
#' 
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal theme labs facet_wrap element_blank element_rect
#' @importFrom reshape2 melt
#' @importFrom dplyr select full_join mutate_if filter arrange all_of
#' @importFrom magrittr %>%
#' @importFrom stats complete.cases as.formula
#' @importFrom utils str
#' 
#' @export
CSR_assign <- function(dAtA,
                       var.name = "Trait",
                       p_adj = "BH",
                       v.equal = FALSE,
                       p.value.cutoff = 0.05,
                       Parallel = TRUE,
                       Verbose = FALSE) {
# ----------------------------------------------------------------------------
# 1. Input Validation and Conversion
# ----------------------------------------------------------------------------
  
  if (missing(dAtA)){
    warning("No data input!")
    return(invisible(NULL))}
  else {
    # Check if the input object is a SummarizedExperiment or phyloseq object.
    # If so, convert it to a MicroEcoTools-compatible data frame.
    if (inherits(dAtA, "SummarizedExperiment") ||
        inherits(dAtA, "phyloseq")) {
      message(
        "Detected SummarizedExperiment/phyloseq input. Converting to MicroEcoTools data frame format..."
      )
      dAtA <- convert_to_microecotools_df(
        dAtA,
        group_col = "Exp.Grp",
        # change if your metadata column name is different
        replicate_col = "Replicate",
        assay_name = "counts"
      )    # adjust if needed for SummarizedExperiment objects
    }
    
    Vis <- TRUE
# ----------------------------------------------------------------------------
# 2. Initialize Variables and Run Welch-ANOVA
# ----------------------------------------------------------------------------
# Initialize the results list and a counter.
    
    CSR_Results <- vector(mode = "list", length = 0)
    j = 1
    DaTa <- dAtA
    
    # Run the overall Welch-ANOVA on the input data.    
    wa <- Welch_ANOVA(
      dAtA = dAtA,
      var.name = var.name,
      p_adj = p_adj,
      Parallel = Parallel
    )

    # Generate the pairwise Welch-ANOVA table.
    if (colnames(dAtA)[2] != "P1") {
      message(
        "  Creating the pairwise Welch-ANOVA table from the input data using pairwise_welch function ..."
      )
      pww <- pairwise_welch(
        dAtA = dAtA,
        var.name = var.name,
        p_adj = p_adj,
        v.equal = v.equal,
        p.value.cutoff = p.value.cutoff,
        Parallel = Parallel
      )
      dAtA <- pww
    }
    
# ----------------------------------------------------------------------------
# 3. Filter the Data for Group Comparisons
# ----------------------------------------------------------------------------
    # Keep rows where the second column equals the first unique group or the third column equals the last unique group.
    
    filtered_dAtA <- dAtA[which(dAtA[, 2] == unique(dAtA[, 2])[1] |
                                  dAtA[, 3] == unique(dAtA[, 3])[length(unique(dAtA[, 3]))]), ]
    
    a <- data.frame(matrix(nrow = length(unique(filtered_dAtA[, 1])), ncol = 3))
    rownames(a) <- unique(filtered_dAtA[, 1])
    a[, 1] <- rownames(a) # First column: trait names.
    a[, 2] <- ""  # Second column: to be filled with CSR assignments.
    a[, 3] <- "-" # Third column: for remarks.
    
    # Get all unique experimental groups and count them.
    Exp.Grp <- unique(c(filtered_dAtA[, 2], filtered_dAtA[, 3]))
    n.Exp.Grp <- length(unique(c(filtered_dAtA[, 2], filtered_dAtA[, 3])))

# ----------------------------------------------------------------------------
# 4. Assign CSR Categories Based on Pairwise Comparisons
# ----------------------------------------------------------------------------
    
    for (i in unique(filtered_dAtA[, 1])) {
      #Assign C, S and R groups
      dAtA_subset <- filtered_dAtA[which(filtered_dAtA[, 1] == i), ]
      
      if (Verbose == "all") {
        print("Data to check category assignments:")
        print(dAtA_subset)
      }
      # For more than two experimental groups, determine CSR categories.
      if (n.Exp.Grp > 2) {
        if (Verbose == TRUE || Verbose == "all") {
          cat("\n=========================================\n")
          cat("\nChecking categories for:\n", dAtA_subset[1, 1], "\n")
          cat("C ",
              all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 11] != "ns") &
                (
                  all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 12] == "medium") |
                    all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 12] == "large")
                ) & all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 4] > 0),
              "\n")
          cat("S ", (
            all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 11] != "ns") &
              (
                all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 12] == "medium") |
                  all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 12] == "large")
              ) &
              all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 4] < 0)
          ), "\n")
        }
        if (all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 11] != "ns") &
            (all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 12] == "medium") |
             all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 12] == "large")) &
            all(dAtA_subset[dAtA_subset[2] == Exp.Grp[1], 4] > 0)) {
          a[i, 2] <- "C"
          if (Verbose) {
            print ("Classified under C category!")
          }
        } else if (all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 11] != "ns") &
                   (all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 12] == "medium") |
                    all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 12] == "large")) &
                   all(dAtA_subset[dAtA_subset[3] == Exp.Grp[length(Exp.Grp)], 4] < 0)) {
          a[i, 2] <- "S"
          if (Verbose == TRUE ||
              Verbose == "all") {
            print ("Classified under S category!")
          }
        } else if (all(dAtA_subset[1:(n.Exp.Grp - 1), 11] == "ns") &
                   all(dAtA_subset[1:(n.Exp.Grp - 1), 12] != "ne")) {
          a[i, 2] <- "CSR"
        } else if (all(dAtA_subset[1:(n.Exp.Grp - 1), 11] == "ns") &
                   any(dAtA_subset[1:(n.Exp.Grp - 1), 12] == "ne")) {
          a[i, 2] <- "NA"
          a[i, 3] <- "Cannot assign CSR probably due to lack of variability in the data. Please inspect the data."
        } else {
          dAtA_subset_CR <- dAtA_subset[which((dAtA_subset[, 2] == Exp.Grp[1] &
                                                 dAtA_subset[, 3] != Exp.Grp[n.Exp.Grp])), ]
          dAtA_subset_CR[, 4] <- dAtA_subset_CR[, 4] * -1
          dAtA_subset_SR <- dAtA_subset[which(dAtA_subset[, 2] != Exp.Grp[1] &
                                                dAtA_subset[, 3] == Exp.Grp[n.Exp.Grp]), ]
          dAtA_subset_R <- rbind(dAtA_subset_CR, dAtA_subset_SR)
          dAtA_subset_R <- dAtA_subset_R[which(
            dAtA_subset_R[, 11] == "*" &
              (dAtA_subset_R[, 12] == "medium" |
                 dAtA_subset_R[, 12] == "large")
          ), ]
          if (Verbose == "all") {
            print ("Data to check if the trait is an R:")
          }
          if (Verbose == "all") {
            if (exists("dAtA_subset_R"))
              print(dAtA_subset_R)
          }
          if (all(dAtA_subset_R[, 4] > 0) &
              any(dAtA_subset_R[, 2] == Exp.Grp[1]) &
              any(dAtA_subset_R[, 3] == Exp.Grp[n.Exp.Grp])) {
            a[i, 2] <- "R"
            if (Verbose == TRUE ||
                Verbose == "all") {
              print ("Classified under R category!")
            }
          }
        }
        if (Verbose == TRUE || Verbose == "all") {
          if (exists("dAtA_subset_R"))
            cat("R ", (
              all(dAtA_subset_R[, 4] > 0) &
                any(dAtA_subset_R[, 2] == Exp.Grp[1]) &
                any(dAtA_subset_R[, 3] == Exp.Grp[n.Exp.Grp])
            ), "\n")
          else{
            cat("R", " FALSE\n")
            a[i, 2] <- "CSR"
            if (Verbose == TRUE ||
                Verbose == "all") {
              if (a[i, 2] == "CSR")
                print ("\nClassified under CSR category!\n")
            }
          }
          
        }
      } else {
# ----------------------------------------------------------------------------
# 5. Special Case: Only Two Experimental Groups
# ----------------------------------------------------------------------------
        message(
          "CSR assignment with only two groups assuming these groups represent no disturbance and press disturbance.\nPlease note that intermediate groups, namely, CS, CR, and SR will be wrongly categorised.\nIf possible please include at least one experimental group representing the intermediate level of disturbance."
        )
        if (all(dAtA_subset[, 11] != "ns") &
            (dAtA_subset[1, 12] == "medium" |
             dAtA_subset[1, 12] == "large") & dAtA_subset[1, 4] > 0) {
          a[i, 2] <- "C"
        } else if (dAtA_subset[1, 11] != "ns" &
                   (dAtA_subset[1, 12] == "medium" |
                    dAtA_subset[1, 12] == "large") & dAtA_subset[1, 4] < 0) {
          a[i, 2] <- "S"
        } else if (dAtA_subset[1, 11] == "ns" &
                   dAtA_subset[1, 12] != "ne")
          a[i, 2] <- "CSR"
        else if (dAtA_subset[1, 11] == "ns" &
                 dAtA_subset[1, 12] == "ne") {
          a[i, 2] <- "NA"
          a[i, 3] <- "Cannot assign CSR probably due to lack of variability in the data. Please inspect the data."
        }
      }
    }
# ----------------------------------------------------------------------------
# 6. Report Traits with No Assignment
# ----------------------------------------------------------------------------
    
    dAtA_remain <- a[, 1][which(a[, 2] == "")]
    if (Verbose == TRUE || Verbose == "all") {
      cat("Traits yet to be assigned with a category:\n",
          dAtA_remain,
          "\n")
    }
# ----------------------------------------------------------------------------
# 7. Process Unassigned Traits for Compound Category Assignment
# ----------------------------------------------------------------------------
    
    for (ii in dAtA_remain) {
      dAtA_subset <- filtered_dAtA[which(filtered_dAtA[, 1] == ii), ]
      dAtA_subset_CR <- data.frame()
      dAtA_subset_CS <- data.frame()
      dAtA_subset_SR <- data.frame()
      
      dAtA_subset_CR <- dAtA_subset[which((dAtA_subset[, 2] == Exp.Grp[1] &
                                             dAtA_subset[, 3] != Exp.Grp[n.Exp.Grp]) &
                                            dAtA_subset[, 11] != "ns"), ]
      dAtA_subset_CS <- dAtA_subset[which(dAtA_subset[, 2] == Exp.Grp[1] &
                                            dAtA_subset[, 3] == Exp.Grp[n.Exp.Grp] &
                                            dAtA_subset[, 11] != "ns"), ]
      dAtA_subset_SR <- dAtA_subset[which(dAtA_subset[, 2] != Exp.Grp[1] &
                                            dAtA_subset[, 3] == Exp.Grp[n.Exp.Grp] &
                                            dAtA_subset[, 11] != "ns"), ]

# ----------------------------------------------------------------------------
# 8. Evaluate and Assign Compound Categories (CR, SR, CS)
# ----------------------------------------------------------------------------
      
      if (Verbose == TRUE || Verbose == "all") {
        cat("=========================================\n")
        cat("Checking compound categories:\n", dAtA_subset[1, 1], "\n")
        if (exists("dAtA_subset_CS") & exists("dAtA_subset_SR")) {
          cat("CR ",
              (
                length(dAtA_subset_CS[, 4]) > 0 &
                  length(dAtA_subset_SR[, 4]) > 0
              ) &
                (
                  all(dAtA_subset_CS[, 4] > 0) &
                    any(dAtA_subset_SR[, 4] > 0) & !any(dAtA_subset_SR[, 4] < 0)
                ),
              "\n")
        } else
          cat("CR ", "  FAlse\n")
      }
      if (Verbose == TRUE || Verbose == "all") {
        if (exists("dAtA_subset_CR") & exists("dAtA_subset_CS")) {
          cat("SR ",
              (
                length(dAtA_subset_CR[, 4]) > 0 &
                  length(dAtA_subset_CS[, 4]) > 0
              ) &
                (
                  any(dAtA_subset_CR[, 4] < 0) &
                    !any(dAtA_subset_CR[, 4] > 0)  & dAtA_subset_CS[, 4] < 0
                ),
              "\n")
        } else
          cat("SR ", "  FAlse\n")
      }
      if (Verbose == TRUE || Verbose == "all") {
        if (exists("dAtA_subset_CR") &
            exists("dAtA_subset_SR") & exists("dAtA_subset_CS")) {
          cat("CS ",
              (
                length(dAtA_subset_CR[, 4]) > 0 &
                  length(dAtA_subset_SR[, 4]) > 0 &
                  length(dAtA_subset_CS[, 11]) == 0
              ) &
                (all(dAtA_subset_CR[, 4] > 0) &
                   all(dAtA_subset_SR[, 4] < 0)),
              "\n")
        } else
          cat("CS ", "  FAlse\n")
      }
      if (Verbose == "all") {
        print ("dAtA_subset_CR")
        if (exists("dAtA_subset_CR"))
          print (dAtA_subset_CR)
        print ("dAtA_subset_SR")
        if (exists("dAtA_subset_SR"))
          print (dAtA_subset_SR)
        print ("dAtA_subset_CS")
        if (exists("dAtA_subset_CS"))
          print (dAtA_subset_CS)
      }
      if (length(dAtA_subset_CR[, 4]) > 0 &
          length(dAtA_subset_CS[, 4]) > 0) {
        if (any(dAtA_subset_CR[, 4] < 0) &
            !any(dAtA_subset_CR[, 4] > 0)  & dAtA_subset_CS[, 4] < 0) {
          a[ii, 2] <- "SR"
          if (Verbose == TRUE ||
              Verbose == "all") {
            print ("Classified under SR category!")
          }
        }
      }

# ----------------------------------------------------------------------------
# 9. Assign Compound CSR Category Based on Conditions
# ----------------------------------------------------------------------------
      
      if (length(dAtA_subset_CS[, 4]) > 0 &
          length(dAtA_subset_SR[, 4]) > 0) {
        if (all(dAtA_subset_CS[, 4] > 0) &
            any(dAtA_subset_SR[, 4] > 0) & !any(dAtA_subset_SR[, 4] < 0)) {
          a[ii, 2] <- "CR"
          if (Verbose == TRUE ||
              Verbose == "all") {
            print ("Classified under CR category!")
          }
        }
      }
      if (length(dAtA_subset_CR[, 4]) > 0 &
          length(dAtA_subset_SR[, 4]) > 0 &
          length(dAtA_subset_CS[, 11]) == 0) {
        if (all(dAtA_subset_CR[, 4] > 0) & all(dAtA_subset_SR[, 4] < 0)) {
          a[ii, 2] <- "CS"
          if (Verbose == TRUE ||
              Verbose == "all") {
            print ("Classified under CS category!")
          }
        }
      }
      if (a[ii, 2] == "") {
        a[ii, 2] <- "CSR"
        if (Verbose == TRUE ||
            Verbose == "all") {
          print ("Classified under CSR category!")
        }
      }
      if (Verbose == TRUE || Verbose == "all")
        cat("\n")
      
# ----------------------------------------------------------------------------
# 10. Report Unassigned Traits (if any)
# ----------------------------------------------------------------------------
      
    }
    if (Verbose == TRUE ||
        Verbose == "all")
      cat("\n=========================================\n")

# ----------------------------------------------------------------------------
# 11. Merge Overall Welch-ANOVA Results with CSR Assignment Table
# ----------------------------------------------------------------------------
    
    colnames(a) <- c(paste(var.name), "CSR assignment", "Remarks_CSR")
    a <- merge(wa[, c(1, 2, 3, 5)], a, all.x = TRUE)

# ----------------------------------------------------------------------------
# 12. Initialize Data Frames for Visualization Data Preparation
# ----------------------------------------------------------------------------
    
    vis_data_C_p <- data.frame()
    vis_data_S_p <- data.frame()
    vis_data_R_p <- data.frame()
    vis_data_CS_p <- data.frame()
    vis_data_CR_p <- data.frame()
    vis_data_SR_p <- data.frame()
    vis_data_C_r <- data.frame()
    vis_data_S_r <- data.frame()
    vis_data_R_r <- data.frame()
    vis_data_CS_r <- data.frame()
    vis_data_CR_r <- data.frame()
    vis_data_SR_r <- data.frame()
    
# ----------------------------------------------------------------------------
# 13. Generate Visualization Data (Using tryCatch to handle errors)
# ----------------------------------------------------------------------------
    
    tryCatch({
      vis_data_C <- a[which(a[, 5] == "C"), ]
      if (length (vis_data_C[, 1]) > 0) {
        j = 1
        for (i in vis_data_C[, 1]) {
          vis_data_C[j, "mean_rel"] <- mean(as.numeric(dAtA[[i]]), na.rm = TRUE)
          j = j + 1
        }
        vis_data_C_p <- vis_data_C[order(vis_data_C$`adjusted Welch-ANOVA p-value`), ]
        vis_data_C_r <- vis_data_C[order(vis_data_C$`mean_rel`, decreasing = TRUE), ]
        if (nrow(vis_data_C) > 5) {
          vis_data_C_p <- vis_data_C_p[1:5, ]
          vis_data_C_r <- vis_data_C_r[1:5, ]
        }
      }
      vis_data_S <- a[which(a[, 5] == "S"), ]
      if (length (vis_data_S[, 1]) > 0) {
        j = 1
        for (i in vis_data_S[, 1]) {
          vis_data_S[j, "mean_rel"] <- mean(as.numeric(dAtA[[i]]), na.rm = TRUE)
          j = j + 1
        }
        vis_data_S_p <- vis_data_S[order(vis_data_S$`adjusted Welch-ANOVA p-value`), ]
        vis_data_S_r <- vis_data_S[order(vis_data_S$`mean_rel`, decreasing = TRUE), ]
        if (nrow(vis_data_S) > 5) {
          vis_data_S_p <- vis_data_S_p[1:5, ]
          vis_data_S_r <- vis_data_S_r[1:5, ]
        }
      }
      vis_data_R <- a[which(a[, 5] == "R"), ]
      if (length (vis_data_R[, 1]) > 0) {
        j = 1
        for (i in vis_data_R[, 1]) {
          vis_data_R[j, "mean_rel"] <- mean(as.numeric(dAtA[[i]]), na.rm = TRUE)
          j = j + 1
        }
        vis_data_R_p <- vis_data_R[order(vis_data_R$`adjusted Welch-ANOVA p-value`), ]
        vis_data_R_r <- vis_data_R[order(vis_data_R$`mean_rel`, decreasing = TRUE), ]
        if (nrow(vis_data_R) > 5) {
          vis_data_R_p <- vis_data_R_p[1:5, ]
          vis_data_R_r <- vis_data_R_r[1:5, ]
        }
      }
      vis_data_CS <- a[which(a[, 5] == "CS"), ]
      if (length (vis_data_CS[, 1]) > 0) {
        j = 1
        for (i in vis_data_CS[, 1]) {
          vis_data_CS[j, "mean_rel"] <- mean(as.numeric(dAtA[[i]]), na.rm = TRUE)
          j = j + 1
        }
        vis_data_CS_p <- vis_data_CS[order(vis_data_CS$`adjusted Welch-ANOVA p-value`), ]
        vis_data_CS_r <- vis_data_CS[order(as.vector(vis_data_CS$`mean_rel`), decreasing = TRUE), ]
        if (nrow(vis_data_CS) > 5) {
          vis_data_CS_p <- vis_data_CS_p[1:5, ]
          vis_data_CS_r <- vis_data_CS_r[1:5, ]
        }
      }
      vis_data_SR <- a[which(a[, 5] == "SR"), ]
      if (length (vis_data_SR[, 1]) > 0) {
        j = 1
        for (i in vis_data_SR[, 1]) {
          vis_data_SR[j, "mean_rel"] <- mean(as.numeric(dAtA[[i]]), na.rm = TRUE)
          j = j + 1
        }
        vis_data_SR_p <- vis_data_SR[order(vis_data_SR$`adjusted Welch-ANOVA p-value`), ]
        vis_data_SR_r <- vis_data_SR[order(vis_data_SR$`mean_rel`, decreasing = TRUE), ]
        if (nrow(vis_data_SR) > 5) {
          vis_data_SR_p <- vis_data_SR_p[1:5, ]
          vis_data_SR_r <- vis_data_SR_r[1:5, ]
        }
      }
      vis_data_CR <- a[which(a[, 5] == "CR"), ]
      if (length (vis_data_CR[, 1]) > 0) {
        j = 1
        for (i in vis_data_CR[, 1]) {
          vis_data_CR[j, "mean_rel"] <- mean(as.numeric(dAtA[[i]]), na.rm = TRUE)
          j = j + 1
        }
        vis_data_CR_p <- vis_data_CR[order(vis_data_CR$`adjusted Welch-ANOVA p-value`), ]
        vis_data_CR_r <- vis_data_CR[order(vis_data_CR$`mean_rel`, decreasing = TRUE), ]
        if (nrow(vis_data_CR) > 5) {
          vis_data_CR_p <- vis_data_CR_p[1:5, ]
          vis_data_CR_r <- vis_data_CR_r[1:5, ]
        }
      }
      vis_data_list_p <- rbind(
        vis_data_C_p,
        vis_data_S_p,
        vis_data_R_p,
        vis_data_CS_p,
        vis_data_CR_p,
        vis_data_SR_p
      )[, 1]
      vis_data_list_r <- rbind(
        vis_data_C_r,
        vis_data_S_r,
        vis_data_R_r,
        vis_data_CS_r,
        vis_data_CR_r,
        vis_data_SR_r
      )[, 1]
      
      CSR_plot_sorted_abundance <- list(
        DaTa = DaTa,
        CSR_cat = a,
        CSR_vis_list = vis_data_list_r,
        var.name = var.name,
        sort.var = "abundance"
      )
      CSR_plot_sorted_p_value <- list(
        DaTa = DaTa,
        CSR_cat = a,
        CSR_vis_list = vis_data_list_p,
        var.name = var.name,
        sort.var = "adjusted P-value"
      )
      
    }, error = function(e) {
      message(
        "Warning: visualisation cannot be performed due to insufficient number of traits categorised as C, S or R "
      )
    })
    # ----------------------------------------------------------------------------
    # 14. Return the Final Results
    # ----------------------------------------------------------------------------
    
    if (sys.nframe() > 1) {
      return(a)
    } else{
      CSR_Results <- list(wa,
                          pww,
                          a,
                          CSR_plot_sorted_abundance,
                          CSR_plot_sorted_p_value)
      names(CSR_Results) <- c(
        "Welch_ANOVA",
        "Pairwise_Welch_ANOVA",
        "CSR_assignments",
        "CSR_plot_sorted_abundance",
        "CSR_plot_sorted_p_value"
      )
      class(CSR_Results$CSR_assignments) <- "CSR_assignments"
      return(CSR_Results)
    }
  }
}



#' Summary method for CSR_assignments objects
#' @param object An object of class "CSR_assignments" returned as part of CSR_Results by CSR_assign.
#' @examples
#' \dontrun{
#' csr_assignment <- CSR_assign(dAtA = CSR_IP2G_data, var.name = "Trait", p_adj = "BH", Parallel = TRUE)
#' summary(csr_assignment)
#' }
#' @return Invisibly returns a table summarizing the CSR assignments.
#' @export
summary.CSR_assignments <- function(object) {
  cat("Summary of CSR Assignments:\n")
  # Simply print the table 'a' which holds the assignments
  print(as.data.frame(unclass(object$CSR_assignments))[, c(1, 2, 5)])
  invisible(object)
}



#' Generate a CSR Plot from Plotting Parameters in the result object generated by the CSR_assign
#'
#' @description This function generates a ggplot object for visualizing CSR (Competitor–Stress tolerant–Ruderal)
#' assignments. It accepts the csr_result object generated by the CSR_assign function and a character argument \code{plot_type}
#' to determine whether to use the parameters for a sorted abundance plot or a sorted adjusted P-value plot.
#'
#' @param CSR_results A list containing two sub-lists:
#'   \itemize{
#'     \item \code{CSR_plot_sorted_abundance}: A list of plotting parameters for a sorted abundance plot.
#'     \item \code{CSR_plot_sorted_p_value}: A list of plotting parameters for a sorted adjusted P-value plot.
#'   }
#'   Each sub-list must include:
#'   \itemize{
#'     \item \code{DaTa}: A data frame containing the data (including metadata and counts).
#'     \item \code{CSR_cat}: A data frame containing CSR assignment results.
#'     \item \code{CSR_vis_list}: A vector of column names (character) to be visualized.
#'     \item \code{var.name}: A character string specifying the variable name (e.g., "Taxa" or "Trait").
#'     \item \code{sort.var}: A character string specifying the sorting criterion (e.g., "abundance" or "adjusted P-value").
#'   }
#' @param plot_type A character string to determine which plot to generate.
#'   Valid options are \code{"abundance"} (default) for the sorted abundance plot, or \code{"p_value"}
#'   for the sorted adjusted P-value plot.
#'
#' @return A ggplot object displaying the CSR plot.
#'
#' @examples
#' \dontrun{
#' 
#'   # Generate and display the sorted abundance plot:
#'   CSR_plot(csr_assignment, plot_type = "abundance")
#'
#'   # Generate and display the sorted adjusted P-value plot:
#'   CSR_plot(csr_assignment, plot_type = "p_value")
#' }
#' @export
CSR_plot <- function(CSR_results,
                     plot_type = c("abundance", "p_value")) {
  # Ensure that plot_type is one of the allowed values
  plot_type <- match.arg(plot_type)
  
  # Select the appropriate plotting parameters from CSR_results
  if (plot_type == "abundance") {
    params <- CSR_results$CSR_plot_sorted_abundance
  } else {
    params <- CSR_results$CSR_plot_sorted_p_value
  }
  
  # Extract required parameters
  DaTa         <- params$DaTa          # Main data frame (with metadata and counts)
  CSR_cat      <- params$CSR_cat       # CSR assignment results data frame
  CSR_vis_list <- params$CSR_vis_list  # Vector of column names to be visualized
  var.name     <- params$var.name      # Name for the variable (e.g., "Taxa" or "Trait")
  sort.var     <- params$sort.var      # Sorting criterion (e.g., "abundance" or "adjusted P-value")
  
  # Use var.name as the plotting variable label
  CSR_plot_variable <- var.name
  
  # Extract the metadata (first two columns, e.g., Exp.Grp and Replicate)
  vis_data_m <- DaTa[, c(1, 2)]
  
  # Select the columns specified in CSR_vis_list from the main data
  vis_data <- dplyr::select(DaTa, dplyr::all_of(CSR_vis_list))
  
  # Combine the metadata with the selected columns
  vis_data <- cbind(vis_data_m, vis_data)
  
  # Filter CSR assignment results (here we keep columns 1 and 5)
  CSR_cat_filtered <- CSR_cat[, c(1, 5)]
  
  # Reshape the data into long format so that each row represents a measurement for a given trait/taxa
  vis_data_0 <- reshape2::melt(
    vis_data,
    id.vars = c(colnames(vis_data)[1], colnames(vis_data)[2]),
    variable.name = CSR_plot_variable,
    value.name = "Abundance"
  )
  
  # Merge the long-format data with the CSR assignment results
  vis_data <- merge(vis_data_0, CSR_cat_filtered, all.x = TRUE)
  
  # Rename columns to standard names
  names(vis_data) <- c(var.name, "Exp.Grp", "Rep", "Abundance", "CSR_categories")
  
  # Convert Exp.Grp to a factor for proper grouping in the plot
  vis_data$Exp.Grp <- as.factor(vis_data$Exp.Grp)
  
  # Create and return the ggplot object
  plot <- ggplot2::ggplot(vis_data,
                          ggplot2::aes(x = Exp.Grp, y = Abundance, color = CSR_categories)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "transparent")
    ) +
    ggplot2::labs(
      title = paste(
        "Top 5",
        CSR_plot_variable,
        "from C, S, R and intermediate categories sorted based on",
        sort.var
      )
    ) +
    ggplot2::facet_wrap(as.formula(paste("~", CSR_plot_variable)), scales = "free") + ggplot2::scale_color_brewer(palette = "Set2")
  
  return(plot)
}



#' @title Theoretical Microbial Ecology Tools (MicroEcoTools)
#' 
#' @description Probability prediction for the CSR assignments.
#' Assigning CSR categories from Grime's framework to the input variable and calculates the probability of assigning these categories using a simulation.
#' This function generates simulated datasets using the input dataset that have a similar multinomial distribution.
#' It is advised to use at least 1000 simulations to get more reliable results (recomended NSim = 10000).
#' After performing the simulation it passes the data to the CSR_assign function to assign CSR categories.
#'
#' This function assigns CSR categories to the input variables using the output table from the \code{\link[MicroEcoTools]{pairwise_welch}} function.
#' To do so, it looks at the pairwise comparisons and assigns C to the ones abundant in no disturbance group, S to the ones abundant in high disturbance group, and R to the ones abundant in intermediate disturbance groups.
#' It also assigns intermediate categories, CR, CS, and SR to groups that are abundant in no-disturbance and intermediate-disturbance, no-disturbance and high-disturbance, and intermediate-disturbance and high-disturbance, respectively.
#' Traits that are not differentially abundant in any groups will be assigned to the CSR category.
#' Please note that if the variability is lacking the group will be categorised as NA.
#' 
#' @param DaTa Dataframe that contains experimental group, replicates, and parameters.
#' @note
#' The input data can be relative abundance or count data. The user must ensure the total count across groups is similar.
#' @param NSim Number of simulations minimum 1000, recomended 10000. Defaults to \code{1000}.
#' @param var.name Variable name used for generating output table, for example, taxa, trait, function.  Defaults to \code{"Trait"}.
#' @param p_adj P-value adjustment method. All methods mentioned in \code{\link[stats]{p.adjust}} function can be used ("BH", "BY", "holm", "hommel" or "none").  Defaults to \code{"BH"}.
#' @param v.equal Assumption of equal variance can be forced using this parameter. By default it is set to FALSE where it cannot be set to FALSE it will switch to TRUE automatically; a note will be shown in the remarks column in those cases.   Defaults to \code{FALSE}.
#' @param p.value.cutoff A cut-off value for calling the P-value significant can be set using this parameter. The default value is 0.05.   Defaults to \code{0.05}.
#' @param Parallel Using this parameter you can make use of more CPU cores to decrease run time for the calculation. By default it is set to TRUE and uses half of the available CPU cores to perform the calculations. Defaults to \code{TRUE}.
#' @param Keep_data You can store the entire calculation results including the simulated communities in a list named CSR_Sim by setting this parameter to TRUE. By default it deletes the simulated data to save space. Defaults to \code{FALSE}.
#' @param NuLl.test Using this parameter you can compare the assigned categories in your observed communities with a community generated with null hypothesis. (Beta version!).   Defaults to \code{FALSE}.
#' @return This function returns a table containing the pairwise statistical comparison results. The output table can be fed into the CSR_assign function to assign CSR categories.
#'
#'@details
#'The input data sample:
#'
#' | Exp.Grp | Replicate | TAXA1 | TAXA2 | TAXA3 |
#'|-----------|:-----------:|-----------:|:-----------:|-----------:|
#'  | X0 | 1 | 10 | 91 | 68 |
#'  | X0 | 2 | 11 | 86 | 70 |
#'  | X0 | 3 | 12 | 84 | 70 |
#'  | X1 | 1 | 15 | 30 | 3452 |
#'  | X1 | 2 | 19 | 45 | 3274 |
#'  | X1 | 3 | 13 | 25 | 3601 |
#'@md
#'  
#' @examples
#' # Generally, the number of replicates in perturbation experiments cannot provide enough statistical power to reveal the effects.
#' # CSR_Simulation generates a number of simulated communities with similar multinomial distributions as the observed communities.
#' # This example uses abundance data from a perturbation experiment described in Santillan et al. (2019).
#' # The dataset \code{\link[MicroEcoTools]{CSR_IP2G_data}} is included in the MicroEcoTools package.
#' # Running the command below generates simulated communities and performs \code{\link[MicroEcoTools]{CSR_assign}} on the observed and simulated communities.
#' # This command automatically assigns C, S, and R traits to traits abundant in undisturbed, highly diturbed and intermediately disturbed experimental groups.
#' # In addition, it assigns intermediate groups CS, CR, and SR to groups prevalent in two categories, e.g., abundant in undisturbed and intermediately disturbed in compare to highly distrubed --> CR. 
#' # At the end, it generates a table showing the percentage of times a trait was assigned to different categories.
#' \dontrun{
#' csr_simulation <- CSR_Simulation(DaTa = CSR_IP2G_data, NSim = 1000, p_adj = "BH", var.name = "Trait", NuLl.test = FALSE, Keep_data = FALSE)
#' }
#' 
#' @importFrom dplyr mutate_if full_join select filter arrange
#' @importFrom magrittr %>%
#' @importFrom stringr str_count
#' @importFrom purrr reduce
#' @importFrom stats rmultinom
#' 
#' @export
CSR_Simulation <- function(DaTa,
                           NSim = 1000,
                           p_adj = "BH",
                           p.value.cutoff = 0.05,
                           Parallel = TRUE,
                           var.name = "Trait",
                           v.equal = FALSE,
                           NuLl.test = FALSE,
                           Keep_data = FALSE) {
  # ----------------------------------------------------------------------------
  # 1. Input Data Validation and Conversion
  # ----------------------------------------------------------------------------
  
  if (missing(DaTa))
    {warning("No data input!")
     return(invisible(NULL))
    }
  else {
    # Check if the input object is a SummarizedExperiment or phyloseq object.
    # If so, convert it to a MicroEcoTools-compatible data frame.
    if (inherits(DaTa, "SummarizedExperiment") ||
        inherits(DaTa, "phyloseq")) {
      message(
        "Detected SummarizedExperiment/phyloseq input. Converting to MicroEcoTools data frame format..."
      )
      DaTa <- convert_to_microecotools_df(
        DaTa,
        group_col = "Exp.Grp",
        # change if your metadata column name is different
        replicate_col = "Replicate",
        assay_name = "counts"
      )    # adjust if needed for SummarizedExperiment objects
    }
    
    if (rowSums(DaTa[1, c(-1, -2)] <= 1) |
        rowSums(DaTa[1, c(-1, -2)] <= 100))
      DaTa[, c(-1, -2)] <- DaTa[, c(-1, -2)] * 100000
    
    message(
      paste(
        "  CSR assignment with",
        NSim,
        "simulated datasets and",
        p_adj,
        "as P-value correction method\n"
      )
    )
    # ----------------------------------------------------------------------------
    # 2. Initialize Result Containers
    # ----------------------------------------------------------------------------
    # Create a list to hold all simulation results.
    
    CSR_Sim <- vector(mode = "list", length = 0)
    CSR_Sim[["data"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["merged_data"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["CSR"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["Verdict"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["Final_Verdict"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["Summary"]] <- vector(mode = "list", length = 0)

    # ----------------------------------------------------------------------------
    # 3. Data Preprocessing for Simulation
    # ----------------------------------------------------------------------------
    # Remove incomplete cases from the input data.
    
    DaTa <- DaTa[complete.cases(DaTa), ]
    InPuT <- DaTa
    InPuT <- as.data.frame(t(InPuT))
    InPuT <- InPuT[c(-1, -2), ]
    expected <- InPuT
    expected <- apply(expected, 2, function(x)
      as.numeric(x))
    expected <- as.data.frame(expected)
    expected <- ceiling(expected)
    rownames(expected) <- rownames(InPuT)
    expected <- expected[which(rowSums(expected) > 0), ]
    
    # ----------------------------------------------------------------------------
    # 4. Run Simulations
    # ----------------------------------------------------------------------------
    message("Starting the simulation\n")
    if (!NuLl.test) {
      CSR_Sim[["data"]] <- lapply(seq(NSim), function(i, expected) {
        # create a list of expected cell counts (list element = row of expected)
        .list <- lapply(apply(expected, 1, list), unlist)
        # sample from these expected cell counts and recombine into a data.frame
        as.data.frame(do.call(rbind, lapply(.list, function(.x)
          t(
            rmultinom(
              n = 1,
              prob = .x,
              size = sum(.x)
            )
          ))))
      }, expected = expected)
      CSR_Sim[["data"]][["Obs"]] <- expected
    }
    
    if (NuLl.test) {
      CSR_Sim[["data"]] <- lapply(seq(NSim), function(i, expected) {
        # create a list of expected cell counts (list element = row of expected)
        .list <- lapply(apply(expected, 1, list), unlist)
        # sample from these expected cell counts and recombine into a data.frame
        as.data.frame(do.call(rbind, lapply(.list, function(.x)
          t(
            rmultinom(
              n = 1,
              prob = rep(sum(.x) / 12, 12),
              size = sum(.x)
            )
          ))))
      }, expected = expected)
    }
    message("  Simulation            ", "\u2714", "\n")
    # ----------------------------------------------------------------------------
    # 5. Perform CSR Assignment for Each Simulation
    # ----------------------------------------------------------------------------

    message("CSR assignment step, please be patient...")
    for (i in 1:length(CSR_Sim[["data"]])) {
      message(paste("\n  Step", i, "out of", length(CSR_Sim[["data"]])))
      rownames(CSR_Sim[[1]][[i]]) <- rownames(expected)
      CSR_Sim[[2]][[i]] <- as.data.frame(t(CSR_Sim[[1]][[i]]))
      CSR_Sim[[2]][[i]] <- cbind(DaTa[1], DaTa[2], CSR_Sim[[2]][[i]])
      CSR_Sim[[3]][[i]] <- suppressMessages(CSR_assign(
        dAtA = CSR_Sim[[2]][[i]],
        var.name = var.name,
        p_adj = p_adj,
        p.value.cutoff = p.value.cutoff,
        v.equal = v.equal,
        Parallel = Parallel
      ))[c(1, 5)]
      if (i == length(CSR_Sim[["data"]]))
        message("\nCSR assignment        ", "\u2713")
    }
    message("CSR assignment            ", "\u2714", "\n")
    
    # ----------------------------------------------------------------------------
    # 6. Aggregate CSR Assignments Across Simulations
    # ----------------------------------------------------------------------------
    # Combine all CSR assignment results from each simulation using a full join.
    CSR_Sim[["Verdict"]] <- CSR_Sim[[3]] %>% reduce(full_join, by = var.name)
    
    # ----------------------------------------------------------------------------
    # 7. Count CSR Categories for Each Trait
    # ----------------------------------------------------------------------------
    message(paste("Counting CSR categories for each ", var.name, "\n"))
    for (i in 1:length(CSR_Sim[["Verdict"]][, 1])) {
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 2] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])C\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 3] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])R\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 4] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])S\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 5] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])CR\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 6] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])CS\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 7] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])SR\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 8] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])CSR\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]]) + 9] <- sum(str_count(CSR_Sim[["Verdict"]][i, -1][1:(NSim +
                                                                                                               1)], "(?<![\\S])NA\\b"))
      if (i == length(CSR_Sim[["Verdict"]][, 1]))
        message("CSR category count    ", "\u2714", "\n")
    }
    
    # ----------------------------------------------------------------------------
    # 8. Compute Final Verdict Table (Percentage Assignment)
    # ----------------------------------------------------------------------------
    message("Counting CSR categories        ", "\u2714", "\n")
    CSR_Sim[["Final_Verdict"]] <- CSR_Sim[["Verdict"]][c(1, (length(CSR_Sim[["data"]]) +
                                                               2):(length(CSR_Sim[["data"]]) + 9))]
    colnames(CSR_Sim[["Final_Verdict"]]) <- c(var.name, "C", "R", "S", "CR", "CS", "SR", "CSR", "NA")
    CSR_Sim[["Final_Verdict"]] <- mutate_if(CSR_Sim[["Final_Verdict"]], is.numeric, ~ . /
                                              (NSim + 1) * 100)
    
    # ----------------------------------------------------------------------------
    # 9. Merge Final Verdict with Original CSR Assignment Table
    # ----------------------------------------------------------------------------
    CSR_Sim[["Final_Verdict"]] <- do.call(cbind, merge(suppressMessages(CSR_assign(DaTa)), CSR_Sim[["Final_Verdict"]], all.x = TRUE))
    colnames(CSR_Sim[["Final_Verdict"]]) <- c(var.name, "Welch-ANOVA p-value", "Method", "Adjusted Welch-ANOVA p-value", "CSR assignment (observed)", "Remarks", "C", "R", "S", "CR", "CS", "SR", "CSR", "NA")
    # ----------------------------------------------------------------------------
    # 10. Return the Final Simulation Results
    # ----------------------------------------------------------------------------
    CSR_Sim[["Summary"]] <- CSR_Sim[["Final_Verdict"]]
    class(CSR_Sim[["Summary"]]) <- "Final_verdict"
    return(CSR_Sim)
  }
}

#' Summary method for CSR_assignments objects
#' @param object An object of class "CSR_assignments" returned as part of CSR_Results by CSR_assign.
#' @examples 
#' \dontrun{
#' csr_simulation <- CSR_Simulation(DaTa = CSR_IP2G_data, NSim = 1000, p_adj = "BH", var.name = "Trait", NuLl.test = FALSE, Keep_data = FALSE)
#' summary(csr_simulation)
#' }
#' @return Invisibly returns a table summarizing the CSR assignments.
#' @export
summary.CSR_simulation <- function(object) {
  cat(
    "Summary of CSR Assignments using simulation method. CSR assignment column shows the assigned group for the observed community while the rest of columns show the percentage each category was assigned to the simulated community.:\n"
  )
  # Simply print the table 'a' which holds the assignments
  print(as.data.frame(unclass(object$Final_Verdict))[, c(1, 5, 7, 8, 9, 10, 11, 12, 13)])
  invisible(object)
}