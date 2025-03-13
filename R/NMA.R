#' @title Theoretical Microbial Ecology Tools (MicroEcoTools)
#'
#' @description The tools implemented in MicroEcoTools help researchers and practitioners apply theoretical frameworks to community data. 
#' The current version covers theories connecting assembly mechanisms to diversity and functions. In particular, the NMA function implements 
#' a null model approach (based on a multinomial distribution) to assess the relative contribution of stochastic versus deterministic assembly 
#' mechanisms.
#'
#' @author Soheil A. Neshat, Ezequiel Santillan, Stefan Wuertz
#' @author (The MicroEcoTools package was developed in Singapore Centre for Environmental Life Sciences Engineering, 
#' Nanyang Technological University, Singapore)
#'
#' @keywords Trait-based, CSR, Grime's triangle, Intermediate stochasticity hypothesis, Null model analysis, community ecology, Theoretical microbial ecology
#'
#' @section Null model analysis (NMA):
#' Using a null model approach, this function tests the relative contribution of assembly mechanisms in a given community. 
#' Based on the observed communities, it simulates random communities using a multinomial distribution. A diversity metric is then used 
#' to calculate the standard effect size between the observed communities and the simulated (null) communities.
#'
#' @param DaTa A data frame containing environmental and community data in long format with column names: 
#'   Exp.Grp, Replicate, Time_point, TAXA, and Count. At least two replicates per experimental group are required.
#' @param NSim A number specifying the number of simulations (values below 1000 are not recommended). Defaults to \code{1000}.
#' @param InDeX A diversity index from the vegan package ("shannon", "simpson", or "invsimpson"). Defaults to \code{"invsimpson"}.
#' @param ObSsIm Logical. If TRUE, a dataset with a multinomial distribution using the observed data as a probability matrix is generated. Defaults to \code{FALSE}.
#' @param p.adj A p-value correction method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none"). Defaults to \code{"BH"}.
#' @param Plot_Level A vector with experimental group titles. Defaults to the unique values of the Exp.Grp column in the input data.
#' @param Keep_Data Logical. If TRUE, simulated data is kept; otherwise, it is removed to save space. Defaults to \code{FALSE}.
#' @param Deb Logical. If TRUE, the function runs in verbose mode. Defaults to \code{FALSE}.
#'
#' @returns A list (NMA_Results) containing:
#'   \itemize{
#'     \item The original data (DaTa)
#'     \item Observed and simulated community data
#'     \item Diversity metrics and test statistics (including p-values and effect sizes)
#'     \item Processed data for plotting (e.g., SES, SI, Kraft models, Cohen's d, Hill diversity)
#'     \item A set of ggplot objects for visualization (if generated)
#'     \item The NMA parameters used in the analysis
#'   }
#'
#' @examples
#' \dontrun{
#'   # Perform null model analysis on experimental data from a perturbation experiment.
#'   nma_results <- NMA(DaTa = NMA_data, NSim = 1000)
#' }
#'
#' @note All environmental variables should be converted to factors.
#' @note Accepted environmental variables: Time_point, Exp.Grp, and Replicate.
#'
#' @references Santillan, Ezequiel, et al. "Frequency of disturbance alters diversity, function, and underlying assembly mechanisms of complex bacterial communities." npj Biofilms and Microbiomes 5.1 (2019): 1-9.
#' @references Kraft, Nathan JB, et al. "Disentangling the drivers of β diversity along latitudinal and elevational gradients." Science 333.6050 (2011): 1755-1758.
#'
#' @importFrom reshape2 dcast melt
#' @importFrom stats p.adjust t.test
#' @importFrom utils ?
#' @importFrom vegan specnumber diversity
#' @importFrom effectsize cohens_d
#' @importFrom ggplot2 aes geom_bar geom_boxplot geom_jitter ggtitle xlab ylab facet_wrap element_text element_blank element_rect theme
#' @importFrom progress progress_bar
#'
#' @export
NMA <- function(DaTa, NSim = 1000, InDeX = "invsimpson", ObSsIm = FALSE, p.adj = "BH", Plot_Level, Keep_Data = FALSE, Deb = FALSE) {
  
  # ---------------------------------------------------------------------------
  # 1. Data Preparation and Validation
  # ---------------------------------------------------------------------------
  # If no data provided, use the demo dataset
  if (missing(DaTa)) {
    DaTa <- NMA_data
    print("NMA analysis using a demo dataset from Santillan et al. 2019.")
  }
  
  # Ensure that DaTa is a data frame
  if (!is.data.frame(DaTa)) {
    ?NMA
    warning("Invalid data format. Please refer to the manual using ?NMA")
    return(invisible(NULL))
  }
  
  # Convert input if it is a SummarizedExperiment or phyloseq object
  if (inherits(DaTa, "SummarizedExperiment") || inherits(DaTa, "phyloseq")) {
    message("Detected SummarizedExperiment/phyloseq input. Converting to MicroEcoTools data frame format...")
    DaTa <- convert_to_microecotools_df(
      DaTa,
      group_col = "Exp.Grp",        # Adjust if needed
      replicate_col = "Replicate",  # Adjust if needed
      assay_name = "counts"         # Adjust if needed
    )
  }
  
  # Validate that NSim is an integer
  if (NSim %% 1 != 0) {
    return("Invalid number of simulations entered. Please try again with an integer number as the number of simulations.")
  }
  
  # Set Plot_Level if missing: use unique experimental groups
  if (missing(Plot_Level)) {
    Plot_Level <- unique(DaTa$Exp.Grp)
  }
  
  # Load pipe operator from magrittr
  `%>%` <- magrittr::`%>%`
  
  # ---------------------------------------------------------------------------
  # 2. Initialize Placeholders for Simulation Results and Statistics
  # ---------------------------------------------------------------------------
  P_SN_NMA_Sim_diversity <- vector(mode = "list", length = 0)
  P_SN_NMA_Sim_data       <- vector(mode = "list", length = 0)
  P_SN_NMA_Observed_data  <- vector(mode = "list", length = 0)
  P_SN_NMA_Obs_diversity  <- vector(mode = "list", length = 0)
  NMA_stat                <- vector(mode = "list", length = 0)
  NMA_kraft               <- vector(mode = "list", length = 0)
  NMA_output              <- vector(mode = "list", length = 0)
  NMA_output_SES          <- vector(mode = "list", length = 0)
  NMA_output_SI           <- vector(mode = "list", length = 0)
  NMA_Results             <- vector(mode = "list", length = 0)
  kraft                   <- vector(mode = "list", length = 0)
  NSimO                   <- 0
  NMA_stat_p              <- vector(length = 0)
  
  # ---------------------------------------------------------------------------
  # 3. Data Format Conversion
  # ---------------------------------------------------------------------------
  # If the "Time_point" column is missing, set it to 1
  if (!("Time_point" %in% colnames(DaTa))) {
    DaTa$Time_point <- 1
  }
  
  # Check if data is in long format (i.e., contains "Count"); if not, melt the data
  if ("Count" %in% colnames(DaTa)) {
    cat("Input data is in long format. Proceeding with the analysis...")
  } else {
    cat("Input data is in wide format. Converting the data into long format...")
    DaTa <- reshape2::melt(
      DaTa,
      id.vars = c("Exp.Grp", "Replicate", "Time_point"),
      variable.name = "TAXA",
      value.name = "Count"
    )
  }
  
  # Convert columns to appropriate types
  DaTa$Exp.Grp   <- as.factor(DaTa$Exp.Grp)
  DaTa$Time_point <- as.factor(DaTa$Time_point)
  DaTa$Replicate <- as.factor(DaTa$Replicate)
  DaTa$TAXA      <- as.character(DaTa$TAXA)
  
  # ---------------------------------------------------------------------------
  # 4. Simulation Setup: Progress Bar
  # ---------------------------------------------------------------------------
  pb <- progress::progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = length(unique(DaTa$Exp.Grp)),
    complete = "=",
    incomplete = "-",
    current = ">",
    clear = FALSE,
    width = 150
  )
  
  # ---------------------------------------------------------------------------
  # 5. Main Simulation Loop: Loop over Experimental Groups and Time Points
  # ---------------------------------------------------------------------------
  tryCatch({
    for (i in unique(DaTa$Exp.Grp)) {
      pb$tick()
      P_SN_NMA_tmp1 <- DaTa %>% dplyr::filter(Exp.Grp == paste(i))
      
      for (j in unique(P_SN_NMA_tmp1$Time_point)) {
        P_SN_NMA_tmp2 <- P_SN_NMA_tmp1 %>% dplyr::filter(Time_point == paste(j))
        NRep <- length(unique(P_SN_NMA_tmp2$Replicate))
        
        if (NRep < 2) {
          print(paste("There is only one replicate for", P_SN_NMA_tmp1[1, 1],
                      "at time point", P_SN_NMA_tmp1[1, 3],
                      ". Please consider removing this experimental group - time point."))
        }
        
        # Reshape data from wide to long format for simulation:
        P_SN_NMA_tmp2 <- reshape2::dcast(
          data.table::as.data.table(P_SN_NMA_tmp2),
          Time_point + Exp.Grp + TAXA ~ Replicate,
          value.var = "Count",
          fun.aggregate = sum
        )
        
        # Extract numeric counts (expected values) and remove rows with zero totals:
        expected <- P_SN_NMA_tmp2 %>% dplyr::select_if(is.numeric)
        expected <- expected[rowSums(expected) != 0, ]
        
        # Store labels for p-value extraction:
        NMA_stat[["p_value"]][["Exp.Grp"]] <- c(NMA_stat[["p_value"]][["Exp.Grp"]], paste(i))
        NMA_stat[["p_value"]][["Time_point"]] <- c(NMA_stat[["p_value"]][["Time_point"]], paste(j))
        
        # ---------------------------------------------------------------------
        # Simulation for each taxon:
        # For each taxon (or trait), compute the total count across replicates,
        # then simulate a new community sample from a multinomial distribution with:
        #   - size: total count for that taxon.
        #   - prob: rep(sum(.x)/NRep, NRep) which normalizes to rep(1/NRep, NRep).
        # This enforces the null hypothesis of uniform allocation across replicates.
        # ---------------------------------------------------------------------
        P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]] <- lapply(seq(NSim), function(i, expected) {
          .list <- lapply(apply(expected, 1, list), unlist)
          as.data.frame(do.call(rbind, lapply(.list, function(.x)
            t(rmultinom(n = 1, prob = rep(sum(.x)/NRep, NRep), size = sum(.x)))
          )))
        }, expected = expected)
        
        # If ObSsIm is TRUE, generate additional observed simulated data using observed probabilities:
        if (ObSsIm) {
          NSimO <- NSim - 1
          P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]] <- lapply(seq(NSimO), function(i, expected) {
            .list <- lapply(apply(expected, 1, list), unlist)
            as.data.frame(do.call(rbind, lapply(.list, function(.x)
              t(rmultinom(n = 1, prob = .x, size = sum(.x)))
            )))
          }, expected = expected)
        }
        
        # Save the observed data for this experimental group and time point
        P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]] <- expected
        
        # Calculate diversity metrics for simulated and observed data
        for (k in 1:NRep) {
          for (l in 1:NSim) {
            P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]] <- c(
              P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]],
              1 - (vegan::diversity(P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]][[l]][[k]], index = InDeX) /
                     vegan::diversity(rowSums(P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]][[l]]), index = InDeX))
            )
            P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]][[paste(j)]] <- c(
              P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]][[paste(j)]],
              mean(as.numeric(lapply(P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]][[l]], 
                                     function(x) vegan::specnumber(x)))) /
                vegan::specnumber(rowSums(as.data.frame(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]])))
            )
          }
          for (m in 1:(NSimO + 1)) {
            P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]] <- c(
              P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]],
              1 - (vegan::diversity(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][[m]][[k]], index = InDeX) /
                     vegan::diversity(rowSums(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][[m]]), index = InDeX))
            )
          }
        }
        
        # Calculate beta-diversity for observed data
        P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]] <- c(
          P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]],
          mean(as.numeric(lapply(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]],
                                 function(x) vegan::specnumber(x)))) /
            vegan::specnumber(rowSums(as.data.frame(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]])))
        )
        
        # Calculate Cohen's d for the difference between observed and simulated diversity
        NMA_stat[["p_value"]][["Cohens_d"]] <- c(
          NMA_stat[["p_value"]][["Cohens_d"]],
          abs(effectsize::cohens_d(P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]],
                                   P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]])$Cohens_d)
        )
        
        # If ObSsIm is TRUE, perform a t-test between observed and simulated diversity values and store the p-value
        if (ObSsIm) {
          NMA_stat[["t-test"]][[paste(i)]][[paste(j)]] <- t.test(
            P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]],
            P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]]
          )
          NMA_stat[["p_value"]][["p_value"]] <- c(
            NMA_stat[["p_value"]][["p_value"]],
            NMA_stat[["t-test"]][[paste(i)]][[paste(j)]]$p.value
          )
        }
      }
      
      # Compute means and standard deviations for simulated and observed diversity indices
      NMA_output[[paste(i)]][["Sim_mean"]]      <- lapply(P_SN_NMA_Sim_diversity[[paste(i)]], mean)
      NMA_output[[paste(i)]][["Sim_sd"]]        <- lapply(P_SN_NMA_Sim_diversity[[paste(i)]], sd)
      NMA_output[[paste(i)]][["Obs_mean"]]        <- lapply(P_SN_NMA_Obs_diversity[[paste(i)]], mean)
      NMA_output[[paste(i)]][["Obs_sd"]]          <- lapply(P_SN_NMA_Obs_diversity[[paste(i)]], sd)
      NMA_output[[paste(i)]][["Sim_mean_beta"]]   <- lapply(P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]], mean)
      NMA_output[[paste(i)]][["Sim_sd_beta"]]     <- lapply(P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]], sd)
      
      # Calculate standardized effect sizes (SES) and stochastic intensities (SI)
      for (j in unique(names(NMA_output[[paste(i)]][["Obs_mean"]]))) {
        NMA_output_SES[[paste(i)]][[paste(j)]] <- abs((NMA_output[[paste(i)]][["Obs_mean"]][[paste(j)]] -
                                                         NMA_output[[paste(i)]][["Sim_mean"]][[paste(j)]]) /
                                                        NMA_output[[paste(i)]][["Sim_sd"]][[paste(j)]])
        NMA_output_SI[[paste(i)]][[paste(j)]] <- 100 - ((abs(NMA_output[[paste(i)]][["Obs_mean"]][[paste(j)]] -
                                                               NMA_output[[paste(i)]][["Sim_mean"]][[paste(j)]]) /
                                                           NMA_output[[paste(i)]][["Obs_mean"]][[paste(j)]]) * 100)
        
        NMA_kraft[["SES"]][[paste(i)]][[paste(j)]] <- abs((P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]] -
                                                             NMA_output[[paste(i)]][["Sim_mean_beta"]][[paste(j)]]) /
                                                            NMA_output[[paste(i)]][["Sim_sd_beta"]][[paste(j)]])
        NMA_kraft[["SI"]][[paste(i)]][[paste(j)]] <- 100 - ((abs(P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]] -
                                                                   NMA_output[[paste(i)]][["Sim_mean_beta"]][[paste(j)]]) /
                                                               P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]]) * 100)
      }
    }
  },
  error = function(e) {
    on.exit(options(show.error.messages = FALSE))
    if (NRep == 1) {
      return(cat("\nError 01: At least one experimental group - time point has fewer than two replicates. Please remove these groups.\n"))
    } else {
      return(cat("\nError 02: Duplicate rows for an experimental group - time point. Please remove or relabel duplicates.\n"))
    }
  }
  )
  
  # Adjust p-values if ObSsIm is TRUE
  if (ObSsIm) {
    NMA_stat[["p_value"]][["Adjusted"]] <- p.adjust(
      NMA_stat[["p_value"]][["p_value"]],
      method = p.adj,
      n = length(NMA_stat[["p_value"]][["p_value"]])
    )
  }
  
  # Remove simulated data to save space if Keep_Data is FALSE
  if (!Keep_Data) {
    P_SN_NMA_Sim_data      <- vector(mode = "list", length = 0)
    P_SN_NMA_Observed_data <- vector(mode = "list", length = 0)
  }
  
  # Aggregate summary data using plyr and reshape2
  tryCatch({
    out_put <- plyr::ldply(NMA_output_SES, data.frame)
    colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
    SES_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")
    
    out_put <- plyr::ldply(NMA_output_SI, data.frame)
    colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
    SI_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")
    
    out_put <- plyr::ldply(NMA_kraft[["SES"]], data.frame)
    colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
    kraft_SES_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")
    
    out_put <- plyr::ldply(NMA_kraft[["SI"]], data.frame)
    colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
    kraft_SI_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")
    
    kraft <- list(kraft_SI_data, kraft_SES_data, NMA_kraft)
    DaTa_casted <- dcast(DaTa, as.formula(paste0(colnames(DaTa[1]),"+",colnames(DaTa[2]),"+",colnames(DaTa[3]),"~",colnames(DaTa[4]))), value.var = colnames(DaTa[5]))
    NMA_diversity_d <- as.data.frame(cbind(DaTa_casted[1:3], apply(DaTa_casted[4:length(DaTa_casted[1,])],1, function(x) vegan::specnumber(x)), apply(DaTa_casted[4:length(DaTa_casted[1,])],1, function(x) exp(vegan::diversity(x, index = "shannon"))), apply(DaTa_casted[4:length(DaTa_casted[1,])],1, function(x) vegan::diversity(x, index = "invsimpson"))))
    colnames(NMA_diversity_d) <- c(colnames(DaTa_casted[1:3]), "H0", "H1", "H2")
  }, error = function(e) {
    on.exit(options(show.error.messages = FALSE))
    return(cat("\nError 03: Unequal number of experimental groups for time points. Please try running the NMA function for each time point individually.\n"))
  })
  
  # Generate final NMA_Results list containig input data, results of calculations, standard effect size and stochastic intensity and diversity data
  tryCatch({
    NMA_parameters <- list(NSim, InDeX, ObSsIm, p.adj, Plot_Level, unique(DaTa$Time_point))
    
    NMA_Results <- list(
      DaTa,
      P_SN_NMA_Observed_data,
      P_SN_NMA_Sim_data,
      P_SN_NMA_Sim_diversity,
      P_SN_NMA_Obs_diversity,
      NMA_stat,
      SES_data,
      SI_data,
      kraft,
      NMA_parameters,
      NMA_diversity_d)
    
    names(NMA_Results) <- c(
      "Input_data",
      "Simulated_observed_data",
      "Simulated_data",
      "Simulated_diversity",
      "Observed_diversity",
      "Statistical_significance",
      "Standard_effect_size_data",
      "Stochastic_intensity_data",
      "SI_SES_species_number_data",
      "NMA_parameters",
      "Hill_table")
    
    return(NMA_Results)
  }, error = function(e) {
    on.exit(options(show.error.messages = Deb))
    if (Deb) {
      message(conditionMessage(e))
    } else {
      return(cat("\nFinished with error!\n"))
    }
  })
}

#' Summary method for NMA objects
#' @param object An object of class "NMA".
#'
#' @return Invisibly returns a table summarizing the Null model analysis results.
#' @examples
#' \dontrun{
#'   nma_results <- nma_results <- NMA(DaTa = NMA_data, NSim = 1000)
#'   summary.NMA(nma_results)
#' }
#'
#' @export
summary.NMA <- function(object) {
  cat("Summary of NMA Results:\nHigher standard effect size values show higher contribution of deterministic assembly mechanisms.\n")
  # Prints the summary of NMA analysis
  summarized_object <- object$Standard_effect_size_data[,-2]
  colnames(summarized_object) <- c("Experimental group", "Standard effect size")
  print(summarized_object)
  invisible(object)
}

#' @title Plotting function for the Null Model Analysis (MicroEcoTools)
#'
#' @description This function takes in the NMA_Results list and plots the requested plot by the user.
#'
#' @param NMA_Results A list returned by the NMA function.
#' @param plot_type A character string specifying the desired plot type. Valid options are:
#'   \itemize{
#'     \item \code{"SES"}: Standard effect size plot using the specified diversity index.
#'     \item \code{"SI"}: Stochastic intensity plot using the specified diversity index.
#'     \item \code{"SES_species_number"}: Standard effect size plot using species richness as diversity measure.
#'     \item \code{"SI_species_number"}: Stochastic intensity plot using species richness as diversity measure.
#'     \item \code{"Cohens_d"}: Standard effect size plot based on Cohen's d.
#'     \item \code{"Hill"}: A list of plots for zero-, first-, and second-order Hill diversity.
#'   }
#'   Defaults to \code{"SES"}.
#'
#' @examples
#' \dontrun{
#'   # Perform null model analysis on experimental data from a perturbation experiment.
#'   nma_results <- NMA(DaTa = NMA_data, NSim = 1000)
#' 
#'   # Generate and display the sorted abundance plot:
#'   NMA_plot(nma_results, plot_type = "SES")
#'
#'   # Generate and display the sorted adjusted P-value plot:
#'   NMA_plot(nma_results, plot_type = "SES_species_number")
#' }
#' @return A ggplot object corresponding to the chosen plot type. If \code{plot_type = "Hill"}, a list of ggplot objects is returned.
#'
#' @export
NMA_plot <- function(NMA_Results, plot_type = c("SES", "SI", "SES_species_number", "SI_species_number", "Cohens_d", "Hill")) {
  # Check that NMA_Results is a list and contains required elements
  if (!is.list(NMA_Results)) {
    warning("NMA_Results must be a list returned by the NMA function.")
    return(invisible(NULL))
  }
  
  required_elements <- c("NMA_parameters", "Standard_effect_size_data", "Stochastic_intensity_data",
                         "SI_SES_species_number_data", "Input_data", "Hill_table")
  missing_elements <- setdiff(required_elements, names(NMA_Results))
  if (length(missing_elements) > 0) {
    warning("The NMA_Results list is missing the following required elements: ", paste(missing_elements, collapse = ", "))
    return(invisible(NULL))
  }
  
  # Define allowed plot types
  allowed_plot_types <- c("SES", "SI", "SES_species_number", "SI_species_number", "Cohens_d", "Hill")
  
  
  # Validate plot_type input
  if (!(plot_type %in% allowed_plot_types)) {
    warning("Invalid plot_type. Please choose one of: ", paste(allowed_plot_types, collapse = ", "))
    return(invisible(NULL))
  }
  
  plot_type <- match.arg(plot_type)
  
  NSim       <- NMA_Results$NMA_parameters[[1]]
  InDeX      <- NMA_Results$NMA_parameters[[2]]
  Plot_Level <- NMA_Results$NMA_parameters[[5]]
  
  SES_data         <- NMA_Results$Standard_effect_size_data
  SI_data          <- NMA_Results$Stochastic_intensity_data
  species_number_SES_data   <- NMA_Results$SI_SES_species_number_data[[2]]
  species_number_SI_data    <- NMA_Results$SI_SES_species_number_data[[1]]
  DaTa             <- NMA_Results$DaTa
  NMA_diversity_d  <- NMA_Results$Hill_table
  NMA_stat <- NMA_Results$Statistical_significance
  
  # Construct the SES plot
  SES_plot <- ggplot2::ggplot(SES_data, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = value)) +
    ggplot2::geom_bar(ggplot2::aes(fill = variable), position = "dodge", stat = "identity", width = 0.5) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Time point")) +
    ggplot2::ggtitle(paste0("Standard effect size plot with ", NSim, " iterations and ", InDeX, " used as diversity index")) +
    ggplot2::xlab("Experimental group") +
    ggplot2::ylab("Standard effect size") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  # Construct the SI plot
  SI_plot <- ggplot2::ggplot(SI_data, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = value)) +
    ggplot2::geom_bar(ggplot2::aes(fill = variable), position = "dodge", stat = "identity", width = 0.5) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Time point")) +
    ggplot2::ggtitle(paste0("Stochastic intensity plot with ", NSim, " iterations and ", InDeX, " used as diversity index")) +
    ggplot2::xlab("Experimental group") +
    ggplot2::ylab("Stochastic intensity (%)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  # Construct the Kraft SES plot (species number)
  species_number_SES_plot <- ggplot2::ggplot(species_number_SES_data, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = value)) +
    ggplot2::geom_bar(ggplot2::aes(fill = variable), position = "dodge", stat = "identity", width = 0.5) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Time point")) +
    ggplot2::ggtitle(paste0("Standard effect size plot with ", NSim, " iterations and species richness as measure of diversity")) +
    ggplot2::xlab("Experimental group") +
    ggplot2::ylab("Standard effect size") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  # Construct the Kraft SI plot (species number)
  species_number_SI_plot <- ggplot2::ggplot(species_number_SI_data, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = value)) +
    ggplot2::geom_bar(ggplot2::aes(fill = variable), position = "dodge", stat = "identity", width = 0.5) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Time point")) +
    ggplot2::ggtitle(paste0("Stochastic intensity plot with ", NSim, " iterations and species richness as measure of diversity")) +
    ggplot2::xlab("Experimental group") +
    ggplot2::ylab("Stochastic intensity (%)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  # Construct the Cohen's d plot
  a <- as.data.frame(NMA_stat[[1]])
  Cohens_d_plot <- ggplot2::ggplot(data = a, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = Cohens_d)) +
    ggplot2::geom_bar(ggplot2::aes(fill = Time_point), position = "dodge", stat = "identity") +
    ggplot2::xlab("Experimental Group") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Time point")) +
    ggplot2::ggtitle(paste0("Standard effect size calculated based on Cohen's d method using ", InDeX, " as diversity index with ", NSim, " iterations")) +
    ggplot2::ylab("Cohen's d") +
    ggplot2::geom_hline(yintercept = c(-0.5, -0.2, 0.2)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  # Construct the Hill diversity plots
  H0 <- melt(NMA_diversity_d[1:4], id.vars = colnames(NMA_diversity_d)[1:3])
  H1 <- melt(NMA_diversity_d[c(-4, -6)], id.vars = colnames(NMA_diversity_d)[1:3])
  H2 <- melt(NMA_diversity_d[c(-4, -5)], id.vars = colnames(NMA_diversity_d)[1:3])
  
  NMA_diversity_H0 <- ggplot2::ggplot(H0, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = value)) +
    ggplot2::geom_boxplot(ggplot2::aes(color = Time_point)) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.4) +
    ggplot2::ggtitle("Alpha diversity: Zero-order Hill number") +
    ggplot2::xlab("Experimental group") +
    ggplot2::ylab("Zero-order Hill number") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Time point"))
  
  NMA_diversity_H1 <- ggplot2::ggplot(H1, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = value)) +
    ggplot2::geom_boxplot(ggplot2::aes(color = Time_point)) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.4) +
    ggplot2::ggtitle("Alpha diversity: First-order Hill number") +
    ggplot2::xlab("Experimental group") +
    ggplot2::ylab("First-order Hill number") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  NMA_diversity_H2 <- ggplot2::ggplot(H2, ggplot2::aes(x = factor(Exp.Grp, levels = Plot_Level), y = value)) +
    ggplot2::geom_boxplot(ggplot2::aes(color = Time_point)) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.4) +
    ggplot2::ggtitle("Alpha diversity: Second-order Hill number") +
    ggplot2::xlab("Experimental group") +
    ggplot2::ylab("Second-order Hill number") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  # Apply facet wrapping if needed (if the grouping variable has multiple levels)
  if (length(unique(SES_data$variable)) > 1) {
    SES_plot <- SES_plot + ggplot2::facet_wrap(~SES_data$variable) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
    SI_plot <- SI_plot + ggplot2::facet_wrap(~SI_data$variable) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
    species_number_SES_plot <- species_number_SES_plot + ggplot2::facet_wrap(~species_number_SES_data$variable) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
    species_number_SI_plot <- species_number_SI_plot + ggplot2::facet_wrap(~species_number_SI_data$variable) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
    Cohens_d_plot <- Cohens_d_plot + ggplot2::facet_wrap(~Time_point) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
    NMA_diversity_H0 <- NMA_diversity_H0 + ggplot2::facet_wrap(~Time_point) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
    NMA_diversity_H1 <- NMA_diversity_H1 + ggplot2::facet_wrap(~Time_point) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
    NMA_diversity_H2 <- NMA_diversity_H2 + ggplot2::facet_wrap(~Time_point) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
  }
  
  # Return the chosen plot based on the plot_type parameter
  if (plot_type == "SES") {
    return(SES_plot)
  } else if (plot_type == "SI") {
    return(SI_plot)
  } else if (plot_type == "SES_species_number") {
    return(species_number_SES_plot)
  } else if (plot_type == "SI_species_number") {
    return(species_number_SI_plot)
  } else if (plot_type == "Cohens_d") {
    return(Cohens_d_plot)
  } else if (plot_type == "Hill") {
    return(list(H0 = NMA_diversity_H0, H1 = NMA_diversity_H1, H2 = NMA_diversity_H2))
  }else return("Invalid choice of plot!")
}

