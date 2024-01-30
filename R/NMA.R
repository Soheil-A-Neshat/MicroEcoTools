#'@title  Theoretical Microbial Ecology Tools (MicroEcoTools)
#'
#'@description Using a null model approach, this package tests the relative contribution of the assembly mechanisms (stochastic and deterministic) in a given community. Based on the observed communities, assuming that the community structure was randomly shaped, using a multinomial distribution, this package simulates random communities. A diversity metric will then be used to calculate the standard effect size between the observed communities and randomly generated communities.
#'
#'@author Soheil A. Neshat, Ezequiel Santillan, Stefan Wuertz
#'@author (The NMA package was developed in Singapore Centre for Environmental Life Sciences Engineering, Nanyang Technological University, Singapore)
#'
#'@keywords Null model analysis, community ecology
#'
#'@param DaTa A data frame containing the environmental and community data in long format with column names of Exp.Grp, Reactor, Time_point, TAXA, and Count Please note that you need at least two replicates for each experimental group.
#'@details
#'The input data sample:
#'
#' | Exp.Grp | Reactor | Time_point | TAXA | Count |
#'|-----------|:-----------:|-----------:|:-----------:|-----------:|
#'  | I | 1 | 1 | rk1 | 68 |
#'  | I | 2 | 1 | rk1 | 70 |
#'  | X0 | 1 | 1 | rk1 | 3452 |
#'  | X0 | 3 | 1 | rk1 | 3274 |
#'  | X0 | 3 | 1 | rk1 | 3601 |
#'@md
#'
#'@param NSim A number to specify the number of simulations (values bellow 1000 are not recomended). Default to 1000 simulations.
#'@param InDeX A diversity index from vegan package ("shannon", "simpson" or "invsimpson"). Default to "invsimpson.
#'@param ObSsIm A dataset with multinomial distribution using the observed dataset as probability matrix will be generated.
#'@param p.adj A p-value correction method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Default to "Benjamini & Hochberg".
#'@param Plot_level A vector with sites/experimental group titles. Default to the unique values in experimental group column in your input data frame.
#'@param Keep_data A TRUE/FALSE variable to specify if you want to keep/remove simulated data. Default to FALSE.
#'
#'@returns A list named NMA_Results containing the data (observed and simulated communities), results of the statistical analysis on the diversity metric between observed and simulated communities, calculated standard effect sizes and stochastic intensities, standard effect size and stochastic intensity plots using the user specified diversity metric (InDeX parameter), standard effect size and stochastic intensity plots using richness (based on Kraft et al. and Santillan et al. papers), and standard effect size plot using Cohen's d method.
#'
#'@examples
#'NMA(DaTa = NMA_data)
#'
#'NMA(DaTa = NMA_data, NSim = 100, InDeX = "invsimpson", ObSsIm = FALSE, p.adj = "BH", Plot_Level = c("L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7"), Keep_Data = FALSE)
#'
#'@note All environmental variables should be converted to factor
#'@note Environmental variables accepted are Time_point, Exp.Grp, and Reactor
#'
#'
#'@references Santillan, Ezequiel, et al. "Frequency of disturbance alters diversity, function, and underlying assembly mechanisms of complex bacterial communities." npj Biofilms and Microbiomes 5.1 (2019): 1-9.
#'@references Kraft, Nathan JB, et al. "Disentangling the drivers of β diversity along latitudinal and elevational gradients." Science 333.6050 (2011): 1755-1758.
#'
#' @export NMA
#'
NMA <- function(DaTa, NSim, InDeX, ObSsIm, p.adj, Plot_Level, Keep_Data){
  if(missing(DaTa)){
    DaTa <- NMA::NMA_data
    print("NMA analysis using a demo dataset from Santillan et al. 2019.")
  }
  if(!is.data.frame(DaTa)){
    ?NMA
    on.exit(options(show.error.messages = FALSE))
    return(cat("Invalid data format. Please refer to the manual using ?NMA"))
  }

  if(dim(DaTa)[2] != 5){
    ?NMA
    return(cat("Invalid data format. Please refer to the manual using ?NMA"))
  }

  if(missing(NSim)){
    NSim <- 1000
    print("Default iteration number will be used (1000 iterations)")
    }
  if(NSim%%1 == 0){ print(paste0("Simulation with ", NSim, " iteration(s)"))
  } else {
  return("Invalid number of simulations entered. Please try again with an integer number as the number of simulations.")
  }
  if(missing(InDeX)){
    InDeX <- "invsimpson"
    print(paste0("The diversity index is set to inverse Simpson"))
  } else print(paste0("The diversity index is set to ", InDeX))
  if(missing(ObSsIm)){
    ObSsIm <- FALSE
  }
  if(ObSsIm == TRUE){
    print("A dataset with multinomial distribution using the observed dataset as probability matrix will be generated.")
  }
  if(missing(p.adj)){
    p.adj <- "BH"
  }
  if(missing(Plot_Level)){
    Plot_Level <- unique(DaTa$Exp.Grp)
    print("The default p-value correction method (Benjamini & Hochberg) will be used.")
  }
  if(missing(Keep_Data)){
    Keep_Data <- FALSE
  }

  `%>%` <- magrittr::`%>%`

  P_SN_NMA_Sim_diversity <- vector(mode = "list", length = 0)
  P_SN_NMA_Sim_data <- vector(mode = "list", length = 0)
  P_SN_NMA_Observed_data <- vector(mode = "list", length = 0)
  P_SN_NMA_Obs_diversity <- vector(mode = "list", length = 0)
  NMA_stat <- vector(mode = "list", length = 0)
  NMA_kraft <- vector(mode = "list", length = 0)
  NMA_output <- vector(mode = "list", length = 0)
  NMA_output_SES <- vector(mode = "list", length = 0)
  NMA_output_SI <- vector(mode = "list", length = 0)
  NMA_Results <<- vector(mode = "list", length = 0)
  kraft <- vector(mode = "list", length = 0)
  NSimO <- 0
  NMA_stat_p <- vector(length =0)

  if ("Time_point" %in% colnames(NMA_data)){cat()}else{DaTa$Time_point <- 1}


  DaTa$Exp.Grp <- as.factor(DaTa$Exp.Grp)
  DaTa$Time_point <- as.factor(DaTa$Time_point)
  DaTa$Reactor <- as.factor(DaTa$Reactor)
  DaTa$TAXA <- as.character(DaTa$TAXA)

  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = length(unique(DaTa$Exp.Grp)),
                         complete = "=",
                         incomplete = "-",
                         current = ">",
                         clear = FALSE,
                         width = 150)

  tImE <- system.time(tryCatch(
    {for (i in unique(DaTa$Exp.Grp)){

    pb$tick()

    P_SN_NMA_tmp1 <- DaTa %>% dplyr::filter(Exp.Grp == paste(i))

    for (j in unique(P_SN_NMA_tmp1$Time_point)){
      P_SN_NMA_tmp2 <- P_SN_NMA_tmp1 %>% dplyr::filter(Time_point == paste(j))

      NRep <- length(unique(P_SN_NMA_tmp2$Reactor))
      if (NRep < 2){
        print(paste("There is only one replicate for ", P_SN_NMA_tmp1[1,1], " experimental group at time point ", P_SN_NMA_tmp1[1,3], ". Please consider removing this experimental group - time point and try again."))
      }
      P_SN_NMA_tmp2 <- reshape2::dcast(data.table::as.data.table(P_SN_NMA_tmp2), Time_point+Exp.Grp+TAXA~Reactor, value.var = "Count", fun.aggregate = sum)
      expected <- P_SN_NMA_tmp2 %>% dplyr::select_if(is.numeric)
      expected <- expected[rowSums(expected) != 0,]

      NMA_stat[["p_value"]][["Exp.Grp"]] <- c(NMA_stat[["p_value"]][["Exp.Grp"]], paste(i))
      NMA_stat[["p_value"]][["Time_point"]] <- c(NMA_stat[["p_value"]][["Time_point"]], paste(j))

      P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]] <- lapply(seq(NSim), function(i, expected){
        .list <- lapply(apply(expected,1,list),unlist)
        as.data.frame(do.call(rbind,lapply(.list, function(.x) t(rmultinom(n = 1, prob = rep(sum(.x)/NRep, NRep),  size = sum(.x) )))))
      }, expected = expected)

      if (ObSsIm == TRUE){
        NSimO <- NSim - 1
        P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]] <- lapply(seq(NSimO), function(i, expected){
          .list <- lapply(apply(expected,1,list),unlist)
          as.data.frame(do.call(rbind,lapply(.list, function(.x) t(rmultinom(n = 1, prob = .x,  size = sum(.x) )))))
        }, expected = expected)}

      P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]] <- expected

      for (k in 1:NRep){
        for (l in 1:NSim){
          P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]] <- c(P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]],
                                                              1-(vegan::diversity(P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]][[l]][[k]], index = InDeX)/vegan::diversity(rowSums(P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]][[l]]), index = InDeX)))
          P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]][[paste(j)]] <- c(P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]][[paste(j)]], mean(as.numeric(lapply(P_SN_NMA_Sim_data[[paste(i)]][[paste(j)]][[l]],
                                                                                                                                                         function(x) vegan::specnumber(x))))/vegan::specnumber(rowSums(as.data.frame(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]]))))


        }
        for (m in 1:(NSimO+1)){
          P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]] <- c(P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]],
                                                              1-(vegan::diversity(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][[m]][[k]], index = InDeX)/vegan::diversity(rowSums(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][[m]]), index = InDeX)))
        }
      }


      P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]] <- c(P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]], mean(as.numeric(lapply(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]],
                                                                                                                                                     function(x) vegan::specnumber(x))))/vegan::specnumber(rowSums(as.data.frame(P_SN_NMA_Observed_data[[paste(i)]][[paste(j)]][["Observed"]]))))
      NMA_stat[["p_value"]][["Cohens_d"]] <- c(NMA_stat[["p_value"]][["Cohens_d"]], abs(effectsize::cohens_d(P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]],P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]])$Cohens_d))
      if (ObSsIm == TRUE){

        NMA_stat[["t-test"]][[paste(i)]][[paste(j)]] <- t.test(P_SN_NMA_Obs_diversity[[paste(i)]][[paste(j)]],P_SN_NMA_Sim_diversity[[paste(i)]][[paste(j)]])
        NMA_stat[["p_value"]][["p_value"]] <- c(NMA_stat[["p_value"]][["p_value"]], NMA_stat[["t-test"]][[paste(i)]][[paste(j)]]$p.value)

      }
    }
    NMA_output[[paste(i)]][["Sim_mean"]] <- lapply(P_SN_NMA_Sim_diversity[[paste(i)]] , mean)
    NMA_output[[paste(i)]][["Sim_sd"]] <- lapply(P_SN_NMA_Sim_diversity[[paste(i)]] , sd)
    NMA_output[[paste(i)]][["Obs_mean"]] <- lapply(P_SN_NMA_Obs_diversity[[paste(i)]] , mean)
    NMA_output[[paste(i)]][["Obs_sd"]] <- lapply(P_SN_NMA_Obs_diversity[[paste(i)]] , sd)
    NMA_output[[paste(i)]][["Sim_mean_beta"]] <- lapply(P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]] , mean)
    NMA_output[[paste(i)]][["Sim_sd_beta"]] <- lapply(P_SN_NMA_Sim_diversity[["beta"]][[paste(i)]] , sd)

    for (j in unique(names(NMA_output[[paste(i)]][["Obs_mean"]]))){
      NMA_output_SES[[paste(i)]][[paste(j)]] <- abs((NMA_output[[paste(i)]][["Obs_mean"]][[paste(j)]]-NMA_output[[paste(i)]][["Sim_mean"]][[paste(j)]])/NMA_output[[paste(i)]][["Sim_sd"]][[paste(j)]])
      NMA_output_SI[[paste(i)]][[paste(j)]] <- 100 - ((abs(NMA_output[[paste(i)]][["Obs_mean"]][[paste(j)]]-NMA_output[[paste(i)]][["Sim_mean"]][[paste(j)]]))/NMA_output[[paste(i)]][["Obs_mean"]][[paste(j)]] * 100)

      NMA_kraft[["SES"]][[paste(i)]][[paste(j)]] <- abs((P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]]-NMA_output[[paste(i)]][["Sim_mean_beta"]][[paste(j)]])/NMA_output[[paste(i)]][["Sim_sd_beta"]][[paste(j)]])
      NMA_kraft[["SI"]][[paste(i)]][[paste(j)]] <- 100 - ((abs(P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]]-NMA_output[[paste(i)]][["Sim_mean_beta"]][[paste(j)]]))/P_SN_NMA_Obs_diversity[["beta"]][[paste(i)]][[paste(j)]] * 100)

    }

  }}
  , error = function(e){
    on.exit(options(show.error.messages = FALSE))
    if(NRep == 1){
      return(cat("\n Error 01: At least one experimental group - time point with less than two replicates present in the dataset. \n Please consider removing experimental group - time point(s) with less than 2 replicates."))
    }else{
    return(cat("\n Error 02: Duplicate rows for an experimental group - time point. \n Please consider removing/relabling the duplicates or remove the experimental group - time points with less than 2 replicates."))
    }}))

  cat("\n")
  print(tImE)


  if (ObSsIm == TRUE){
    NMA_stat[["p_value"]][["Adjusted"]] <- p.adjust(NMA_stat[["p_value"]][["p_value"]], method = p.adj, n = length(NMA_stat[["p_value"]][["p_value"]]))
  }

  if (Keep_Data == FALSE){
    P_SN_NMA_Sim_data <- vector(mode = "list", length = 0)
    P_SN_NMA_Observed_data <- vector(mode = "list", length = 0)
  }
  tryCatch(
    {
  out_put <- plyr::ldply(NMA_output_SES,data.frame)
  colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
  SES_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")

  out_put <- plyr::ldply(NMA_output_SI,data.frame)
  colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
  SI_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")

  out_put <- plyr::ldply(NMA_kraft[["SES"]],data.frame)
  colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
  kraft_SES_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")

  out_put <- plyr::ldply(NMA_kraft[["SI"]],data.frame)
  colnames(out_put) <- c("Exp.Grp", unique(as.character(DaTa$Time_point)))
  kraft_SI_data <- reshape2::melt(data.table::as.data.table(out_put), id.vars = "Exp.Grp")

  kraft <- list(kraft_SI_data, kraft_SES_data, NMA_kraft)
    },error = function(e){
      on.exit(options(show.error.messages = FALSE))
      return(cat("\n Error 03: Unequal number of experimental groups for time points. \nPlease try running the NMA function for each time point individually."))
    }
  )

  tryCatch(
    {

  SES_plot <- ggplot2::ggplot(SES_data, ggplot2::aes(x= factor(Exp.Grp, levels = Plot_Level), y=value)) + ggplot2::geom_bar(ggplot2::aes(fill= variable), position = "dodge", stat = "identity", width = 0.5) + ggplot2::guides(fill = ggplot2::guide_legend(title="Time point")) +
    ggplot2::ggtitle(paste0("Standard effect size plot with ", NSim, " iterations and ", InDeX, " used as diversity index")) +  ggplot2::xlab("Experimental group") + ggplot2::ylab("Standard effect size") + ggplot2::theme(axis.title.x = ggplot2::element_text(size=18),
                                                                                                                                                                                         axis.title.y = ggplot2::element_text(size=18),
                                                                                                                                                                                         panel.background = ggplot2::element_blank(),
                                                                                                                                                                                         panel.grid.major = ggplot2::element_blank(),
                                                                                                                                                                                         panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                                         plot.background = ggplot2::element_blank())

  SI_plot <- ggplot2::ggplot(SI_data, ggplot2::aes(x= factor(Exp.Grp, levels = Plot_Level), y=value)) + ggplot2::geom_bar(ggplot2::aes(fill= variable), position = "dodge", stat = "identity", width = 0.5) + ggplot2::guides(fill = ggplot2::guide_legend(title="Time point")) +
    ggplot2::ggtitle(paste0("Stochastic intensity plot with ", NSim, " iterations and ", InDeX, " used as diversity index")) +  ggplot2::xlab("Experimental group") + ggplot2::ylab("Stochastic intensity (%)") + ggplot2::theme(axis.title.x = ggplot2::element_text(size=18),
                                                                                                                                                                                             axis.title.y = ggplot2::element_text(size=18),
                                                                                                                                                                                             panel.background = ggplot2::element_blank(),
                                                                                                                                                                                             panel.grid.major = ggplot2::element_blank(),
                                                                                                                                                                                             panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                                             plot.background = ggplot2::element_blank())

  Kraft_SES_plot <- ggplot2::ggplot(kraft_SES_data, ggplot2::aes(x= factor(Exp.Grp, levels = Plot_Level), y=value)) + ggplot2::geom_bar(ggplot2::aes(fill= variable), position = "dodge", stat = "identity", width = 0.5) + ggplot2::guides(fill = ggplot2::guide_legend(title="Time point")) +
    ggplot2::ggtitle(paste0("Standard effect size plot with ", NSim, " iterations and and species richness as measure of diversity")) +  ggplot2::xlab("Experimental group") + ggplot2::ylab("Standard effect size") + ggplot2::theme(axis.title.x = ggplot2::element_text(size=18),
                                                                                                                                                                    axis.title.y = ggplot2::element_text(size=18),
                                                                                                                                                                    panel.background = ggplot2::element_blank(),
                                                                                                                                                                    panel.grid.major = ggplot2::element_blank(),
                                                                                                                                                                    panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                    plot.background = ggplot2::element_blank())

  Kraft_SI_plot <- ggplot2::ggplot(kraft_SI_data, ggplot2::aes(x= factor(Exp.Grp, levels = Plot_Level), y=value)) + ggplot2::geom_bar(ggplot2::aes(fill= variable), position = "dodge", stat = "identity", width = 0.5) + ggplot2::guides(fill = ggplot2::guide_legend(title="Time point")) +
    ggplot2::ggtitle(paste0("Stochastic intensity plot with ", NSim, " iterations and species richness as measure of diversity")) +  ggplot2::xlab("Experimental group") + ggplot2::ylab("Stochastic intensity (%)") + ggplot2::theme(axis.title.x = ggplot2::element_text(size=18),
                                                                                                                                                                        axis.title.y = ggplot2::element_text(size=18),
                                                                                                                                                                        panel.background = ggplot2::element_blank(),
                                                                                                                                                                        panel.grid.major = ggplot2::element_blank(),
                                                                                                                                                                        panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                                                        plot.background = ggplot2::element_blank())

  a <- as.data.frame(NMA_stat[[1]])
  Cohens_d_plot <- ggplot2::ggplot(data = a, ggplot2::aes(x= factor(Exp.Grp, levels = Plot_Level), y=Cohens_d)) + ggplot2::geom_bar(ggplot2::aes(fill = Time_point), position = "dodge", stat = "identity") + ggplot2::xlab("Experimental Group") + ggplot2::guides(fill = ggplot2::guide_legend(title="Time point")) +
    ggplot2::ggtitle(paste0("Standard effect size calculated based on Cohen's d method using ",InDeX," as diversity index with ", NSim, " iterations")) + ggplot2::ylab("Cohen's d")  + ggplot2::geom_hline(yintercept = c(-0.5, -0.2, 0.2)) + ggplot2::theme(axis.title.x = ggplot2::element_text(size=18),
                                                                                                                                            axis.title.y = ggplot2::element_text(size=18),
                                                                                                                                            panel.background = ggplot2::element_blank(),
                                                                                                                                            panel.grid.major = ggplot2::element_blank(),
                                                                                                                                            panel.grid.minor = ggplot2::element_blank(),
                                                                                                                                            plot.background = ggplot2::element_blank())

  if (length(unique(SES_data$variable)) > 1){
    SES_plot <- SES_plot + ggplot2::facet_wrap(~SES_data$variable) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=0.95,vjust=0.2))
    SI_plot <- SI_plot + ggplot2::facet_wrap(~SI_data$variable) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=0.95,vjust=0.2))
    Kraft_SES_plot <- Kraft_SES_plot + ggplot2::facet_wrap(~kraft_SES_data$variable) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=0.95,vjust=0.2))
    Kraft_SI_plot <- Kraft_SI_plot + ggplot2::facet_wrap(~kraft_SES_data$variable) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=0.95,vjust=0.2))
    Cohens_d_plot <- Cohens_d_plot + ggplot2::facet_wrap(~Time_point) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=0.95,vjust=0.2))
  }
  NMA_parameters <- list(NSim, InDeX, ObSsIm, p.adj, Plot_Level, unique(DaTa$Time_point))

  NMA_Results <<- list(P_SN_NMA_Observed_data, P_SN_NMA_Sim_data, P_SN_NMA_Sim_diversity, P_SN_NMA_Obs_diversity, NMA_stat, SES_data, SI_data, kraft, SES_plot, SI_plot, Kraft_SES_plot, Kraft_SI_plot, Cohens_d_plot, NMA_parameters)
  names(NMA_Results) <<-c("Simulated_observed_data", "Simulated_data", "Simulated_diversity", "Observed_diversity", "Statistical_significance", "Standard_effect_size_data", "Stochastic_intensity_data", "Kraft_model_data", "SES_plot", "SI_plot", "Kraft_SES_plot", "Kraft_SI_plot", "Cohens_d_plot", "NMA_parameters")

  if(exists("NMA_Results")){
    print("The results are stored in NMA_Results list. You can see the plots using NMA_Results$SES_plot, NMA_Results$Kraft_SES_plot, and NMA_Results$Cohens_d_plot.")
  }
    },error = function(e){
      on.exit(options(show.error.messages = FALSE))
      return(cat("\nFinished with error! \n"))
}
)
  }


