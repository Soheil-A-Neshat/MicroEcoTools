#' @title  Welch-ANOVA comparison (MicroEcoTools)
#'
#' @description Statistical comparison for traits/taxa between experimental groups using Welch-ANOVA or ANOVA. This function uses oneway test assuming non-equal variance. In cases where the assumptions are not valid it chooses equal variance by automatically setting var.equal parameter to TRUE.
#'
#' @param dAtA Dataframe that contains experimental group, replicates, and parameters.
#' @param var.name Variable name used for generating output table, for example, taxa, trait, function.
#' @param p_adj P-value adjustment method. All methods mentioned in ?p.adjust function can be used ("BH", "BY", "holm", "hommel" or "none").
#' @param Parallel Using this parameter you can make use of more CPU cores to decrease run time for the calculation. By default it is set to TRUE and uses half of the available CPU cores to perform the calculations.
#' @return This function returns a table containing the statistical comparison results.
#'@details
#'The input data sample:
#'
#' | Exp.Grp | Reactor | TAXA1 | TAXA2 | TAXA3 |
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
#' Welch_ANOVA(dAtA = CSR_IP2G_data, var.name = "taxa", p_adj = "BH", Parallel)
#' Welch_ANOVA(dAtA = CSR_TAXA_data, var.name = "taxa", p_adj = "BH", Parallel)
#' @export
Welch_ANOVA <- function(dAtA, var.name, p_adj, Parallel){
  if(missing(dAtA)) print("No data input!") else {
    a <- data.frame()
    if(missing(var.name)) var.name <- "Variable"
    if(missing(p_adj)) p_adj <- "BH"
    if(missing(Parallel)) Parallel <- TRUE
    total_iterations <- length(dAtA[1,]) - 2
    
    if(Parallel == TRUE){
      ncores <- parallel::detectCores() - 2
      if (ncores < 1) ncores <- 1 
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      message(paste0("  Parallel Welch-ANOVA comparisons with ", ncores, " CPU cores. Please be patient..."))
      
      a_result <- foreach::foreach (i = 1:total_iterations, .combine = rbind) %dopar% {
        f <- as.formula(noquote(paste0("`",names(dAtA)[i+2],"`","~",paste(names(dAtA)[1],collapse = " + "))))
        a[i,1] <- names(dAtA)[i+2]
        if (is.na(oneway.test(f, data = dAtA, var.equal = FALSE)$p.value)){
          a[i,2] <- oneway.test(f, data = dAtA, var.equal = TRUE)$p.value
          a[i,3] <- oneway.test(f, data = dAtA, var.equal = TRUE)$method
          a[i,4] <- "*"}else {
            a[i,2] <- oneway.test(f, data = dAtA, var.equal = FALSE)$p.value
            a[i,3] <- oneway.test(f, data = dAtA, var.equal = FALSE)$method
            a[i,4] <- "-"}
        Sys.sleep(1 / total_iterations)
        return(a[i, , drop = FALSE])
      }
      
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
      a <- a_result
      
    }else{
    pb <- progress::progress_bar$new(format = "  Welch ANOVA [:bar] :percent in :elapsed", clear = FALSE, total = total_iterations)

    for (i in 1:total_iterations){
      f <- as.formula(noquote(paste0("`",names(dAtA)[i+2],"`","~",paste(names(dAtA)[1],collapse = " + "))))
      a[i,1] <- names(dAtA)[i+2]
      if (is.na(oneway.test(f, data = dAtA, var.equal = FALSE)$p.value)){
        a[i,2] <- oneway.test(f, data = dAtA, var.equal = TRUE)$p.value
        a[i,3] <- oneway.test(f, data = dAtA, var.equal = TRUE)$method
        a[i,4] <- "*"}else {
          a[i,2] <- oneway.test(f, data = dAtA, var.equal = FALSE)$p.value
          a[i,3] <- oneway.test(f, data = dAtA, var.equal = FALSE)$method
          a[i,4] <- "-"}
      Sys.sleep(1 / total_iterations)
      pb$tick()
    }

    pb$terminate()
    }
    colnames(a) <- c(var.name, "Welch-ANOVA p-value", "Method", "Remarks_WA")
    a$`adjusted Welch-ANOVA p-value`<-p.adjust(a$`Welch-ANOVA p-value`, method = p_adj)
    return(a)
  }}

#' @title Pairwise Welch-ANOVA (MicroEcoTools)
#'
#' @description Pairwise comparison between experimental groups using Welch-ANOVA or ANOVA. This function uses oneway test assuming non-equal variance to perform pairwise comparison.
#' In cases where the assumptions are not valid it chooses equal variance by automatically setting var.equal parameter to TRUE it can be forced to use equal/not equal variance by setting v.equal to TRUE or FALSE.
#'
#' @param dAtA Dataframe that contains experimental group, replicates, and parameters.
#' @param var.name Variable name used for generating output table, for example, taxa, trait, function.
#' @param p_adj P-value adjustment method. All methods mentioned in ?p.adjust function can be used ("BH", "BY", "holm", "hommel" or "none").
#' @param v.equal Assumption of equal variance can be forced using this parameter. By default it is set to FALSE where it cannot be set to FALSE it will switch to TRUE automatically; a note will be shown in the remarks column in those cases.
#' @param p.value.cutoff A cut-off value for calling the P-value significant can be set using this parameter. The default value is 0.05.
#' @param Parallel Using this parameter you can make use of more CPU cores to decrease run time for the calculation. By default it is set to TRUE and uses half of the available CPU cores to perform the calculations.
#' @return This function returns a table containing the pairwise statistical comparison results. The output table can be fed into the CSR_assign function to assign CSR categories.
#' @details
#'The input data sample:
#'
#' | Exp.Grp | Reactor | TAXA1 | TAXA2 | TAXA3 |
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
#' pairwise_welch(dAtA = CSR_IP2G_data)
#' pairwise_welch(dAtA = CSR_IP2G_data, var.name = "taxa", p_adj = "BH")
#' pairwise_welch(dAtA = CSR_IP2G_data, var.name = "taxa", p_adj = "BH", v.equal = FALSE, p.value.cutoff = 0.05, Parallel = TRUE)
#' pairwise_welch(dAtA = CSR_TAXA_data)
#' pairwise_welch(dAtA = CSR_TAXA_data, var.name = "taxa", p_adj = "BH")
#' pairwise_welch(dAtA = CSR_TAXA_data, var.name = "taxa", p_adj = "BH", v.equal = FALSE, p.value.cutoff = 0.05, Parallel = TRUE)
#' @export
pairwise_welch <- function(dAtA, var.name, p_adj, v.equal, p.value.cutoff, Parallel){
  if(missing(dAtA)) print("No data input!") else {
    if(missing(var.name)) var.name <- "Variable"
    if(missing(p_adj)) p_adj <- "BH"
    if(missing(v.equal)) v.equal <- FALSE
    if(missing(p.value.cutoff)) p.value.cutoff <- 0.05
    if(missing(Parallel)) Parallel <- TRUE
    library(progress)
    print(Sys.time())
    if(rowSums(dAtA[1,c(-1,-2)] <= 1) | rowSums(dAtA[1,c(-1,-2)] <= 100)) {
      message("  Relative abundance data detected. Converting to count data by multiplying the relative abundances by 100,000 ...")
      dAtA[,c(-1,-2)] <- dAtA[,c(-1,-2)]*100000}
    dAtA <- dAtA[,-2]
    n <- (factorial(length(unique(dAtA[,1]))))/(2*factorial(length(unique(dAtA[,1]))-2))
    c <- t(combn(unique(dAtA[,1]),2))
    a <- data.frame()
    a <- rep(names(dAtA)[2:length(dAtA[1,])], each = n)
    a <- as.data.frame(a)
    a <- cbind(a, unlist(rep(as.character(c[,1],length(names(dAtA)[2:length(dAtA[1,])])))))
    a <- cbind(a, unlist(rep(as.character(c[,2],length(names(dAtA)[2:length(dAtA[1,])])))))


    if(Parallel == TRUE){
    message(paste("  Pairwise comparison of variables with variable name", var.name, ", P-value cutoff of" , p.value.cutoff, "in parallel mode."))
    ncores <- parallel::detectCores() - 2
    if (ncores < 1) ncores <- 1 
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    # Create a progress bar

    message("  Parallel processing of pairwise comparisons step 1/4. Please be patient...")

    a_result <- foreach::foreach (i = 1:length(a[,1]), .combine = rbind) %dopar% {

      if (v.equal == FALSE) {
        if (sd(dAtA[dAtA[1] == a[i,2], a[i,1]]) == 0 | sd(dAtA[dAtA[1] == a[i,3], a[i,1]]) == 0) v.equal <- TRUE
      }

      a[i,4] <- mean(dAtA[dAtA[1] == a[i,2], a[i,1]]) - mean(dAtA[dAtA[1] == a[i,3], a[i,1]])
      f <- as.formula(noquote(paste0("`", a[i,1], "`", "~", paste(names(dAtA)[1], collapse = " + "))))

      a[i,5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = v.equal)$statistic
      a[i,6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = v.equal)$parameter[1]
      a[i,7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = v.equal)$parameter[2]
      a[i,8] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = v.equal)$p.value
      a[i,9] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = v.equal)$method

      if (is.na(a[i,5])) {
        a[i,5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = TRUE)$statistic
        a[i,6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = TRUE)$parameter[1]
        a[i,7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = TRUE)$parameter[2]
        a[i,8] <- NA
        a[i,9] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA)[1], a[i,1])], var.equal = TRUE)$method
      }
      a[i,10] <- a[i,8]

      return(a[i, , drop = FALSE])
    }

    parallel::stopCluster(cl)
    foreach::registerDoSEQ()

    a <- a_result

    message("  Pairwise Welch ANOVA P-value adjustment step 2/4")
    for (i in unlist(unique(a[1]))){
      a[a[1] == i ,10] <- p.adjust(a[a[1] == i , 8], method = p_adj)
    }


    total_iterations <- length(a[,1])
    pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA effect size calculation step 3/4 [:bar] :percent in :elapsed", clear = FALSE, total = total_iterations)
    for (i in 1:length(a[,1])){
      if(!is.na(a[i,8])){
        if(a[i,10] < p.value.cutoff) a[i,11] <- "*" else a[i,11] <- "ns"} else a[i,11] <- "ns"
        pb$tick()
    }
    pb$terminate()

    message("  Pairwise Welch ANOVA effect size interpretation step 4/4. Please be patient...")

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    a <- foreach::foreach (i = 1:length(a[,1]), .combine = rbind) %dopar% {

      tryCatch({a[i,12] <- effectsize::interpret_hedges_g(effectsize::hedges_g(dAtA[dAtA[1] == a[i,2],  a[i,1]], dAtA[dAtA[1] == a[i,3],  a[i,1]], pooled_sd = FALSE)$Hedges_g)},
               error = function(e){a[i,12] <- "ne"})
      a[i,12][is.na(a[i,12])] <- "ne"

      return (a[i, , drop = FALSE])
    }

    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    print(Sys.time())

    } else {
      message(paste("Pairwise comparison of variables with variable name", var.name, ", P-value cutoff of" , p.value.cutoff, "using a single core. Please consider using parallel mode by setting Parallel parameter to TRUE: Parallel=TRUE"))
      total_iterations <- length(a[,1])
      pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA 1/4 [:bar] :percent in :elapsed", clear = FALSE, total = total_iterations)

      for (i in 1:length(a[,1])){
        if(v.equal == FALSE){
          if(sd(dAtA[dAtA[1] == a[i,2],  a[i,1]]) == 0 | sd(dAtA[dAtA[1] == a[i,3],  a[i,1]]) == 0) v.equal <- TRUE
        }
        a[i,4] <- mean(dAtA[dAtA[1] == a[i,2],  a[i,1]]) - mean(dAtA[dAtA[1] == a[i,3],  a[i,1]])
        f <- as.formula(noquote(paste0("`",a[i,1],"`","~",paste(names(dAtA)[1],collapse = " + "))))

        a[i,5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = v.equal)$statistic
        a[i,6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = v.equal)$parameter[1]
        a[i,7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = v.equal)$parameter[2]
        a[i,8] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = v.equal)$p.value
        a[i,9] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = v.equal)$method

        if(is.na(a[i,5])){
          a[i,5] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = TRUE)$statistic
          a[i,6] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = TRUE)$parameter[1]
          a[i,7] <- oneway.test(f, data = dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = TRUE)$parameter[2]
          a[i,8] <- NA
          a[i,9] <- oneway.test(f, data =dAtA[dAtA[1] == a[i,2] | dAtA[1] == a[i,3], c(names(dAtA[1]),a[i,1])], var.equal = TRUE)$method
        }

        pb$tick()
      }
      pb$terminate()

      a[,10] <- a[,8]


      message("  Pairwise Welch ANOVA P-value adjustment 2/4")
      for (i in unlist(unique(a[1]))){
        a[a[1] == i ,10] <- p.adjust(a[a[1] == i , 8], method = p_adj)
      }


      total_iterations <- length(a[,1])
      pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA effect size calculation 3/4 [:bar] :percent in :elapsed", clear = FALSE, total = total_iterations)
      for (i in 1:length(a[,1])){
        if(!is.na(a[i,8])){
          if(a[i,10] < 0.05) a[i,11] <- "*" else a[i,11] <- "ns"} else a[i,11] <- "ns"
          pb$tick()
      }
      pb$terminate()

      total_iterations <- length(a[,1])
      pb <- progress::progress_bar$new(format = "  Pairwise Welch ANOVA effect size interpretation 4/4 [:bar] :percent in :elapsed", clear = FALSE, total = total_iterations)
      for (i in 1:length(a[,1])){
        tryCatch({a[i,12] <- effectsize::interpret_hedges_g(effectsize::hedges_g(dAtA[dAtA[1] == a[i,2],  a[i,1]], dAtA[dAtA[1] == a[i,3],  a[i,1]], pooled_sd = FALSE)$Hedges_g)},
                 error = function(e){a[i,12] <- "ne"})
        a[i,12][is.na(a[i,12])] <- "ne"
        pb$tick()
      }
      pb$terminate()
      colnames(a) <- c(var.name, "P1", "P2", "Distance", "Statistic", "dF1", "dF2", "P-value", "Method", "Adj.P-value", "Significance", "Effect_size")
      print(Sys.time())
      return(a)
      }

    colnames(a) <- c(var.name, "P1", "P2", "Distance", "Statistic", "dF1", "dF2", "P-value", "Method", "Adj.P-value", "Significance", "Effect_size")
    return(a)
  }  }

#' @title Competitor - Stress tolerant - Ruderal assignment (MicroEcoTools)
#' @description Assigning CSR categories from Grime's framework to the input variable.
#' This function assigns CSR categories to the input variables using the output table from the ?pairwise_welch function.
#' To do so, it looks at the pairwise comparisons and assigns C to the ones abundant in no disturbance group, S to the ones abundant in high disturbance group, and R to the ones abundant in intermediate disturbance groups.
#' It also assigns intermediate categories, CR, CS, and SR to groups that are abundant in no-disturbance and intermediate-disturbance, no-disturbance and high-disturbance, and intermediate-disturbance and high-disturbance, respectively.
#' Traits that are not differentially abundant in any groups will be assigned to the CSR category.
#' Please note that if the variability is lacking the group will be categorised as NA.
#' @param dAtA Dataframe that contains experimental group, replicates, and parameters.
#' @param var.name Variable name used for generating output table, for example, taxa, trait, function.
#' @param p_adj P-value adjustment method. All methods mentioned in ?p.adjust function can be used ("BH", "BY", "holm", "hommel" or "none").
#' @param v.equal Assumption of equal variance can be forced using this parameter. By default it is set to FALSE where it cannot be set to FALSE it will switch to TRUE automatically; a note will be shown in the remarks column in those cases.
#' @param p.value.cutoff A cut-off value for calling the P-value significant can be set using this parameter. The default value is 0.05.
#' @param Parallel Using this parameter you can make use of more CPU cores to decrease run time for the calculation. By default it is set to TRUE and uses half of the available CPU cores to perform the calculations.
#' @param Vis Using this parameter you can make visualize the top 5 taxa/trait sorted based on P-value or abundance. The graphs will be saved as CSR_plot_sorted_p_value and CSR_plot_sorted_abundance.
#' @return This function returns a table containing the pairwise statistical comparison results. The output table can be fed into the CSR_assign function to assign CSR categories.
#' 
#'@details
#'The input data sample:
#'
#' | Exp.Grp | Reactor | TAXA1 | TAXA2 | TAXA3 |
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
#' CSR_assign(dAtA = CSR_IP2G_data)
#' CSR_assign(dAtA = CSR_IP2G_data, var.name = "taxa", p_adj = "BH")
#' CSR_assign(dAtA = CSR_IP2G_data, var.name = "taxa", p_adj = "BH", v.equal = FALSE, p.value.cutoff = 0.05, Parallel = TRUE)
#' CSR_assign(dAtA = CSR_TAXA_data)
#' CSR_assign(dAtA = CSR_TAXA_data, var.name = "taxa", p_adj = "BH")
#' CSR_assign(dAtA = CSR_TAXA_data, var.name = "taxa", p_adj = "BH", v.equal = FALSE, p.value.cutoff = 0.05, Parallel = TRUE)
#' @export
CSR_assign <- function(dAtA, var.name, p_adj, p.value.cutoff, Parallel, Vis) {
  if(missing(dAtA)) print("No data input!") else {
    if(missing(var.name)) var.name <- "Trait"
    if(missing(p_adj)) p_adj <- "BH"
    if(missing(p.value.cutoff)) p.value.cutoff <- 0.05
    if(missing(Parallel)) Parallel <- TRUE
    if(missing(Vis)) Vis <- TRUE
    j=1
    DaTa <- dAtA
    wa <- Welch_ANOVA(dAtA = dAtA, var.name = var.name, p_adj = "BH")
    if (colnames(dAtA[2]) != "P1"){
      message("  Creating the pairwise Welch-ANOVA table from the input data using pairwise_welch function ...")
      pww <- pairwise_welch(dAtA = dAtA, var.name = var.name, p_adj = p_adj, p.value.cutoff = p.value.cutoff, Parallel = Parallel)
      dAtA <- pww
    }
    
    
    filtered_dAtA <- dAtA[which(dAtA[,2]==unique(dAtA[,2])[1] | dAtA[,3]==unique(dAtA[,3])[length(unique(dAtA[,3]))]),]
    
    a <- data.frame(matrix(nrow = length(unique(filtered_dAtA[,1])), ncol = 3))
    rownames(a) <- unique(filtered_dAtA[,1])
    a[,1] <- rownames(a)
    a[,2] <- ""
    a[,3] <- "-"
    Exp.Grp <- unique(c(filtered_dAtA[,2],filtered_dAtA[,3]))
    n.Exp.Grp <- length(unique(c(filtered_dAtA[,2],filtered_dAtA[,3])))

    for (i in unique(filtered_dAtA[,1])){
      #Assign C, S and R groups
      dAtA_subset <- filtered_dAtA[which(filtered_dAtA[,1] == i),]

      if (n.Exp.Grp > 2) {
        if(all(dAtA_subset[dAtA_subset[2]==Exp.Grp[1],11] != "ns") & (all(dAtA_subset[dAtA_subset[2]==Exp.Grp[1],12] == "medium") | all(dAtA_subset[dAtA_subset[2]==Exp.Grp[1],12] == "large")) & all(dAtA_subset[dAtA_subset[2]==Exp.Grp[1],4] > 0)) {
          a[i,2] <- "C"} else if (all(dAtA_subset[dAtA_subset[3]==Exp.Grp[length(Exp.Grp)],11] != "ns") & (all(dAtA_subset[dAtA_subset[3]==Exp.Grp[length(Exp.Grp)],12] == "medium") | all(dAtA_subset[dAtA_subset[3]==Exp.Grp[length(Exp.Grp)],12] == "large")) & all(dAtA_subset[dAtA_subset[3]==Exp.Grp[length(Exp.Grp)],4] < 0)) {
            a[i,2] <- "S" } else if(all(dAtA_subset[1:(n.Exp.Grp-1),11] == "ns") & all(dAtA_subset[1:(n.Exp.Grp-1),12] != "ne"))a[i,2] <- "CSR" else if(all(dAtA_subset[1:(n.Exp.Grp-1),11] == "ns") & any(dAtA_subset[1:(n.Exp.Grp-1),12] == "ne")){a[i,2] <- "NA"
            a[i,3] <- "Cannot assign CSR probably due to lack of variability in the data. Please inspect the data."} else {
              dAtA_subset_CR <- dAtA_subset[which((dAtA_subset[,2] == Exp.Grp[1] & dAtA_subset[,3] != Exp.Grp[n.Exp.Grp])),]
              dAtA_subset_CR[,4] <- dAtA_subset_CR[,4] * -1
              dAtA_subset_SR <- dAtA_subset[which(dAtA_subset[,2] != Exp.Grp[1] & dAtA_subset[,3] == Exp.Grp[n.Exp.Grp]),]
              dAtA_subset_R <- rbind(dAtA_subset_CR, dAtA_subset_SR)
              dAtA_subset_R <- dAtA_subset_R[which(dAtA_subset_R[,11] == "*" & (dAtA_subset_R[,12] == "medium" | dAtA_subset_R[,12] == "large")),]
              if(any(dAtA_subset_R[,4] > 0) & !any(dAtA_subset_R[,4] < 0)) {
                a[i,2] <- "R"}
            }
      }else {
        message("CSR assignment with only two groups assuming these groups represent no disturbance and press disturbance.\nPlease note that intermediate groups, namely, CS, CR, and SR will be wrongly categorised.\nIf possible please include at least one experimental group representing the intermediate level of disturbance.")
        if(all(dAtA_subset[,11] != "ns") & (dAtA_subset[1,12] == "medium" | dAtA_subset[1,12] == "large") & dAtA_subset[1,4] > 0) {
          a[i,2] <- "C"} else if (dAtA_subset[1,11] != "ns" & (dAtA_subset[1,12] == "medium" | dAtA_subset[1,12] == "large") & dAtA_subset[1,4] < 0) {
            a[i,2] <- "S" } else if(dAtA_subset[1,11] == "ns" & dAtA_subset[1,12] != "ne") a[i,2] <- "CSR" else if(dAtA_subset[1,11] == "ns" & dAtA_subset[1,12] == "ne"){a[i,2] <- "NA"
            a[i,3] <- "Cannot assign CSR probably due to lack of variability in the data. Please inspect the data."}
        }}

    dAtA_remain <- a[,1][which(a[,2] == "")]

    for (ii in dAtA_remain){
      dAtA_subset <- filtered_dAtA[which(filtered_dAtA[,1] == ii),]
      dAtA_subset_CR <- dAtA_subset[which((dAtA_subset[,2] == Exp.Grp[1] & dAtA_subset[,3] != Exp.Grp[n.Exp.Grp]) & dAtA_subset[,11] != "ns"),]
      dAtA_subset_SR <- dAtA_subset[which(dAtA_subset[,2] != Exp.Grp[1] & dAtA_subset[,3] == Exp.Grp[n.Exp.Grp] & dAtA_subset[,11] != "ns"),]
      dAtA_subset_CS <- dAtA_subset[which(dAtA_subset[,2] == Exp.Grp[1] & dAtA_subset[,3] == Exp.Grp[n.Exp.Grp] & dAtA_subset[,11] != "ns"),]
      if (length(dAtA_subset_CR[,4]) > 0 & length(dAtA_subset_CS[,4]) > 0){
      if(any(dAtA_subset_CR[,4] < 0) & !any(dAtA_subset_CR[,4] > 0)  & dAtA_subset_CS[,4] < 0) {
        a[ii,2] <- "SR"} }
      if(length(dAtA_subset_CS[,4]) > 0 & length(dAtA_subset_SR[,4]) > 0) {
      if(all(dAtA_subset_CS[,4] > 0) & any(dAtA_subset_SR[,4] > 0) & !any(dAtA_subset_SR[,4] < 0)) {
          a[ii,2] <- "CR"} }
      if(length(dAtA_subset_CR[,4]) > 0 & length(dAtA_subset_SR[,4]) > 0 & length(dAtA_subset_CS[,11] == 0)) {
      if(all(dAtA_subset_CR[,4] > 0) & all(dAtA_subset_SR[,4] < 0)) {
            a[ii,2] <- "CS"} }
      if (a[ii,2] == "") {
        a[ii,2] <- "CSR"}


            }

    colnames(a) <- c(paste(var.name), "CSR categories", "Remarks_CSR")
    a <- merge(wa[, c(1, 2, 3, 5)], a, all.x = TRUE)
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
    if(Vis == TRUE){
      vis_data_C <- a[which(a[,5] == "C" ),]
      if(length (vis_data_C[,1]) > 0){
      j = 1
      for (i in vis_data_C[,1]){
        vis_data_C[j,"mean_rel"] <- mean(dAtA[[i]], na.rm = TRUE)
        j=j+1}
      vis_data_C_p <- vis_data_C[order(vis_data_C$`adjusted Welch-ANOVA p-value`),]
      vis_data_C_r <- vis_data_C[order(vis_data_C$`mean_rel`, decreasing = TRUE),]
      if (nrow(vis_data_C) > 5) {
        vis_data_C_p <- vis_data_C_p[1:5, ]
        vis_data_C_r <- vis_data_C_r[1:5, ]}
      }
      vis_data_S <- a[which(a[,5] == "S" ),]
      if(length (vis_data_S[,1]) > 0){
      j = 1
      for (i in vis_data_S[,1]){
        vis_data_S[j,"mean_rel"] <- mean(dAtA[[i]], na.rm = TRUE)
        j=j+1}
      vis_data_S_p <- vis_data_S[order(vis_data_S$`adjusted Welch-ANOVA p-value`),]
      vis_data_S_r <- vis_data_S[order(vis_data_S$`mean_rel`, decreasing = TRUE),]
      if (nrow(vis_data_S) > 5) {
        vis_data_S_p <- vis_data_S_p[1:5, ]
        vis_data_S_r <- vis_data_S_r[1:5, ]}
      }
      vis_data_R <- a[which(a[,5] == "R" ),]
      if(length (vis_data_R[,1]) > 0){
      j = 1
      for (i in vis_data_R[,1]){
        vis_data_R[j,"mean_rel"] <- mean(dAtA[[i]], na.rm = TRUE)
        j=j+1}
      vis_data_R_p <- vis_data_R[order(vis_data_R$`adjusted Welch-ANOVA p-value`),]
      vis_data_R_r <- vis_data_R[order(vis_data_R$`mean_rel`, decreasing = TRUE),]
      if (nrow(vis_data_R) > 5) {
        vis_data_R_p <- vis_data_R_p[1:5, ]
        vis_data_R_r <- vis_data_R_r[1:5, ]}
      }
      vis_data_CS <- a[which(a[,5] == "CS" ),]
      if(length (vis_data_CS[,1]) > 0){
      j = 1
      for (i in vis_data_CS[,1]){
        vis_data_CS[j,"mean_rel"] <- mean(dAtA[[i]], na.rm = TRUE)
        j=j+1}
      vis_data_CS_p <- vis_data_CS[order(vis_data_CS$`adjusted Welch-ANOVA p-value`),]
      vis_data_CS_r <- vis_data_CS[order(as.vector(vis_data_CS$`mean_rel`), decreasing = TRUE),]
      if (nrow(vis_data_CS) > 5) {
        vis_data_CS_p <- vis_data_CS_p[1:5, ]
        vis_data_CS_r <- vis_data_CS_r[1:5, ]}
      }
      vis_data_SR <- a[which(a[,5] == "SR" ),]
      if(length (vis_data_SR[,1]) > 0){
      j = 1
      for (i in vis_data_SR[,1]){
        vis_data_SR[j,"mean_rel"] <- mean(dAtA[[i]], na.rm = TRUE)
        j=j+1}
      vis_data_SR_p <- vis_data_SR[order(vis_data_SR$`adjusted Welch-ANOVA p-value`),]
      vis_data_SR_r <- vis_data_SR[order(vis_data_SR$`mean_rel`, decreasing = TRUE),]
      if (nrow(vis_data_SR) > 5) {
        vis_data_SR_p <- vis_data_SR_p[1:5, ]
        vis_data_SR_r <- vis_data_SR_r[1:5, ]}
      }
      vis_data_CR <- a[which(a[,5] == "CR" ),]
      if(length (vis_data_CR[,1]) > 0){
      j = 1
      for (i in vis_data_CR[,1]){
        vis_data_CR[j,"mean_rel"] <- mean(dAtA[[i]], na.rm = TRUE)
        j=j+1}
      vis_data_CR_p <- vis_data_CR[order(vis_data_CR$`adjusted Welch-ANOVA p-value`),]
      vis_data_CR_r <- vis_data_CR[order(vis_data_CR$`mean_rel`, decreasing = TRUE),]
      if (nrow(vis_data_CR) > 5) {
        vis_data_CR_p <- vis_data_CR_p[1:5, ]
        vis_data_CR_r <- vis_data_CR_r[1:5, ]}
      }
      vis_data_list_p <- rbind(vis_data_C_p, vis_data_S_p, vis_data_R_p, vis_data_CS_p, vis_data_CR_p, vis_data_SR_p)[,1]
      vis_data_list_r <- rbind(vis_data_C_r, vis_data_S_r, vis_data_R_r, vis_data_CS_r, vis_data_CR_r, vis_data_SR_r)[,1]
      
      CSR_Plot <- function(DaTa, CSR_vis_list, CSR_cat, var.name, sort.var){
        CSR_plot_variable <- var.name
        vis_data_m <- DaTa[,c(1,2)]
        
        vis_data <- DaTa %>% select(all_of(CSR_vis_list))
        vis_data <- cbind(vis_data_m, vis_data)
        CSR_cat_filtered <- CSR_cat[,c(1,5)]
        vis_data_0 <- reshape2::melt(vis_data, id.vars = c(colnames(vis_data[1]), colnames(vis_data[2])), variable.name = CSR_plot_variable, value.name = "Abundance")
        str(vis_data_0)
        str(CSR_cat_filtered)
        vis_data <- merge(vis_data_0, CSR_cat_filtered, all.x = TRUE)
        str(vis_data)
        names(vis_data) <- c(var.name, "Exp.Grp", "Rep", "Abundance", "CSR_categories")
        vis_data$Exp.Grp <- as.factor(vis_data$Exp.Grp)
        vis_data <-vis_data
        
        CSR_plot <- ggplot2::ggplot(vis_data, aes(x = Exp.Grp, y = Abundance, color = CSR_categories)) +
          geom_point() + theme_minimal() + theme(
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "transparent")) +
          labs(title = paste( "Top 5", CSR_plot_variable, "from C, S, R and intermediate categories sorted based on", sort.var)) + facet_wrap(as.formula(paste("~", CSR_plot_variable)), scales = "free")
        CSR_plot
      }
      
      CSR_plot_sorted_abundance <<- CSR_Plot(DaTa = DaTa, CSR_cat = a, CSR_vis_list = vis_data_list_r, var.name = var.name,  sort.var = "abundance")
      CSR_plot_sorted_p_value <<- CSR_Plot(DaTa = DaTa, CSR_cat = a, CSR_vis_list = vis_data_list_p, var.name = var.name, sort.var = "adjusted P-value")
    }

    return(a)
  }}



#' @title Theoretical Microbial Ecology Tools (MicroEcoTools)
#' 
#' @description Probability prediction for the CSR assignments.
#' Assigning CSR categories from Grime's framework to the input variable and calculates the probability of assigning these categories using a simulation.
#' This function generates simulated datasets using the input dataset that have a similar multinomial distribution.
#' It is advised to use at least 1000 simulations to get more reliable results (recomended NSim = 10000).
#' After performing the simulation it passes the data to the CSR_assign function to assign CSR categories.
#'
#' This function assigns CSR categories to the input variables using the output table from the ?pairwise_welch function.
#' To do so, it looks at the pairwise comparisons and assigns C to the ones abundant in no disturbance group, S to the ones abundant in high disturbance group, and R to the ones abundant in intermediate disturbance groups.
#' It also assigns intermediate categories, CR, CS, and SR to groups that are abundant in no-disturbance and intermediate-disturbance, no-disturbance and high-disturbance, and intermediate-disturbance and high-disturbance, respectively.
#' Traits that are not differentially abundant in any groups will be assigned to the CSR category.
#' Please note that if the variability is lacking the group will be categorised as NA.
#' 
#' @param DaTa Dataframe that contains experimental group, replicates, and parameters.
#' @param NSim Number of simulations minimum 1000, recomended 10000.
#' @param var.name Variable name used for generating output table, for example, taxa, trait, function.
#' @param p_adj P-value adjustment method. All methods mentioned in ?p.adjust function can be used ("BH", "BY", "holm", "hommel" or "none").
#' @param v.equal Assumption of equal variance can be forced using this parameter. By default it is set to FALSE where it cannot be set to FALSE it will switch to TRUE automatically; a note will be shown in the remarks column in those cases.
#' @param p.value.cutoff A cut-off value for calling the P-value significant can be set using this parameter. The default value is 0.05.
#' @param Parallel Using this parameter you can make use of more CPU cores to decrease run time for the calculation. By default it is set to TRUE and uses half of the available CPU cores to perform the calculations.
#' @param Keep_data You can store the entire calculation results including the simulated communities in a list named CSR_Sim by setting this parameter to TRUE. By default it deletes the simulated data to save space.
#' @return This function returns a table containing the pairwise statistical comparison results. The output table can be fed into the CSR_assign function to assign CSR categories.
#'
#'@details
#'The input data sample:
#'
#' | Exp.Grp | Reactor | TAXA1 | TAXA2 | TAXA3 |
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
#' CSR_Simulation(DaTa = CSR_IP2G_data)
#' CSR_Simulation(DaTa = CSR_IP2G_data, NSim = 10000, p_adj = "BH")
#' CSR_Simulation(DaTa = CSR_IP2G_data, NSim = 10000, p_adj = "BH", var.name = "TAXA", NuLl.test = FALSE, Keep_data = FALSE)
#' CSR_Simulation(DaTa = CSR_TAXA_data)
#' CSR_Simulation(DaTa = CSR_TAXA_data, NSim = 10000, p_adj = "BH")
#' CSR_Simulation(DaTa = CSR_TAXA_data, NSim = 10000, p_adj = "BH", var.name = "TAXA", NuLl.test = FALSE, Keep_data = FALSE)

#' @export
CSR_Simulation <- function(DaTa, NSim, p_adj, p.value.cutoff, Parallel, var.name, NuLl.test, Keep_data){

  if(!"dplyr" %in% (.packages())){
    require("dplyr")
  }
  if(!"magrittr" %in% (.packages())){
    require("magrittr")
  }
  if(!"tidyverse" %in% (.packages())){
    require("tidyverse")
  }
  if(missing(DaTa)) print("No data input!") else {
    message(paste("CSR assignment with:"))
    if(missing(NSim)) NSim <- 1000
    message(paste("  number of simulation =", NSim))
    if(missing(var.name)) var.name <- "Trait"
    message(paste("  variable name =",  var.name))
    if(missing(p_adj)) p_adj <- "BH"
    message(paste("  P-value correction method =",  p_adj))
    if(missing(NuLl.test)) NuLl.test <- FALSE
    if(missing(Keep_data)) Keep_data <- FALSE
    if(missing(NuLl.test)) v.equal <- FALSE
    if(missing(Keep_data)) p.value.cutoff <- 0.05
    if(missing(Parallel)) Parallel <- TRUE
    if(rowSums(DaTa[1,c(-1,-2)] <= 1) | rowSums(DaTa[1,c(-1,-2)] <= 100)) DaTa[,c(-1,-2)] <- DaTa[,c(-1,-2)]*100000

    message(paste("  CSR assignment with", NSim, "simulated datasets and", p_adj, "as P-value correction method\n"))
    CSR_Sim <- vector(mode = "list", length = 0)
    CSR_Sim[["data"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["merged_data"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["CSR"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["Verdict"]] <- vector(mode = "list", length = 0)
    CSR_Sim[["Final_Verdict"]] <- vector(mode = "list", length = 0)
    
    DaTa <- DaTa[complete.cases(DaTa),]
    InPuT <- DaTa
    InPuT <- as.data.frame(t(InPuT))
    InPuT <- InPuT[c(-1,-2),]
    expected <- InPuT
    expected <- apply(expected,2, function(x) as.numeric(x))
    expected <- as.data.frame(expected)
    expected <- ceiling(expected)
    rownames(expected) <- rownames(InPuT)
    expected <- expected[which(rowSums(expected) > 0),]

    message("Starting the simulation\n")
    if (NuLl.test == FALSE){
      CSR_Sim[["data"]] <- lapply(seq(NSim), function(i, expected){
        # create a list of expected cell counts (list element = row of expected)
        .list <- lapply(apply(expected,1,list),unlist)
        # sample from these expected cell counts and recombine into a data.frame
        as.data.frame(do.call(rbind,lapply(.list, function(.x) t(rmultinom(n = 1, prob = .x,  size = sum(.x) )))))
      }, expected = expected)
      CSR_Sim[["data"]][["Obs"]] <- expected}

    if (NuLl.test == TRUE){
      CSR_Sim[["data"]] <- lapply(seq(NSim), function(i, expected){
        # create a list of expected cell counts (list element = row of expected)
        .list <- lapply(apply(expected,1,list),unlist)
        # sample from these expected cell counts and recombine into a data.frame
        as.data.frame(do.call(rbind,lapply(.list, function(.x) t(rmultinom(n = 1, prob = rep(sum(.x)/12, 12),  size = sum(.x) )))))
      }, expected = expected)}

    message("  Simulation            ","\u2714","\n")
    message("CSR assignment step, please be patient...")
    for (i in 1:length(CSR_Sim[["data"]])){
      message(paste("\n  Step", i, "out of", length(CSR_Sim[["data"]])))
      rownames(CSR_Sim[[1]][[i]]) <- rownames(expected)
      CSR_Sim[[2]][[i]] <- as.data.frame(t(CSR_Sim[[1]][[i]]))
      CSR_Sim[[2]][[i]] <- cbind(DaTa[1],DaTa[2],CSR_Sim[[2]][[i]])
      CSR_Sim[[3]][[i]] <- CSR_assign(dAtA = CSR_Sim[[2]][[i]], var.name = var.name, p_adj = p_adj ,p.value.cutoff = p.value.cutoff,Parallel = Parallel, Vis = FALSE)[c(1,5)]
      if (i == length(CSR_Sim[["data"]])) message("\nCSR assignment        ", "\u2713")
    }
    message("CSR assignment            ","\u2714","\n")
    CSR_Sim[["Verdict"]] <- CSR_Sim[[3]] %>% reduce(full_join, by=var.name)
    message(paste("Counting CSR categories for each ",var.name,"\n"))
    for (i in 1:length(CSR_Sim[["Verdict"]][,1])){
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+2] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])C\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+3] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])R\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+4] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])S\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+5] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])CR\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+6] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])CS\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+7] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])SR\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+8] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])CSR\\b"))
      CSR_Sim[["Verdict"]][i, length(CSR_Sim[["data"]])+9] <- sum(str_count(CSR_Sim[["Verdict"]][i,-1][1:(NSim+1)], "(?<![\\S])NA\\b"))
      if (i == length(CSR_Sim[["Verdict"]][,1])) message("CSR category count    ","\u2714","\n")
    }
    message("Counting CSR categories        ","\u2714","\n")
    CSR_Sim[["Final_Verdict"]] <- CSR_Sim[["Verdict"]][c(1,(length(CSR_Sim[["data"]])+2):(length(CSR_Sim[["data"]])+9))]
    colnames(CSR_Sim[["Final_Verdict"]]) <- c(var.name, "C", "R", "S", "CR", "CS", "SR", "CSR", "NA")
    CSR_Sim[["Final_Verdict"]] <- mutate_if(CSR_Sim[["Final_Verdict"]], is.numeric, ~ . /(NSim+1)* 100)

    CSR_Sim[["Final_Verdict"]] <- merge(CSR_assign(DaTa), CSR_Sim[["Final_Verdict"]], all.x = TRUE)

    if(Keep_data == TRUE){
      CSR_Sim <<- CSR_Sim
    }
    return(CSR_Sim[["Final_Verdict"]])
  }}
