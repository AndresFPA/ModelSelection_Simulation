# 21-11-2024
# Analysis without including the model K=1
# Load libraries
library(dplyr)

# Source evaluation
setwd("C:/Users/perezalo/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/All (50+50)")
load("Results/FinalResults.Rdata")
load("Results/design.Rdata")
source("evaluation.R")

# Set working directory
setwd("C:/Users/perezalo/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/All (50+50)/Fit")

# Load and re-analyze the data considering only models from 2 to 6
# Prepare new load function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Given that the name of the files can have "-1" or "-2" at the end, use regex to load them
ik <- 0
for(i in 1:144){
  for(k in 1:100){
    ik <- ik + 1
    # browser()
    file         <- list.files(pattern = paste0("^FitRow", i, "Rep", k, "-"), all.files = T)
    results      <- loadRData(file)
    results[1, ] <- NA # Remove the first row corresponding to K = 1
    res          <- vector(mode = "list", length = 1); res[[1]] <- results; names(res) <- "Overview" # This line used only so evaluation has correct input (needs a list)
    eval         <- evaluation(res = res, clus = design[i, "nclus"])
    Results_final[ik, 1:13] <- unlist(eval)
  }
}

colnames(Results_final) <- c("entropyR2",
                             "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
                             "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac",
                             "Replication", "Condition")

# setwd("C:/Users/perezalo/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/All (50+50)/Results")
# save(Results_final, file = "FinalResults2to6.Rdata")
# load("FinalResults2to6.Rdata")

# Merge datasets
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")
col_order <- c("Condition", "Replication", "nclus", "ngroups", "coeff", "N_g", "balance", "sd", 
               "entropyR2", "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
               "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac")
Results_final <- Results_final[, col_order]
rm(col_order)

# REMEMBER THAT THE RESULTS MUST TYPE "FACTOR"
measures <- Results_final %>% dplyr::select(Chull:ICL_fac) %>% as.matrix() %>% as.data.frame()
measures <- lapply(X = measures, FUN = factor, levels = c("-1", "0", "1"), labels = c("-1", "0", "1")) %>% as.data.frame()
Results_final[, 10:21] <- measures

####################################################################################################
############################ TABLES - CLUSTER AND PARAMETER RECOVERY ###############################
####################################################################################################
# Check mean results per simulation factor
# Create a function for this
count_results <- function(data, by, type = "count"){
  # Extract necessary columns
  reduced <- data %>% dplyr::select(Chull:ICL_fac)
  
  #Initialize objects to store
  counted        <- vector(mode = "list", length = ncol(reduced))
  names(counted) <- colnames(reduced)
  final <- c()
  # browser()
  # Count per column
  for(i in 1:ncol(reduced)){
    if(length(by) == 1 && by == "total"){
      counted[[i]] <- data %>% count(get(colnames(reduced)[i]), .drop = F) %>% filter(!is.na(`get(colnames(reduced)[i])`)) # Count and remove NA
      
      if(type == "relative"){
        counted[[i]][, ncol(counted[[i]])] <- round(counted[[i]][, ncol(counted[[i]]), drop = F]/sum(counted[[i]][, ncol(counted[[i]]), drop = F]), 3)
      }
      # browser()
      colnames(counted[[i]]) <- c("result", colnames(reduced)[i]) # change colnames
      ifelse(test = i == 1, yes = final <- counted[[i]][, 1, drop = F], no = final <- final) # for the first iteration, keep the group variable and results column
      final <- cbind(final, counted[[i]][, ncol(counted[[i]]), drop = F]) # Add the results for each measure
      
    } else {
      # browser()
      counted[[i]] <- data %>% group_by(across(all_of(by))) %>% count(get(colnames(reduced)[i]), .drop = F) %>% filter(!is.na(`get(colnames(reduced)[i])`)) # count per measure
      # browser()
      if(type == "relative"){
        counted[[i]][, ncol(counted[[i]])] <- round(counted[[i]][, ncol(counted[[i]])]/sum(counted[[i]][1:3, ncol(counted[[i]])]), 3)
      }
      
      colnames(counted[[i]]) <- c(by, "result", colnames(reduced)[i]) # change colnames
      ifelse(test = i == 1, yes = final <- counted[[i]][, c(seq_len(length(by)), (length(by) + 1))], no = final <- final) # for the first iteration, keep the group variable and results column
      final <- cbind(final, counted[[i]][, ncol(counted[[i]]), drop = F]) # Add the results for each measure
    }
    
  }
  return(final)
}

# Pre-check to know if there are NAs
colSums(apply(Results_final, 2, is.na))

# Main effects
count_results(data = Results_final, by = "total", type = "relative")

K_res   <- count_results(data = Results_final, by = c("nclus"), type = "relative")
N_res   <- count_results(data = Results_final, by = c("N_g"), type = "relative")
G_res   <- count_results(data = Results_final, by = c("ngroups"), type = "relative")
B_res   <- count_results(data = Results_final, by = c("coeff"), type = "relative")
Bal_res <- count_results(data = Results_final, by = c("balance"), type = "relative")
sd_res  <- count_results(data = Results_final, by = c("sd"), type = "relative")
