library(dplyr)

# Source evaluation
source("C:/Users/User/Documents/GitHub/ModelSelection_Simulation/evaluation.R")
load("C:/Users/User/Documents/GitHub/ModelSelection_Simulation/Results/FinalResults.Rdata")
load("~/GitHub/ModelSelection_Simulation/Results/design.Rdata")

# Set working directory
setwd("C:/Users/User/Documents/GitHub/ModelSelection_Simulation/Results/Fit")

# Load and re-analyze the data considering only models from 2 to 6
# Given that the name of the files can have "-1" or "-2" at the end, use regex to load them
ik <- 0
for(i in 1:144){
  for(k in 1:50){
    ik <- ik + 1
    # browser()
    file <- list.files(pattern = paste0("^FitRow", i, "Rep", k, "-"), all.files = T)
    load(file)
    Overview[1, ] <- NA
    res  <- vector(mode = "list", length = 1); res[[1]] <- Overview; names(res) <- "Overview" # This line used only so evaluation has correct input
    eval <- evaluation(res = res, clus = design[i, "nclus"])
    Results_final[ik, 1:13] <- unlist(eval)
  }
}

colnames(Results_final) <- c("entropyR2",
                             "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
                             "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac",
                             "Replication", "Condition")

save(Results_final, file = "FinalResults2ndBest.Rdata")
load("~/GitHub/ModelSelection_Simulation/Results/FinalResults2ndBest.Rdata")

# Merge datasets
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")
col_order <- c("Condition", "Replication", "nclus", "ngroups", "coeff", "N_g", "balance", "sd", 
               "entropyR2", "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
               "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac")
Results_final <- Results_final[, col_order]
rm(col_order)

# Make the results factor
measures <- Results_final %>% dplyr::select(Chull:ICL_fac) %>% as.matrix() %>% as.data.frame()

measures <- lapply(X = measures, FUN = factor, levels = c("-1", "0", "1"), labels = c("Under", "Correct", "Over")) %>% as.data.frame()

Results_final[, 10:21] <- measures
Results_final[, "entropyR2"] <- as.numeric(Results_final[, "entropyR2"])

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
