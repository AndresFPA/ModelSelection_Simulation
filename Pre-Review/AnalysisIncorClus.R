library(dplyr)

# Source evaluation
source("C:/Users/User/Documents/GitHub/ModelSelection_Simulation/evaluationIncorClus.R")
load("C:/Users/User/Documents/GitHub/ModelSelection_Simulation/Results/FinalResults.Rdata")
load("~/GitHub/ModelSelection_Simulation/Results/design.Rdata")

# Set working directory
setwd("C:/Users/User/Documents/GitHub/ModelSelection_Simulation/Results/Fit")

# Load and re-analyze the data searching the second best model
# Given that the name of the files can have "-1" or "-2" at the end, use regex to load them
ik <- 0
for(i in 1:144){
  for(k in 1:50){
    ik <- ik + 1
    # browser()
    file <- list.files(pattern = paste0("^FitRow", i, "Rep", k, "-"), all.files = T)
    load(file)
    res  <- vector(mode = "list", length = 1); res[[1]] <- results_ov; names(res) <- "Overview" # This line used only so evaluation has correct input
    eval <- evaluation(res = res, clus = design[i, "nclus"])
    Results_final[ik, 1:13] <- unlist(eval)
  }
}

colnames(Results_final) <- c("entropyR2",
                             "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
                             "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac",
                             "Replication", "Condition")

save(Results_final, file = "FinalResultsIncorClus.Rdata")
load("~/GitHub/ModelSelection_Simulation/Results/FinalResultsIncorClus.Rdata")

# Merge datasets
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")
col_order <- c("Condition", "Replication", "nclus", "ngroups", "coeff", "N_g", "balance", "sd", 
               "entropyR2", "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
               "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac")
Results_final <- Results_final[, col_order]
rm(col_order)

# # Make the results factor
# measures <- Results_final %>% dplyr::select(Chull:ICL_fac) %>% as.matrix() %>% as.data.frame()
# 
# measures <- lapply(X = measures, FUN = factor, levels = c("-1", "0", "1"), labels = c("Under", "Correct", "Over")) %>% as.data.frame()
# 
# Results_final[, 10:21] <- measures
# Results_final[, "entropyR2"] <- as.numeric(Results_final[, "entropyR2"])

####################################################################################################
############################ TABLES - CLUSTER AND PARAMETER RECOVERY ###############################
####################################################################################################
# Get mean results for under- and over-selection
# Over-selection
Results_final %>% filter(sign(Chull) == 1) %>% dplyr::select(Chull) %>% colMeans()
Results_final %>% filter(sign(AIC)   == 1) %>% dplyr::select(AIC)   %>% colMeans()
Results_final %>% filter(sign(AIC3)  == 1) %>% dplyr::select(AIC3)  %>% colMeans()
Results_final %>% filter(sign(BIC_G) == 1) %>% dplyr::select(BIC_G) %>% colMeans()
Results_final %>% filter(sign(BIC_N) == 1) %>% dplyr::select(BIC_N) %>% colMeans()
Results_final %>% filter(sign(ICL)   == 1) %>% dplyr::select(ICL)   %>% colMeans()

# Under-selection
Results_final %>% filter(sign(Chull) == -1) %>% dplyr::select(Chull) %>% colMeans()
Results_final %>% filter(sign(AIC)   == -1) %>% dplyr::select(AIC)   %>% colMeans()
Results_final %>% filter(sign(AIC3)  == -1) %>% dplyr::select(AIC3)  %>% colMeans()
Results_final %>% filter(sign(BIC_G) == -1) %>% dplyr::select(BIC_G) %>% colMeans()
Results_final %>% filter(sign(BIC_N) == -1) %>% dplyr::select(BIC_N) %>% colMeans()
Results_final %>% filter(sign(ICL)   == -1) %>% dplyr::select(ICL)   %>% colMeans()
