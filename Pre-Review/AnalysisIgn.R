library(lavaan)
library(qwraps2)
library(fpp3)
library(dplyr)
library(xtable)
library(ggpubr)
library(ggplot2)
library(ggthemes)
# library(Cairo)
# CairoWin()

# Set wd
setwd("~/GitHub/ModelSelection_Simulation/Results")

# Load empty final results matrix
load("FinalResults.Rdata")
colnames(Results_final)[1:13] <- c("entropyR2", "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
                                   "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac")
load("design.Rdata")

# Merge datasets
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")
col_order <- c("Condition", "Replication", "nclus", "ngroups", "coeff", "N_g", "balance", "sd", 
               "entropyR2", "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL", 
               "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac")
Results_final <- Results_final[, col_order]
rm(col_order)

# Fill in the matrix with all results
ncond <- unique(Results_final$Condition) # How many conditions?
K <- length(unique(Results_final$Replication)) # How many replications?

for (i in ncond) {
  test <- NA
  test <- try(load(paste0("ResultIgnRow", i, ".Rdata")))
  if(!c(class(test) == "try-error")){
    Results_final[(K*(i-1)+1):(i*K), 9:21] <- ResultsRow
  }
}

# remove uncomplete entries
Results_final <- Results_final[!is.na(Results_final$BIC_G), ]

# Turn NAs from Chull into FALSE input (Chull was not able to select any model)
# apply(X = apply(X = Results_final, MARGIN = 2, FUN = is.na), MARGIN = 2, FUN = sum)
# Results_final$`Chull Scree` <- ifelse(test = is.na(Results_final$`Chull Scree`), yes = FALSE, no = Results_final$`Chull Scree`)

# Transform to factor
# changed <- Results_final %>% dplyr::select(Chull:ICL_fac) %>% as.matrix() %>% as.data.frame() %>% 
#   mutate(
#     Chull     = recode(Chull, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     BIC_G     = recode(BIC_G, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     BIC_N     = recode(BIC_N, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     AIC       = recode(AIC, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     AIC3      = recode(AIC3, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     ICL       = recode(ICL, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     Chull_fac = recode(Chull_fac, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     BIC_G_fac = recode(BIC_G_fac, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     BIC_N_fac = recode(BIC_N_fac, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     AIC_fac   = recode(AIC_fac, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     AIC3_fac  = recode(AIC3_fac, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1"),
#     ICL_fac   = recode(ICL_fac, "1" = "0", "TRUE" = "0", "over" = "1", "under" = "-1")
#   )

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

count_results(data = Results_final, by = c("sd", "N_g"), type = "relative")

mean(Results_final$entropyR2)
View(Results_final %>% group_by(nclus, N_g, ngroups, coeff, balance, sd) %>% summarise(across(entropyR2, mean)))

Results_final %>% group_by(nclus) %>% summarise(across(entropyR2, mean))
Results_final %>% group_by(N_g) %>% summarise(across(entropyR2, mean))
Results_final %>% group_by(ngroups) %>% summarise(across(entropyR2, mean))
Results_final %>% group_by(coeff) %>% summarise(across(entropyR2, mean))
Results_final %>% group_by(balance) %>% summarise(across(entropyR2, mean))
Results_final %>% group_by(sd) %>% summarise(across(entropyR2, mean))