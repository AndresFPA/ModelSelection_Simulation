# 12-11-2024
# Analysis CHull - Simulation K = 1
# Load packages
library(tidyr)

# Define relevant objects
n_cond <- 18 # Number of conditions
K      <- 100 # Number of replications

# Create a matrix to store the maximum scree ratio
results_matrix <- matrix(data = NA, 
                         nrow = (n_cond * K),
                         ncol = 5)

colnames(results_matrix)   <- c("Condition", "Replication", "Scree", "Selected_K", "Scree_fac")
results_matrix             <- as.data.frame(results_matrix) 
results_matrix$Condition   <- rep(x = 1:18, each = 100)
results_matrix$Replication <- rep(x = 1:100, times = 18)

# Load results and input into the matrix
# Prepare new load function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#setwd("C:/Users/perezalo/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/K1/Fit")
setwd("C:/Users/User/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/K1/Fit")

ik <- 0
for(i in 1:n_cond){
  for(k in 1:K){
    ik <- ik + 1
    file <- list.files(pattern = paste0("^FitRow", i, "Rep", k, "-"), all.files = T)
    results_row <- loadRData(file)
    Chull_fac   <- results_row$Overview$Chull_fac
    Chull       <- results_row$Overview$Chull
    results_matrix[ik, "Scree"]      <- ifelse(test = all(is.na(Chull)), yes = NA, no = max(Chull, na.rm = T))
    results_matrix[ik, "Selected_K"] <- ifelse(test = all(is.na(Chull)), yes = NA, no = which.max(Chull))
    results_matrix[ik, "Scree_fac"]  <- ifelse(test = all(is.na(Chull_fac)), yes = NA, no = max(Chull_fac, na.rm = T))
  }
}

mean(x = results_matrix$Scree[is.finite(results_matrix$Scree)])
table(x = results_matrix$Selected_K[is.finite(results_matrix$Selected_K)])#/sum(table(x = results_matrix$Selected_K[is.finite(results_matrix$Selected_K)]))

mean(x = results_matrix$Scree_fac[is.finite(results_matrix$Scree_fac)])
