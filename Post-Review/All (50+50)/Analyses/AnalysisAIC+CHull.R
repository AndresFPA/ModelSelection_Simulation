# 14-11-2024
# Analysis AIC + CHull
# Load libraries
library(tidyr)

# Load design matrix
# Set wd
# setwd("C:/Users/perezalo/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/All (50+50)/Results")
setwd("C:/Users/User/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/All (50+50)/Results")
load("design.Rdata")

# # Define relevant objects
ncond <- nrow(design) # How many conditions?
K     <- 100 # How many replications?

# Create a matrix to store the maximum scree ratio
results_matrix <- matrix(data = NA, 
                         nrow = (ncond * K),
                         ncol = 5)

colnames(results_matrix)   <- c("Condition", "Replication", "AIC", "Chull", "True")
results_matrix             <- as.data.frame(results_matrix) 
results_matrix$AIC         <- I(vector(mode = "list", length = 2))
results_matrix$Chull       <- I(vector(mode = "list", length = 2))
results_matrix$Chull_first <- NA
results_matrix$AIC_first   <- NA
results_matrix$Condition   <- rep(x = 1:144, each = 100)
results_matrix$Replication <- rep(x = 1:100, times = 144)
results_matrix$True        <- rep(x = design$nclus, each = 100) 

# Load results and input into the matrix
# Prepare new load function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# setwd("C:/Users/perezalo/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/All (50+50)/Fit")
setwd("C:/Users/User/OneDrive - Tilburg University/1. Papers/Paper 2/Methodology (Journal)/Review/Simulation results/ModelSelection_Simulation/Post-Review/All (50+50)/Fit")

# Get best two results for AIC and Chull
ik <- 0
for(i in 1:ncond){
  for(k in 1:100){
    # browser()
    ik <- ik + 1
    file <- list.files(pattern = paste0("^FitRow", i, "Rep", k, "-"), all.files = T)
    results_row <- loadRData(file)
    AIC_res     <- list(c(which.min(results_row$AIC),
                          which(results_row$AIC %in% sort(results_row$AIC, decreasing = F)[2]))) # Best two results
    Chull_res   <- list(c(which.max(results_row$Chull),
                          which(results_row$Chull %in% sort(results_row$Chull, decreasing = T)[2]))) # Best two results
    results_matrix[ik, "AIC"][[1]]   <- AIC_res
    results_matrix[ik, "Chull"][[1]] <- Chull_res
  }
}

# Quick evaluation, do model selection considering the results that repeats the most for both measures
results_matrix$AIC_Chull <- NA
ik <- 0
for(i in 1:ncond){
  for(k in 1:100){
    ik <- ik + 1
    tmp_row  <- c(results_matrix$AIC[ik][[1]], results_matrix$Chull[ik][[1]]) # Extract and put into a vector the results of both measures
    freq_row <- table(tmp_row) # frequency table
    test     <- sum(max(freq_row) == freq_row) # Which is the winner model? 
    if(test > 1){ # Is there no clear winner?
      results_matrix$AIC_Chull[ik] <- results_matrix$Chull[ik][[1]][1] # Use the first choice from CHull (more consistent)
    } else {
      results_matrix$AIC_Chull[ik] <- as.numeric(names(which.max(freq_row))) # If there is a clear winner, select it
    }
    
    # Just in case also save result CHull and AIC
    results_matrix[ik, "Chull_first"] <- results_matrix$Chull[ik][[1]][1]
    results_matrix[ik, "AIC_first"] <- results_matrix$AIC[ik][[1]][1]
  }
}

# FINALLY! Check if we are under- over- or correctly selecting the clusters
results_matrix$Performance <- NA
results_matrix$Performance <- ifelse(test = results_matrix$True == results_matrix$AIC_Chull, yes = 0, 
                                     no = ifelse(test = results_matrix$True < results_matrix$AIC_Chull, yes = 1, 
                                                 no = -1))

# CHull alone performance
results_matrix$Chull_per <- NA
results_matrix$Chull_per <- ifelse(test = results_matrix$True == results_matrix$Chull_first, yes = 0, 
                                   no = ifelse(test = results_matrix$True < results_matrix$Chull_first, yes = 1, 
                                               no = -1))

# AIC alone performance
results_matrix$AIC_per <- NA
results_matrix$AIC_per <- ifelse(test = results_matrix$True == results_matrix$AIC_first, yes = 0, 
                                 no = ifelse(test = results_matrix$True < results_matrix$AIC_first, yes = 1, 
                                             no = -1))

# Check individual percentages
sum(results_matrix$AIC_per == 0 & results_matrix$Chull_per == 0)
sum(results_matrix$AIC_per != 0 & results_matrix$Chull_per == 0)
sum(results_matrix$AIC_per == 0 & results_matrix$Chull_per != 0)

# Which models were selected by CHull?
table(results_matrix$Chull_first)
