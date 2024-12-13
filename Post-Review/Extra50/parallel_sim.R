# Set the working directory
setwd("~/researchdrive/TSB0031 mmgsem (Projectfolder)/Extra50")

# Source the relevant functions
source("MMG-SEM.R")
source("E_Step.R")
source("ModelSelection.R")
source("DataGeneration.R")
source("evaluation.R")
source("evaluationARI.R")

# Simulation Design
# Which factors are going to be tested? For now:
nclus   <- c(2, 4)         # Number of clusters
ngroups <- c(24, 48)       # Number of groups
coeff   <- c(0.3, 0.4)     # Initial regression parameters
N_g     <- c(50, 100, 200) # Sample size per groups
balance <- c("bal", "unb") # Cluster size
sd      <- c(0, 0.05, 0.1) # Differences within a cluster (in betas)

model <- '
    # factor loadings
    F1 =~ x1 + x2 + x3 + x4 + x5
    F2 =~ z1 + z2 + z3 + z4 + z5
    F3 =~ m1 + m2 + m3 + m4 + m5
    F4 =~ y1 + y2 + y3 + y4 + y5
    
    # Regression parameters
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'
S1 <- '
    # factor loadings
    F1 =~ x1 + x2 + x3 + x4 + x5
    F2 =~ z1 + z2 + z3 + z4 + z5
    F3 =~ m1 + m2 + m3 + m4 + m5
    F4 =~ y1 + y2 + y3 + y4 + y5
'

S2 <- '
    # Regression parameters
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'

# Get design matrix
design <- expand.grid(nclus, ngroups, coeff, N_g, balance, sd, model) 
colnames(design) <- c("nclus", "ngroups", "coeff", "N_g", "balance", "sd", "model")
rownames(design) <- NULL
rm(balance, coeff, N_g, nclus, sd, ngroups) 

# Functions for the simulation
# First, to avoid stopping due to errors, create a function with data generation and MMGSEM
# Errors come from non positive definite cov matrices. This code allows the re-sample
genDat_analysis <- function(seed, RowDesign, k, NonInv){
  # Set seed per design condition (row) and replication (K)
  set.seed(seed)
  # browser()
  # Generate data
  #SimData <- do.call(what = DataGeneration, args = Design[RowDesign, ])$SimData
  SimData <- DataGeneration(model     = design[RowDesign, "model"], 
                            nclus     = design[RowDesign, "nclus"], 
                            ngroups   = design[RowDesign, "ngroups"], 
                            reg_coeff = design[RowDesign, "coeff"], 
                            N_g       = design[RowDesign, "N_g"], 
                            balance   = design[RowDesign, "balance"], 
                            sd        = design[RowDesign, "sd"])
  
  # Save data
  save(SimData, file = paste("Data/Data", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  save(seed, file = paste("Seed/Seed", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  # Run model selection from 1 to 6 clusters
  # 1. BOTH RES AND LOAD NON-INV ARE INCLUDED
  results <- ModelSelection(dat = SimData$SimData, 
                            S1 = S1, 
                            S2 = S2,
                            group = "group",
                            clusters = c(1, 6),
                            nstarts = 25,
                            seed = seed, 
                            group.equal = "loadings",
                            group.partial = NonInv)
}

# Create SAFE function (in case of errors)
genDat_analysis <- purrr::safely(.f = genDat_analysis, otherwise = NULL)

# Main simulation function
do_sim <- function(RowDesign){
  # Create the original clustering matrix for comparison below
  original <- create_original(balance = design[RowDesign, "balance"], 
                              ngroups = design[RowDesign, "ngroups"], 
                              nclus   = design[RowDesign, "nclus"])
  
  # Create matrix to store results
  # 12 columns for: BIC_G, BIC_N, AIC, AIC3, Chull, ICL 
  # There are 2 columns for each model selection measure (for factors and observed LL)
  ResultsRow <- matrix(data = NA, nrow = (K), ncol = 13)
  
  # Create second matrix for the ARI
  ResultsRowARI <- matrix(data = NA, nrow = (K), ncol = 2)
  
  # Define non-invariances
  NonInv <- c("F1 =~ x2", "F1 =~ x3",
              "F2 =~ z2", "F2 =~ z3",
              "F3 =~ m2", "F3 =~ m3",
              "F4 =~ y2", "F4 =~ y3")
  
  for(k in 1:K){
    k <- (k + 50)
    print(paste("Replication", k, "out of", K))
    
    # Code to re-sample in case the covariance matrix is non positive definite
    attempts <- 10
    for(j in 1:attempts){
      # Seed will change if there is an error
      # Seed includes a +50 to reflect that these are the extra 50 replications required by reviewers
      ctimes <- system.time(
        results <- genDat_analysis(seed = (RowDesign * k * j), RowDesign = RowDesign, k = k, NonInv = NonInv)
      ) 
      test <- results$result
      if(!is.null(test)){
        # If there was no error, break the loop and continue
        break
      }
    }
    
    # Save computation times
    save(ctimes, file = paste("Times/Time", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))

    # Save results if necessary
    results <- test
    save(results, file = paste("Fit/Fit", "Row", RowDesign, "Rep", k, "-", j, ".Rdata" , sep = ""))
    
    # Evaluate the results
    Evaluated            <- evaluation(res = results, clus = design[RowDesign, "nclus"])
    
    EvaluatedARI         <- evaluationARI(z_gks    = results$Models[[design[RowDesign, "nclus"]]]$posteriors,
                                          original = original,
                                          nclus    = design[RowDesign, "nclus"])
    # Store the results
    colnames(ResultsRow) <- colnames(Evaluated) 
    ResultsRow[(k - 50), ]      <- unlist(Evaluated) 
    
    colnames(ResultsRowARI) <- c("ARI", "CC")
    ResultsRowARI[(k - 50), ] <- unlist(EvaluatedARI)
  }
  
  # Save the results for each row
  save(ResultsRow, file = paste("Results/Result", "Row", RowDesign,".Rdata" , sep =""))
  
  save(ResultsRowARI, file = paste("Results/Result", "Row", "ARI", RowDesign,".Rdata" , sep =""))
  
  # Return the final results
  return(results)
}

# Set number of replications per condition
K <- 50

# Prepare object to save
Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K, ncol = 13))
Results_final$Replication <- rep(x = 1:K, times = nrow(design))
Results_final$Condition <- rep(x = 1:nrow(design), each = K)

save(Results_final, file = "FinalResults.Rdata")
save(design, file = "design.Rdata")

# ###################################################################### #
# ######################## START PARALLELIZATION ####################### #
# ###################################################################### #

library(foreach)
library(parallel)
library(doParallel)

# Get object names
obj <- objects()

# Setup parallel cluster
cl <- makeCluster(2)
registerDoParallel(cl)

# Divide rows per cluster
# rows_divided <- split(1:2, 1:2)

# Export the necessary function and variables to the cluster
clusterEvalQ(cl, {
  library(lavaan)
  library(MASS)
  library(combinat)
})
clusterEvalQ(cl, setwd("~/researchdrive/TSB0031 mmgsem (Projectfolder)/Extra50"))
parallel::clusterExport(cl, varlist = obj)

# Parallel execution over RowDesign
# results <- parLapply(cl = cl, X = 1:2, fun = do_sim)
results <- foreach(RowDesign = 1:2) %dopar% {
  do_sim(RowDesign)
}

stopCluster(cl)





