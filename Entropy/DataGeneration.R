library(MASS)
library(matrixcalc)

#' DataGeneration
#'
#' Generates simulated data according to the manipulated conditions. 
#'
#' INPUT: Arguments required by the function
#' @param model: full model in lavaan syntax (string).
#' @param nclus: number of clusters of the data. 
#' @param ngroups: number of groups for the generated data.
#' @param N_g: number of observations per group.
#' @param reg_coeff: regression coefficients of cluster 1. Must a numeric vector of length 3.
#' @param balance: can take values "balanced" and "unbalanced". Determines if the clusters will have the same
#'        number of groups.
#' @param reliability: reliability per item (manipulates ratio between loading and unique variance).
#' @param NonInvSize: size of the loading non-invariance (for the simulation: 0.2 and 0.4).
#' @param NonInvItems: number of items affected by the non-invariance (for now, only 2).
#' @param NonInvG: proportion of groups affected by the non-invariance (e.g., 0.25, 0.5).
#' @param randomVarX: define the variance of the exogenous variable. If TRUE, then random values will be defined
#'        for all groups. If FALSE, the variance will be 1 for all groups.
#'
#' OUTPUT
#' @return SimData: generated data
#' @return NonInvIdx: Index of the groups that presented the non-invariant loadings 
#' @return psi_g: generated psi matrix
#' @return OrTheta: generated theta matrix
#' @return cov_eta: generated cov_eta matrix (phi in the paper)

# Note:
# In the code, we use exog and endog to name the latent variables, while in the paper we use F1, F2, ..., etc.
# The corresponding matches are:
# F1 = exog2
# F2 = exog1
# F3 = endog1
# F4 = endog2

# Note 2:
# To ease the reading of the code, the names of the regression parameters can be seen below:
#   # exog1 -> endog2 = B1
#   # exog1 -> endog1 = B2
#   # endog1 -> endog2 = B3
#   # endog1 -> endog2 = B4


PopR2Entropy <- function(model, step1model, nclus, ngroups, N_g,
                         reg_coeff, balance, sd,
                         reliability = "high", NonInvSize = 0.4, # The factors below are fixed in this simulation
                         NonInvItems = 2, NonInvG = 0.5, NonInvType = "random",
                         ResRange = 0.2, randomVarX = T){
  
  # Get number of variables
  par_table <- lavaanify(model) # Create parameter table
  lat_var <- lavNames(par_table, "lv")
  obs_var <- lavNames(par_table, "ov")
  m <- length(lat_var) # How many latent variables?
  p <- length(obs_var) # How many observed variables?
  
  # Identify type of latent variable
  endog1 <- lat_var[(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog2 <- lat_var[!c(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]
  
  # Reorder latent variables
  lat_var <- c(exog, endog)
  
  # Get a cluster label for each group (GperK)
  # E.g., If balanced: G = 6, K = 2. GperK = 111222
  if (balance == "bal"){
    GperK <- rep(x = 1:nclus, each = (ngroups/nclus))
  } else if (balance == "unb"){
    largest <- ngroups*.75; smaller <- ngroups - largest # largest and smaller clusters
    GperK <- c(rep(x = 1, times = largest), rep(x = 2:nclus, each = smaller/(nclus - 1)))
  }
  
  # Give names to the regression parameters (only for better understanding):
  # Get regression parameters for each CORRECT cluster
  # Cluster 1
  B1 <- rnorm(n = sum(GperK == 1), mean = 0, sd = sd)
  B2 <- rnorm(n = sum(GperK == 1), mean = reg_coeff, sd = sd)
  B3 <- rnorm(n = sum(GperK == 1), mean = reg_coeff, sd = sd)
  B4 <- rnorm(n = sum(GperK == 1), mean = reg_coeff, sd = sd)
  
  # Cluster 2
  B1 <- c(B1, rnorm(n = sum(GperK == 2), mean = reg_coeff, sd = sd))
  B2 <- c(B2, rnorm(n = sum(GperK == 2), mean = 0, sd = sd))
  B3 <- c(B3, rnorm(n = sum(GperK == 2), mean = reg_coeff, sd = sd))
  B4 <- c(B4, rnorm(n = sum(GperK == 2), mean = reg_coeff, sd = sd))
  
  if (nclus == 4){
    # Cluster 3
    B1 <- c(B1, rnorm(n = sum(GperK == 3), mean = reg_coeff, sd = sd))
    B2 <- c(B2, rnorm(n = sum(GperK == 3), mean = reg_coeff, sd = sd))
    B3 <- c(B3, rnorm(n = sum(GperK == 3), mean = 0, sd = sd))
    B4 <- c(B4, rnorm(n = sum(GperK == 3), mean = reg_coeff, sd = sd))
    
    # Cluster 4
    B1 <- c(B1, rnorm(n = sum(GperK == 4), mean = reg_coeff, sd = sd))
    B2 <- c(B2, rnorm(n = sum(GperK == 4), mean = reg_coeff, sd = sd))
    B3 <- c(B3, rnorm(n = sum(GperK == 4), mean = reg_coeff, sd = sd))
    B4 <- c(B4, rnorm(n = sum(GperK == 4), mean = 0, sd = sd))
  }
  
  # beta is group-specific (small difference within-cluster).
  beta <- array(data = 0, dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))
  colnames(beta) <- rownames(beta) <- lat_var
  
  # Fill in correct regression parameters
  for(g in 1:ngroups){
    # beta
    beta[endog2, exog[1], g] <- B1[g]
    beta[endog1, exog[1], g] <- B2[g]
    beta[endog1, exog[2], g] <- B3[g]
    beta[endog2, endog1, g] <- B4[g]
  }
  
  # Get the POPULATION regressions
  B1_pop <- numeric(nclus)
  B2_pop <- numeric(nclus)
  B3_pop <- numeric(nclus)
  B4_pop <- numeric(nclus)
  
  # Cluster 1
  B1_pop[1] <- 0
  B2_pop[1] <- reg_coeff
  B3_pop[1] <- reg_coeff
  B4_pop[1] <- reg_coeff
  
  # Cluster 2
  B1_pop[2] <- reg_coeff
  B2_pop[2] <- 0
  B3_pop[2] <- reg_coeff
  B4_pop[2] <- reg_coeff
  
  if (nclus == 4){
    # Cluster 3
    B1_pop[3] <- reg_coeff
    B2_pop[3] <- reg_coeff
    B3_pop[3] <- 0
    B4_pop[3] <- reg_coeff
    
    # Cluster 4
    B1_pop[4] <- reg_coeff
    B2_pop[4] <- reg_coeff
    B3_pop[4] <- reg_coeff
    B4_pop[4] <- 0
  }
  
  # Generate the POPULATION structural parameters
  # BETA - Regression parameters
  beta_pop <- array(data = 0, dim = c(m, m, nclus), dimnames = list(lat_var, lat_var))
  colnames(beta_pop) <- rownames(beta_pop) <- lat_var
  
  # Initialize
  for(k in 1:nclus){
    # beta
    beta_pop[endog2, exog[1], k] <- B1_pop[k]
    beta_pop[endog1, exog[1], k] <- B2_pop[k]
    beta_pop[endog1, exog[2], k] <- B3_pop[k]
    beta_pop[endog2, endog1, k] <- B4_pop[k]
  }
  
  # psi is group- and cluster-specific. ngroups matrices are needed.
  # Generate enough psi matrices
  psi_g <- array(data = diag(m), dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))
  
  # Generate group-specific values for psi sampling from uniform distributions
  if (randomVarX == T){
    exog_var1 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_var2 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_cov <- runif(n = ngroups, min = -0.30, max = 0.30)
    endo_var1 <- runif(n = ngroups, min = 0.75, max = 1.25) # Total endogenous variance
    endo_var2 <- runif(n = ngroups, min = 0.75, max = 1.25) # Total endogenous variance
  } else { # DEPRECATED
    exog_var1 <- rep(1, times = ngroups) 
    exog_var2 <- rep(1, times = ngroups) 
    exog_cov <- rep(1, times = ngroups) 
  }
  
  # Insert the corresponding group- and cluster-specific parts of psi
  for(g in 1:ngroups){
    # Insert group-specific parts
    psi_g[exog[1], exog[1], g] <- exog_var1[g]
    psi_g[exog[2], exog[2], g] <- exog_var2[g]
    psi_g[exog[1], exog[2], g] <- exog_cov[g]
    psi_g[exog[2], exog[1], g] <- exog_cov[g]

    # Insert the group-and-cluster-specific parts
    # For the endogenous variances, start from the total var (endog_var) and subtract the explained variance by the regression
    psi_g[endog1, endog1, g] <- endo_var1[g] - ((B2[g]^2 * exog_var1[g]) +
                                                  (B3[g]^2 * exog_var2[g]) +
                                                  (2 * B2[g] * B3[g] * exog_cov[g]))

    psi_g[endog2, endog2, g] <- endo_var2[g] - ((B1[g]^2 * exog_var1[g]) +
                                                  (B4[g]^2 * endo_var1[g]) +
                                                  (2 * B1[g] * B4[g] * ((B2[g] * exog_var1[g]) + (B3[g] * exog_cov[g]))))
  }
  
  # browser()
  
  # The following double-loop is only done for the entropy computation
  # Create the group-cluster parameters of psi_gk
  psi_gk <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
  for(k in 1:nclus){
    for(g in 1:ngroups){
      # Insert empty matrix in the corresponding part
      psi_gk[[g, k]] <- matrix(data = diag(m), nrow = m, ncol = m, dimnames = list(lat_var, lat_var))
      
      # Insert group-specific parts
      psi_gk[[g, k]][exog[1], exog[1]] <- exog_var1[g]
      psi_gk[[g, k]][exog[2], exog[2]] <- exog_var2[g]
      psi_gk[[g, k]][exog[1], exog[2]] <- exog_cov[g]
      psi_gk[[g, k]][exog[2], exog[1]] <- exog_cov[g]
      
      # Insert the group-and-cluster-specific parts
      # For the endogenous variances, start from the total var (endog_var) and subtract the explained variance by the regression
      psi_gk[[g, k]][endog1, endog1] <- endo_var1[g] - ((B2_pop[k]^2 * exog_var1[g]) +
                                                       (B3_pop[k]^2 * exog_var2[g]) +
                                                       (2 * B2_pop[k] * B3_pop[k] * exog_cov[g]))
      
      psi_gk[[g, k]][endog2, endog2] <- endo_var2[g] - ((B1_pop[k]^2 * exog_var1[g]) +
                                                       (B4_pop[k]^2 * endo_var1[g]) +
                                                       (2 * B1_pop[k] * B4_pop[k] * ((B2_pop[k] * exog_var1[g]) + (B3_pop[k] * exog_cov[g]))))
    }
  }
  
  # Create the covariance matrix of the factors (phi in the paper)
  I <- diag(m) # Identity matrix
  cov_eta <- array(data = 0, dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))

  for(g in 1:ngroups){
    cov_eta[, , g] <- solve(I - beta[, , g]) %*% psi_g[, , g] %*% solve(t(I - beta[, , g]))
  }
  
  # Generate the measurement model parameters
  # Lambda (depends on reliability lvl)
  if (reliability == "low"){load <- .4} else if (reliability == "high"){load <- .6}
  loadings <- sqrt(load)
  
  # Create Invariant Lambda
  Lambda <- matrix(data = rep(x = c(1, rep(loadings, 4), rep(0, p)), times = m)[1:(p*m)], nrow = p, ncol = m)
  
  # Create non-invariant Lambda (depending on the type of non-invariance)
  if(NonInvType == "fixed"){
    LambdaNonInv <- matrix(data = rep(x = c(1, 
                                            (loadings - NonInvSize), 
                                            (loadings + NonInvSize), 
                                            rep(loadings, (4 - NonInvItems)), 
                                            rep(0, p)), times = m)[1:(p*m)], nrow = p, ncol = m)
  } else if (NonInvType == "random"){
    # Sample random non-invariances
    # NonInvariantLoadings <- sample(x = c(runif(100, min = (loadings - NonInvSize) - .1, max = (loadings - NonInvSize) + .1), 
    #                                      runif(100, min = (loadings + NonInvSize) - .1, max = (loadings + NonInvSize) + .1)),
    #                                size = NonInvItems*4)
    
    NonInvariantLoadings <- sample(x = c(runif(100, min = (loadings - NonInvSize), max = (loadings - NonInvSize)), 
                                         runif(100, min = (loadings + NonInvSize), max = (loadings + NonInvSize))),
                                   size = NonInvItems*4)
    
    # Create a non-invariant lambda matrix
    LambdaNonInv <- matrix(data = c(1, NonInvariantLoadings[1:NonInvItems], rep(loadings, (4 - NonInvItems)), rep(0, p),
                                    1, NonInvariantLoadings[(NonInvItems + 1):(NonInvItems*2)], rep(loadings, (4 - NonInvItems)), rep(0, p),
                                    1, NonInvariantLoadings[((NonInvItems*2) + 1):(NonInvItems*3)], rep(loadings, (4 - NonInvItems)), rep(0, p),
                                    1, NonInvariantLoadings[((NonInvItems*3) + 1):(NonInvItems*4)], rep(loadings, (4 - NonInvItems))),
                           nrow = p, ncol = m)
  }
  
  # Theta
  Theta <- array(data = 0, dim = c(p, p, ngroups))
  
  # Sample group-specific (diagonal of) theta values from an uniform distribution
  for (g in 1:ngroups){
    Theta[, , g] <- diag(runif(n = p, min = ((1 - load) - (ResRange/2)), max = ((1 - load) + (ResRange/2))))
  }
  
  # Generate sample covariance matrix (sigma) per groups
  Sigma <- array(data = 0, dim = c(p, p, ngroups))
  
  # Include non-invariance in the loadings according to NonInvG
  # Which groups are non-invariant? (random sampling)
  NonInvIdx <- sample(x = 1:ngroups, size = NonInvG*ngroups, replace = F)
  for(g in 1:ngroups){
    # Non-invariance - How many groups?
    if (g %in% NonInvIdx){
      Sigma[, , g] <- LambdaNonInv %*% cov_eta[, , g] %*% t(LambdaNonInv) + Theta[, , g]
    } else if (!c(g %in% NonInvIdx)){
      Sigma[, , g] <- Lambda %*% cov_eta[, , g] %*% t(Lambda) + Theta[, , g]
    }
  }

  
  # Data Generation final
  # For now, mu would be 0 as we are only interested in centered variables
  # SimData <- c()
  # for(g in 1:ngroups){
  #   tmp <- mvrnorm(n = N_g, mu = rep(0, p), Sigma = Sigma[, , g], empirical = T)
  #   SimData <- rbind(SimData, tmp)
  # }
  # 
  # # Add the final labels
  # group <- rep(x = c(1:ngroups), each = N_g) # Group variable
  # SimData <- cbind(SimData, group)
  # colnames(SimData) <- c(obs_var, "group")
  
  ####################################### DATA IS GENERATED #########################################
  ###################################################################################################
  ######################################### ESTIMATE MODEL ##########################################
  # browser()
  # "STEP 1"
  # Run CFA on the generated data with true values as starting values
  # fake_partable <- lavaan::parTable(cfa(
  #   model = step1model, data = SimData, group = "group",
  #   estimator = "ML", group.equal = "loadings",
  #   se = NULL, test = "none", baseline = FALSE, h1 = FALSE,
  #   implied = FALSE, loglik = FALSE,
  #   meanstructure = FALSE, do.fit = F
  # ))
  # browser()
  # S1output <- lavaan::cfa(
  #   model = step1model, data = SimData, group = "group",
  #   estimator = "ML", group.equal = "loadings",
  #   se = NULL, test = "none", baseline = FALSE, h1 = FALSE,
  #   implied = FALSE, loglik = FALSE,
  #   meanstructure = FALSE, do.fit = F
  # )
  
  # Define some important objects
  # How many groups?
  # ngroups <- lavaan::lavInspect(S1output, "ngroups")
  # N_gs    <- lavaan::lavInspect(S1output, "nobs") # nobs per group
  # 
  # # all estimated model matrices, per group
  # EST       <- lavaan::lavInspect(S1output, "est", add.class = FALSE, add.labels = TRUE)
  # theta_gs  <- lapply(EST, "[[", "theta")
  # lambda_gs <- lapply(EST, "[[", "lambda")
  # cov_eta_es   <- lapply(EST, "[[", "psi") # cov_eta name refers to Variance of eta (eta being the latent variables)
  
  ## "Run" the "second step"
  # Create random partition (NOT SURE IF THIS IS OKAY)
  # cl <- 0
  # while (cl < 1) { # "while loop" to make sure all clusters get at least one group
  #   z_gks <- t(replicate(ngroups, sample(x = c(rep(0, (nclus - 1)), 1))))
  #   cl <- min(colSums(z_gks))
  # }
  # 
  # pi_ks <- colMeans(z_gks)
  
  # Create correct partition
  z_gks <- create_original(balance = balance, 
                           ngroups = ngroups, 
                           nclus = nclus)
  
  pi_ks <- colMeans(z_gks)
  
  # Estimate loglikelihood using true values as starting values
  # Initialize matrices to store loglikelihoods
  loglik_gks  <- matrix(data = 0, nrow = ngroups, ncol = nclus)
  loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)
  
  # Initialize Sigma
  Sigma <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
  I <- diag(length(lat_var)) # Identity matrix based on number of latent variables. Used later
  
  # Estimate LL
  for(k in 1:nclus){
    for(g in 1:ngroups){
      # Estimate Sigma (factor covariance matrix of step 2)
      Sigma[[g, k]] <- solve(I - beta_pop[, , k]) %*% psi_gk[[g, k]] %*% solve(t(I - beta_pop[, , k]))
      Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]])) # Force to be symmetric
      
      # Estimate the loglikelihood
      loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
        sample.mean = rep(0, nrow(cov_eta[, , g])),
        sample.nobs = N_g, # Use original sample size to get the correct loglikelihood
        # sample.nobs = N_gks[g, k],
        sample.cov  = cov_eta[, , g], #cov_eta_es[[g]], # Factor covariance matrix from step 1
        Mu          = rep(0, nrow(cov_eta[, , g])),
        Sigma       = Sigma[[g, k]] # Factor covariance matrix from step 2
      )
      
      loglik_gks[g, k] <- loglik_gk
      loglik_gksw[g, k] <- log(pi_ks[k]) + loglik_gk # weighted loglik
    }
  }
  
  # Get total loglikelihood
  # First, deal with arithmetic underflow by subtracting the maximum value per group
  max_gs <- apply(loglik_gksw, 1, max) # Get max value per row
  minus_max <- sweep(x = loglik_gksw, MARGIN = 1, STATS = max_gs, FUN = "-") # Subtract the max per row
  exp_loglik <- exp(minus_max) # Exp before summing for total loglikelihood
  loglik_gsw <- log(apply(exp_loglik, 1, sum)) # Sum exp_loglik per row and then take the log again
  LL <- sum((loglik_gsw + max_gs)) # Add the maximum again and then sum them all for total loglikelihood
  
  # Define E-step for the final z_gks computation
  EStep <- function(pi_ks, ngroup, nclus, loglik){
    
    max_g <-rep(0,ngroup)
    z_gks <- matrix(NA,nrow = ngroup,ncol = nclus)
    
    for(g in 1:ngroup){
      for(k in 1:nclus){
        z_gks[g,k] <- log(pi_ks[k])+loglik[g,k]
      }
      max_g[g] <- max(z_gks[g,]) # prevent arithmetic underflow 
      z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclus))
    }
    
    # divide by the rowwise sum of the above calculated part 
    z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
    # z_gks <- round(z_gks, digits = 16)
    # z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
    
    return(z_gks)
  }
  
  z_gks <- EStep(
    pi_ks = pi_ks, ngroup = ngroups,
    nclus = nclus, loglik = loglik_gks
  )
  
  # DO THE ENTROPY COMPUTATION
  # Entropy
  # Code from github user daob (Oberski, 2019): https://gist.github.com/daob/c2b6d83815ddd57cde3cebfdc2c267b3
  # p is the prior or posterior probabilities
  entropy <- function(p) {
    p <- p[p > sqrt(.Machine$double.eps)] # since Lim_{p->0} p log(p) = 0
    sum(-p * log(p))
  }
  # sum_entropy <- sum(apply(z_gks, 1, entropy)) # Total entropy
  
  # Entropy R2
  entropy.R2 <- function(pi, post) { 
    error_prior <- entropy(pi) # Class proportions
    error_post <- mean(apply(post, 1, entropy))
    R2_entropy <- (error_prior - error_post) / error_prior
    R2_entropy
  }
  
  R2_entropy <- entropy.R2(pi = pi_ks, post = z_gks)
  
  # Return data
  return(R2_entropy)
}


# Function to create the original cluster matrix
create_original <- function(balance, ngroups, nclus){
  if (balance == "unb"){
    unb <- c(rep(0, ngroups), rep(1, (ngroups*.25)/(nclus - 1)))
    original <- matrix(data = c(rep(1, ngroups*.75), rep(unb, nclus - 1)), nrow = ngroups, ncol = nclus)
  } else {
    data <- c(rep(x = c(rep(1, (ngroups/nclus)), rep(0, ngroups)), times = nclus))
    data <- data[-c((length(data)-ngroups+1):length(data))]
    original <- matrix(data = data, nrow = ngroups, ncol = nclus)
  }
  return(original)
}
