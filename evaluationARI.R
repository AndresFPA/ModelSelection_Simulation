evaluationARI <- function(z_gks, original, nclus){
  # Transform posterior classification probabilities to hard clustering for the sake of the evaluation
  # clusResult <- ifelse(test = z_gks > .5, yes = 1, no = 0)
  clusResult <- t(apply(z_gks, 1, function(x) as.numeric(x == max(x))))
  colnames(clusResult) <- paste("Cluster", seq_len(nclus))
  
  # Organize columns in the same way as the original order (NOT USED NOW)
  # i.e. first cluster has the first x groups, second cluster the second x groups, etc.
  # clusResult <- clusResult[, names(sort(apply(X = clusResult, MARGIN = 2, FUN = which.max)))]
  
  # Get "true" cluster labels
  true_or <- original
  for (i in 1:ncol(original)){
    original[original[, i] != 0, i] <- i
  }
  or_vec <- c(original)[c(original) != 0]
  or_vec <- factor(x = or_vec, levels = 1:nclus) # Just in case there is an empty cluster
  
  # Evaluation 1 - Cluster Recovery ----------------------------------------------------------------
  # Cluster Recovery 1. Misclassification Error Rate
  
  # Check all possible permutations of cluster results order
  perm <- permn(x = ncol(clusResult))
  n_perm <- length(perm)
  
  # Initialize necessary objects
  permutedClusters <- vector(mode = "list", length = n_perm)
  clusVec <- vector(mode = "list", length = n_perm)
  MisClassError <- numeric(n_perm)
  
  # Check every permutation individually
  for(i in 1:n_perm){
    permutedClusters[[i]] <- clusResult[, perm[[i]]]
    
    # Label the groups depending on the cluster 
    # (i.e., groups in cluster 1 will have a 1, groups in cluster 2 will have a 2, etc.)
    for (j in 1:ncol(clusResult)){
      permutedClusters[[i]][permutedClusters[[i]][, j] != 0, j] <- j
    }
    
    # Re order columns by order of groups (NOT USED NOW)
    # permutedClusters[[i]] <- permutedClusters[[i]][, names(sort(apply(X = permutedClusters[[i]], MARGIN = 2, FUN = which.max)))]
    
    # Get a vector of cluster labels
    clusVec[[i]] <- rowSums(permutedClusters[[i]])
    
    # Transform into a factor to add all the necessary levels (in case some cluster does not have any groups)
    clusVec[[i]] <- factor(x = clusVec[[i]], levels = 1:nclus)
    
    # Create a confusion matrix and calculate misclassification error rate
    conf_mat <- table(clusVec[[i]], or_vec)
    MisClassError[i] <- 1 - (sum(diag(conf_mat))/sum(conf_mat))
  }
  
  # Get the result from the best permutation
  MisClassErrorOut <- min(MisClassError)
  
  # Which permutation is the best one? (will be used later)
  bestPerm <- permutedClusters[[which.min(MisClassError)]]
  bestPermVec <- clusVec[[which.min(MisClassError)]]
  #browser()
  # Cluster Recovery 2. Adjusted Rand Index (ARI)
  ARI_res <- adjrandindex(part1 = bestPermVec, part2 = or_vec)
  
  # Cluster Recovery 3. Correct clustering?
  CorrectClus <- all(bestPermVec == or_vec)
  
  # Evaluation 3 - Proportion of correct results ----------------------------------------------------
  # Proportion of global maxima
  # Round the loglikelihoods to avoid differences of floating points
  # ProporGlobalMax <- mean(round(runs_LL, 2) == round(global_max, 2)) #/length(runs_LL)
  # 
  # # EXTRA: Average exogenous variable variance ------ 
  # exo_mean <- mean(unlist(lapply(X = psi_gks[, 1], FUN = '[[', 1)))
  # cov_mean <- mean(unlist(lapply(X = psi_gks[, 1], FUN = '[[', 2)))
  
  # Return results ---------------------------------------------------------------------------------
  return(list(ClusRecovery = list(ARI = ARI_res, CorrectClus = CorrectClus)))
}

# Function to create the original cluster matrix
create_original <- function(balance, ngroups, nclus){
  if (balance == "unbalanced"){
    unb <- c(rep(0, ngroups), rep(1, (ngroups*.25)/(nclus - 1)))
    original <- matrix(data = c(rep(1, ngroups*.75), rep(unb, nclus - 1)), nrow = ngroups, ncol = nclus)
  } else {
    data <- c(rep(x = c(rep(1, (ngroups/nclus)), rep(0, ngroups)), times = nclus))
    data <- data[-c((length(data)-ngroups+1):length(data))]
    original <- matrix(data = data, nrow = ngroups, ncol = nclus)
  }
  return(original)
}



# computation of adjusted rand index
adjrandindex <- function(part1,part2){
  part1 <- as.numeric(part1)
  part2 <- as.numeric(part2)
  
  IM1 <- diag(max(part1))
  IM2 <- diag(max(part2))
  A <- IM1[part1,]
  B <- IM2[part2,]
  
  T <- t(A)%*%B
  N <- sum(T)
  Tc <- apply(T,2,sum)
  Tr <- apply(T,1,sum)
  a <- (sum(T^2) - N)/2
  b <- (sum(Tr^2) - sum(T^2))/2
  c <- (sum(Tc^2) - sum(T^2))/2
  d <- (sum(T^2) + N^2 - sum(Tr^2) - sum(Tc^2))/2
  ARI <- (choose(N,2)*(a + d) - ((a+b)*(a+c)+(c+d)*(b+d)))/(choose(N,2)^2 - ((a+b)*(a+c)+(c+d)*(b+d)))
  
  return(ARI)
}

