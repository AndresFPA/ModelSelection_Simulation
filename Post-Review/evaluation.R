#' Evaluation
#'
#' Evaluation function to assess the performance of the model according to the simulation results. 
#'
#' INPUT: Arguments required by the function
#' @param beta 
#' @param z_gks
#' @param original
#' @param global_max
#' @param runs_LL
#' @param nclus
#' @param reg_coef
#' @param reg_diff

library("combinat")
library("FARI")

evaluation <- function(res, clus){
  # browser()
  # Create the matrix to be returned
  res           <- res$Overview
  rmv.idx       <- which(colnames(res) %in% c("Clusters", "LL", "nrpar", "LL_fac", "nrpar_fac"))
  evaluated     <- res[1, -rmv.idx, drop = F]
  evaluated[, ] <- NA
  evaluated     <- as.data.frame(evaluated)
  
  res <- as.data.frame(res)
  
  # Fill in the matrix: Does the selected model correspond with the true model?
  # TRUE  = We matched the correct clustering
  # Under = We under selected
  # Over  = We over selected
  
  evaluated$R2entropy <- res$R2entropy[[clus]]
  
  evaluated$Chull <- ifelse(all(is.na(res$Chull)), NA, ifelse(which.max(res$Chull) == clus, 0, ifelse(which.max(res$Chull) < clus, -1, 1)))
  evaluated$BIC_G <- ifelse(which.min(res$BIC_G) == clus, 0, ifelse(which.min(res$BIC_G) < clus, -1, 1))
  evaluated$BIC_N <- ifelse(which.min(res$BIC_N) == clus, 0, ifelse(which.min(res$BIC_N) < clus, -1, 1))
  evaluated$AIC   <- ifelse(which.min(res$AIC)   == clus, 0, ifelse(which.min(res$AIC)   < clus, -1, 1))
  evaluated$AIC3  <- ifelse(which.min(res$AIC3)  == clus, 0, ifelse(which.min(res$AIC3)  < clus, -1, 1))
  evaluated$ICL   <- ifelse(which.min(res$ICL)   == clus, 0, ifelse(which.min(res$ICL)   < clus, -1, 1))
  
  evaluated$Chull_fac <- ifelse(all(is.na(res$Chull_fac)), NA, ifelse(which.max(res$Chull_fac) == clus, 0, ifelse(which.min(res$Chull_fac) < clus, -1, 1)))
  evaluated$BIC_G_fac <- ifelse(which.min(res$BIC_G_fac) == clus, 0, ifelse(which.min(res$BIC_G_fac) < clus, -1, 1))
  evaluated$BIC_N_fac <- ifelse(which.min(res$BIC_N_fac) == clus, 0, ifelse(which.min(res$BIC_N_fac) < clus, -1, 1))
  evaluated$AIC_fac   <- ifelse(which.min(res$AIC_fac)   == clus, 0, ifelse(which.min(res$AIC_fac)   < clus, -1, 1))
  evaluated$AIC3_fac  <- ifelse(which.min(res$AIC3_fac)  == clus, 0, ifelse(which.min(res$AIC3_fac)  < clus, -1, 1))
  evaluated$ICL_fac   <- ifelse(which.min(res$ICL_fac)   == clus, 0, ifelse(which.min(res$ICL_fac)   < clus, -1, 1))
  
  return(evaluated)
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










