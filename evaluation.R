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

evaluation <- function(res, true_clus){
  
  # Create the matrix to be returned
  evaluated <- matrix(data = NA, nrow = 1, ncol = ncol(res[, -c(1:3, 8:9)]))
  colnames(evaluated) <- colnames(res[, -c(1:3, 8:9)])
  evaluated <- as.data.frame(evaluated)
  
  res <- as.data.frame(res)
  
  # Fill in the matrix: Does the selected model correspond with the true model?
  evaluated$`Chull Scree` <- ifelse(all(is.na(res$`Chull Scree`)), yes = NA, no = which.max(res$`Chull Scree`) == true_clus)
  evaluated$BIC_G         <- which.min(res$BIC_G) == true_clus
  evaluated$BIC_N         <- which.min(res$BIC_N) == true_clus
  evaluated$AIC           <- which.min(res$AIC) == true_clus
  
  evaluated$`Chull Scree_fac` <- which.max(res$`Chull Scree_fac`) == true_clus
  evaluated$BIC_G_fac         <- which.min(res$BIC_G_fac) == true_clus
  evaluated$BIC_N_fac         <- which.min(res$BIC_N_fac) == true_clus
  evaluated$AIC_fac           <- which.min(res$AIC_fac) == true_clus
  
  return(evaluated)
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










