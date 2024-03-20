library(lavaan)

# Simulation Design
# Which factors are going to be tested? For now:
nclus   <- c(2, 4)         # Number of clusters
ngroups <- c(24, 48)       # Number of groups
coeff   <- c(0.3, 0.4)     # Initial regression parameters
N_g     <- c(50, 100, 200) # Sample size per groups
balance <- c("bal", "unb") # Cluster size
sd      <- c(0, 0.05, 0.1) # Differences within a cluster (in betas)

# reliability <- c("low")
# NonInvSize <- c(0.6)
# ResRange <- 0.2
# NonInvItems <- 2
# NonInvG <- c(0.50)
# NonInvType <- c("fixed")

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
design <- expand.grid(nclus, ngroups, coeff, N_g, balance, sd, model) # , reliability, NonInvSize, ResRange,
# NonInvItems, NonInvG, NonInvType)
colnames(design) <- c("nclus", "ngroups", "coeff", "N_g", "balance", "sd", "model")
# "reliability", "NonInvSize", "ResRange", "NonInvItems", "NonInvG", "NonInvType")

rownames(design) <- NULL
rm(balance, coeff, N_g, nclus, ngroups, sd) #, NonInvG, NonInvItems, NonInvSize, reliability, ResRange)

# Remove unnecessary rows
design <- design[design$ngroups == 24, ]

# Rename rownames
rownames(design) <- seq_len(nrow(design))

# Run mini "sim"
R2 <- vector(mode = "list", length = 72)
K <- 300

set.seed(1)
for(i in 1:72){
  print(i)
  R2[[i]] <- numeric(K)
  for(k in 1:K){
    print(k)
    R2[[i]][k] <- PopR2Entropy(model      = model, 
                               step1model = S1, 
                               nclus      = design[i, "nclus"], 
                               ngroups    = 192, 
                               reg_coeff  = design[i, "coeff"], 
                               N_g        = design[i, "N_g"], 
                               balance    = design[i, "balance"], 
                               sd         = design[i, "sd"])
  }
}

mean(R2)

design$R2 <- unlist(lapply(X = R2, FUN = mean))

summary(unlist(lapply(X = R2, FUN = mean)))


















