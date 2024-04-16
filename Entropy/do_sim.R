library(lavaan)
library(fpp3)
library(ggplot2)
library(scales)

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

save(R2, file = "popR2.Rdata")

design$R2 <- unlist(lapply(X = R2, FUN = mean))

summary(unlist(lapply(X = R2, FUN = mean)))

res_nclus   <- design %>% group_by(nclus)   %>% summarise(across(R2, mean)) %>% add_column(Factor = "nclus", .before = "nclus") %>% rename("Level" = "nclus")
res_sd      <- design %>% group_by(sd)      %>% summarise(across(R2, mean)) %>% add_column(Factor = "sd", .before = "sd") %>% rename("Level" = "sd")
res_N_g     <- design %>% group_by(N_g)     %>% summarise(across(R2, mean)) %>% add_column(Factor = "N_g", .before = "N_g") %>% rename("Level" = "N_g")
res_coeff   <- design %>% group_by(coeff)   %>% summarise(across(R2, mean)) %>% add_column(Factor = "coeff", .before = "coeff") %>% rename("Level" = "coeff")
res_balance <- design %>% group_by(balance) %>% summarise(across(R2, mean)) %>% add_column(Factor = "balance", .before = "balance") %>% rename("Level" = "balance")

# Table
res_r2 <- rbind(res_N_g, res_sd, res_coeff) %>% select(Level, R2) %>% pivot_wider(., names_from = "Level", values_from = "R2")

# Plot
res_r2_long <- design %>% group_by(sd, coeff) %>% summarise(across(R2, mean))
ggplot(data = res_r2_long, mapping = aes(x = sd, y = R2)) + 
  facet_wrap(~coeff) + labs(x = "Within-cluster Differences", y = "R2 Entropy") +
  geom_point(color = "green4") + 
  geom_line(linewidth = 0.5, color = "green4") + 
  scale_x_continuous(labels = scales::label_number(scale = 1, accuracy = 0.01),
                     sec.axis = sec_axis(~ . , name = "Size of Regression Parameter", breaks = NULL, labels = NULL))




















