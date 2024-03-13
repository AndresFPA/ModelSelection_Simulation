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
  test <- try(load(paste0("ResultRow", i, ".Rdata")))
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

# BAR GRAPHS
# By number of clusters
a <- count_results(data = Results_final, by = "nclus", type = "relative") %>% filter(result == "Correct")

a1 <- a %>% dplyr::select(nclus, Chull:ICL) %>% pivot_longer(cols = Chull:ICL, names_to = "Measure", values_to = "Value")
a2 <- a %>% dplyr::select(nclus, Chull_fac:ICL_fac) %>% pivot_longer(cols = Chull_fac:ICL_fac, names_to = "Measure", values_to = "Value")

plot1 <- ggplot(data = a1, aes(x = Measure, y = Value)) + facet_grid(~nclus) +  
  geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 4)

plot2 <- ggplot(data = a2, aes(x = Measure, y = Value)) + facet_grid(~nclus) +  
  geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 4)

ggarrange(plotlist = list(plot1, plot2), common.legend = T, legend = "bottom", nrow = 1)

# By sd
a <- count_results(data = Results_final, by = "sd", type = "relative") %>% filter(result == "Correct")

a1 <- a %>% dplyr::select(sd, Chull:ICL) %>% pivot_longer(cols = Chull:ICL, names_to = "Measure", values_to = "Value")
a2 <- a %>% dplyr::select(sd, Chull_fac:ICL_fac) %>% pivot_longer(cols = Chull_fac:ICL_fac, names_to = "Measure", values_to = "Value")

plot1 <- ggplot(data = a1, aes(x = Measure, y = Value)) + facet_grid(~sd) +  
  geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 3)

plot2 <- ggplot(data = a2, aes(x = Measure, y = Value)) + facet_grid(~sd) +  
  geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 3)

ggarrange(plotlist = list(plot1, plot2), common.legend = T, legend = "bottom", nrow = 1)

# HEATMAP
a <- count_results(data = Results_final, by = c("sd", "N_g", "nclus"), type = "relative") %>% filter(result == "Correct")

a1 <- a %>% dplyr::select(N_g, sd, nclus, Chull:ICL) %>% pivot_longer(cols = Chull:ICL, names_to = "Measure", values_to = "Value")

plot <- ggplot(data = a1, aes(x = sd, y = Measure)) + facet_grid(nclus~N_g) +
  geom_tile(aes(fill = Value)) + geom_text(aes(label = Value)) + 
  scale_fill_gradient(low = "yellow", high = "red") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Sample Size", breaks = NULL, labels = NULL)) # +
  # scale_y_discrete(sec.axis = sec_axis(~ . , name = "Number of clusters", breaks = NULL, labels = NULL))

plot




# TABLES
K_res   <- count_results(data = Results_final, by = c("nclus"), type = "relative")   %>% select(nclus:ICL)
N_res   <- count_results(data = Results_final, by = c("N_g"), type = "relative")     %>% select(N_g:ICL)
G_res   <- count_results(data = Results_final, by = c("ngroups"), type = "relative") %>% select(ngroups:ICL)
B_res   <- count_results(data = Results_final, by = c("coeff"), type = "relative")   %>% select(coeff:ICL)
Bal_res <- count_results(data = Results_final, by = c("balance"), type = "relative") %>% select(balance:ICL)
sd_res  <- count_results(data = Results_final, by = c("sd"), type = "relative")      %>% select(sd:ICL)

K_res   <- data.table::transpose(K_res, keep.names = "rn")
N_res   <- data.table::transpose(N_res, keep.names = "rn")
G_res   <- data.table::transpose(G_res, keep.names = "rn")
B_res   <- data.table::transpose(B_res, keep.names = "rn")
Bal_res <- data.table::transpose(Bal_res, keep.names = "rn")
sd_res  <- data.table::transpose(sd_res, keep.names = "rn")

list2 <- list(a, b, c, d, e, f, g, h)
current <- c()

for(i in 1:length(list2)){
  tmp <- list2[[i]]
  tmp$Factor <- colnames(tmp)[1]
  colnames(tmp)[1] <- c("Level")
  tmp$Level <- as.factor(tmp$Level)
  tmp <- tmp[, c(1, ncol(tmp), c(2:(ncol(tmp) - 1)))]
  current <- rbind(current, tmp)
}

rm(a, b, c, d, e, f, g, h)
current <- current[, c("Factor", "Level", "ARI", "CorrectClus", "fARI", 
                       "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4")]
# for_paper <- current %>% pivot_wider(names_from = `Non-inv Included`, values_from = c(ARI:RMSE_C)) %>% 
#   dplyr::select(Factor, Level, ARI_yes, CorrectClus_yes, fARI_yes, RMSE_A_yes, RMSE_B_yes, RMSE_C_yes)
for_paper <- current
colnames(for_paper) <- c("Factor", "Level", "ARI", "CorrectClus", "fARI", "beta1", "beta2", "beta3", "beta4")

# Re-organize table
# Get row indices of each factor
factors <- unique(for_paper$Factor)
for(i in 1:length(factors)){ 
  assign(x = paste0("rn_", factors[i]), value = which(for_paper$Factor == factors[i]))
}

for_paper <- for_paper[c(rn_coeff, rn_ngroups, rn_N_g, rn_nclus,
                         rn_balance, rn_reliability, rn_NonInvG,
                         rn_NonInvSize), ]

rm(rn_coeff, rn_ngroups, rn_N_g, rn_nclus, rn_balance, rn_reliability, rn_NonInvG, rn_NonInvSize)

# Add total - Cluster
yes_tot <- t(apply(Results_final[, c("ARI", "CorrectClus", "fARI")], 2, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""
yes_tot <- yes_tot[, c("Factor", "Level", "ARI", "CorrectClus", "fARI")]

#Final - Cluster
for_paper_clus <- for_paper[, c("Factor", "Level", "ARI", "CorrectClus", "fARI")]
for_paper_clus <- rbind(for_paper_clus, yes_tot)
print(xtable(for_paper_clus, digits = 3), include.rownames = F)

rm(yes_tot, tmp)

# Add total - Parameter
yes_tot <- t(apply(Results_final[, c("RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4")], 2, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""
yes_tot <- yes_tot[, c("Factor", "Level", "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4")]
colnames(yes_tot) <- c("Factor", "Level", "beta1", "beta2", "beta3", "beta4")

#Final - Parameter
for_paper_par <- for_paper[, c("Factor", "Level", "beta1", "beta2", "beta3", "beta4")]
for_paper_par <- rbind(for_paper_par, yes_tot)
# for_paper_par[, c("beta1", "beta2", "beta3", "beta4")] <- round(for_paper_par[, c("beta1", "beta2", "beta3", "beta4")], digits = 3)

rm(yes_tot)

print(xtable(for_paper_par, digits = 3), include.rownames = F)

####################################################################################################
##################################### TABLE 3 - GLOBAL MAXIMA ######################################
####################################################################################################
Global <- Results_final
Global$`Global Max %` <- ifelse(test = Global$`Global Max %` == 0, yes = 0, no = 1)

# Check mean results per simulation factor
a <- Global %>% group_by(NonInvIncl, nclus) %>% summarise(across(MisClass:Exo_var, mean))
b <- Global %>% group_by(NonInvIncl, ngroups) %>% summarise(across(MisClass:Exo_var, mean))
c <- Global %>% group_by(NonInvIncl, N_g) %>% summarise(across(MisClass:Exo_var, mean))
d <- Global %>% group_by(NonInvIncl, coeff) %>% summarise(across(MisClass:Exo_var, mean))
e <- Global %>% group_by(NonInvIncl, balance) %>% summarise(across(MisClass:Exo_var, mean))
f <- Global %>% group_by(NonInvIncl, reliability) %>% summarise(across(MisClass:Exo_var, mean))
g <- Global %>% group_by(NonInvIncl, NonInvSize) %>% summarise(across(MisClass:Exo_var, mean))
h <- Global %>% group_by(NonInvIncl, NonInvG) %>% summarise(across(MisClass:Exo_var, mean))

list2 <- list(a, b, c, d, e, f, g, h)
current <- c()

for(i in 1:length(list2)){
  tmp <- list2[[i]]
  tmp$Factor <- colnames(tmp)[2]
  colnames(tmp)[1:2] <- c("Non-inv Included", "Level")
  tmp$Level <- as.factor(tmp$Level)
  tmp <- tmp[, c(1, ncol(tmp), c(2:(ncol(tmp) - 1)))]
  current <- rbind(current, tmp)
}

rm(a, b, c, d, e, f, g, h)

Global2 <- current[, c("Factor", "Level", "Global Max %")]
1 - mean(Global2$`Global Max %`)
mean(Global$`Global Max %`)

Local <- Global %>% filter(`Global Max %` == 0) %>% select(Condition:NonInvG, `Global Max %`)
Local %>% count(nclus)
Local %>% count(ngroups)
Local %>% count(N_g)
Local %>% count(coeff)
Local %>% count(balance)
Local %>% count(reliability)
Local %>% count(NonInvSize)
Local %>% count(NonInvG)

####################################################################################################
######################################## COR fARI - RMSE ###########################################
####################################################################################################

# Total correlations
Results_final %>% select(fARI, RMSE_B1:RMSE_B4) %>% cor()

# Correlation depending on the number of clusters
Results_final %>% filter(nclus == 2) %>% select(fARI, RMSE_B1:RMSE_B4) %>% cor()
Results_final %>% filter(nclus == 4) %>% select(fARI, RMSE_B1:RMSE_B4) %>% cor()









