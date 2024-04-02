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

# HEATMAP
a <- count_results(data = Results_final, by = c("sd", "nclus"), type = "relative") %>% filter(result == "Correct")

a1 <- a %>% dplyr::select(sd, nclus, Chull:ICL) %>% pivot_longer(cols = Chull:ICL, names_to = "Measure", values_to = "Proportion")

a1$Measure <- factor(a1$Measure, levels = c("BIC_N", "ICL", "BIC_G", "AIC3", "AIC", "Chull"))

plot <- ggplot(data = a1, aes(x = sd, y = Measure)) + facet_grid(~nclus) +
  geom_tile(aes(fill = Proportion)) + geom_text(aes(label = Proportion), size = 3.2) + 
  scale_fill_gradient(low = "yellow", high = "green4") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of clusters", breaks = NULL, labels = NULL)) +
  labs(x = expression("Within-cluster differences (" * sigma[beta] * ")"),  # Combines text with Greek letter, no space
       y = expression("Model Selection measure")) +
  scale_y_discrete(labels = c("Chull" = "CHull",
                              "AIC" = "AIC",
                              "AIC3" = expression(AIC[3]),
                              "BIC_G" = expression(BIC[G]),
                              "ICL" = "ICL",
                              "BIC_N" = expression(BIC[N])))
  # +
  # scale_y_discrete(sec.axis = sec_axis(~ . , name = "Number of clusters", breaks = NULL, labels = NULL))

plot

# TABLES
K_res   <- count_results(data = Results_final, by = c("nclus"), type = "relative")   %>% select(nclus:ICL)
N_res   <- count_results(data = Results_final, by = c("N_g"), type = "relative")     %>% select(N_g:ICL)
G_res   <- count_results(data = Results_final, by = c("ngroups"), type = "relative") %>% select(ngroups:ICL)
B_res   <- count_results(data = Results_final, by = c("coeff"), type = "relative")   %>% select(coeff:ICL)
Bal_res <- count_results(data = Results_final, by = c("balance"), type = "relative") %>% select(balance:ICL)
sd_res  <- count_results(data = Results_final, by = c("sd"), type = "relative")      %>% select(sd:ICL)
to_res  <- count_results(data = Results_final, by = "total", type = "relative")      %>% select(result:ICL)

K_res   <- K_res %>% 
  pivot_longer(cols = -c(nclus, result), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = nclus, values_from = value, names_prefix = "nclus") %>% 
  relocate(metric) %>% arrange(metric)
N_res   <- N_res %>% 
  pivot_longer(cols = -c(N_g, result), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = N_g, values_from = value, names_prefix = "N_g") %>% 
  relocate(metric) %>% arrange(metric) %>% select(contains("N_g"))
G_res   <- G_res %>% 
  pivot_longer(cols = -c(ngroups, result), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = ngroups, values_from = value, names_prefix = "ngroups") %>% 
  relocate(metric) %>% arrange(metric) %>% select(contains("ngroups"))
B_res   <- B_res %>% 
  pivot_longer(cols = -c(coeff, result), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = coeff, values_from = value, names_prefix = "coeff") %>% 
  relocate(metric) %>% arrange(metric) %>% select(contains("coeff"))
Bal_res <- Bal_res %>% 
  pivot_longer(cols = -c(balance, result), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = balance, values_from = value, names_prefix = "balance") %>% 
  relocate(metric) %>% arrange(metric) %>% select(contains("balance"))
sd_res  <- sd_res %>% 
  pivot_longer(cols = -c(sd, result), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = sd, values_from = value, names_prefix = "sd") %>% 
  relocate(metric) %>% arrange(metric) %>% select(contains("sd"))
to_res  <- to_res %>% 
  pivot_longer(cols = -c(result), names_to = "metric", values_to = "value") %>%
  relocate(metric) %>% arrange(metric) %>% select(value)

final <- cbind(K_res, N_res, G_res, B_res, Bal_res, sd_res, to_res)

# Add totals per column
total_col <- cbind("total", final %>% group_by(result) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE)))
colnames(total_col)[1] <- "metric"; colnames(total_col)[ncol(total_col)] <- "value"
final <- rbind(final, total_col)
print(xtable(final), include.rownames=FALSE)

####################################################################################################
############################################### ARI ################################################
####################################################################################################
# Load a second Results_final
ARI_res <- Results_final %>% dplyr::select(Condition:sd)

# Modify for the ARI results
ARI_res$ARI <- NA; ARI_res$CC <- NA

# Fill in the matrix with all results
ncond <- unique(ARI_res$Condition) # How many conditions?
K <- length(unique(ARI_res$Replication)) # How many replications?

for (i in ncond) {
  test <- NA
  test <- try(load(paste0("ResultRowARI", i, ".Rdata")))
  if(!c(class(test) == "try-error")){
    ARI_res[(K*(i-1)+1):(i*K), 9:10] <- ResultsRowARI
  }
}

Results_final$ARI <- ARI_res$ARI

final_ARI <- cbind(
  Results_final %>% group_by(Chull) %>% summarise(across(ARI, mean)) %>% rename(Result = Chull, Chull = ARI) %>% drop_na(),
  Results_final %>% group_by(AIC)   %>% summarise(across(ARI, mean)) %>% rename(Result = AIC,   AIC   = ARI) %>% select(AIC),
  Results_final %>% group_by(AIC3)  %>% summarise(across(ARI, mean)) %>% rename(Result = AIC3,  AIC3  = ARI) %>% select(AIC3),
  Results_final %>% group_by(BIC_G) %>% summarise(across(ARI, mean)) %>% rename(Result = BIC_G, BIC_G = ARI) %>% select(BIC_G),
  Results_final %>% group_by(BIC_N) %>% summarise(across(ARI, mean)) %>% rename(Result = BIC_N, BIC_N = ARI) %>% select(BIC_N),
  Results_final %>% group_by(ICL)   %>% summarise(across(ARI, mean)) %>% rename(Result = ICL,   ICL   = ARI) %>% select(ICL)  
)

final_ARI %>% select(Result, AIC, AIC3, BIC_G, BIC_N, Chull, ICL) %>% xtable() %>% print(., include.rownames=FALSE)







# # BAR GRAPHS
# # By number of clusters
# a <- count_results(data = Results_final, by = "nclus", type = "relative") %>% filter(result == "Correct")
# 
# a1 <- a %>% dplyr::select(nclus, Chull:ICL) %>% pivot_longer(cols = Chull:ICL, names_to = "Measure", values_to = "Value")
# a2 <- a %>% dplyr::select(nclus, Chull_fac:ICL_fac) %>% pivot_longer(cols = Chull_fac:ICL_fac, names_to = "Measure", values_to = "Value")
# 
# plot1 <- ggplot(data = a1, aes(x = Measure, y = Value)) + facet_grid(~nclus) +  
#   geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 4)
# 
# plot2 <- ggplot(data = a2, aes(x = Measure, y = Value)) + facet_grid(~nclus) +  
#   geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 4)
# 
# ggarrange(plotlist = list(plot1, plot2), common.legend = T, legend = "bottom", nrow = 1)
# 
# # By sd
# a <- count_results(data = Results_final, by = "sd", type = "relative") %>% filter(result == "Correct")
# 
# a1 <- a %>% dplyr::select(sd, Chull:ICL) %>% pivot_longer(cols = Chull:ICL, names_to = "Measure", values_to = "Value")
# a2 <- a %>% dplyr::select(sd, Chull_fac:ICL_fac) %>% pivot_longer(cols = Chull_fac:ICL_fac, names_to = "Measure", values_to = "Value")
# 
# plot1 <- ggplot(data = a1, aes(x = Measure, y = Value)) + facet_grid(~sd) +  
#   geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 3)
# 
# plot2 <- ggplot(data = a2, aes(x = Measure, y = Value)) + facet_grid(~sd) +  
#   geom_col(aes(fill = Measure)) + scale_y_continuous(limits = c(0,1)) + geom_text(aes(label = Value), vjust = -0.5, size = 3)
# 
# ggarrange(plotlist = list(plot1, plot2), common.legend = T, legend = "bottom", nrow = 1)



