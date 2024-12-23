rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(ape) # general phylogenetic analysis
library(caper) # PGLS
library(geiger) # compare taxa in data and tree
# remotes::install_github("uyedaj/treeplyr")
library(treeplyr) # general phylogenetic analysis
library(ggtree) # phylogenetic tree visualization
library(aplot) # combine subplots
library(ggheatmap) #devtools::install_github("XiaoLuo-boy/ggheatmap")
library(pheatmap)
# library(rr2) # calculate r2 for PGLS models


# ------------------------------------------------------------------------------
# Run PGLS using the constructed phylotree
# ------------------------------------------------------------------------------
# Ref: 
# https://nhcooper123.github.io/macro-module-2020/phylogenetic-generalised-least-squares-pgls-in-r.html
load("traitDataFujian-spavg-phylo-step2.RData")

# phylotree_result[[1]] is the resulting phylogenetic tree (type: phylo).
phylotree_touse = phylotree_result[[1]]

# preliminary check, make sure all these statements are true
# is.binary(phylotree_touse)
# is.rooted(phylotree_touse)
# is.ultrametric(phylotree_touse)

# ------------------------------------------------------------------------------
# Step1: combine raw data and tree to a uniform format (using z-transformed and log-transformed values)
# ------------------------------------------------------------------------------
traitDataIndv_normalized$speciesFullNameConverted = ""
traitDataIndv_normalized = traitDataIndv_normalized %>%
  mutate(speciesFullNameConverted = gsub(" ", "_", traitDataIndv_normalized$speciesFullName))

# check that the species names match up in the tree and the data
phylotree_check_result = name.check(
  phy = phylotree_touse, data = traitDataIndv_normalized, data.names = traitDataIndv_normalized$speciesFullNameConverted)

# combine tree and data to exclude species that are not in both
# result tree: 40 tips and 39 internal nodes
traitPhyloDataSH = make.treedata(tree = phylotree_touse, data = traitDataIndv_normalized,
                                 name_column = "speciesFullNameConverted")
traitPhyloDataSH$dat$tiplabel = traitPhyloDataSH$phy$tip.label # restore the removed name column

# Save data and tree separately, facilitating downstream analyses.
traitPhyloDataSH_data = as.data.frame(traitPhyloDataSH$dat)
traitPhyloDataSH_tree = traitPhyloDataSH$phy

# Separate genus and species for each binominal
traitPhyloDataSH_data_sep = traitPhyloDataSH_data %>% separate(speciesFullName, into = c("genus_name", "sp_name"), sep=" ") %>%
  dplyr::select(genus_name, sp_name)

traitPhyloDataSH_data = cbind(traitPhyloDataSH_data, traitPhyloDataSH_data_sep)
rm(traitPhyloDataSH_data_sep) # clear intermediate dataframe

# set tree node label to NULL to avoid Error: Labels duplicated between tips and nodes in phylogeny
# Ref: https://stackoverflow.com/questions/51261388/duplicate-tips-label-removal-importing-phylogenetic-tree-in-r-for-comparison
traitPhyloDataSH_tree$node.label = NULL

# ------------------------------------------------------------------------------
# Step2: combine raw data and tree to a uniform format (using only log-transformed values)
# TODO: this part is modified 7/11-24, in response to Sandy's advice.
# ------------------------------------------------------------------------------

traitDataIndv_spavg_log$speciesFullNameConverted = ""
traitDataIndv_spavg_log = traitDataIndv_spavg_log %>%
  mutate(speciesFullNameConverted = gsub(" ", "_", traitDataIndv_spavg_log$speciesFullName))

# check that the species names match up in the tree and the data
phylotree_check_result = name.check(
  phy = phylotree_touse, data = traitDataIndv_spavg_log, data.names = traitDataIndv_spavg_log$speciesFullNameConverted)

# combine tree and data to exclude species that are not in both
# result tree: 40 tips and 39 internal nodes
traitPhyloDataSH_logonly = make.treedata(tree = phylotree_touse, data = traitDataIndv_spavg_log,
                                 name_column = "speciesFullNameConverted")
traitPhyloDataSH_logonly$dat$tiplabel = traitPhyloDataSH_logonly$phy$tip.label # restore the removed name column

# Save data and tree separately, facilitating downstream analyses.
traitPhyloDataSH_data_logonly = as.data.frame(traitPhyloDataSH_logonly$dat)
traitPhyloDataSH_tree_logonly = traitPhyloDataSH_logonly$phy

# Separate genus and species for each binominal
traitPhyloDataSH_data_sep = traitPhyloDataSH_data_logonly %>% separate(speciesFullName, into = c("genus_name", "sp_name"), sep=" ") %>%
  dplyr::select(genus_name, sp_name)

traitPhyloDataSH_data_logonly = cbind(traitPhyloDataSH_data_logonly, traitPhyloDataSH_data_sep)
rm(traitPhyloDataSH_data_sep) # clear intermediate dataframe

# set tree node label to NULL to avoid Error: Labels duplicated between tips and nodes in phylogeny
# Ref: https://stackoverflow.com/questions/51261388/duplicate-tips-label-removal-importing-phylogenetic-tree-in-r-for-comparison
traitPhyloDataSH_tree_logonly$node.label = NULL

# SAVE CONVERTED TREE FOR DOWNSTREAM ANALYSIS ----
save(traitPhyloDataSH, traitPhyloDataSH_tree, traitPhyloDataSH_data, 
     traitPhyloDataSH_logonly, traitPhyloDataSH_tree_logonly, traitPhyloDataSH_data_logonly, 
     traitDataIndv_normalized, traitDataIndv_spavg, traitDataIndv_spavg_log, traitDataIndv,
     file = "traitDataFujian-fmtdata-step3.RData")


# ------------------------------------------------------------------------------
# Step2: run PGLS analysis
# ------------------------------------------------------------------------------
rm(list=ls())
load("traitDataFujian-fmtdata-step3.RData")

# exploratory analysis: plot RTD against SLA
# Plot eyesize against body mass, coloured by family
# ggplot(traitPhyloDataSH_data, aes(x = `SLA(cm2/g)`, y = `RTD (g/cm3)`, colour = genus_name)) +
#   geom_point(na.rm = T) + theme_bw()

PGLSdata = comparative.data(phy = traitPhyloDataSH_tree, data = traitPhyloDataSH_data,
                            names.col = tiplabel, vcv = T, na.omit = F, warn.dropped = T)

# fit PGLS model: LNC & LMA
model.pgls.LNC_LMA <- pgls(LNC ~ LMA, data = PGLSdata, lambda = "ML")
par(mfrow = c(2, 2)) # Validation 1
plot(model.pgls.LNC_LMA) # Validation 1
anova.pgls(model.pgls.LNC_LMA)
plt.pgls.LNC_LMA = ggplot(traitPhyloDataSH_data, aes(x = LNC, y = LMA)) + 
  geom_point(color="#3f72af", alpha=0.8) +
  geom_abline(slope = coefficients(model.pgls.LNC_LMA)[2], intercept = coefficients(model.pgls.LNC_LMA)[1]) +
  geom_text(aes(label = paste0("PGLS: y=", round(coefficients(model.pgls.LNC_LMA)[2],3), "x", intercept = round(coefficients(model.pgls.LNC_LMA)[1],3)),x = 1,y = 1)) +
  labs(x="normalized(log(LNC))", y="normalized(log(LMA))") + 
  theme_bw()
ggsave(filename = "plt_pglsCorr_LNC_LMA.pdf", plot = plt.pgls.LNC_LMA, width = 4, height = 4)
model.pgls.LNC_LMA.summary = summary(model.pgls.LNC_LMA) # Adjusted R square


# fit PGLS model: RD & SRL
model.pgls.RD_SRL <- pgls(RD ~ SRL, data = PGLSdata, lambda = "ML")
par(mfrow = c(2, 2)) # Validation 1
plot(model.pgls.RD_SRL) # Validation 1
anova.pgls(model.pgls.RD_SRL)
plt.pgls.RD_SRL = ggplot(traitPhyloDataSH_data, aes(x = RD, y = SRL)) + 
  geom_point(color="#3f72af", alpha=0.8) +
  geom_abline(slope = coefficients(model.pgls.RD_SRL)[2], intercept = coefficients(model.pgls.RD_SRL)[1]) +
  geom_text(aes(label = paste0("PGLS: y=", round(coefficients(model.pgls.RD_SRL)[2],3), "x", intercept = round(coefficients(model.pgls.RD_SRL)[1],3)),x = 1,y = 1)) +
  labs(x="normalized(log(RD))", y="normalized(log(SRL))") + 
  theme_bw()
ggsave(filename = "plt_pglsCorr_RD_SRL.pdf", plot = plt.pgls.RD_SRL, width = 4, height = 4)
summary(model.pgls.RD_SRL) # Adjusted R square


# fit PGLS model: RTD & RNC
# remove outliers
PGLSdata_RTD_RNC_dataonly = traitPhyloDataSH_data %>%
  dplyr::filter(RNC < 3) %>%
  dplyr::filter(RTD > -2)
PGLSdata_RTD_RNC = comparative.data(phy = traitPhyloDataSH_tree, data = PGLSdata_RTD_RNC_dataonly,
                            names.col = tiplabel, vcv = T, na.omit = F, warn.dropped = T)

model.pgls.RTD_RNC <- pgls(RTD ~ RNC, data = PGLSdata_RTD_RNC, lambda = "ML")
par(mfrow = c(2, 2)) # Validation 1 
plot(model.pgls.RTD_RNC) # Validation 1
anova.pgls(model.pgls.RTD_RNC)
plt.pgls.RTD_RNC = ggplot(traitPhyloDataSH_data, aes(x = RTD, y = RNC)) + 
  geom_point(color="#3f72af", alpha=0.8) +
  geom_abline(slope = coefficients(model.pgls.RTD_RNC)[2], intercept = coefficients(model.pgls.RTD_RNC)[1]) +
  geom_text(aes(label = paste0("PGLS: y=", round(coefficients(model.pgls.RTD_RNC)[2],3), "x", intercept = round(coefficients(model.pgls.RTD_RNC)[1],3)),x = 0,y = 1)) +
  labs(x="normalized(log(RTD))", y="normalized(log(RNC))") + 
  theme_bw()
ggsave(filename = "plt_pglsCorr_RTD_RNC.pdf", plot = plt.pgls.RTD_RNC, width = 4, height = 4)
summary(model.pgls.RTD_RNC) # Adjusted R square


# fit PGLS model: SRR25 & Rdark25(P)
# remove outliers
PGLSdata_SRR25_Rdark25_dataonly = traitPhyloDataSH_data %>%
  # dplyr::filter(SRR25 < 2) %>%
  dplyr::filter(Rdark25 < 2)
PGLSdata_SRR25_Rdark25 = comparative.data(phy = traitPhyloDataSH_tree, data = PGLSdata_SRR25_Rdark25_dataonly,
                                    names.col = tiplabel, vcv = T, na.omit = F, warn.dropped = T)

model.pgls.SRR25_Rdark25 <- pgls(SRR25 ~ Rdark25P, data = PGLSdata, lambda = "ML")
par(mfrow = c(2, 2)) # Validation 1
plot(model.pgls.SRR25_Rdark25) # Validation 1
plot(pgls.profile(model.pgls.SRR25_Rdark25, "lambda")) # Validation 2
anova.pgls(model.pgls.SRR25_Rdark25) # ANOVA
plt.pgls.SRR25_Rdark25 = ggplot(traitPhyloDataSH_data, aes(x = SRR25, y = Rdark25)) + 
  geom_point(color="#3f72af", alpha=0.8) +
  geom_abline(slope = coefficients(model.pgls.SRR25_Rdark25)[2], intercept = coefficients(model.pgls.SRR25_Rdark25)[1]) +
  geom_text(aes(label = paste0("PGLS: y=", round(coefficients(model.pgls.SRR25_Rdark25)[2],3), "x", intercept = round(coefficients(model.pgls.SRR25_Rdark25)[1],3)),x = 1,y = 1)) +
  labs(x="normalized(log(SRR25))", y="normalized(log(Rdark25))") + 
  theme_bw()
ggsave(filename = "plt_pglsCorr_SRR25_Rdark25.pdf", plot = plt.pgls.SRR25_Rdark25, width = 4, height = 4)
summary(model.pgls.SRR25_Rdark25) # Adjusted R square


# fit PGLS model: SRR25 & rs25
model.pgls.SRR25_rs25 <- pgls(SRR25 ~ rs25, data = PGLSdata, lambda = 0.86)
# model.pgls.SRR25_rs25 <- pgls(SRR25 ~ rs25, data = PGLSdata, lambda = "ML")
par(mfrow = c(2, 2)) # Validation 1
plot(model.pgls.SRR25_rs25) # Validation 1
plot(pgls.profile(model.pgls.SRR25_rs25, "lambda")) # Validation 2
anova.pgls(model.pgls.SRR25_rs25) # ANOVA
plt.pgls.SRR25_rs25 = ggplot(traitPhyloDataSH_data, aes(x = SRR25, y = rs25)) + 
  geom_point(color="#3f72af", alpha=0.8) +
  geom_abline(slope = coefficients(model.pgls.SRR25_rs25)[2], intercept = coefficients(model.pgls.SRR25_rs25)[1]) +
  geom_text(aes(label = paste0("PGLS: y=", round(coefficients(model.pgls.SRR25_rs25)[2],3), "x", intercept = round(coefficients(model.pgls.SRR25_Rdark25)[1],3)),x = 1,y = 1)) +
  labs(x="normalized(log(SRR25))", y="normalized(log(rs25))") + 
  theme_bw()
ggsave(filename = "plt_pglsCorr_SRR25_rs25.pdf", plot = plt.pgls.SRR25_rs25, width = 4, height = 4)
summary(model.pgls.SRR25_rs25) # Adjusted R square


# # fit PGLS model: SRR25 & rs25
# model.pgls.Vcmax25_LNC <- pgls(Vcmax25 ~ LNC, data = PGLSdata, lambda = "ML")
# par(mfrow = c(2, 2)) # Validation 1
# plot(model.pgls.Vcmax25_LNC) # Validation 1
# plot(pgls.profile(model.pgls.Vcmax25_LNC, "lambda")) # Validation 2
# anova.pgls(model.pgls.Vcmax25_LNC) # ANOVA
# plt.pgls.Vcmax25_LNC = ggplot(traitPhyloDataSH_data, aes(x = Vcmax25, y = LNC)) + 
#   geom_point(color="#3f72af", alpha=0.8) +
#   geom_abline(slope = coefficients(model.pgls.Vcmax25_LNC)[2], intercept = coefficients(model.pgls.Vcmax25_LNC)[1]) +
#   geom_text(aes(label = paste0("PGLS: y=", round(coefficients(model.pgls.Vcmax25_LNC)[2],3), "x", intercept = round(coefficients(model.pgls.Vcmax25_LNC)[1],3)),x = 1,y = 1)) +
#   labs(x="normalized(log(Vcmax25))", y="normalized(log(LNC))") + 
#   theme_bw()
# ggsave(filename = "plt_pglsCorr_SRR25_rs25.pdf", plot = plt.pgls.SRR25_rs25, width = 4, height = 4)


# ------------------------------------------------------------------------------
# Batch PGLS implementation to test variable correlation
# ------------------------------------------------------------------------------
# Define the list of variables
variables <- c("RTD","SRL","RD","RNC","RCC","RPC","SRR25","RDMC","LMA","LNC","LCC","LPC","Rdark25P", "Rlight25P","rs25","Vcmax25","Asat") # original

variables <- c("vH","WD","Ks","TLP","SWCleaf","SWCbranch","H","DBH", "RTD","SRL","RD","RNC","RCC","RPC","SRR25","RDMC") # hydraulic

# Initialize an empty matrix to store the Pr values
pr_matrix <- matrix(NA, nrow = length(variables), ncol = length(variables), 
                    dimnames = list(variables, variables))

# Loop through each pair of variables
for (i in 1:length(variables)) {
  for (j in 1:length(variables)) {
    if (i != j) {
      # Fit the PGLS model for the variable pair
      formula <- as.formula(paste(variables[i], "~", variables[j]))
      model.pgls <- tryCatch({
        pgls(formula, data = PGLSdata, lambda = "ML")
      }, error = function(e) return(NULL)) # Error handling
      
      # Perform ANOVA and extract Pr value if the model is fit
      if (!is.null(model.pgls)) {
        anova_result <- anova.pgls(model.pgls)
        pr_value <- anova_result[1, 5]  # Extract Pr(>F)
        pr_matrix[i, j] <- pr_value
      }
    }
  }
}

# Convert the Pr matrix into a heatmap-friendly format
# pr_matrix[is.na(pr_matrix)] <- 1  # Set NA values to 1 for the heatmap
pr_matrix_rmNS = pr_matrix
pr_matrix_rmNS = round(pr_matrix_rmNS, 2)
pr_matrix_rmNS = ifelse(pr_matrix>0.1,0.1,pr_matrix_rmNS)



# Plot the heatmap using pheatmap
#pdf("plt_PGLS_Prhtmap03.pdf", width = 7, height = 6.5) # original
pdf("plt_PGLS_Prhtmap04_hydraulicRES.pdf", width = 8, height = 7) # hydraulic
pheatmap(pr_matrix_rmNS,
         # color = colorRampPalette(c("#BF3131", "#FFFFFF"))(100),
         color = colorRampPalette(c("#006363", "#FFFFFF"))(100),
         display_numbers = TRUE, number_format = "%.2f",
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "PGLS: Pr(>F)"
         )
dev.off()

