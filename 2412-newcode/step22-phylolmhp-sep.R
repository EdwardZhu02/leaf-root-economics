rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(ape) # general phylogenetic analysis
library(treeplyr) # general phylogenetic analysis
library(phylolm.hp) # devtools::install_github('laijiangshan/phylolm.hp',build_vignettes = TRUE)

# Load: traitDataIndv, traitDataIndv_spavg, traitDataIndv_spavg_log, traitDataIndv_spavg_log_ztransform, spTaxaNames_withfam, phylotree_result
load("traitDataFujian-spavg-phylo-step2.RData")

# Part 1: Data preparation -----------------------------------------------------
# 
### SET DATA TO USE ###
traitData_touse = traitDataIndv_spavg_log

# Convert species Latin name `A B` to `A_B` to match with tree `tip_label`
traitData_touse = traitData_touse %>% dplyr::mutate(tipLabelMatched=gsub(" ", "_", speciesFullName)) %>%
  dplyr::select(tipLabelMatched, everything())

# Extract phylo-tree constructed in step2
# phylotree_result[[1]] is the resulting phylogenetic tree (type: phylo).
phylotree_touse = phylotree_result[[1]] 
rm(phylotree_result)


# Part 2: Curate trait data for analysis ---------------------------------------
# 
phylolm_data = traitData_touse %>% dplyr::select(LMA,LNC,LPC)
rownames(phylolm_data) = traitData_touse$tipLabelMatched

phylolm_fit = phylolm(LMA ~ LNC + LPC,
                      data=phylolm_data, phy=phylotree_touse, model="lambda")
# Extract chi
phylolm_result_chi = phylolm_fit[["X"]]

phylolm_fithp =  phyloglm.hp(phylolm_fit, commonality=TRUE)
plot(phylolm_fithp,commonality=TRUE)

