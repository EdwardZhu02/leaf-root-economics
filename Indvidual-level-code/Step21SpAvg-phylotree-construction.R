rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
# devtools::install_github("jinyizju/V.PhyloMaker2")
library(V.PhyloMaker2) # phylo-tree construction

load(file = "2502-indv-level-code/traitDataFujian-SpAvg-step1.RData")

# DEPRECATED: USE ORIGINAL CODE INSTEAD.


# # Step4: Build phylogenetic tree using PhyloMaker2 -----------------------------
# # TODO: Update 20/11-24: all species are now bound to the tree after the addition of
# # family names from genus names
# spTaxaNames_tobuild = spTaxaNames_withfam %>%
#   dplyr::rename(species = speciesFullName, genus = genus_name, family = family_name) %>%
#   dplyr::select(-sp_name) %>%
#   dplyr::select(species,genus,family) # re-order
# 
# phylotree_result = phylo.maker(sp.list = spTaxaNames_tobuild, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL)
# # phylotree_result[[1]] is the resulting phylogenetic tree.
# # summary(phylotree_result[[1]])
