# ------------------------------------------------------------------------------
# Plot phylogeny tree, generate trait value heat map besides the tree
# ------------------------------------------------------------------------------
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

load("traitDataFujian-spavg-phylo-step2.RData")

### SET DATA TO USE ###
traitData_touse = traitDataIndv_spavg_log

# Perform data normalization:
# not log and z-score transformation, but to normalize value in each column to range 0-1

# Function to normalize columns in a dataframe to range 0-1, based on a column range
func_normalize_0to1_bycol <- function(x) {
  # validate if a column is numeric, if yes, do calculation, if not, return the original value
  if (!is.numeric(x)) {
    return(x)
  } else {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
}

traitData_touse = as.data.frame(lapply(traitData_touse, func_normalize_0to1_bycol))

# Start phylo data preparation ---------------------------------------------------

# Convert species Latin name `A B` to `A_B` to match with tree `tip_label`
traitData_touse = traitData_touse %>% dplyr::mutate(tipLabelMatched=gsub(" ", "_", speciesFullName)) %>%
  dplyr::select(tipLabelMatched, everything())

# Extract phylo-tree constructed in step2
# phylotree_result[[1]] is the resulting phylogenetic tree (type: phylo).
phylotree_touse = phylotree_result[[1]] 
rm(phylotree_result)

# Combine tree and data
traitData_touse_PhyloDataCombined = make.treedata(
  tree = phylotree_touse, data = traitData_touse, name_column = "tipLabelMatched")

# Save data and tree separately, facilitating downstream analyses.
traitData_touse_PhyloDataCombined_data = as.data.frame(traitData_touse_PhyloDataCombined$dat) %>% 
  dplyr::mutate(tipLabelMatched=gsub(" ", "_", speciesFullName)) %>%
  dplyr::select(tipLabelMatched, everything())
traitData_touse_PhyloDataCombined_tree = traitData_touse_PhyloDataCombined$phy

# End phylo data preparation ---------------------------------------------------

# Plot phylogenetic tree
plt.phylo.tree = ggtree(traitData_touse_PhyloDataCombined_tree, color="darkblue", layout = "rectangular") # + geom_tiplab(size=2, color="black")

pltdata.phylo.heatmap = traitData_touse_PhyloDataCombined_data %>% dplyr::select(-speciesFullName, family_name, -genus_name, -sp_name, speciesFullName_ori, GrowthForm)
pltdata.phylo.heatmap[is.na(pltdata.phylo.heatmap)] = NA
pltdata.phylo.heatmap_tiplabel = pltdata.phylo.heatmap$tipLabelMatched

# Specify traits to draw plot ------------------------------------------------
pltdata.phylo.heatmap = pltdata.phylo.heatmap %>% dplyr::select(-tipLabelMatched) %>%
  dplyr::select(
    RTD,SRL,RD,RNC,SRA,RPC,RCC,RDMC,SRR25, # FR traits
    LMA,LNC,LCC,LPC,Rdark25P,Rlight25P,Vcmax25,Asat, # leaf traits
    rs25, # stem traits
    vH,WD,Ks,TLP,SWCleaf,SWCbranch,H,DBH # hydraulic traits
  ) %>% # TRAITS
  dplyr::rename(
    R_RTD = RTD, R_SRL = SRL, R_RD = RD, R_RNC = RNC, R_SRA = SRA, R_RPC = RPC, R_RCC = RCC, R_RDMC = RDMC, R_SRR25 = SRR25,
    L_LMA = LMA, L_LNC = LNC, L_LCC = LCC, L_LPC = LPC,
    L_Rdark25 = Rdark25P, # need rename for shorter names
    L_Rlight25 = Rlight25P, # need rename for shorter names
    L_Vcmax25 = Vcmax25, L_Asat = Asat,
    S_rs25 = rs25,
    H_vH = vH, H_WD = WD, H_Ks = Ks, H_TLP = TLP, H_SWCleaf = SWCleaf, H_SWCbrch = SWCbranch, H_H = H, H_DBH = DBH
  )

# pltdata.phylo.heatmap[pltdata.phylo.heatmap>3] = 3 # normalize values that is too large (Added 10/10-24)
rownames(pltdata.phylo.heatmap) = pltdata.phylo.heatmap_tiplabel

# Curate annotation rows (trait kind)
pltdata.phylo.heatmap.annot_traitkind = data.frame(
  TraitKind = rep(c("root (R_)","leaf (L_)","stem (S_)","hydraulic (H_)"), times = c(9,8,1,8))
)
rownames(pltdata.phylo.heatmap.annot_traitkind) = colnames(pltdata.phylo.heatmap)

# Curate annotation columns (life form)
pltdata.phylo.heatmap.annot_lifeform = traitData_touse_PhyloDataCombined_data %>% dplyr::select(GrowthForm) 
pltdata.phylo.heatmap.annot_lifeform = pltdata.phylo.heatmap.annot_lifeform %>%
  mutate(GrowthForm = as.factor(pltdata.phylo.heatmap.annot_lifeform$GrowthForm))
rownames(pltdata.phylo.heatmap.annot_lifeform) = pltdata.phylo.heatmap_tiplabel

# Curate color
growthFormCol <- c("#F9E400","#D4352D","#FFAF00","#CC976B")
names(growthFormCol) <- c("shrub","tree", "tree+shrub", "liana")
traitKindCol <- c("#FFE1FF","#E4B1F0","#7E60BF","#433878")
names(traitKindCol) <- c("roo
                         t (R_)","leaf (L_)","stem (S_)","hydraulic (H_)")
color_list = list(GrowthForm=growthFormCol, TraitKind=traitKindCol)

# Start plotting
plt.phylo.heatmap = ggheatmap(pltdata.phylo.heatmap, scale = "none",
                              annotation_rows = pltdata.phylo.heatmap.annot_lifeform,
                              annotation_cols = pltdata.phylo.heatmap.annot_traitkind,
                              annotation_color = color_list,
                              border = "grey",
                              cluster_rows = F, cluster_cols = F,
                              color = colorRampPalette(c("#91CAE8","#FFFFFF","#F48892"))(100),
                              legendName = "log(TraitValue)\n(normalized to 0-1)")

plt.phylo.heatmap = plt.phylo.heatmap %>% ggheatmap_theme(1, # the heatmap
                theme =list(
                  theme(axis.text.x = element_text(angle = 90, hjust=1, face = "bold"),
                        axis.text.y = element_text(colour = "black",face = "italic")),
                  theme(legend.title = element_text(face = "bold")),
                  theme(legend.title = element_text(face = "bold")),
                  theme(legend.title = element_text(face = "bold")),
                  theme(legend.title = element_text(face = "bold"))
                ))
# ggheatmap_plotlist(plt.phylo.heatmap)


# Plot phylogenentic tree besides heatmap to show genus-level trait conservation patterns
plt.phylo.total = plt.phylo.heatmap %>% insert_left(plt.phylo.tree, width = 0.25)
ggsave(filename = "plt_phylotree_heatmap_all.pdf", plot = plt.phylo.total, width = 8.2, height = 6)
