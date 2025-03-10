rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(ape) # general phylogenetic analysis
library(geiger) # compare taxa in data and tree
library(treeplyr) # general phylogenetic analysis
library(phytools) # general phylogenetic analysis, phyl.pca
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(ggplotify)
library(ggrepel)
library(ggbiplot) # devtools::install_github("vqv/ggbiplot")

library(patchwork) # plot merging
library(cowplot) # plot merging

library(reshape2) # scree plot, stacked barplot

# Run phylogenetically-informed PCA using the constructed phylotree
# Phylo-PCA need the data to be scaled to unit variance, so a z-score transformation is needed prior to analyses.
# 
# TODO: Updated 21/11-24, use original and perform binding in the latter script
# Load: traitDataIndv, traitDataIndv_spavg, traitDataIndv_spavg_log, traitDataIndv_spavg_log_ztransform, spTaxaNames_withfam, phylotree_result
load("traitDataFujian-spavg-phylo-step2.RData") 

# Part 1: combine raw data and tree --------------------------------------------
# 
### SET DATA TO USE ###
traitData_touse = traitDataIndv_spavg_log # perform scaling later (line 71, 7/12-24)
### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC")


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

# set tree node label to NULL to avoid Error: Labels duplicated between tips and nodes in phylogeny
# Ref: https://stackoverflow.com/questions/51261388/duplicate-tips-label-removal-importing-phylogenetic-tree-in-r-for-comparison
traitData_touse_PhyloDataCombined_tree$node.label = NULL


# Part 2: Curate data for PCA analysis (normalized data) -----------------------
phyloPCAData_numonly = traitData_touse %>% ungroup() %>% 
  dplyr::select(tipLabelMatched,GrowthForm,all_of(traitName_touse)) %>% 
  na.omit()
phyloPCAData_numonly_rownames = phyloPCAData_numonly$tipLabelMatched

# used for growth form annotation in PCA biplot
phyloPCAData_meta_rmna = phyloPCAData_numonly
phyloPCAData_numonly = phyloPCAData_numonly %>% dplyr::select(-tipLabelMatched, -GrowthForm)

# scale traits based on columns (mean=0, SD=1, same as z-transform) 
phyloPCAData_numonly[, traitName_touse] = apply(phyloPCAData_numonly[, traitName_touse], 2, scale)
rownames(phyloPCAData_numonly) = phyloPCAData_numonly_rownames

# Step 3: Perform phylo-informed PCA -------------------------------------------
phyloPCA.results = phyl.pca(traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly, mode="corr", method='lambda')
phyloPCA.results.prcomp = as.prcomp(phyloPCA.results)


# Step 4: Visualization --------------------------------------------------------
# 2-D PCA biplot: GrowthForm as grouping factor
# using ggbiplot version 0.6.2 allow custom ellipse fill color and themes
# installed on black TBk laptop
plt_phyloPCA_biplot_ax12 = ggbiplot( 
  phyloPCA.results.prcomp, choices = 1:2, 
  scale=1, # covariance biplot, scale all to unit variance (1)
  groups = phyloPCAData_meta_rmna$GrowthForm, ellipse = F, circle = T) +
  scale_color_manual(values = c("#383838", "#5CA7C7", "#D4352D", "#FBCE6A")) +
  theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())

plt_phyloPCA_contribplot_ax12 = fviz_pca_var(
  phyloPCA.results.prcomp, col.var = "contrib",
  axes=c(1,2),
  gradient.cols = c("#7998AD", "#EECA40", "#F07673")) +
  theme_classic() +
  theme(legend.direction = 'horizontal', legend.position = 'none') + labs(title = "")


plt_phyloPCA_biplot_ax23 = ggbiplot( 
  phyloPCA.results.prcomp, choices = 2:3, 
  scale=1, # covariance biplot, scale all to unit variance (1)
  groups = phyloPCAData_meta_rmna$GrowthForm, ellipse = F, circle = T) +
  scale_color_manual(values = c("#383838", "#5CA7C7", "#D4352D", "#FBCE6A")) +
  theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())

plt_phyloPCA_contribplot_ax23 = fviz_pca_var(
  phyloPCA.results.prcomp, col.var = "contrib",
  axes=c(2,3),
  gradient.cols = c("#7998AD", "#EECA40", "#F07673")) +
  theme_classic() +
  theme(legend.direction = 'horizontal', legend.position = 'none') + labs(title = "")


# Save plot after merging
plt_ax12 = plt_phyloPCA_biplot_ax12 + plt_phyloPCA_contribplot_ax12 + 
  #plot_annotation(title = "Phylo-PCA, PC1-2") +
  plot_layout(guides = "collect") & theme(legend.position='bottom')

plt_ax23 = plt_phyloPCA_biplot_ax23 + plt_phyloPCA_contribplot_ax23 + 
  #plot_annotation(title = "Phylo-PCA, PC2-3") +
  plot_layout(guides = "collect") & theme(legend.position='bottom')

# Combine the two plots side by side
plt_ax12_ax23_combined = plot_grid(plt_ax12, plt_ax23, ncol = 1)
ggsave(plot = plt_ax12_ax23_combined, filename = "plt_phyloPCA_ax123Combined_RLESLPCRDMC.pdf",
       width = 6, height = 8)
ggsave(plot = plt_ax12_ax23_combined, filename = "plt_phyloPCA_ax123Combined_RLESLPCRDMC.png",
       width = 6, height = 8, dpi=300)

# Scree plot to visualize contributions of different axes
phyloPCA.screeplot = fviz_screeplot(phyloPCA.results.prcomp, addlabels = TRUE, ylim = c(0, 40),
                                    barfill="#C47070", barcolor="#C47070", linecolor="black") +
  labs(title="RES+LES(+LPC,RDMC)") + theme_classic()
ggsave(filename = "plt_phyloPCA_screeplot-RESLESLPCRDMC.pdf", plot = phyloPCA.screeplot, width = 3, height = 3)
ggsave(filename = "plt_phyloPCA_screeplot-RESLESLPCRDMC.png", plot = phyloPCA.screeplot, width = 3, height = 3, dpi=300)


# Step 5: Compare with non-phylo PCA, grouped barplot --------------------------
screeData_compare_RLES = data.frame(
  Dimensions = c("1","2","3","4","5","6"),
  NonPhylo = c(34.4,27.6,21.5,10,6.3,0.2),
  Phylo = c(35,27,22.1,9.5,6.3,0.2))
screeData_compare_RLESRDMC = data.frame(
  Dimensions = c("1","2","3","4","5","6","7"),
  NonPhylo = c(35.5,26.6,21.3,9.3,5.4,1.8,0.1),
  Phylo = c(34.3,27.6,21.7,8.8,5.4,2,0.1))
screeData_compare_RLESLPC = data.frame(
  Dimensions = c("1","2","3","4","5","6","7"),
  NonPhylo = c(32,27.6,21.2,9,7.2,2.8,0.1),
  Phylo = c(32.7,27.2,21.1,9.1,6.8,3,0.1))
screeData_compare_RLESLPCRDMC = data.frame(
  Dimensions = c("1","2","3","4","5","6","7","8"),
  NonPhylo = c(32.5,25.7,22.5,8.8,6.3,2.5,1.6,0.1),
  Phylo = c(32.2,25,23.4,8.9,6.1,2.6,1.7,0.1))

screeData_compare_touse = screeData_compare_RLES
#screeData_compare_touse = screeData_compare_RLESLPC
#screeData_compare_touse = screeData_compare_RLESRDMC
#screeData_compare_touse = screeData_compare_RLESLPCRDMC
screeData_compare_touse = melt(screeData_compare_touse,variable.name="Category",value.name = "ExplainedVariance")

# Grouped bar plot
plt_scree_compare = ggplot(screeData_compare_touse) +
  geom_bar(aes(x=Dimensions, y=ExplainedVariance, fill=Category), color="black", # border
           stat="identity", width=0.8, position="dodge") + 
  
  geom_point(aes(x=Dimensions, y=ExplainedVariance, color=Category), size=2, position=position_dodge(width = 0.6)) +
  geom_line(aes(x=Dimensions, y=ExplainedVariance, color=Category, group = Category), linewidth=1, position=position_dodge(width = 0.6)) +
  
  scale_fill_manual(values = c("#7998AD","#C47070")) + # bar coloring
  scale_color_manual(values = c("#000080","#B22222")) + # point and line coloring
  
  geom_text(aes(x = Dimensions, y = ExplainedVariance, group = Category, label = paste0(ExplainedVariance, "%")),
            position = position_dodge(width = 0.6),size = 3.5,vjust = 0.8,hjust=-0.7,angle=90) + 
  labs(title="RES+LES", y="Variance Explained (%)") + 
  theme_classic() +
  scale_y_continuous(limits = c(0, 45), expand = c(0, 0)) + # Removes the blank space
  theme(legend.title = element_blank(), 
        legend.position = c(1, 1),
        legend.justification = c(1, 1)
        # legend.background = element_rect(fill = "white", color = "black")
        )

ggsave(filename = "plt_PCAcompare_scree-RESLES.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)
#ggsave(filename = "plt_PCAcompare_scree-RESLESLPC.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)
#ggsave(filename = "plt_PCAcompare_scree-RESLESRDMC.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)
#ggsave(filename = "plt_PCAcompare_scree-RESLESLPCRDMC.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)

