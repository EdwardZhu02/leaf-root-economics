# ****************************************************************
# TODO: Modified for use with individual-level analysis workflow
# in response to HW's advice, 2025.2
#
# Keeping the full 87 coupled aboveground and belowground sampled
# individuals for analysis
# ****************************************************************

rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(patchwork) # plot merging
library(cowplot) # plot merging

load("2502-indv-level-code/traitDataFujian-Ind-step1.RData")

### SET DATA TO USE ### 
# TODO: perform scaling later (line 36, 27/12-24)
traitDataPCA_touse = traitDataIndv_SelectedTraits_log # perform scaling later (line 36, 27/12-24)

### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC")

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ 
# because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",SpeciesFullName), SiteID, SampleID, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier_spgf,SiteID,all_of(traitName_touse)) %>% 
  dplyr::mutate(SiteID = as.factor(SiteID)) %>%
  na.omit()

# TODO: Added 27/12-24, concordant with angle calculation
# scale traits based on columns (mean=0, SD=1, same as z-transform) 
nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)
# fix rownames
rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = nonphyloPCAData_numonly
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% dplyr::select(-PCAIdentifier_spgf, -SiteID)

#-----------------------------------------------------------
# Perform PCA using FactoMineR (original result)
nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = TRUE, graph = FALSE)

#-----------------------------------------------------------
# Varimax Rotation Procedure
# Apply varimax rotation on the variable loadings from the PCA result
varimax_result <- varimax(nonphyloPCAresult$var$coord)
rotated_loadings <- as.matrix(varimax_result$loadings)

# Compute rotated individual scores by applying the rotation matrix
rotated_scores <- nonphyloPCAresult$ind$coord %*% varimax_result$rotmat

# Create a new PCA object by updating the original PCA result with rotated values.
nonphyloPCAresult_varimax_rotated <- nonphyloPCAresult
nonphyloPCAresult_varimax_rotated$var$coord <- rotated_loadings
nonphyloPCAresult_varimax_rotated$ind$coord <- rotated_scores  

# ------------------------------------------------------------------------------
# Scree plot: contributions of axes
# ------------------------------------------------------------------------------
nonphyloPCA.screeplot = fviz_screeplot(
  nonphyloPCAresult, addlabels = TRUE, ylim = c(0, 40),
  barfill="#7998AD", barcolor="#7998AD", linecolor="black") +
  labs(title="R-LES") + theme_classic()

nonphyloPCA.varimax.screeplot = fviz_screeplot(
  nonphyloPCAresult_varimax_rotated, addlabels = TRUE, ylim = c(0, 40),
  barfill="#7998AD", barcolor="#7998AD", linecolor="black") +
  labs(title="R-LES, rotated") + theme_classic()


#ggsave(filename = "2502-indv-level-code/Ind_nonphyPCA_scree-RLES.pdf", plot = nonphyloPCA.screeplot, width = 3, height = 3)
#ggsave(filename = "2502-indv-level-code/Ind_nonphyPCA_scree-RLES_varimax.pdf", plot = nonphyloPCA.varimax.screeplot, width = 3, height = 3) # no difference to the original one

# ------------------------------------------------------------------------------
# 2-D PCA biplot: SiteID as grouping factor
# ------------------------------------------------------------------------------
# 
# TODO: Added 27/12-24
# SET GLOBAL LIMIT FOR contribution levels ('contrib') FOR LEGEND MERGING
contrib_limits = c(5, 30)
#
# 
plt_nonphyloPCA_biplot_ax12 = fviz_pca_biplot(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#94c6cd", "#cb9475")) + # hilltop,valley
  annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())


plt_nonphyloPCA_biplot_ax23 = fviz_pca_biplot(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID, 
  axes=c(2,3),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#94c6cd", "#cb9475")) + # hilltop, valley
  annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())


plt_nonphyloPCA_biplot_ax12.varimax = fviz_pca_biplot(
  nonphyloPCAresult_varimax_rotated, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#94c6cd", "#cb9475")) + # hilltop,valley
  annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())


plt_nonphyloPCA_biplot_ax23.varimax = fviz_pca_biplot(
  nonphyloPCAresult_varimax_rotated, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID, 
  axes=c(2,3),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#94c6cd", "#cb9475")) + # hilltop, valley
  annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())


# Save plot after merging ------------------------------------------------------
plt_bi_ax123 = plt_nonphyloPCA_biplot_ax12 / plt_nonphyloPCA_biplot_ax23 + 
  #plot_annotation(title = "PC2-3, biplot and contribution plot") +
  plot_layout(guides = "collect") & theme(legend.position='bottom')

# ggsave(plot = plt_bi_ax123, filename = "2502-indv-level-code/IndSite_nonphyPCA_ax123_RLESRDMCLPC.pdf",
#        width = 2.5, height = 5.5)

plt_bi_ax123.varimax = plt_nonphyloPCA_biplot_ax12.varimax / plt_nonphyloPCA_biplot_ax23.varimax + 
  #plot_annotation(title = "PC2-3, biplot and contribution plot") +
  plot_layout(guides = "collect") & theme(legend.position='bottom')

ggsave(plot = plt_bi_ax123.varimax, filename = "2502-indv-level-code/IndSite_nonphyPCA_ax123_RLES.varimax.pdf",
       width = 2.5, height = 5.5)
