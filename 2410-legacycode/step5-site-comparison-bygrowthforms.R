rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggsignif)  # For adding significance markers
library(broom)     # For tidying the t-test results
library(ggpubr)
library(plotly)
library(patchwork)

# Perform across-site comparison using community-mean trait values
load("traitDataFujian-spavg-phylo-step2.RData")

# use data frame 'traitDataIndv', the original, not species averaged or normalized data for analysis
traitDataIndv_touse = traitDataIndv # backup original df, while copy it for this script

# --------------------------------------------------------------------------
# STEP1: add SiteLoc column, and filter observations into different DFs ####
# --------------------------------------------------------------------------
traitDataIndv_touse$SiteLoc = NA
traitDataIndv_touse$SiteLoc = ifelse(startsWith(traitDataIndv_touse$SpeciesID, "SH"), "hilltop", 
                                     ifelse(startsWith(traitDataIndv_touse$SpeciesID, "WH"), "valley", traitDataIndv_touse$SiteLoc))

# traitDataIndv_shrub = traitDataIndv_touse %>% filter(GrowthForm == "shrub")
# traitDataIndv_tree = traitDataIndv_touse %>% filter(GrowthForm == "tree")

# --------------------------------------------------------------------------
# STEP2: reshape data for plotting (box plot, long data) ####
# --------------------------------------------------------------------------
traitDataIndv_byLocGF_mean_box = traitDataIndv_touse %>%
  dplyr::select(where(is.numeric), GrowthForm, SiteLoc, SampleID) %>%   # Select numeric columns and SiteLoc
  pivot_longer(cols = -c(GrowthForm, SiteLoc, SampleID), names_to = "trait", values_to = "value") %>% 
  group_by(trait) %>%
  filter(!is.na(value))


traitlist_touse = c("SRR25","rs25","Asat","Vcmax25","RCC","RNC","RPC","RDMC","LMA","LPC")
signif_comparison = list(c("hilltop", "valley"))
tomerge_plot_list_tree = list() # empty list for plot merging
tomerge_plot_list_shrub = list() # empty list for plot merging

# Exclude liana entries (mostly NA values)
for (i in 1:length(traitlist_touse)) {
  
  traitDataIndv_filtered_tree = traitDataIndv_byLocGF_mean_box %>% 
    dplyr::filter(trait == traitlist_touse[i]) %>%
    dplyr::filter(GrowthForm == "tree")
  
  traitDataIndv_filtered_shrub = traitDataIndv_byLocGF_mean_box %>% 
    dplyr::filter(trait == traitlist_touse[i]) %>%
    dplyr::filter(GrowthForm == "shrub")
  
  
  boxplt_SiteComp_tree = ggplot(
    traitDataIndv_filtered_tree, aes(x = SiteLoc, y = value)) +
    geom_boxplot(alpha=0.1, fill="#E54B4B") +
    geom_jitter(size=0.1) +
    stat_compare_means(method="t.test", label="p.signif") +
    labs(subtitle = traitlist_touse[i]) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  
  boxplt_SiteComp_shrub = ggplot(
    traitDataIndv_filtered_shrub, aes(x = SiteLoc, y = value)) +
    geom_boxplot(alpha=0.1, fill="#167C80") +
    geom_jitter(size=0.1) +
    stat_compare_means(method="t.test", label="p.signif") +
    labs(subtitle = traitlist_touse[i]) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  
  # add to the total plot list
  tomerge_plot_list_tree[[i]] = boxplt_SiteComp_tree
  tomerge_plot_list_shrub[[i]] = boxplt_SiteComp_shrub
  
}

# Combine plots with patchwork
boxplt_GrowthFormCombined_tree = tomerge_plot_list_tree[[1]]
for (i in 2:length(traitlist_touse)) {
  boxplt_GrowthFormCombined_tree = boxplt_GrowthFormCombined_tree + tomerge_plot_list_tree[[i]]
}
boxplt_GrowthFormCombined_shrub = tomerge_plot_list_shrub[[1]]
for (i in 2:length(traitlist_touse)) {
  boxplt_GrowthFormCombined_shrub = boxplt_GrowthFormCombined_shrub + tomerge_plot_list_shrub[[i]]
}

boxplt_GrowthFormCombined_tree = boxplt_GrowthFormCombined_tree + 
  plot_annotation(title = "site-comparison: tree") +
  plot_layout(guides = "collect", ncol=5) & theme(legend.position='none')

boxplt_GrowthFormCombined_shrub = boxplt_GrowthFormCombined_shrub + 
  plot_annotation(title = "site-comparison: shrub") +
  plot_layout(guides = "collect", ncol=5) & theme(legend.position='none')


ggsave(plot = boxplt_GrowthFormCombined_tree, filename = "boxplt_SiteComp-byGF-tree.pdf",
       width = 5.5, height = 7)
# Save config: part 02
ggsave(plot = boxplt_GrowthFormCombined_shrub, filename = "boxplt_SiteComp-byGF-shrub.pdf",
       width = 5.5, height = 7)
