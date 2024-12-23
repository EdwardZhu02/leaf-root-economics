rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggsignif)  # For adding significance markers
library(broom)     # For tidying the t-test results
library(ggpubr)
library(patchwork) # for merging plots
library(plotly)

# Perform across-site comparison using community-mean trait values
load("traitDataFujian-spavg-phylo-step2.RData")

# use data frame 'traitDataIndv', the original, not species averaged or normalized data for analysis
traitDataIndv_touse = traitDataIndv # backup original df, while copy it for this script

# reshape data for plotting (box plot, long data)
traitDataIndv_commean_box = traitDataIndv_touse %>%
  dplyr::select(where(is.numeric), GrowthForm, SampleID) %>%   # Select numeric columns and SiteLoc
  pivot_longer(cols = -c(GrowthForm, SampleID), names_to = "trait", values_to = "value") %>% 
  group_by(trait) %>%
  filter(!is.na(value))

# Filter only data of interest to generate the plot
traitlist_touse = c(
  "RTD","SRL","RD","RNC","SRA","RPC","RCC","RDMC","SRR25",
  "rs25","LMA","LNC","LPC","Rdark25P","Rlight25P","Vcmax25","Asat",
  "vH","SWCleaf","SWCbranch","WD","Ks","TLP","H","DBH" # hydraulic traits
  )


# Define significant comparison pairs
signif_comparison = list(c("liana", "shrub"), c("shrub", "tree"),c("liana", "tree"))
tomerge_plot_list = list() # empty list for plot merging

for (i in 1:length(traitlist_touse)) {
  # current_traitname = traitlist_touse[i]
  traitDataIndv_commean_box_filtered = traitDataIndv_commean_box %>% 
    dplyr::filter(trait == traitlist_touse[i])
  
  # Exclude lianas in FR traits ------------------------------------------------
  if (traitlist_touse[i] %in% c("RTD","SRL","RD","RNC","SRA","RPC","RCC","SRR","SRR25","RDMC")) {
    traitDataIndv_commean_box_filtered = traitDataIndv_commean_box_filtered %>% 
      dplyr::filter(GrowthForm != "liana")
  }
  # ----------------------------------------------------------------------------

  GrowthForm_distinct = traitDataIndv_commean_box_filtered[,1] %>% dplyr::distinct()
  
  if ("liana" %in% GrowthForm_distinct$GrowthForm){
    boxplt_GrowthFormIndv = ggplot(
      traitDataIndv_commean_box_filtered, aes(x = GrowthForm, y = value, fill = GrowthForm, color = GrowthForm)) +
      geom_boxplot(alpha=0.1) +
      geom_jitter(size=0.1) +
      stat_compare_means(comparisons=signif_comparison, method="wilcox.test", label="p.signif") +
      labs(subtitle = traitlist_touse[i]) +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
    
  } else {
    boxplt_GrowthFormIndv = ggplot(
      traitDataIndv_commean_box_filtered, aes(x = GrowthForm, y = value, fill = GrowthForm, color = GrowthForm)) +
      geom_boxplot(alpha=0.1) +
      geom_jitter(size=0.1) +
      #stat_compare_means(method="t.test", label="p.signif") + # Response to Runxi's comment, 25/11-24
      stat_compare_means(method="wilcox.test", label="p.signif") +
      labs(subtitle = traitlist_touse[i]) +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  
  # add to the total plot list
  tomerge_plot_list[[i]] = boxplt_GrowthFormIndv

}

# Combine plots with patchwork
boxplt_GrowthFormCombined = tomerge_plot_list[[1]]
for (i in 2:length(traitlist_touse)) {
  boxplt_GrowthFormCombined = boxplt_GrowthFormCombined + tomerge_plot_list[[i]]
}

boxplt_GrowthFormCombined = boxplt_GrowthFormCombined + 
  # plot_annotation(title = "Comparison: growth forms") +
  plot_layout(guides = "collect", ncol=4) & 
  theme(legend.position='none')

boxplt_GrowthFormCombined
# Save config: part 01
ggsave(plot = boxplt_GrowthFormCombined, filename = "plts-growthform-comparison/boxplt_GFcomp01.pdf",
       width = 8.3, height = 8.1)
# Save config: part 02
ggsave(plot = boxplt_GrowthFormCombined, filename = "plts-growthform-comparison/boxplt_GFcomp02.pdf",
       width = 7, height = 6.8)
