rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggsignif)  # For adding significance markers
library(broom)     # For tidying the t-test results
library(ggpubr)
library(plotly)

# Perform across-site comparison using community-mean trait values
load("traitDataFujian-fmtdata-step3.RData")

# use data frame 'traitDataIndv', the original, not species averaged or normalized data for analysis
traitDataIndv_touse = traitDataIndv # backup original df, while copy it for this script


# --------------------------------------------------------------------------
# STEP1: add SiteLoc column, and filter observations into different DFs ####
# --------------------------------------------------------------------------
traitDataIndv_touse$SiteLoc = NA
traitDataIndv_touse$SiteLoc = ifelse(startsWith(traitDataIndv_touse$SpeciesID, "SH"), "hilltop", 
                  ifelse(startsWith(traitDataIndv_touse$SpeciesID, "WH"), "valley", traitDataIndv_touse$SiteLoc))

traitDataIndv_top = traitDataIndv_touse %>% filter(SiteLoc == "hilltop")
traitDataIndv_bottom = traitDataIndv_touse %>% filter(SiteLoc == "valley")

# --------------------------------------------------------------------------
# STEP2: calculate species mean trait values and SD, TOP and BOTTOM site, sep. ####
# --------------------------------------------------------------------------
traitDataIndv_spavg_withsite = traitDataIndv_touse %>%
  group_by(`speciesFullName`, `SiteLoc`) %>%
  dplyr::summarize(
    GrowthForm = first(GrowthForm),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)) # Calculate mean for all numeric columns
  )

traitDataIndv_commean_withsite = traitDataIndv_spavg_withsite %>%
  group_by(`SiteLoc`) %>%
  dplyr::summarize(
    across(where(is.numeric), list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)))
  ) %>% ungroup()


# reshape data for plotting (W -> L)
traitDataIndv_commean_long = traitDataIndv_commean_withsite %>%
  pivot_longer(cols = -SiteLoc, 
               names_to = c("trait", ".value"), 
               names_pattern = "(.+)_(mean|sd)",
               values_drop_na = TRUE)

# --------------------------------------------------------------------------
# STEP3: Perform t-tests for each trait comparing the two sites ####
# --------------------------------------------------------------------------
t_test_rawdata = traitDataIndv_touse %>%
  dplyr::select(where(is.numeric), SiteLoc, SampleID) %>%   # Select numeric columns and SiteLoc
  pivot_longer(cols = -c(SiteLoc, SampleID), names_to = "trait", values_to = "value") %>% 
  group_by(trait) %>%
  filter(!is.na(value))

t_test_results = traitDataIndv_touse %>%
  dplyr::select(where(is.numeric), SiteLoc) %>%   # Select numeric columns and SiteLoc
  pivot_longer(cols = -SiteLoc, names_to = "trait", values_to = "value") %>% 
  group_by(trait) %>% 
  do(tidy(t.test(value ~ SiteLoc, data = .))) %>%   # Perform t-test for each trait
  dplyr::select(trait, p.value) # Keep the p-value for significance

# Combine t-test results with the data for plotting
traitDataIndv_commean_long <- traitDataIndv_commean_long %>%
  left_join(t_test_results, by = "trait")


# --------------------------------------------------------------------------
# STEP4: Filter only data of interest to generate the plot
#        So that significant marks are displayed in a readable manner
# --------------------------------------------------------------------------
traitDataIndv_commean_long_BACKUP = traitDataIndv_commean_long # make backup

traitlist_touse = c("RTD","SRL","RD","RNC","SRA","RPC","RCC","RDMC","SRR25","rs25","LMA","LNC","LPC","Rdark25P","Rlight25P","Vcmax25","Asat") # Part 1

traitlist_touse = c("vH","SWCleaf","SWCbranch","WD","Ks","TLP","H","DBH") # Part 2, hydraulic traits


traitDataIndv_commean_long = traitDataIndv_commean_long %>%
  dplyr::filter(trait %in% traitlist_touse)


# --------------------------------------------------------------------------
# STEP5: data visualization ####
# --------------------------------------------------------------------------
facet_labeller_pvalue_func = function(oristr){
  
  for (idx in 1:length(oristr)) {
    pvalue_trait = traitDataIndv_commean_long %>% dplyr::filter(trait == oristr[idx])
    # pvalue_fmt = format(pvalue_trait$p.value[1], digits = 3, scientific = F)
    pvalue_fmt = signif(pvalue_trait$p.value[1], digits = 3)
    signif_indicator = ifelse(pvalue_fmt < 0.01, "**", ifelse(pvalue_fmt < 0.05, "*", "NS"))
    
    output_label = paste0(oristr[idx], "(", signif_indicator, ", p=", pvalue_fmt, ")")
    oristr[idx] = output_label
  }
  return(oristr)
}

barplt_sitecomp01 = ggplot(traitDataIndv_commean_long, aes(x = SiteLoc, y = mean, fill = SiteLoc)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7, alpha=0.9) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                position = position_dodge(width = 0.9), 
                width = 0.2) +
  facet_wrap(~trait, scales = "free_y", labeller = labeller(trait=facet_labeller_pvalue_func)) +
  theme_minimal() +
  labs(x = "Traits", y = "Mean Trait Value", fill = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("#EB8317", "#10375C"))

# Save config: part 01
# ggsave(plot = barplt_sitecomp01, filename = "plts-site-comparison/barplt_sitecomp01.pdf",
#        width = 10, height = 6)

# Save config: part 02
ggsave(plot = barplt_sitecomp01, filename = "plts-site-comparison/barplt_sitecomp02.pdf",
      width = 6.7, height = 4)
  

# ------------------------------------------------------------------------------
# AUX plotting
# ------------------------------------------------------------------------------
# Plot individual trait pairs, finding outliers
traitPair_IndvComp = t_test_rawdata %>% dplyr::filter(trait == "SRR25")
boxplt_traitpaircomp01 = ggplot(traitPair_IndvComp, aes(x = SiteLoc, y = value, fill = SiteLoc)) +
  geom_boxplot() +
  theme_minimal()
  #geom_signif(comparisons = list(c("top", "bottom")), map_signif_level = TRUE, y_position = 0.5)

ggplotly(boxplt_traitpaircomp01)
# # +
#   
#   # Add significance markers using p-value thresholds
#   geom_signif(comparisons = list(c("top", "bottom")), 
#               aes(annotations = ifelse(p.value < 0.05, "*", "")),
#               y_position = max(traitDataIndv_commean_long$mean + traitDataIndv_commean_long$sd, na.rm = TRUE) + 0.1,
#               tip_length = 0.03,
#               vjust = 0.5)
