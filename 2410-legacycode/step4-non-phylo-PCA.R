rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(rgl) # 3-D PCA plot

# ------------------------------------------------------------------------------
# Perform PCA analysis
# ------------------------------------------------------------------------------
load("traitDataFujian-fmtdata-step3.RData")

# TODO: remove LIANA entries in the data due to its minority in PCA analysis
# traitDataIndv_normalized = traitDataIndv_normalized %>% filter(GrowthForm != "liana") # (41->37 species, but only 1 contain FR trait measurements)

# # FILTER OUTLIERS (for tentative analysis, 13/10-24)
# traitDataIndv_normalized = traitDataIndv_normalized %>% 
#   filter(RTD > -2)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataIndv_normalized %>% dplyr::select(RTD,SRL,RD,RNC,SRR25,Rdark25P)
rownames(nonphyloPCAData_numonly) = traitDataIndv_normalized$speciesFullName
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% na.omit() # remove lines /w NA values

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = traitDataIndv_normalized %>% 
  dplyr::select(speciesFullName,GrowthForm,RTD,SRL,RD,RNC,SRR25,Rdark25P) %>% na.omit()

nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = T, graph = T) 

#phyloPCA.results = phyl.pca(traitPhyloDataSH_tree, phyloPCAData_numonly, method = 'lambda')
#phyloPCA.results.prcomp = as.prcomp(phyloPCA.results)

# ------------------------------------------------------------------------------
# 2-D PCA biplot: GrowthForm as grouping factor
# ------------------------------------------------------------------------------
plt_nonphyloPCA_biplot01 = fviz_pca_biplot(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$GrowthForm, 
  addEllipses=T, palette = c("#E7B800", "#a8e6cf", "#fecea8")) +
  labs(title = "") + theme_bw() + theme(legend.direction = 'horizontal', legend.position = 'top')

ggsave(plot = plt_nonphyloPCA_biplot01, filename = "plt_nonphyPCA_biplot_RESResp.pdf", width = 4.5, height = 4.7)

# ------------------------------------------------------------------------------
# Try to plot 3-D PCA graph showing PC1-PC3 (GrowthForm as grouping factor)
# Ref: https://zhuanlan.zhihu.com/p/435342487
# ------------------------------------------------------------------------------
nonphyloPCA.ind = data.frame(get_pca_ind(nonphyloPCAresult)$coord[,1:3]) %>%
  dplyr::mutate(
    speciesFullName = nonphyloPCAData_meta_rmna$speciesFullName,
    GrowthForm = nonphyloPCAData_meta_rmna$GrowthForm) # individuals
colnames(nonphyloPCA.ind)= c("PC1","PC2","PC3","speciesFullName","GrowthForm")

nonphyloPCA.var = data.frame(get_pca_var(nonphyloPCAresult)$coord[,1:3]) # variables
colnames(nonphyloPCA.var)= c("PC1","PC2","PC3")

# START! =======================================================================
# Specify colors --
color_list = c("#00adb5", "#303841", "#ff5722")
nonphyloPCA.ind$GrowthForm = factor(nonphyloPCA.ind$GrowthForm,levels = c("shrub","tree","liana"))
#color_list = c("#00adb5", "#303841")
#nonphyloPCA.ind$GrowthForm = factor(nonphyloPCA.ind$GrowthForm,levels = c("shrub","tree"))
labels.col = color_list[as.numeric(nonphyloPCA.ind$GrowthForm)]


open3d()
par3d(family="serif", cex=1.2, font=1) # set params
plot3d(nonphyloPCA.ind[,1:3],
       type="p", 
       col=labels.col,
       lwd = 1,
       size=5,
       ticktype = "detailed",
)

# Prepare arrow origin (0,0,0) DF
origin = data.frame(PC1= rep(0,nrow(nonphyloPCA.var)),
                    PC2= rep(0,nrow(nonphyloPCA.var)),
                    PC3= rep(0,nrow(nonphyloPCA.var)),
                    row.names = rownames(nonphyloPCA.var))
# Draw arrow
# Env coord *3 for better ornamental purposes
for (i in 1:nrow(nonphyloPCA.var[1:nrow(nonphyloPCA.var),])) {
  arrow3d(p0 = origin[i,], 
          p1 = nonphyloPCA.var[i,]*3, 
          type = "lines",
          width = 3,
          thickness = 0.618 * width, 
          theta = 0.45,
          s = 0.2,
          n = 2,
          col = "red",
          plot = TRUE)
}
# Add variable label
text3d(nonphyloPCA.var[1:nrow(nonphyloPCA.var),]*3, 
       texts=rownames(nonphyloPCA.var[1:nrow(nonphyloPCA.var),]),
       offset = 0.5, pos = 3, col="black", family = "serif", cex=1.5)

# Add legend (rerun after adjusting for plot size to maintain resolution)
legend3d("right", levels(nonphyloPCA.ind$GrowthForm), pch=16, col=color_list)

# After adjusting for optimal angle, run this line to save plot.
rgl.postscript("plt_nonphyloPCA_3D_RESLES.pdf", fmt = "pdf", drawText = T)
close3d()



# --------------------------------------------------------------------------
# Extra: added 14/10-24
# 2-D PCA biplot: SiteLoc as grouping factor (redo species-level average)
# --------------------------------------------------------------------------
traitDataIndv_touse = traitDataIndv # backup original df, while copy it for this script

traitDataIndv_touse$SiteLoc = NA
traitDataIndv_touse$SiteLoc = ifelse(startsWith(traitDataIndv_touse$SpeciesID, "SH"), "hilltop", 
                                     ifelse(startsWith(traitDataIndv_touse$SpeciesID, "WH"), "valley", traitDataIndv_touse$SiteLoc))

traitDataIndv_top = traitDataIndv_touse %>% filter(SiteLoc == "hilltop")
traitDataIndv_bottom = traitDataIndv_touse %>% filter(SiteLoc == "valley")

# curate individual ID for dplyr summarize
traitDataIndv_touse$SpSiteLoc = paste(traitDataIndv_touse$SpeciesID, traitDataIndv_touse$SiteLoc, sep = "_")

# calculate the species mean values for each trait, discarding NA values.
traitDataIndv_splocavg = traitDataIndv_touse %>%
  group_by(`SpSiteLoc`) %>%
  dplyr::summarize(
    GrowthForm = first(GrowthForm),
    SiteLoc = first(SiteLoc),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)) # Calculate mean for all numeric columns
  )

# function for data normalization
func_log_z_transform <- function(obj_df) {
  
  numeric_cols <- sapply(obj_df, is.numeric)
  numeric_data <- obj_df[, numeric_cols]
  
  # Log-transform each trait (signed log transformation, dealing with negative values)
  log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
  # Z-transform the log-transformed data
  z_transformed <- scale(log_transformed, center = TRUE, scale = TRUE)
  
  obj_df[, numeric_cols] <- z_transformed
  return(obj_df)        
}

traitDataIndv_splocavg_normalized = func_log_z_transform(traitDataIndv_splocavg)


# Curate data for PCA analysis
nonphyloPCAData_numonly2 = traitDataIndv_splocavg_normalized %>% dplyr::select(RTD,SRL,RD,RNC,LMA,LNC,SRR25,Rdark25P)
rownames(nonphyloPCAData_numonly2) = traitDataIndv_splocavg_normalized$SpSiteLoc
nonphyloPCAData_numonly2 = nonphyloPCAData_numonly2 %>% na.omit() # remove lines /w NA values

# used for annotation in PCA biplot
nonphyloPCAData_meta2_growthform = traitDataIndv_splocavg_normalized %>% 
  dplyr::select(SpSiteLoc,GrowthForm,RTD,SRL,RD,RNC,LMA,LNC,SRR25,Rdark25P) %>% na.omit()
nonphyloPCAData_meta2_siteloc = traitDataIndv_splocavg_normalized %>% 
  dplyr::select(SpSiteLoc,SiteLoc,RTD,SRL,RD,RNC,LMA,LNC,SRR25,Rdark25P) %>% na.omit()

nonphyloPCAData_meta2_siteloc$SiteLoc = as.factor(nonphyloPCAData_meta2_siteloc$SiteLoc)

nonphyloPCAresult2 = PCA(nonphyloPCAData_numonly2, scale.unit = T, graph = T) 

# Visualization
plt_nonphyloPCA_biplot02 = fviz_pca_biplot(
  nonphyloPCAresult2, label = "var", habillage=nonphyloPCAData_meta2_siteloc$SiteLoc, 
  addEllipses=T, palette = c("#95e8d7", "#de95ba")) +
  labs(title = "") + theme_bw() + theme(legend.direction = 'horizontal', legend.position = 'top')

ggsave(plot = plt_nonphyloPCA_biplot02, filename = "plt_nonphyPCA_biplot_RESLESResp02.pdf", width = 4.5, height = 4.7)

