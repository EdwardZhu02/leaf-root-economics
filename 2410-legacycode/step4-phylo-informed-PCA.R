rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(ape) # general phylogenetic analysis
library(phytools) # phyl.pca
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(ggplotify)
library(ggrepel)
library(ggbiplot) # devtools::install_github("vqv/ggbiplot")
library(rgl) # 3-D PCA plot

# ------------------------------------------------------------------------------
# Run phylogenetically-informed PCA using the constructed phylotree
# ------------------------------------------------------------------------------
# Ref: https://rdrr.io/cran/phytools/man/phyl.pca.html
par(mfrow = c(1, 1)) # set plotting param

load("traitDataFujian-fmtdata-step3.RData")

# `traitPhyloDataSH_data` was constructed using the log-transformed, averaged values
# so there is no need to transform again.

# Curate data for PCA analysis (normalized data)
phyloPCAData_numonly = traitPhyloDataSH_data %>% dplyr::select(RTD,SRL,RD,RNC,LMA,LNC)

rownames(phyloPCAData_numonly) = traitPhyloDataSH_data$tiplabel
phyloPCAData_numonly = phyloPCAData_numonly %>% na.omit() # remove lines /w NA values

# used for growth form annotation in PCA biplot (normalized data)
phyloPCAData_meta_rmna = traitPhyloDataSH_data %>% 
  dplyr::select(speciesFullName,GrowthForm,RTD,SRL,RD,RNC,LMA,LNC) %>% na.omit()


phyloPCA.results = phyl.pca(traitPhyloDataSH_tree, phyloPCAData_numonly, method = 'lambda')
phyloPCA.results.prcomp = as.prcomp(phyloPCA.results)


# ------------------------------------------------------------------------------
# Scree plot to visualize contributions of different axes
# ------------------------------------------------------------------------------
phyloPCA.screeplot = fviz_screeplot(phyloPCA.results.prcomp, addlabels = TRUE, ylim = c(0, 80)) +
  labs(title="phylo-PCA, RES+Resp (FR+L)")
# Save if needed.
ggsave(filename = "plt_phyloPCA_screeplot-RESResp.pdf", plot = phyloPCA.screeplot, width = 4, height = 4)

# ------------------------------------------------------------------------------
# Plot a default, 2-D PCA graph with PC1-PC2
# ------------------------------------------------------------------------------
phyloPCA.biplot = ggbiplot(phyloPCA.results.prcomp, obs.scale = 1, var.scale = 1,
         groups = phyloPCAData_meta_rmna$GrowthForm, ellipse = T, circle = F, ellipse.fill = F) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + theme(legend.direction = 'horizontal', legend.position = 'top')
# Save if needed.
ggsave(filename = "plt_phyloPCA_biplot-RESLESResp.pdf", plot = phyloPCA.biplot, width = 4.5, height = 3.8)

# ------------------------------------------------------------------------------
# Try to plot 3-D PCA graph showing PC1-PC3
# Ref: https://zhuanlan.zhihu.com/p/435342487
# ------------------------------------------------------------------------------
phyloPCA.prcomp.ind = data.frame(phyloPCA.results.prcomp$x[,1:3]) %>%
  dplyr::mutate(
    speciesFullName = phyloPCAData_meta_rmna$speciesFullName,
    GrowthForm = phyloPCAData_meta_rmna$GrowthForm) # individuals

phyloPCA.prcomp.var = data.frame(phyloPCA.results.prcomp$rotation[,1:3]) # variables
#phyloPCA.prcomp.var = data.frame(get_pca_var(phyloPCA.results.prcomp)$coord[,1:3]) # variables
#colnames(phyloPCA.prcomp.var)= c("PC1","PC2","PC3")

# START! =======================================================================
# Specify colors --
color_list = c("#00adb5", "#303841", "#ff5722")
phyloPCA.prcomp.ind$GrowthForm = factor(phyloPCA.prcomp.ind$GrowthForm,levels = c("shrub","tree","liana")) # 设置分组的因子水平
labels.col = color_list[as.numeric(phyloPCA.prcomp.ind$GrowthForm)]


open3d()
par3d(family="serif", cex=1.2, font=1) # set params
plot3d(phyloPCA.prcomp.ind[,1:3],
       type="p", 
       col=labels.col,
       lwd = 1,
       size=5,
       ticktype = "detailed",
)

# Prepare arrow origin (0,0,0) DF
origin = data.frame(PC1= rep(0,nrow(phyloPCA.prcomp.var)),
                    PC2= rep(0,nrow(phyloPCA.prcomp.var)),
                    PC3= rep(0,nrow(phyloPCA.prcomp.var)),
                    row.names = rownames(phyloPCA.prcomp.var))
# Draw arrow
# Env coord *3 for better ornamental purposes
for (i in 1:nrow(phyloPCA.prcomp.var[1:nrow(phyloPCA.prcomp.var),])) {
  arrow3d(p0 = origin[i,], 
          p1 = phyloPCA.prcomp.var[i,]*3, 
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
text3d(phyloPCA.prcomp.var[1:nrow(phyloPCA.prcomp.var),]*3, 
       texts=rownames(phyloPCA.prcomp.var[1:nrow(phyloPCA.prcomp.var),]),
       offset = 0.5, pos = 3, col="black", family = "serif", cex=1.5)

# Add legend (rerun after adjusting for plot size to maintain resolution)
legend3d("right", levels(phyloPCA.prcomp.ind$GrowthForm), pch=16, col=color_list)

# After adjusting for optimal angle, run this line to save plot.
rgl.postscript("plt_phyloPCA_3D_RESLES.pdf", fmt = "pdf", drawText = T)
close3d()
