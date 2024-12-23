rm(list=ls())
# devtools::install_github('laijiangshan/phylolm.hp',build_vignettes = TRUE)
library(phylolm.hp)
library(tidyverse)
library(dplyr)

# Run phylolm.hp with test data
# Ref1: https://mp.weixin.qq.com/s/L-y6C_oSBcvpaPcqF0PvkA
# Ref2: https://mp.weixin.qq.com/s/QbFlhcGHOxZwnZzVZgV8Lg

#连续性状数据
set.seed(231)
tre <- rcoal(60) # simulate a tree
taxa <- sort(tre$tip.label) # extract taxa names from tree
b0 <- 0    
b1 <- 0.3    
b2 <- 0.5
b3 <- 0.4

x <- rTrait(n=1, phy=tre, model="lambda", parameters=list(ancestral.state=0, sigma2=15, lambda=0.9))          
x2 <- rTrait(n=1, phy=tre, model="lambda", parameters=list(ancestral.state=0, sigma2=10, lambda=0.9))  
x3 <- rTrait(n=1, phy=tre, model="lambda", parameters=list(ancestral.state=0, sigma2=13, lambda=0.9))       

y <- b0 + b1 * x + b2 * x2 + b3 * x3 + rTrait(
  n=1, phy=tre, model="lambda", parameters=list(ancestral.state=0, sigma2=5, lambda=0.9))            

dat <- data.frame(trait=y[taxa], pred=x[taxa], pred2=x2[taxa],pred3=x3[taxa])

fit <- phylolm(trait ~ pred + pred2 + pred3, data=dat, phy=tre, model="lambda")

phyloglm.hp(fit,commonality=TRUE)

iv=list(env1="pred",env2=c("pred2","pred3"))

phyloglm.hp(fit,iv)

plot(phyloglm.hp(fit,iv))

plot(phyloglm.hp(fit,iv,commonality=T),commonality=T)


#二元性状数据

set.seed(123456)

tre <- rtree(50)

x1 <- rTrait(n=1, phy=tre)  

x2 <- rTrait(n=1, phy=tre)

x3 <- rTrait(n=1, phy=tre)

X <- cbind(rep(1, 50), x1, x2, x3)

y <- rbinTrait(n=1, phy=tre, beta=c(-1, 0.9, 0.9, 0.5), alpha=1, X=X)

dat <- data.frame(trait01=y, predictor1=x1, predictor2=x2, predictor3=x3)

fit <- phyloglm(trait01 ~ predictor1 + predictor2 + predictor3, phy=tre, data=dat)

phyloglm.hp(fit)

iv=list(env1="predictor1",env2=c("predictor2","predictor3"))

phyloglm.hp(fit,iv)

plot(phyloglm.hp(fit,iv))