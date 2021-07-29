#Fit models to traits, RGR, and their interactions.
rm(list=ls())
library(boot) #for logit transforms.

#set output path.----
output.path <- 'data/surv_models_fitted.rds'

#load data.----
d <- readRDS('data/phase1_for_analysis.rds')
dat <- d$all.dat[,c('spp','phase1_surv_yr2','phase1_RGR_yr2','bg.PC1','bg.PC2','bg.PC3','ag.PC1','ag.PC2','ag.PC3','all.PC1')]
colnames(dat)[2:3] <- c('survival','RGR')
#get logit transformed survival.
dat$surv.logit <- logit(dat$survival)
dat <- dat[complete.cases(dat),]

#fit regression models.----
m.rgr    <- lm(surv.logit ~ RGR        , data = dat)
m.ag     <- lm(surv.logit ~ ag.PC1     , data = dat)
m.bg     <- lm(surv.logit ~ bg.PC1     , data = dat)
m.wp     <- lm(surv.logit ~ all.PC1    , data = dat)
m.ag_rgr <- lm(surv.logit ~ RGR*ag.PC1 , data = dat)
m.bg_rgr <- lm(surv.logit ~ RGR*bg.PC1 , data = dat)
m.wp_rgr <- lm(surv.logit ~ RGR*all.PC1, data = dat)


#save model list.
mod.list <- list(m.ag, m.bg, m.wp, m.rgr, m.ag_rgr, m.bg_rgr, m.wp_rgr)
mod.lab  <-    c('m.ag','m.bg','m.wp','m.rgr','m.rgr_ag','m.rgr_bg','m.rgr_wp')
names(mod.list) <- mod.lab

#get rsq and rsq.adj dataframe.----
rsq     <- list()
rsq.adj <- list()
for(i in 1:length(mod.list)){
  rsq    [[i]] <- summary(mod.list[[i]])$r.squared
  rsq.adj[[i]] <- summary(mod.list[[i]])$adj.r.squared
}
rsq     <- unlist(rsq)
rsq.adj <- unlist(rsq.adj)
rsq.sum <- data.frame(mod.lab, rsq, rsq.adj)

#save output.----
output <- list(mod.list, rsq.sum, dat)
names(output) <- c('models','rsq.sum','data')
saveRDS(output, output.path)
