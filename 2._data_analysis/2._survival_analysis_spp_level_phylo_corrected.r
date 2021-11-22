#Fit models to traits, RGR, and their interactions.
rm(list=ls())
library(boot) #for logit transforms.
library(MCMCglmm)
library(phytools)
library(caper)
library(betareg)
library(performance)
library(dplyr)

#set output path.----
output.path <- 'data/phylo_surv_models_fitted.rds'

#load data.----
d <- readRDS('data/phase1_for_analysis.rds')
phy <- read.tree('data/PRTSBtree.tre')
dat <- d$all.dat[,c('spp','phase1_surv_yr2','phase1_RGR_yr2','bg.PC1','bg.PC2','bg.PC3','ag.PC1','ag.PC2','ag.PC3','all.PC1','ag.sub.PC1','ag.sub.PC2','ag.sub.PC3')]
colnames(dat)[2:3] <- c('survival','RGR')

#get logit transformed survival.
dat$surv.logit <- logit(dat$survival)
dat <- dat[complete.cases(dat),]
dat$spp <- recode_factor(dat$spp, SAMSAM = "ALBSAM")

#clean up the tree.
phy$node.label <- NULL

#phylo-analysis priors.
phy <- force.ultrametric(phy)
rnd <- inverseA(phy)$Ainv
#Will Pearse- priors on residual phylogenetic autocorrelation terms.
priors <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))

#clean up species names so we can actually match phylogeny to dataframe.
spp <- as.character(dat$spp)
spp <- spp[order(spp)]
species <- phy$tip.label[order(phy$tip.label)]
ref <- data.frame(spp,species)
dat <- merge(dat,ref, all.x=T)
dat <- dat[complete.cases(dat),]


#fit regression models w/o phylogeny.----
m.rgr    <- lm(surv.logit ~ RGR        , data = dat)
m.ag     <- lm(surv.logit ~ ag.PC1     , data = dat)
m.bg     <- lm(surv.logit ~ bg.PC1     , data = dat)
m.wp     <- lm(surv.logit ~ all.PC1    , data = dat)
m.ag_rgr <- lm(surv.logit ~ RGR*ag.PC1 , data = dat)
m.bg_rgr <- lm(surv.logit ~ RGR*bg.PC1 , data = dat)
m.wp_rgr <- lm(surv.logit ~ RGR*all.PC1, data = dat)

#save model list.
mod.list.null <- list(m.ag, m.bg, m.wp, m.rgr, m.ag_rgr, m.bg_rgr, m.wp_rgr)
mod.lab  <-    c('m.ag','m.bg','m.wp','m.rgr','m.rgr_ag','m.rgr_bg','m.rgr_wp')
names(mod.list.null) <- mod.lab

#get rsq and rsq.adj dataframe.
rsq     <- list()
rsq.adj <- list()
for(i in 1:length(mod.list.null)){
  rsq    [[i]] <- summary(mod.list.null[[i]])$r.squared
  rsq.adj[[i]] <- summary(mod.list.null[[i]])$adj.r.squared
}
rsq     <- unlist(rsq)
rsq.adj <- unlist(rsq.adj)
rsq.sum <- data.frame(mod.lab, rsq, rsq.adj)

#wrap output.
output.null <- list(mod.list.null, rsq.sum, dat)
names(output.null) <- c('models','rsq.sum','data')


#fit regression models WITH phylogeny.----
m.rgr        <- MCMCglmm(surv.logit ~ RGR            , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.ag         <- MCMCglmm(surv.logit ~  ag.PC1        , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.ag.sub     <- MCMCglmm(surv.logit ~  ag.sub.PC1    , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.bg         <- MCMCglmm(surv.logit ~  bg.PC1        , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.wp         <- MCMCglmm(surv.logit ~ all.PC1        , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.ag_rgr     <- MCMCglmm(surv.logit ~ RGR* ag.PC1    , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.bg_rgr     <- MCMCglmm(surv.logit ~ RGR* bg.PC1    , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.wp_rgr     <- MCMCglmm(surv.logit ~ RGR*all.PC1    , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.ag.sub_rgr <- MCMCglmm(surv.logit ~ RGR* ag.sub.PC1, random=~species, ginverse=list(species=rnd), data=dat, prior = priors)

#save model list.
mod.list.phylo <- list(m.ag, m.bg, m.wp, m.ag.sub, m.rgr, m.ag_rgr, m.bg_rgr, m.wp_rgr, m.ag.sub_rgr)
mod.lab  <-    c('m.ag','m.bg','m.wp', 'm.ag.sub','m.rgr','m.rgr_ag','m.rgr_bg','m.rgr_wp','m.ag.sub_rgr')
names(mod.list.phylo) <- mod.lab

#add fitted values and residuals to model objects.
for(i in 1:length(mod.list.phylo)){
  mod.list.phylo[[i]]$fitted.values <- predict(mod.list.phylo[[i]])
  mod.list.phylo[[i]]$observed      <- dat$surv.logit
  mod.list.phylo[[i]]$residuals     <- mod.list.phylo[[i]]$observed - mod.list.phylo[[i]]$fitted.values
}

#get rsq and rsq.adj dataframe.
rsq     <- list()
rsq.adj <- list()
for(i in 1:length(mod.list.phylo)){
  y <- dat[complete.cases(dat),]$surv.logit
  model <- mod.list.phylo[[i]]
  x <- model$fitted.values
  n.pred <- model$Fixed$nfl - 1
  fit <- lm(y~x)
  rsq    [[i]] <- summary(fit)$r.squared
  rsq.adj.est <- 1 - (1-rsq[[i]])*(nrow(dat) - 1) / (nrow(dat) - n.pred - 1)
  rsq.adj[[i]] <- rsq.adj.est
}
rsq     <- unlist(rsq)
rsq.adj <- unlist(rsq.adj)
rsq.sum <- data.frame(mod.lab, rsq, rsq.adj)

#wrap output.
output.phylo <- list(mod.list.phylo, rsq.sum, dat)
names(output.phylo) <- c('models','rsq.sum','data')

#save output.----
output <- list(output.null, output.phylo)
names(output) <- c('raw.models','phylo.models')
saveRDS(output, output.path)

##----

# compare results of beta regression to MCMCglmm (Reviewer 3 R2 response)

# fit regression models WITH beta regression (no phylogeny) ----
# convergence issues w/ species level random effect so removed 
m.rgr.beta       <- betareg(survival ~ RGR, data=dat)
m.ag.beta        <- betareg(survival ~  ag.PC1, data=dat)
m.ag.sub.beta     <- betareg(survival ~  ag.sub.PC1, data=dat)
m.bg.beta        <- betareg(survival ~  bg.PC1, data=dat)
m.wp.beta        <- betareg(survival ~ all.PC1, data=dat)
m.ag_rgr.beta     <- betareg(survival ~ RGR* ag.PC1, data=dat)
m.bg_rgr.beta   <- betareg(survival ~ RGR* bg.PC1, data=dat)
m.wp_rgr.beta     <- betareg(survival ~ RGR*all.PC1, data=dat)
m.ag.sub_rgr.beta <- betareg(survival ~ RGR* ag.sub.PC1, data=dat)

#save model list.
mod.list.beta<- list(m.ag.beta, m.bg.beta, m.wp.beta, m.ag.sub.beta, m.rgr.beta, m.ag_rgr.beta, 
                     m.bg_rgr.beta, m.wp_rgr.beta, m.ag.sub_rgr.beta)
mod.lab.beta <-    c('m.ag','m.bg','m.wp', 'm.ag.sub','m.rgr','m.rgr_ag','m.rgr_bg','m.rgr_wp','m.ag.sub_rgr')
names(mod.list.beta) <- mod.lab

#add fitted values and residuals to model objects.
for(i in 1:length(mod.list.beta)){
  mod.list.beta[[i]]$fitted.values <- predict(mod.list.beta[[i]])
  mod.list.beta[[i]]$observed      <- dat$survival
  mod.list.beta[[i]]$residuals     <- mod.list.beta[[i]]$observed - mod.list.beta[[i]]$fitted.values
}

#get rsq and rsq.adj dataframe.
rsq.pseudo     <- list()
for(i in 1:length(mod.list.beta)){
  rsq.pseudo    [[i]] <- summary(mod.list.beta[[i]])$pseudo.r.squared
}
rsq.pseudo     <- unlist(rsq.pseudo)
rsq.sum <- data.frame(mod.lab.beta, rsq.pseudo)

#wrap output.
output.beta <- list(mod.list.beta, rsq.sum, dat)
names(output.beta) <- c('models','rsq.sum','data')

# compare R2
output.beta$rsq.sum
output.phylo$rsq.sum

# plot fitted vs. fitted (phylogenetic vs. beta regression)m - best model 
plot(scale(mod.list.beta$m.rgr_bg$fitted.values), scale(mod.list.phylo$m.rgr_bg$fitted.values), 
     xlab = "beta regression (scaled fitted values)", 
     ylab = "phylogenetic regression (scaled fitted values)", 
     main = "best model (Figure 1) - scaled fitted values")
abline(0, 1)


plot(scale(mod.list.beta$m.rgr_bg$residuals), scale(mod.list.phylo$m.rgr_bg$residuals),
     xlab = "beta regression (scaled residuals)", 
     ylab = "phylogenetic regression (scaled residuals)", 
     main = "best model (Figure 1) - scaled residuals")
abline(0, 1)

# plot residuals



