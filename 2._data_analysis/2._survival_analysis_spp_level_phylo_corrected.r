#Fit models to traits, RGR, and their interactions.
rm(list=ls())
library(boot) #for logit transforms.
library(MCMCglmm)
library(phytools)
library(caper)

#set output path.----
output.path <- 'data/phylo_surv_models_fitted.rds'

#load data.----
d <- readRDS('data/phase1_for_analysis.rds')
phy <- read.tree('data/PRTSBtree.tre')
dat <- d$all.dat[,c('spp','phase1_surv_yr2','phase1_RGR_yr2','bg.PC1','bg.PC2','bg.PC3','ag.PC1','ag.PC2','ag.PC3','all.PC1')]
colnames(dat)[2:3] <- c('survival','RGR')
#get logit transformed survival.
dat$surv.logit <- logit(dat$survival)
dat <- dat[complete.cases(dat),]

#clean up the tree.
phy$node.label <- NULL
phy <- force.ultrametric(phy)
rnd <- inverseA(phy)$Ainv
#Will Pearse- priors on residual phylogenetic autocorrelation terms.
priors <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))


#clean up species names so we can actually match phylogeny to dataframe.
spp <- as.character(dat$spp)
spp <- ifelse(spp == 'SAMSAM','ALBSAM',spp)
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
m.rgr     <- MCMCglmm(surv.logit ~ RGR        , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.ag      <- MCMCglmm(surv.logit ~  ag.PC1    , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.bg      <- MCMCglmm(surv.logit ~  bg.PC1    , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.wp      <- MCMCglmm(surv.logit ~ all.PC1    , random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.ag_rgr  <- MCMCglmm(surv.logit ~ RGR* ag.PC1, random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.bg_rgr  <- MCMCglmm(surv.logit ~ RGR* bg.PC1, random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
m.wp_rgr  <- MCMCglmm(surv.logit ~ RGR*all.PC1, random=~species, ginverse=list(species=rnd), data=dat, prior = priors)

#save model list.
mod.list.phylo <- list(m.ag, m.bg, m.wp, m.rgr, m.ag_rgr, m.bg_rgr, m.wp_rgr)
mod.lab  <-    c('m.ag','m.bg','m.wp','m.rgr','m.rgr_ag','m.rgr_bg','m.rgr_wp')
names(mod.list.phylo) <- mod.lab

#add fitted values to model objects.
for(i in 1:length(mod.list.phylo)){
  mod.list.phylo[[i]]$fitted.values <- predict(mod.list.phylo[[i]])
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
