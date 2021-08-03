#Merge phase 1 growth, survival and traits.
rm(list=ls())
source('functions/crib_fun.r')

#set output path.----
output.path <- 'data/phase1_for_analysis.rds'

#load data.----
wp <- read.csv('data/phase1_whole_plant_traits.csv')  #'whole plant' traits
lf <- read.csv('data/phase1_leaf_traits.csv')         #leaf traits
lf$leafDens_g_cm3 <- as.numeric(lf$leafDens_g_cm3)
rt <- read.csv('data/phase1_root_traits.csv')         #root traits.
ps <- read.csv('data/VertTraitCoord_PsDataFinal.csv') #photosynthetic traits.
#growth survival at species level.
gs <- read.csv('data/VertTraitCoord_survivalGrowthTraits_speciesLevel.csv')

#Clean and average wp traits.----
wp$cn.foliar <- wp$C_ug / wp$N_ug

#subset to q-traits / responses of interest.
wp.sub <- wp[,c('spp','crownRadius_cm','rootLateral_cm','rootDepth_cm',
                'SMF','LMF','RMF','area_allLeaves_cm2','perimeter_allLeaves_cm',
                'cn.foliar','d13C','lwpDelta_mean')]
wp.avg <- aggregate(. ~ spp, data = wp.sub, FUN = mean, na.rm=TRUE, na.action=NULL)
wp.avg[is.na(wp.avg)] <- NA #convert NaN entries to NA.

#Clean and average foliar traits.----
#make unique plant IDs.
lf$plantID <- paste0(lf$spp,"_",lf$plantID)

#aggregate by plant.
lf.sub <- lf[,!(colnames(lf) %in% c('experiment','soilType','spp','leafRep'))]
lf.sub <- aggregate(. ~ plantID, data = lf.sub, FUN=mean, na.rm=TRUE, na.action=NULL)

#Drop SDMC and LDMC, better captured by SMF and LMF.
lf.sub$SDMC_mg_g <- NULL
lf.sub$LDMC_mg_g <- NULL
#drop leaf tougness because one of our species had leaves too small to measure it.
lf.sub$toughness_g <- NULL

#Drop some other aboveground traits that are correlated with other traits and hard to justify. Mostly image analysis auto-output.----
lf.sub$perimeter_wholeLeaf_cm            <- NULL
lf.sub$minorAxisLength_largestLeaflet_cm <- NULL
lf.sub$majorAxisLength_largestLeaflet_cm <- NULL
lf.sub$elongationRatio_largestLeaflet_dimless <- NULL
lf.sub$shapeFactor_largestLeaflet_dimless <- NULL
lf.sub$feretDiameterRatio_largestLeaflet_dimless <- NULL
lf.sub$area_largestLeaflet_cm2 <- NULL
lf.sub$perimeter_largestLeaflet_cm <- NULL

#aggregate by aboveground traits by species.
lf.sub$spp <- substr(lf.sub$plantID, 1, 6)
lf.sub$plantID <- NULL
lf.avg <- aggregate(. ~ spp, data = lf.sub, FUN=mean, na.rm=TRUE, na.action=NULL)

#clean and average root traits.----
#make unique plant and leaf IDs.
rt$plantID <- paste0(rt$spp,"_",rt$plantID)

#aggregate by plant.
rt.sub <- rt[,!(colnames(rt) %in% c('experiment','soilType','spp','imageID'))]
rt.sub <- aggregate(. ~ plantID, data = rt.sub, FUN=mean, na.rm=TRUE, na.action=NULL)

#Drop kimura root length, almost perfectly correlated with SRL.
#Drip RDMC, its better capture as RMF.
rt.sub$SRL_kimura_m_g <- NULL
rt.sub$RDMC_lessThan2mm_mg_g <- NULL

#aggregate by species.
rt.sub$spp <- substr(rt.sub$plantID, 1, 6)
rt.sub$plantID <- NULL
rt.avg <- aggregate(. ~ spp, data = rt.sub, FUN=mean, na.rm=TRUE, na.action=NULL)

#Clean and average photosynthetic (ps) traits.----
ps$plantID <- paste0(ps$spp,'_',ps$plantID)
ps.sub <- ps[,c('plantID','E','A','WUE')]
ps.sub <- aggregate(. ~ plantID, data = ps.sub, FUN=mean, na.rm=TRUE, na.action=NULL)
ps.sub$spp <- substr(ps.sub$plantID, 1, 6)
ps.sub$plantID <- NULL
ps.avg <- aggregate(. ~ spp, data = ps.sub, FUN=mean, na.rm=TRUE, na.action=NULL)
colnames(ps.avg) <- c('spp','E','Amax','WUE')

#re-assign wp traits to above (foliar) and belowground (root) trait dataframes.----
lf.avg <- cbind(lf.avg, wp.avg[,c('crownRadius_cm','SMF','LMF','area_allLeaves_cm2','perimeter_allLeaves_cm','cn.foliar','d13C')])
rt.avg <- cbind(rt.avg, wp.avg[,c('rootLateral_cm','rootDepth_cm','RMF')])


#clean growth-survival, other wp traits.----
gs.sub <- gs[,c('species','spp','leafHabit','leafHabitSimple','leafType','dispersalType',
                'phase1_surv_yr1','phase1_surv_yr2','phase1_RGR_yr1','phase1_RGR_yr2',
                'seedMass','woodDensity','deltaLWP')]
#subset to species actually in phase 1.
gs.sub <- gs.sub[gs.sub$spp %in% rt.avg$spp,]

#Make dispersal type into whether or not it is wind dispersed.
gs.sub$dispersal.wind <- ifelse(gs.sub$dispersalType == 'wind',1, 0)

#make leaf habit simple 0-1 for whetehr or not deciduous.
gs.sub$deciduous <- ifelse(gs.sub$leafHabitSimple == 'deciduous', 1, 0)
gs.sub$leaf.compound <- ifelse(gs.sub$leafType == 'compound',1, 0)

#merge more aboveground traits into lf.avg.
lf.avg$spp <- as.character(lf.avg$spp)
lf.avg <- cbind(lf.avg, gs.sub[,c('deciduous','leaf.compound','deltaLWP',
                                  'dispersal.wind','seedMass','woodDensity')])

#grab growth and survival data of interest, name as gs.avg.
gs.avg <- gs.sub[,c('spp','phase1_surv_yr2','phase1_RGR_yr2')]

#Convert survival rates to cribari-neto transformed proportions.
gs.avg$phase1_surv_yr2 <- crib_fun(gs.avg$phase1_surv_yr2 / 100)

#drop photosyntehtic traits into the lf.avg frame.----
ps.avg <- ps.avg[ps.avg$spp %in% lf.avg$spp,]
ps.avg <- ps.avg[order(match(ps.avg$spp, lf.avg$spp)),]
lf.avg <- merge(lf.avg, ps.avg)

#Drop aboveground traits for a reviewer.----
lf.sub <- lf.avg
lf.sub$E <- NULL
lf.sub$d13C <- NULL
lf.sub$perimeter_allLeaves_cm <- NULL
lf.sub$leaf.compound <- NULL
lf.sub$woodDensity <- NULL
lf.sub$WUE <- NULL
lf.sub$petioleLength_mm <- NULL
lf.sub$LMF <- NULL

#Generate PC1-3 vectors for above and belowground traits.----
#only use species we have RGR for.
gs.sub <- as.character(gs.avg[complete.cases(gs.avg),]$spp)

#belowground PCA
rt.pca <- rt.avg
rt.pca <- rt.pca[rt.pca$spp %in% gs.sub,]
lab <- rt.pca$spp
rt.pca$spp <- NULL
rownames(rt.pca) <- lab
colnames(rt.pca) <- c('root diameter','SRL','RTD','root lateral extent','root depth','RMF')
rt.pca.mod <- prcomp(rt.pca, center = T, scale = T)
rt.pca.out <- data.frame(rt.pca.mod$x[,1:3])
colnames(rt.pca.out) <- c('bg.PC1','bg.PC2','bg.PC3')
rownames(rt.pca.out) <- lab
rt.pca.out$spp <- lab

#aboveground PCA
lf.pca <- lf.avg
lf.pca <- lf.pca[lf.pca$spp %in% gs.sub,]
lab <- lf.pca$spp
rownames(lf.pca) <- lab
colnames(lf.pca) <- c('spp','petiole length', 'thickness', 'leaf area', 'SLA', 'leaf density', 'crown radius', 'SMF', 'LMF',
                     'total leaf area', 'total leaf perimeter', 'C:N', 'd13C', 'leaf habit', 'leaf compoundness', 'deltaLWP',
                    'dispersal mode', 'seed mass', 'wood density', 'E', 'Amax', 'WUE')
# try to fix delta 13 C symbol
#colnames(lf.pca) <- c('spp','petiole length', 'thickness', 'leaf area', 'SLA', 'leaf density', 'crown radius', 'SMF', 'LMF',
 #                     'total leaf area', 'total leaf perimeter', 'C:N', expression(paste(delta^{13}, "C")), 'leaf habit', 'leaf compoundness', 'deltaLWP',
  #                    'dispersal mode', 'seed mass', 'wood density', 'E', 'Amax', 'WUE')
lf.pca$spp <- NULL
lf.pca.mod <- prcomp(lf.pca, center = T, scale = T)
lf.pca.out <- data.frame(lf.pca.mod$x[,1:3])
colnames(lf.pca.out) <- c('ag.PC1','ag.PC2','ag.PC3')
rownames(lf.pca.out) <- lab
lf.pca.out$spp <- lab

#aboveground subset PCA.
lf.sub.pca <- lf.sub
lf.sub.pca <- lf.sub.pca[lf.sub.pca$spp %in% gs.sub,]
lab <- lf.sub.pca$spp
rownames(lf.sub.pca) <- lab
colnames(lf.sub.pca) <- c('spp','thickness', 'leaf area', 'SLA', 'leaf density', 'crown radius', 'SMF',
                      'total leaf area', 'C:N','leaf habit', 'deltaLWP',
                      'dispersal mode', 'seed mass', 'Amax')
lf.sub.pca$spp <- NULL
lf.sub.pca.mod <- prcomp(lf.sub.pca, center = T, scale = T)
lf.pca.sub.out <- data.frame(lf.sub.pca.mod$x[,1:3])
colnames(lf.pca.sub.out) <- c('ag.sub.PC1','ag.sub.PC2','ag.sub.PC3')
rownames(lf.pca.sub.out) <- lab
lf.pca.sub.out$spp <- lab

#all trait PCA.
rt.pca$spp <- rownames(rt.pca)
lf.pca$spp <- rownames(lf.pca)
all.pca <- merge(rt.pca, lf.pca, all.x = T)
all.pca <- all.pca[complete.cases(all.pca),]
lab <- all.pca$spp
all.pca$spp <- NULL
all.pca.mod <- prcomp(all.pca, center = T, scale = T)
all.pca.out <- data.frame(all.pca.mod$x[,1:3])
colnames(all.pca.out) <- c('all.PC1','all.PC2','all.PC3')
rownames(all.pca.out) <- lab
all.pca.out$spp <- lab

#Make list, do merging.----
all.dat <- merge(gs.avg , lf.avg)
all.dat <- merge(all.dat, rt.avg)
#merge beloground PCA to all.dat
to_merge <- data.frame(rt.pca.out)
to_merge$spp <- rownames(to_merge)
all.dat <- merge(all.dat, to_merge, all.x = T)
#merge aboveground PCA to all.dat
to_merge <- data.frame(lf.pca.out)
to_merge$spp <- rownames(to_merge)
all.dat <- merge(all.dat, to_merge, all.x = T)
#merge whole-plant PCA to all.dat
to_merge <- data.frame(all.pca.out)
to_merge$spp <- rownames(to_merge)
all.dat <- merge(all.dat, to_merge, all.x = T)
#merge in subset of ag PCA vectors.
to_merge <- data.frame(lf.pca.sub.out)
all.dat <- merge(all.dat, to_merge, all.x = T)

#form output list.
output <- list(gs.avg,rt.avg,lf.avg,all.dat,lf.pca.out,rt.pca.out,lf.pca.mod,rt.pca.mod,all.pca.mod,lf.pca.sub.out,lf.sub.pca.mod)
names(output) <- c('gs','bg','ag','all.dat','ag.pca','bg.pca','ag.pca.mod','bg.pca.mod','all.pca.mod','ag.sub.pca','ag.sub.pca.mod')

#save output.----
saveRDS(output, output.path)

#end script.----
