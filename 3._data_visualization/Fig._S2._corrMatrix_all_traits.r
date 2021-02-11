#correlation matrix of all traits
rm(list=ls())
library(corrplot)
library(RColorBrewer)

#set output path.----
output.path <- 'figures/Fig._S2.png'

#load data.----
d <- readRDS('data/phase1_for_analysis.rds')

#subset to species w/ complete RGR data.----
gs <- d$gs[complete.cases(d$gs),]
all <- d$all.dat[d$all.dat$spp %in% gs$spp,]

#pull data for correlation matrix ----
all.pca.dat <- all
rownames(all.pca.dat) <- all.pca.dat$spp
all.pca.dat$spp <- NULL

# only pull continuous variables
corr.dat <- all.pca.dat[c(3:14, 17, 19:29)]

# rename columns to match
colnames(corr.dat) <- c('petiole length', 'thickness', 'leaf area', 'SLA', 'leaf density', 'crown radius', 'SMF', 'LMF',
                        'total leaf area', 'total leaf perimeter', 'C:N', ':paste(delta^{13},"C")', ':paste(psi[diurnal])',
                        'seed mass', 'wood density', 'E', 'Amax', 'WUE', 'root diameter','SRL','RTD','root lateral extent','root depth','RMF')

# generate pearson correlation matrix and get p-values
M <- cor(corr.dat, method = "pearson")
res1 <- cor.mtest(as.matrix(corr.dat), conf.level = .95)

## plot a pretty correlation matrix
#Specify png output.----
png(output.path, width=9, height= 6, units='in', res=300)

# plot the correlation matrix 
corrplot(M, p.mat = res1$p, method = "square", type = "upper", tl.col = "black",
                sig.level = c(.001, .01, 0.05), pch.cex = .9,
                insig = "label_sig", pch.col = "white",  diag = FALSE,
                col = brewer.pal(n = 8, name = "RdBu"))

#end plot.----
dev.off()

###