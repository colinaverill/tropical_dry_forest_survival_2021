#visualizing best fits.
#clear environment, load packages.
rm(list=ls())
library(boot)

#set output path.----
output.path <- 'figures/Reviewer_Fig._1_additional_models.png'

#load data.----
d <- readRDS('data/surv_models_fitted.rds')
d <- readRDS('data/phylo_surv_models_fitted.rds')
d <- d$phylo.models

#png save line.----
png(output.path, width=5, height= 5, units='in', res=300)



#global plot settings.----
line = 1
cex = 1.2
side = 3
adj=-.25

par(mfrow = c(1,1),
    mar = c(6,5,2,2))
olab.cex <- 1.2

#panel 1- rsq adjusted values ranked.----
z <- d$rsq.sum
#z <- z[order(z$rsq.adj),]
z$lab <- c('aboveground\n traits','belowground\n traits','whole-plant\n traits','aboveground\n traits subset',
           'relative growth\n rate (RGR)',
           'aboveground\n traits + RGR',
           'belowground\n traits + RGR',
           'whole-plant\n traits + RGR',
           'aboveground\n traits subset + RGR')

lty.o <- par("lty") #remember current grahpic setting.
par(lty = 0) #set outer lines to zero width.
barplot(z$rsq.adj, col = 'hotpink', xaxt='n', horiz = F, ylim = c(-0.1, 0.8))
ylab <- bquote({R^2}[adj])
mtext(ylab, side = 2, line=2.5, cex = olab.cex)
#mtext("(A)", side=side, line=line, cex=cex, adj=adj)
par(lty = lty.o) #reset graphic setting.
pos <- c(0:(length(z$lab)-1))
N <- 9.6
pos <- seq(0,N, by = N/8)
text(cex=0.8, x= 0.7+pos, y=-0.02, z$lab, xpd=TRUE, srt = 90, adj=1)

#end plot.----
dev.off()
