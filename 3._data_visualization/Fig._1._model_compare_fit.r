#visualizing best fits.
#clear environment, load packages.
rm(list=ls())
library(boot)

#set output path.----
output.path <- 'figures/Fig._1.png'

#load data.----
d <- readRDS('data/surv_models_fitted.rds')

#png save line.----
png(output.path, width=9, height= 5, units='in', res=300)

#global plot settings.----
par(mfrow = c(1,2),
    mar = c(6,5,2,2))
olab.cex <- 1.2

#panel 1- rsq adjusted values ranked.----
z <- d$rsq.sum
#z <- z[order(z$rsq.adj),]
z$lab <- c('aboveground\n traits','belowground\n traits','relative growth\n rate (RGR)','aboveground\n traits + RGR','belowground\n traits + RGR')

lty.o <- par("lty") #remember current grahpic setting.
par(lty = 0) #set outer lines to zero width.
barplot(z$rsq.adj, col = 'hotpink', xaxt='n', horiz = F, ylim = c(-0.1, 0.8))
ylab <- bquote({R^2}[adj])
mtext(ylab, side = 2, line=2.5, cex = olab.cex)
par(lty = lty.o) #reset graphic setting.
pos <- c(0:4)
N <- 5
pos <- seq(0,N, by = N/4)
text(cex=1, x= 0.7+pos, y=-0.02, z$lab, xpd=TRUE, srt = 90, adj=1)
#text(cex=1, x= 0.7+pos, y=-0.02, z$lab, xpd=TRUE, srt = 45, adj=1)
#text(cex=1, x= 0.7+pos, y=1, z$lab, xpd=TRUE, srt = 90, adj=3)


#Panel 2 - visualize best model.----
m <- d$models$m.rgr_bg
plot(inv.logit(m$model[,1]) ~ inv.logit(m$fitted.values), bty = 'l', pch = 16, cex = 2, xlab = NA, ylab = NA)
abline(lm(inv.logit(m$model[,1]) ~ inv.logit(m$fitted.values)), lwd = 2, lty = 1)
#setup Rsq label.
rsq <- round(d$rsq.sum[d$rsq.sum$mod.lab == 'm.rgr_bg',]$rsq.adj, 2)
rsq.lab <- bquote({R^2}[adj] == .(rsq))
#drop rsq label
#mtext('Best model:\n belowground traits + RGR', side = 3, line = -2, adj = 0.05)
mtext('best model', side = 3, line = -1, adj = 0.05)
mtext(rsq.lab, side = 3, line = -2.5, adj = 0.05)
#x and y axis labels.
mtext('observed survival' , line = 2.5, side = 2, cex = olab.cex)
mtext('estimated survival', line = 2.5, side = 1, cex = olab.cex)

#end plot.----
dev.off()
