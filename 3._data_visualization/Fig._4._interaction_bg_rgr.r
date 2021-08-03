#Visualizing interactions between belowground PC1 and RGR.
rm(list=ls())
library(boot)
library(viridis)
library(MCMCglmm)
source('functions/predict_MC.r')

#set output path.----
output.path <- 'figures/Fig._4.png'

#load data.----
#d <- readRDS('data/surv_models_fitted.rds')
#fit <- d$models$m.rgr_bg
#d <- d$data
d <- readRDS('data/phylo_surv_models_fitted.rds')
fit <- d$phylo.models$models$m.rgr_bg
d   <- d$phylo.models$data

#grab low and high RGR values, make low and high RGR dat, detrend.----
rgr.rel <- quantile(d$RGR, probs = c(0.2, 0.8))
lo.dat <- data.frame(rgr.rel[1], d$bg.PC1)
hi.dat <- data.frame(rgr.rel[2], d$bg.PC1)
colnames(lo.dat) <- c('RGR','bg.PC1')
colnames(hi.dat) <- c('RGR','bg.PC1')

#generate detrended values.
lo.detrend <- fit$residuals + predict_MC(fit, newdat = lo.dat)$predicted
hi.detrend <- fit$residuals + predict_MC(fit, newdat = hi.dat)$predicted

#Generate predictions across RGR range w/ different PC1 scores.----
pc1 <- quantile(d$bg.PC1, probs = c(0.05, 0.5, 0.95)) #5%, 50% and 95% quantiles of PC1
rgr <- quantile(d$RGR   , probs = c(0.10, 0.90)) #10% and 90% quantles of RGR.
x <- seq(rgr[1], rgr[2], by = (rgr[2] - rgr[1]) / 100)

#make your dataframes.
dat.list <- list()
for(i in 1:length(pc1)){
  dat.list[[i]] <- data.frame(x, pc1[i])
  colnames(dat.list[[i]]) <- c('RGR','bg.PC1')
}

#Generate predictions for each dataframe.
pred.list <- list()
se.list <- list()
for(i in 1:length(dat.list)){
  output <- predict_MC(fit, newdat = dat.list[[i]])
  preds <- output$predicted
  lo <- inv.logit(output$loSE)
  hi <- inv.logit(output$hiSE)
  se.out <- list(lo, hi)
  names(se.out) <- c('lo','hi')
  pred.list[[i]] <- inv.logit(preds)
    se.list[[i]] <-se.out
}

#global chart options.----
line = 0.6
cex = 0.9
side = 3
adj=-.34

limy <- c(0,1)
magma.cols <- magma(100, begin = 0)
cols <- c(magma.cols[80], magma.cols[60], magma.cols[40])
cols <- c('#0072B2',magma.cols[40],'#E69F00')
trans <- 0.3 #shading transparency

#png save line.----
png(output.path, width=7, height= 5, units='in', res=300)
par(mar = c(4,4,2,2))

#setup plot matrix for layout.----
m <- rbind(c(1, 1, 2), c(1, 1, 3))
layout(m)


#Panel 1 - interaction plot.----
plot(pred.list[[1]] ~ x, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0,1), cex.axis = 1.1)
for(i in 1:length(pred.list)){
  lines(smooth.spline(pred.list[[i]] ~ x), col = cols[i], lwd = 2)
  polygon(c(x, rev(x)),c(se.list[[i]]$hi, rev(se.list[[i]]$lo)), col=adjustcolor(cols[i], trans), lty=0)
}

#axis labels
mtext('estimated survival rate (%)', side = 2, line = 2.5, cex = 0.9)
rgr.lab <- expression(paste("relative growth rate (ln[cm] day"^"-1",")"))
mtext(rgr.lab, side = 1, line = 3, cex = 0.9)
mtext("(A)", side=side, line=line, cex=cex, adj=-0.13)

#legend.
legend(x = 0.0008, y = 1, c('thick and deep','average','thin and shallow'), col = cols, bty = 'n', lty = 1, lwd = 2, cex = 1.2)

#Panel 2 - survival vs. PC1 @low  RGR.----
plot(inv.logit(lo.detrend) ~ d$bg.PC1, pch = 16, cex = 1.5, bty = 'l', ylab = NA, xlab = NA, ylim = c(0,1))
abline(lm(inv.logit(lo.detrend) ~ d$bg.PC1), lwd  = 2, lty = 2)
mtext('survival at low RGR', side = 3, line = -1.0, cex = 0.7)
mtext('survival rate (%)'      , side = 2, line =  2.2, cex = 0.7)
mtext('belowground PC1'    , side = 1, line =  2.2, cex = 0.7)
mtext("(B)", side=side, line=line, cex=cex, adj=adj)

#Panel 3 - survival vs. PC1 @high RGR.----
plot(inv.logit(hi.detrend) ~ d$bg.PC1, pch = 16, cex = 1.5, bty = 'l', ylab = NA, xlab = NA, ylim = c(0,1))
abline(lm(inv.logit(hi.detrend) ~ d$bg.PC1), lwd  = 2, lty = 1)
mtext('survival at high RGR', side = 3, line = -1.0, cex = 0.7)
mtext('survival rate (%)'       , side = 2, line =  2.2, cex = 0.7)
mtext('belowground PC1'     , side = 1, line =  2.2, cex = 0.7)
mtext("(C)", side=side, line=line, cex=cex, adj=adj)

#end plot.----
dev.off()
