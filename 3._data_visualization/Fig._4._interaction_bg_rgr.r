#Visualizing interactions between belowground PC1 and RGR.
rm(list=ls())
library(boot)
library(viridis)

#set output path.----
output.path <- 'figures/Fig._4.png'

#load data.----
d <- readRDS('data/surv_models_fitted.rds')
fit <- d$models$m.rgr_bg
d   <- d$data

#grab low and high RGR values, make low and high RGR dat, detrend.----
rgr.rel <- quantile(d$RGR, probs = c(0.2, 0.8))
lo.dat <- data.frame(rgr.rel[1], d$bg.PC1)
hi.dat <- data.frame(rgr.rel[2], d$bg.PC1)
colnames(lo.dat) <- c('RGR','bg.PC1')
colnames(hi.dat) <- c('RGR','bg.PC1')

#generate detrended values.
lo.detrend <- resid(fit) + predict(fit, newdata = lo.dat)
hi.detrend <- resid(fit) + predict(fit, newdata = hi.dat)

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
  preds <- predict(fit, newdata = dat.list[[i]], se.fit = T)
  lo <- inv.logit(preds$fit - preds$se.fit)
  hi <- inv.logit(preds$fit + preds$se.fit)
  se.out <- list(lo, hi)
  names(se.out) <- c('lo','hi')
  pred.list[[i]] <- inv.logit(preds$fit)
  se.list[[i]] <-se.out
}

#global chart options.----
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
plot(pred.list[[1]] ~ x, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0,1))
for(i in 1:length(pred.list)){
  lines(smooth.spline(pred.list[[i]] ~ x), col = cols[i], lwd = 2)
  polygon(c(x, rev(x)),c(se.list[[i]]$hi, rev(se.list[[i]]$lo)), col=adjustcolor(cols[i], trans), lty=0)
}

#axis labels
mtext('estimated survival rate', side = 2, line = 2.5, cex = 1)
mtext('relative growth rate'   , side = 1, line = 2.5, cex = 1)

#legend.
legend(x = 0.0008, y = 1, c('thick and deep','average','thin and shallow'), col = cols, bty = 'n', lty = 1, lwd = 2, cex = 1.2)

#Panel 2 - survival vs. PC1 @low  RGR.----
plot(inv.logit(lo.detrend) ~ d$bg.PC1, pch = 16, cex = 1.5, bty = 'l', ylab = NA, xlab = NA, ylim = c(0,1))
abline(lm(inv.logit(lo.detrend) ~ d$bg.PC1), lwd  = 2, lty = 2)
mtext('survival at low RGR', side = 3, line = -1.0, cex = 0.7)
mtext('survival rate'      , side = 2, line =  2.2, cex = 0.7)
mtext('belowground PC1'    , side = 1, line =  2.2, cex = 0.7)

#Panel 3 - survival vs. PC1 @high RGR.----
plot(inv.logit(hi.detrend) ~ d$bg.PC1, pch = 16, cex = 1.5, bty = 'l', ylab = NA, xlab = NA, ylim = c(0,1))
abline(lm(inv.logit(hi.detrend) ~ d$bg.PC1), lwd  = 2, lty = 1)
mtext('survival at high RGR', side = 3, line = -1.0, cex = 0.7)
mtext('survival rate'       , side = 2, line =  2.2, cex = 0.7)
mtext('belowground PC1'     , side = 1, line =  2.2, cex = 0.7)

#end plot.----
dev.off()
