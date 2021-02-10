#Visualize PCAs and the traits that best correlate.
rm(list=ls())
library(gridExtra)
library(factoextra)

#set output path.----
output.path <- 'figures/Fig._3.png'

#load data.----
d <- readRDS('data/phase1_for_analysis.rds')

#subset to species w/ complete RGR data.----
gs <- d$gs[complete.cases(d$gs),]
bg <- d$bg[d$bg$spp %in% gs$spp,]
ag <- d$ag[d$ag$spp %in% gs$spp,]

#generate belowground pca and correlation data.----
ag.pca.dat <- ag
rownames(ag.pca.dat) <- ag.pca.dat$spp
ag.pca.dat$spp <- NULL
ag.pca <- d$ag.pca.mod
ag.pca.x <- ag.pca$x

#which traits most correlated with PC1?
y <- ag.pca$x[,1]
ag.pc.rsq <- list()
for(i in 1:ncol(ag.pca.dat)){
  mod <- lm(y ~ ag.pca.dat[,i])
  ag.pc.rsq[[i]] <- summary(mod)$r.squared
}
ag.pc.rsq <- unlist(ag.pc.rsq)
ag.pc.rsq <- data.frame(ag.pc.rsq, colnames(ag.pca.dat))
colnames(ag.pc.rsq) <- c('rsq','trait')
ag.pc.rsq <- ag.pc.rsq[order(ag.pc.rsq$rsq, decreasing = T),]

#specify trait labels.
ag.pca.traitlab
ag.trait.lab <- c('crown radius','total leaf area')
subpanel.lab <- c('B','C')

#Make PCA plot and top 4 regressions w/ PC1 plot.----
p.var.1 <- round(summary(ag.pca)$importance[2,1]*100,1) #grab proportion variance explained.
p.var.2 <- round(summary(ag.pca)$importance[2,2]*100,1) #grab proportion variance explained.
ag.pca.plot <- fviz_pca_biplot(ag.pca, geom='point', repel = T,
                               xlab = paste0('PC1 (',p.var.1,'% variance explained)'),
                               ylab = paste0('PC2 (',p.var.2,'% variance explained)'),
                               title = NULL,
                               select.var = list(contrib = 10)) + #subset to 10 most important vectors
  labs(tag = 'A')

ag.scatter <- list()
for(i in 1:2){
  x <- ag.pca.dat[,colnames(ag.pca.dat) %in% ag.pc.rsq$trait[i]]
  scatter.dat <- data.frame(x,y)
  fit <- lm(y~x, data = scatter.dat)
  ag.scatter[[i]] <- ggplot(scatter.dat, aes(x=x, y=y), size = 3) + 
    geom_point() + 
    theme_bw()  + #drop gray background.
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #drop gridlines
    xlab(ag.trait.lab[i]) +  #x axis label.
    ylab(expression(paste("Aboveground PC1"))) + #y axis label.
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
    scale_y_continuous(expand = expand_scale(mult = c(.01, .02))) + #change where y-axis cuts off.
    scale_x_continuous(expand = expand_scale(mult = c(.01, .01)))  + #change where x-axis cuts off.
    theme(axis.text.x=element_text(size=rel(0.5))) +                  #reduce x-axis text size.
    geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], size = 1) + #add regression line.
    #theme(plot.margin = unit(c(unit.scale,2*unit.scale,2*unit.scale,unit.scale), 'cm')) +
    #geom_text(x=0.50, y = 0.05, label=expression(paste(R^2,'= 0.41'))) +
    labs(tag = subpanel.lab[i])
}


#Specify png output.----
png(output.path, width=9, height= 6, units='in', res=300)

#Drop panels with grid arrange.----
grid.arrange(ag.pca.plot, ag.scatter[[1]], ag.scatter[[2]],
             widths = c(4,2),
             heights = c(4,4),
             layout_matrix = rbind(c(1,2),
                                   c(1,3)
             )
)

#end plot.----
dev.off()
