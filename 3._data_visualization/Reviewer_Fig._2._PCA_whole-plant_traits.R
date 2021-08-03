#Visualize whole plant PCs and the traits that best correlate.
rm(list=ls())
library(gridExtra)
library(factoextra)

#set output path.----
output.path <- 'figures/Reviewer_Fig._2_whole-plant_PCA.png'

#load data.----
d <- readRDS('data/phase1_for_analysis.rds')

#subset to species w/ complete RGR data.----
gs <- d$gs[complete.cases(d$gs),]
bg <- d$bg[d$bg$spp %in% gs$spp,]
ag <- d$ag[d$ag$spp %in% gs$spp,]
wp <- merge(bg,ag)

#generate whole-plant pca and correlation data.----
wp.pca.dat <- wp
rownames(wp.pca.dat) <- wp.pca.dat$spp
wp.pca.dat$spp <- NULL
wp.pca <- d$all.pca.mod
wp.pca.x <- wp.pca$x

#which traits most correlated with PC1?
y <- wp.pca$x[,1]
wp.pc.rsq <- list()
for(i in 1:ncol(wp.pca.dat)){
  mod <- lm(y ~ wp.pca.dat[,i])
  wp.pc.rsq[[i]] <- summary(mod)$r.squared
}
wp.pc.rsq <- unlist(wp.pc.rsq)
wp.pc.rsq <- data.frame(wp.pc.rsq, colnames(wp.pca.dat))
colnames(wp.pc.rsq) <- c('rsq','trait')
wp.pc.rsq <- wp.pc.rsq[order(wp.pc.rsq$rsq, decreasing = T),]

#specify trait labels.
wp.trait.lab <- c('crown radius (cm)',expression(paste("total leaf area (cm"^"2",")")))
subpanel.lab <- c('(B)','(C)')


#Make PCA plot and top 2 regressions w/ PC1 plot.----
p.var.1 <- round(summary(wp.pca)$importance[2,1]*100,1) #grab proportion variance explained.
p.var.2 <- round(summary(wp.pca)$importance[2,2]*100,1) #grab proportion variance explained.
wp.pca.plot <- fviz_pca_biplot(wp.pca, geom='point', repel = T,
                               xlab = paste0('whole-plant PC1 (',p.var.1,'% variance explained)'),
                               ylab = paste0('whole-plant PC2 (',p.var.2,'% variance explained)'),
                               title = NULL,
                               select.var = list(contrib = 10)) + #subset to 10 most important vectors
  labs(tag = '(A)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +    #drop gridlines
  theme(axis.text=element_text(size=10)) +    #increase axis text size.
  theme(axis.line = element_line(colour = "black")) # add axis line

wp.scatter <- list()
for(i in 1:2){
  x <- wp.pca.dat[,colnames(wp.pca.dat) %in% wp.pc.rsq$trait[i]]
  scatter.dat <- data.frame(x,y)
  fit <- lm(y~x, data = scatter.dat)
  wp.scatter[[i]] <- ggplot(scatter.dat, aes(x=x, y=y), size = 3) + 
    geom_point() + 
    theme_bw()  + #drop gray background.
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #drop gridlines
    xlab(wp.trait.lab[i]) +  #x axis label.
    ylab(expression(paste("whole-plant PC1"))) + #y axis label.
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
    scale_y_continuous(expand = expand_scale(mult = c(.02, .02))) + #change where y-axis cuts off.
    scale_x_continuous(expand = expand_scale(mult = c(.02, .03)))  + #change where x-axis cuts off.
    theme(axis.text.x=element_text(size=rel(0.9))) +                  #reduce x-axis text size.
    geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], size = 1) + #add regression line.
    #theme(plot.margin = unit(c(unit.scale,2*unit.scale,2*unit.scale,unit.scale), 'cm')) +
    #geom_text(x=0.50, y = 0.05, label=expression(paste(R^2,'= 0.41'))) +
    labs(tag = subpanel.lab[i])
}


#Specify png output.----
png(output.path, width=9, height= 6, units='in', res=300)

#Drop panels with grid arrange.----
grid.arrange(wp.pca.plot, wp.scatter[[1]], wp.scatter[[2]],
             widths = c(4,2),
             heights = c(4,4),
             layout_matrix = rbind(c(1,2),
                                   c(1,3)
             )
)

#end plot.----
dev.off()
