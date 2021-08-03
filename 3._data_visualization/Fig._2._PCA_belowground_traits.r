#Merge this with figure 3 - aboveground traits.
#Pick two traits, not 4. Root diameter and Root depth. They correlate well and make sense.
#Do the same thing for the aboevground traits.
#Visualize PCAs and the traits that best correlate.
rm(list=ls())
#library(ggbiplot) #I like the factoextra one better, allows you to subset number of vectors.
library(gridExtra)
library(factoextra)


#set output path.----
output.path <- 'figures/Fig._2.png'

#load data.----
d <- readRDS('data/phase1_for_analysis.rds')

#subset to species w/ complete RGR data.----
gs <- d$gs[complete.cases(d$gs),]
bg <- d$bg[d$bg$spp %in% gs$spp,]
ag <- d$ag[d$ag$spp %in% gs$spp,]

#generate belowground pca and correlation data.----
bg.pca.dat <- bg
rownames(bg.pca.dat) <- bg.pca.dat$spp
bg.pca.dat$spp <- NULL
bg.pca <- d$bg.pca.mod
bg.pca.x <- bg.pca$x

#which traits most correlated with PC1?
y <- bg.pca$x[,1]
bg.pc.rsq <- list()
for(i in 1:ncol(bg.pca.dat)){
  mod <- lm(y ~ bg.pca.dat[,i])
  bg.pc.rsq[[i]] <- summary(mod)$r.squared
}
bg.pc.rsq <- unlist(bg.pc.rsq)
bg.pc.rsq <- data.frame(bg.pc.rsq, colnames(bg.pca.dat))
colnames(bg.pc.rsq) <- c('rsq','trait')
bg.pc.rsq <- bg.pc.rsq[order(bg.pc.rsq$rsq, decreasing = T),]

#pick the two traits we like - root diameter and root depth.-----
bg.pc.rsq <- bg.pc.rsq[bg.pc.rsq$trait %in% c('Mean_Dia_C_mm','RMF'),]

#specify trait labels.
bg.trait.lab <- c('root diameter (mm)','RMF')
subpanel.lab <- c('(B)','(C)')
bg.pc.rsq.annotate <- c(expression(~R^2~" = 0.71"), expression(~R^2~" = 0.39")) 

#Make PCA plot and top 2 regressions w/ PC1 plot.----
#bg.pca.plot <- ggbiplot(bg.pca,  varname.adjust = 1.1) + labs(tag = 'A')
#drop biplot.
p.var.1 <- round(summary(bg.pca)$importance[2,1]*100,1) #grab proportion variance explained.
p.var.2 <- round(summary(bg.pca)$importance[2,2]*100,1) #grab proportion variance explained.
bg.pca.plot <- fviz_pca_biplot(bg.pca, geom='point', repel = T,
                               xlab = paste0('belowground PC1 (',p.var.1,'% variance explained)'),
                               ylab = paste0('belowground PC2 (',p.var.2,'% variance explained)'),
                               title = NULL) +
                               labs(tag = '(A)')+
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +    #drop gridlines
               theme(axis.text=element_text(size=10)) +    #increase axis text size.
               theme(axis.line = element_line(colour = "black")) # add axis line
  

bg.scatter <- list()
for(i in 1:2){
  x <- bg.pca.dat[,colnames(bg.pca.dat) %in% bg.pc.rsq$trait[i]]
  scatter.dat <- data.frame(x,y)
  fit <- lm(y~x, data = scatter.dat)
  bg.scatter[[i]] <- ggplot(scatter.dat, aes(x=x, y=y), size = 3) + 
    geom_point() + 
    theme_bw()  + #drop gray background.
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #drop gridlines
    xlab(bg.trait.lab[i]) +  #x axis label.
    ylab(expression(paste("belowground PC1"))) + #y axis label.
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
    scale_y_continuous(expand = expand_scale(mult = c(.02, .02))) + #change where y-axis cuts off.
    scale_x_continuous(expand = expand_scale(mult = c(.02, .02)))  + #change where x-axis cuts off.
    theme(axis.text.x=element_text(size=rel(0.9))) +                  #reduce x-axis text size.
    geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], size = 1) + #add regression line.
    labs(tag = subpanel.lab[i])+
    ggtitle(bg.pc.rsq.annotate[i])+
    theme(plot.title = element_text(size=10))
  
  }

#Specify png output.----
png(output.path, width=9, height= 6, units='in', res=300)

#Drop panels with grid arrange.----
grid.arrange(bg.pca.plot, bg.scatter[[1]], bg.scatter[[2]],
             widths = c(4,2),
             heights = c(4,4),
             layout_matrix = rbind(c(1,2),
                                   c(1,3)
                               )
             )

#end plot.----
dev.off()
