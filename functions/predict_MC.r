#Functino to calcualte prediction mean and 95% CI from MCMCglmm object.
predict_MC <- function(model, newdat, n.sim=1000, q=c(0.025,0.975)){
  #grab parameter posterior.
  post <- as.matrix(model$Sol)
  colnames(post)[1] <- 'intercept'
  #grab sample size.
  N <- length(model$error.term)
  
  #we are assuming this is the special case of RGR:bg.PC1 interaction.
  newdat$intercept <- 1
  newdat$`RGR:bg.PC1` <- newdat$RGR * newdat$bg.PC1
  newdat <- newdat[,order(match(colnames(newdat), colnames(post)))]
  
  #Do monte-carlo to get prediction mean and 95% quantiles
  mc.out <- list()
  for(i in 1:n.sim){
            par <- post[sample(nrow(post),1),] #draw a random parameter combination.
    mc.out[[i]] <- as.matrix(newdat) %*% par   #multiply predictor matrix by paramters.
  }
  #calculate mean and quantiles.
  mc.out <- do.call(cbind, mc.out)
  mu <- rowMeans(mc.out)
  quantile.lo <- apply(mc.out, 1, quantile, probs=q[1])
  quantile.hi <- apply(mc.out, 1, quantile, probs=q[2])
  #calculate standard error.
  se <- (quantile.hi - mu)/sqrt(N)
  mu.loSE <- mu - se
  mu.hiSE <- mu + se
  
  #return output.
  output <- list(mu, quantile.lo, quantile.hi,se,mu.loSE,mu.hiSE)
  names(output) <- c('predicted','lo95','hi95','se','loSE','hiSE')
  return(output)
}