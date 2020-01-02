#' @title Parameters Estimation of Mixture of Two Normal Distributions by EM Algorithms
#' @description Estimate the parameters of mixture of two normal distribution by EM algorithms 
#' @param p the initial of the weights of two normal distributions
#' @param mu the initial of the population means of two normal distributions
#' @param sig the initial of the population deviations of two normal distributions
#' @param data the data which generated from the mixture of two normal distributions
#' @return a list of estimated values of parameters
#' @examples
#' \dontrun{
#' n <- 200
#' u <- runif(n)
#' x <- ifelse(u<.7,rnorm(n),rnorm(3,9))
#' gmm2(c(0.65,0.35),c(-0.1,3.21),c(1.59,7.99),x)
#' }
#' @export
gmm2 <- function(p,mu,sig,data){
  n <- length(data)
  dat <- matrix(rep(data,2),ncol=2,byrow=F)
  prob <- matrix(0,nrow=n,ncol=2)
  
  while(TRUE){
    for(i in 1:2)prob[,i] <- as.matrix(sapply(dat[,i],dnorm,mu[i],sig[i]))
    pmatrix <- matrix(rep(p,n),ncol=2,byrow=T)
    a <- pmatrix*prob
    b <- matrix(rep(rowSums(a),2),ncol=2)
    w <- a/b
    
    oldp <- p
    oldmu <- mu
    oldsig <- sig

    p <- colSums(w)/n
    mu <- colSums(dat*w)/colSums(w) 
    d <- matrix(rep(mu,n),ncol=2,byrow=T)
    sig <- colSums(w*(dat-d)^2)/colSums(w)

    if((sum((p-oldp)^2) < 1e-4) &
       (sum((mu-oldmu)^2) < 1e-4) &
       (sum((sig-oldsig)^2) < 1e-4) )break
  }
  return(list(p=p,mu=mu,sig=sig))
}
