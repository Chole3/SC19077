#' @title Predict the Density Value at certain pointsPoint From the Mixture of Two Normal Distributions by MLE
#' @description Estimate the parameters of mixture of two normal distribution by MLE and then predict the density value at certain points
#' @param x the specific point or points from the mixture of two normal distributions
#' @param initial a vector consisting the initial values of weight, means and deviations of two normal distributions
#' @param data the data which generated from the mixture of two normal distributions
#' @return an estimated value
#' @examples
#' \dontrun{
#' n <- 200
#' u <- runif(n)
#' x <- ifelse(u<.7,rnorm(n),rnorm(n,3,9))
#' mixnormal2(0.8,c(.65,-0.1,1.59,3.21,7.99),x)
#' }
#' @export
mixnormal2 <- function(x,initial,data){
  LL <- function(params, data=data) { 
    t1 <- dnorm(data, params[2], params[3])
    t2 <- dnorm(data, params[4], params[5])
    f <- params[1] * t1 + (1 - params[1]) * t2
    ll <- sum(log(f))
    return(-ll)
  }
  res <- nlminb(initial,LL,
                data=data,
                lower=c(.00001,-Inf,.00001,-Inf,.00001),
                upper=c(0.99999,Inf,Inf,Inf,Inf))
  y <- res$par[1]*dnorm(x,res$par[2],res$par[3])+(1-res$par[1])*dnorm(x,res$par[4],res$par[5])
  return(y)
}

