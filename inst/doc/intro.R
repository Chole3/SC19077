## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
n <- 200
u <- runif(n)
x <- ifelse(u<.7,rnorm(n),rnorm(n,3,9))
mixnormal2(0.8,c(.65,-0.1,1.59,3.21,7.99),x)

## -----------------------------------------------------------------------------
gmm <- function(p,mu,sig,data){
  n <- length(data)
  dat <- matrix(rep(data,2),ncol=2,byrow=F)
  prob <- matrix(0,nrow=n,ncol=2)
  
  while(TRUE){
    for(i in 1:2)
      prob[,i] <- as.matrix(sapply(dat[,i],dnorm,mu[i],sig[i]))
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

## -----------------------------------------------------------------------------
n <- 200
u <- runif(n)
x <- ifelse(u<.7,rnorm(n),rnorm(n,3,9))
gmm(c(0.65,0.35),c(-0.1,3.21),c(1.59,7.99),x)

## -----------------------------------------------------------------------------
library(knitr)
num <- 1:5
Mon <- c(' ',' ','English',' ','Politcs')
Tues <- c(' ',' ','Probability',' ',' ')
Wed <- c('Linear Model',' ',' ',' ',' ')
Thur <- c('Probability',' ',' ',' ',' ')
Fri <- c(' ','Computing','Linear Model',' ',' ')
schedule <- data.frame(num,Mon,Tues,Wed,Thur,Fri)
kable(schedule)

## -----------------------------------------------------------------------------
library(MASS)
attach(iris)
train <- sample.int(150,100)
fit <- lda(Species~.,iris[train,])
plot(fit,main="QDA for iris Dataset")

## -----------------------------------------------------------------------------
sigma <- 1:4 
f1 <- function(x){
x/(sigma[1]^2)*exp(-x^2/(2*sigma[1]^2))}

f2 <- function(x){
x/(sigma[2]^2)*exp(-x^2/(2*sigma[2]^2))}

f3 <- function(x){
x/(sigma[3]^2)*exp(-x^2/(2*sigma[3]^2))}

f4 <- function(x){
x/(sigma[4]^2)*exp(-x^2/(2*sigma[4]^2))}
t <- seq(0,10,0.01)
curve(f1,0,10)
curve(f2,0,10)
curve(f3,0,10)
curve(f4,0,10)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1e2
rn <- matrix(rep(0,n*4),ncol=4)

i <- c(0,0,0,0) 
j <- c(0,0,0,0) 

for(k in 1:4){

  while(i[k]<n){
  u <- runif(1) 
  j[k] <- j[k]+1
  x <- rgamma(1,shape=2,scale=2*sigma[k]^2)
  
   if(u <= exp((x-x^2)/(2*sigma[k]^2))){
    i[k] <- i[k]+1
    rn[i[k],k] <- x
   }
  
  }
  
}
j
i

## -----------------------------------------------------------------------------

hist(rn[,1],prob=T,xlab='x',ylab='f(x)',main=expression(f(x)==x/sigma[1]^2*e^(-x^2/2*sigma[1]^2)))
lines(t,f1(t))

hist(rn[,2],prob=T,xlab='x',ylab='f(x)',main=expression(f(x)==x/sigma[2]^2*e^(-x^2/2*sigma[2]^2)))
lines(t,f2(t))

hist(rn[,3],prob=T,xlab='x',ylab='f(x)',main=expression(f(x)==x/sigma[3]^2*e^(-x^2/2*sigma[3]^2)))
lines(t,f3(t))

hist(rn[,4],prob=T,xlab='x',ylab='f(x)',main=expression(f(x)==x/sigma[4]^2*e^(-x^2/2*sigma[4]^2)))
lines(t,f4(t))


## -----------------------------------------------------------------------------
set.seed(2)
n <- 1e3
u <- runif(n)
rn <- ifelse(u<.75,rnorm(n),rnorm(n,3,1))
t <- seq(-6,6,.01)
hist(rn,prob=T,xlim=c(-6,6),ylim=c(0,0.4),main='Normal Location Mixtrue of N(0,1) and N(3,1) for p1=0.75')
lines(t,0.75*dnorm(t,0,1)+0.25*dnorm(t,3,1))

## -----------------------------------------------------------------------------
p1 <- seq(.05,.95,.05) 
mix <- matrix(rep(0,n*19),ncol=19)

for(i in 1:19){
  u <- runif(n)
  mix[,i] <- ifelse(u<p1[i],rnorm(n),rnorm(n,3,1))
}

x <- seq(-4,8,.01)
for(i in 1:19){
  p <- 0.05*i
  hist(mix[,i],prob=T,xlim=c(-4,8),ylim=c(0,0.5),main='')
  legend(5,0.5,title="p1",p,xjust=0,box.lty = 0)
 lines(x,p1[i]*dnorm(x,0,1)+(1-p1[i])*dnorm(x,3,1))
}

## ----echo=T-------------------------------------------------------------------
rmvn <- function(d,n,Sigma){
  
  T <- matrix(rep(0,d*d),nrow=d)
  for(i in 1:d){
   x <- rchisq(1,n-i+1)
   T[i,i] <- sqrt(x)
  }
  for(i in 2:d){
   for(j in 1:(i-1) ){
      T[i,j] <- rnorm(1)
   }
  }
  A <- T %*% t(T)
  
  ev <- eigen(Sigma)
  lambda <- ev$values
  v <- ev$vectors

  L <- v %*% diag(sqrt(lambda)) %*% t(v)
  
  rn <- L %*% A %*% t(L)
  
}


## ----echo=T-------------------------------------------------------------------
n <- 5
d <- 3
B <- matrix(c(1,1,1,1,2,4,1,3,9),nrow=3)
Sigma <- (B/10) %*% t(B/10)
Sigma
r <- rmvn(d,n,Sigma)
r

## -----------------------------------------------------------------------------
real <- cos(0) - cos(pi/3)
real

set.seed(1)
n <- 1e4
x <- runif(n,0,pi/3)
theta <- round(mean(pi/3*sin(x)),5)
theta <- real 

## -----------------------------------------------------------------------------
set.seed(4)
g <- function(x){pi/3*sin(x)}
f <-function(x){pi/3*x}

n <- 1e4
r <- 1e3
u <- runif(r,0,pi/3)
c <- -cov(f(u),g(u))/var(f(u))
c

MC1 <- MC2 <- numeric(n)
for (i in 1:n){
  v <- runif(r,0,pi/3)
  MC1[i] <- mean(g(v))
  MC2[i] <- mean(g(v)+c*(f(v)-pi/3*pi/6))
}
theta1 <- mean(MC1)
theta2 <- mean(MC2)
theta1
theta2
(var(MC1)-var(MC2))/var(MC1)

## -----------------------------------------------------------------------------
set.seed(2)

g <- function(x){exp(-x)/(1+x^2)}
n <- 1e3
r <- 1e4

MC1 <- MC2 <- numeric(n)
for (i in 1:n){
  u <- runif(r/2)
  v <- runif(r/2)
  MC1[i] <- round(mean(g(c(u,v))),5)
  MC2[i] <- round(mean(g(c(u,1-u))),5)
}
alpha1 <- round(mean(MC1),5)
alpha2 <- round(mean(MC2),5)
alpha1
alpha2
(var(MC1)-var(MC2))/var(MC1)

## -----------------------------------------------------------------------------
set.seed(3)
k <- 5
Finvers <- function(x){-log(1-(1-exp(-1))*x)}
a <- c(rep(0,5),1)
for (i in 2:5){
  a[i] <- Finvers((i-1)/5)
}
a

## -----------------------------------------------------------------------------
n <- 100
r <- 100
g <- function(x){exp(-x)/(1+x^2) * (x>0) * (x<1)}
f <- function(x){exp(-x)/(1-exp(-1))}
p <- function(x){k*f(x)}
Ginvers <- function(a,u){-log(exp(-a)-(1-exp(-1))/5*u)}

theta <- numeric(n)
b <- numeric(k)
for (i in 1:n){
  for (j in 1:k){
    u <- runif(r/k,a[j],a[j+1])
    x <- Ginvers(a[j],u)
    b[j] <- mean( g(x)/p(x) * (x>a[j]) * (x<a[j+1]) )
  }
  theta[i] <- sum(b)
}
round(mean(theta),5)

se <- 0.0970314
(se^2 - var(theta))/se^2

## -----------------------------------------------------------------------------
set.seed(1)
m <- 1e4
n <- 20
alpha <- .05
UCL <- LCL <- numeric(m)
for (i in 1:m){
  x <-  rchisq(n,2)
  UCL[i]=mean(x)+var(x)/sqrt(n)*qt(alpha/2,df=n-1,lower.tail = F)
  LCL[i]=mean(x)-var(x)/sqrt(n)*qt(alpha/2,df=n-1,lower.tail = F)
}
r <- mean(LCL < 2 & 2 < UCL)
r
abs(r - (1-alpha)) - abs(0.956 - (1-alpha))
se <- sqrt(r*(1-r)/m)
se
(0.00689 - se)/0.00689

## ----eval=TRUE----------------------------------------------------------------
set.seed(2) 
m <- 100
a <- 1e3
n <- 50
q <- c(.025,.05,.95,.975)
xq <- qnorm(q,0,sqrt(6/n))

sk <- function(x){
  xbar <- mean(x)
  m3 <- mean((x-xbar)^3)
  m2 <- mean((x-xbar)^2)
  return(m3/m2^1.5)
}
 
b <- numeric(a)
for(j in 1:a){
    x <- rnorm(n)
    b[j] <- sk(x)
  }
M <- quantile(b,q)

for (i in 1:(m-1)){
   for(j in 1:a){
      x <- rnorm(n)
      b[j] <- sk(x)
   }
   r <- quantile(b,q)
   M <- rbind(M,r)
}

xq.est <- apply(M,2,mean)
se.xq.est <- apply(M,2,sd)
t <- dnorm(xq,0,sqrt(6*(n-2)/((n+1)*(n+3))))
se <- sqrt( q*(1-q)/n/t^2 )
knitr::kable(rbind(xq,xq.est), formate = "html",col.names = c(".025",".05","0.95","0.975"))
knitr::kable(rbind(se,se.xq.est), formate = "html",col.names = c(".025",".05","0.95","0.975"))

## -----------------------------------------------------------------------------
set.seed(1)
alpha <- .05
a <- 1
b <- 10
dis <- seq(0,1,.05)
m <- 1e4
n <- 30
cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/(n+1)/(n+3)))
pwr1 <- numeric(length(dis))
sk <- function(x){
  xbar <- mean(x)
  return(mean((x-xbar)^3)/(mean((x-xbar)^2))^1.5)
}

for(i in 1:length(dis)){
  d <- dis[i]
  sktest <- numeric(m)
  for(j in 1:m){
   u <- runif(n)
   x <- ifelse(u<d, rbeta(n,a,a),rbeta(n,b,b))
    sktest[j] <- as.integer(abs(sk(x)) >= cv)
  }
 pwr1[i] <- mean(sktest) 
}
plot(dis,pwr1,type="b",xlab="d",ylab="power",main="Power of skewness test against Beta(a,a)")
abline(h=alpha,lty = 3)

## -----------------------------------------------------------------------------
set.seed(2)
alpha <- .05
t1 <- 2
t2 <- 9#parameter of beta distribution
dis <- seq(0,1,.05)#probabilities to sample
m <- 1e4
n <- 30
cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/(n+1)/(n+3)))
pwr2 <- numeric(length(dis))
sk <- function(x){
  xbar <- mean(x)
  return(mean((x-xbar)^3)/(mean((x-xbar)^2))^1.5)
}

for(i in 1:length(dis)){
  d <- dis[i]
  sktest <- numeric(m)
  for(j in 1:m){
   u <- runif(n)
   x <- ifelse(u<d, rt(n,t1),rt(n,t2))
   sktest[j] <- as.integer(abs(sk(x)) >= cv)
  }
 pwr2[i] <- mean(sktest) 
}
plot(dis,pwr2,type="b",xlab="d",ylab="power",main="Power of skewness test against t(v)")
abline(h=alpha,lty=3)

## -----------------------------------------------------------------------------
set.seed(3)
alpha <- .05
mu0 <- 1
n <- seq(20,500,30)
m <- 100
p.hat <- numeric(length(n))
for(i in 1:length(n)){
  N <- n[i]
  ttest <- numeric(m)
  for(j in 1:m){
    x <- rchisq(N,1)
    ttest[j] <- t.test(x,mu = mu0)$p.value
  }
  p.hat[i] <- mean(ttest <= alpha)
}
se <- sqrt(p.hat*(1-p.hat)/m)

plot(n,p.hat,ylab="tIe",type="b",pch=20)
abline(h=alpha,lty=3)
lines(n, p.hat+se ,lty=3)
lines(n,p.hat-se, lty=3)

## -----------------------------------------------------------------------------
set.seed(4)
alpha <- .05
mu0 <- 1
n <- seq(20,500,30)
m <- 100
p.hat <- numeric(length(n))
for(i in 1:length(n)){
  N <- n[i]
  ttest <- numeric(m)
  for(j in 1:m){
    x <- runif(N,0,2)
    ttest[j] <- t.test(x,mu = mu0)$p.value
  }
  p.hat[i] <- mean(ttest <= alpha)
}
se <- sqrt(p.hat*(1-p.hat)/m)

plot(n,p.hat,ylab="tIe",type="b",pch=20)
abline(h=alpha,lty=3)
lines(n, p.hat+se ,lty=3)
lines(n,p.hat-se, lty=3)

## -----------------------------------------------------------------------------
set.seed(6)
alpha <- .05
mu0 <- 1
n <- seq(20,500,30)
m <- 100
p.hat <- numeric(length(n))
for(i in 1:length(n)){
  N <- n[i]
  ttest <- numeric(m)
  for(j in 1:m){
    x <- rexp(N,1)
    ttest[j] <- t.test(x,mu = mu0)$p.value
  }
  p.hat[i] <- mean(ttest <= alpha)
}
se <- sqrt(p.hat*(1-p.hat)/m)

plot(n,p.hat,ylab="tIe",type="b",pch=20)
abline(h=alpha,lty=3)
lines(n, p.hat+se ,lty=3)
lines(n,p.hat-se, lty=3)

## -----------------------------------------------------------------------------
set.seed(7)
n <- 30
m <- 100
p1 <- 0.651
p2 <- 0.676
x <- rbinom(n,size=m,prob=p1)
y <- rbinom(n,size=m,prob=p2)
ttest <- t.test(x,y,conf.level = alpha)
p <- ttest$p.value
p

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----results = "hide"---------------------------------------------------------
library(bootstrap)
data(scor)
attach(scor)
pairs(scor,labels=c("mec","vec","alg","ana","sta"))
round(cor(scor),4)

## -----------------------------------------------------------------------------
set.seed(1)
library(boot)
n <- nrow(scor)
r <- 100
f12 <- function(x,i)cor(x[i,1],x[i,2])
f34 <- function(x,i)cor(x[i,3],x[i,4])
f35 <- function(x,i)cor(x[i,3],x[i,5])
f45 <- function(x,i)cor(x[i,4],x[i,5])
cor12 <- boot(data=scor, statistic=f12, R=r)
cor34 <- boot(data=scor, statistic=f34, R=r)
cor35 <- boot(data=scor, statistic=f35, R=r)
cor45 <- boot(data=scor, statistic=f45, R=r)
boot.est <- cbind(cor12[1],cor34[1],cor35[1],cor45[1])
boot.est <- round(as.numeric(boot.est),3)
y <- cbind(cor12$t,cor34$t,cor35$t,cor45$t)
se.hat <- round(apply(y,2,sd),4)
f = data.frame(boot.est,se.hat,row.names=c("rho12","rho34","rho35","rho45"))
knitr::kable(f)

## ----results = "hide"---------------------------------------------------------
set.seed(2)
n <- 20
m <- 100
sk1 <- 0

library(boot)
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
b <- function(x,i)sk(x[i])

n.ci.norm <- n.ci.basic <- n.ci.perc <- matrix(NA,m,2)
n.cover <- n.misleft <- n.misright <- numeric(3)

for (i in 1:m){
   x <- rnorm(n)
   n.boot.est <- boot(data=x,statistic=b,R=m)
   n.boot.ci <- boot.ci(n.boot.est,type=c("norm","basic","perc"))
   n.ci.norm[i,] <- n.boot.ci$normal[2:3]
   n.ci.basic[i,] <- n.boot.ci$basic[4:5]
   n.ci.perc[i,] <- n.boot.ci$percent[4:5]
 } 

n.norm.indic <- n.basic.indic <- n.perc.indic <- numeric(m)
n.norm.left <- n.norm.right <- numeric(m)
n.basic.left <- n.basic.right <- numeric(m)
n.perc.left <- n.perc.right <- numeric(m)
for (i in 1:m){
  n.norm.indic[i] <- as.integer(n.ci.norm[i,1] <=  0 & 0 <= n.ci.norm[i,2])
  n.basic.indic[i] <- as.integer(n.ci.basic[i,1] <= 0 & 0 <= n.ci.basic[i,2])
  n.perc.indic[i] <- as.integer(n.ci.perc[i,1] <= 0 & 0 <= n.ci.perc[i,2])
  n.norm.right[i] <- as.integer(n.ci.norm[i,2] < 0)
  n.basic.right[i] <- as.integer(n.ci.basic[i,2] < 0)
  n.perc.right[i] <- as.integer(n.ci.perc[i,2] < 0)
  n.norm.left[i] <- as.integer(0 < n.ci.norm[i,1] )
  n.basic.left[i] <- as.integer(0 < n.ci.basic[i,1])
  n.perc.left[i] <- as.integer(0 < n.ci.perc[i,1])
}
n.cover <- c(mean(n.norm.indic),mean(n.basic.indic),mean(n.perc.indic))
n.misleft <- c(mean(n.norm.left),mean(n.basic.left),mean(n.perc.left))
n.misright <- c(mean(n.norm.right),mean(n.basic.right),mean(n.perc.right))
n.cover
n.misright
n.misleft

## ----results = "hide"---------------------------------------------------------
set.seed(3)
n <- 20
m <- 100

f <- function(x){dchisq(x,5)*(x-5)^3/sqrt(10)^3}
sk2 <- round(integrate(f,0,Inf)$value,3)
sk2

ch.ci.norm <- ch.ci.basic <- ch.ci.perc <- matrix(NA,m,2)
ch.cover <- ch.misleft <- ch.misright <- numeric(3)

for (i in 1:m){
   x <- rchisq(n,5)
   ch.boot.est <- boot(data=x,statistic=b,R=m)
   ch.boot.ci <- boot.ci(ch.boot.est,type=c("norm","basic","perc"))
   ch.ci.norm[i,] <- ch.boot.ci$normal[2:3]
   ch.ci.basic[i,] <- ch.boot.ci$basic[4:5]
   ch.ci.perc[i,] <- ch.boot.ci$percent[4:5]
 } 

ch.norm.indic <- ch.basic.indic <- ch.perc.indic <- numeric(m)
ch.norm.left <- ch.norm.right <- numeric(m)
ch.basic.left <- ch.basic.right <- numeric(m)
ch.perc.left <- ch.perc.right <- numeric(m)
for (i in 1:m){
  ch.norm.indic[i] <- as.integer(ch.ci.norm[i,1] <=  sk2 & sk2 <= ch.ci.norm[i,2])
  ch.basic.indic[i] <- as.integer(ch.ci.basic[i,1] <= sk2 & sk2 <= ch.ci.basic[i,2])
  ch.perc.indic[i] <- as.integer(ch.ci.perc[i,1] <= sk2 & sk2 <= ch.ci.perc[i,2])
  ch.norm.right[i] <- as.integer(ch.ci.norm[i,2] < sk2)
  ch.basic.right[i] <- as.integer(ch.ci.basic[i,2] < sk2)
  ch.perc.right[i] <- as.integer(ch.ci.perc[i,2] < sk2)
  ch.norm.left[i] <- as.integer(sk2 < ch.ci.norm[i,1] )
  ch.basic.left[i] <- as.integer(sk2 < ch.ci.basic[i,1])
  ch.perc.left[i] <- as.integer(sk2 < ch.ci.perc[i,1])
}
ch.cover <- c(mean(ch.norm.indic),mean(ch.basic.indic),mean(ch.perc.indic)) 
ch.misleft <- c(mean(ch.norm.left),mean(ch.basic.left),mean(ch.perc.left)) 
ch.misright <- c(mean(ch.norm.right),mean(ch.basic.right),mean(ch.perc.right)) 
ch.cover #0.719 0.702 0.71
ch.misright #0.26 0.264 0.29
ch.misleft #0.021 0.034    0

## ----echo=FALSE,eval=TRUE-----------------------------------------------------
f <- data.frame(n.cover,ch.cover,row.names=c("Standard","Basic","Percentile"))
knitr::kable(f)

## ----echo=FALSE,eval=TRUE-----------------------------------------------------
f <- data.frame(n.misright,n.misleft,row.names=c("Standard","Basic","Percentile"))
knitr::kable(f)

## ----echo=FALSE,eval=TRUE-----------------------------------------------------
f <- data.frame(ch.misright,ch.misleft,row.names=c("Standard","Basic","Percentile"))
knitr::kable(f)

## -----------------------------------------------------------------------------
library(bootstrap)
n <- nrow(scor)
cov1 <- (n-1)/n*cov(scor)
value1 <- unlist(eigen(cov1, symmetric=T, only.values=T))
theta.hat <- max(value1)/sum(value1)

theta.jack <- numeric(n)
for (i in 1:n){
  cov2 <- (n-1)/n*cov(scor[-i,])
  value2 <- unlist(eigen(cov2, symmetric=T, only.values=T))
  theta.jack[i] <- max(value2)/sum(value2)
}

bias.jack <- (n-1)*(mean(theta.jack) - theta.hat)
se.jack <- sqrt( (n-1) * mean( (theta.jack - mean(theta.jack))^2 ))
bias.jack <- round(bias.jack,5)
se.jack <- round(se.jack,4)
print(c(bias = bias.jack,se = se.jack))

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
modl <- c("Linear","Quadratic","Exponential","CubicPolynomial")

which.best <- function(dat){
  n <- length(dat[,1])
  
  y1 <- y2 <- y3 <- y4 <- numeric(n)
  for (i in 1:n){
    x <- dat[-i,1]
    y <- dat[-i,2]
    x0 <- dat[i,1]
    
    cv1 <- lm( y ~ x )
    a1 <- cv1$coef
    y1[i] <- a1[1] + a1[2]*x0
    
    cv2 <- lm( y ~ x + I(x^2))
    a2 <- cv2$coef
    y2[i] <- a2[1] + a2[2]*x0 + a2[3]*x0^2
    
    cv3 <- lm( log(y) ~ x)
    a3 <- cv3$coef
    y3[i] <- exp( a3[1] + a3[2]*x0 )
   
    cv4 <- lm( y ~ poly(x,degree=3))
    a4 <- cv4$coef
    y4[i] <- a4[1] + a4[2]*x0 + a4[3]*x0^2 + a4[4]*x0^3
  }
  y <- cbind(y1,y2,y3,y4)
  
  Y <- dat[,2]
  e <- numeric(4)
  for (i in 1:4)e[i] <- mean((Y-y[,i])^2)
  
  SST <- sum((Y - mean(Y))^2)
  SSR <- numeric(4)
  for (i in 1:4)SSR[i] <- sum((Y-y[,i])^2)
  k <- c(1,2,1,3)
  R2 <- SSR/SST
  adjR2 <- (1-(1-R2))*(n-1)/(n-k-1)
  
  num <- c(which.min(e),which.max(adjR2))
  return(c(CrossValidation=modl[num[1]],adjustedR2=modl[num[2]]))
}

which.best(ironslag)

## -----------------------------------------------------------------------------
library(GeneralizedHyperbolic)
rw.Metropolis <- function(x0,N,sigma){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N){
    y <- rnorm(1,x[i-1],sigma)
    if (u[i] <= dskewlap(y) / dskewlap(x[i-1]))
      x[i] <- y
    else{
      x[i] <- x[i-1]
      k <- k+1
      
    }
  }
  return(list(x=x,k=k))
}

sigma <- c(.05,.1,seq(.2,1,.4),seq(2,6,1),10)
n <- length(sigma)
N <- 2000
x0 <- rnorm(1)
rw <- matrix(nrow=N,ncol=n)
Num_of_Reject <- numeric(n)
for (i in 1:n){
 rw[,i] <- rw.Metropolis(x0,N,sigma[i])$x
 Num_of_Reject[i] <- rw.Metropolis(x0,N,sigma[i])$k
}

## -----------------------------------------------------------------------------
Accept_Rate <- round(1-Num_of_Reject/N,3)
knitr::kable(rbind(Num_of_Reject,Accept_Rate),formate="html",col.names=c("0.05","0.1","0.2","0.6","1","2","3","4","5","6","10"))

refline <- qskewlap(c(.025,.975))
for (i in 1:n){
  plot(1:N,rw[,i],type="l",xlab=bquote(sigma==.(sigma[i])),ylab="X")
  abline(h=refline)
}


## -----------------------------------------------------------------------------
isTRUE(exp(log(.001))==.001)
isTRUE(log(exp(.001))==.001)
isTRUE(log(exp(.001))==exp(log(.001)))

isTRUE(all.equal(log(exp(.001)),.001))
isTRUE(all.equal(exp(log(.001)),.001))
isTRUE(all.equal(exp(log(.001)),log(exp(.001))))

## -----------------------------------------------------------------------------
#####11.4
K <- c(4:25,100,500,1e3)
UP <- sqrt(K)
n <- length(K)
A <- numeric(n)

for (i in 1:n){
  k <- K[i]
  f <- function(a){
   pt(sqrt(a^2*(k-1)/(k-a^2)),k-1) - pt(sqrt(a^2*k/(k+1-a^2)),k)
  }
  a <- UP[i]/2+1
  out <- uniroot(f,c(0.5,a))
  A[i] <- out$root
}


## -----------------------------------------------------------------------------
### 11.5
K1 <- c(4:25,100)
n1 <- length(K1)
B <- numeric(n1)

ck <- function(k,a)sqrt((k*(a^2))/(k+1-a^2))
integ_f1 <- function(u)(1+(u^2)/(k-1))^(-k/2)
integ_f2 <- function(u)(1+(u^2)/k)^((-1-k)/2)

f <- function(a){
gamma(k/2)/(gamma((k-1)/2)*sqrt(k-1))*integrate(integ_f1,0,ck(k-1,a))$value-gamma((k+1)/2)/(gamma(k/2)*sqrt(k))*integrate(integ_f2,0,ck(k,a))$value
}

for (i in 1:n1){
  k <- K1[i]
  B[i] <- uniroot(f,lower=0.05,upper=(sqrt(k)/2+1))$root
}

## -----------------------------------------------------------------------------
A <- round(A[1:n1],4)
B <- round(B,4)
err <- A-B

f <- data.frame(rbind(A,B,err),row.names=c("11.4","11.5","error"))
knitr::kable(f)

## -----------------------------------------------------------------------------
library(rootSolve)
r <- 1e3
n1 <- 28
n2 <- 24
n3 <- 41
n4 <- 70
p0 <- .1
q0 <- .1
L <- c(.1,.1)
tol <- .Machine$double.eps^0.5
L.old <- L+1
logML <- numeric(r)
for(j in 1:r){
  logML[j] <- 2*p0*n1*log(p0)/(2-p0-2*q0)+2*q0*n2*log(q0)/(2-q0-2*p0)+2*n3*log(1-p0-q0)+n1*(2-2*p0-2*q0)*log(2*p0*(1-p0-q0))/(2-p0-2*q0)+n2*(2-2*p0-2*q0)*log(2*q0*(1-p0-q0))/(2-q0-2*p0)+n4*log(2*p0*q0)
  g <- function(x){
    F1 <- 2*p0*n1/((2-p0-2*q0)*x[1])-2*n3/(1-x[1]-x[2])+n1*(2-2*p0-2*q0)*(1-2*x[1]-x[2])/((2-p0-2*q0)*x[1]*(1-x[1]-x[2]))-n2*(2-2*p0-2*q0)/((2-q0-2*p0)*(1-x[1]-x[2]))+n4/x[1]
    F2 <- 2*q0*n2/((2-q0-2*p0)*x[2])-2*n3/(1-x[1]-x[2])-n1*(2-2*p0-2*q0)/((2-p0-2*q0)*(1-x[1]-x[2]))+n2*(2-2*p0-2*q0)*(1-2*x[2]-x[1])/((2-q0-2*p0)*x[2]*(1-x[1]-x[2]))+n4/x[2]
    c(F1=F1,F2=F2)
  }
  ss=multiroot(f=g,star=c(.1,.1))
  L=ss$root
  if (sum(abs(L-L.old)/L.old)<tol) break
  L.old <- L
}
L.old
plot(logML,type = "l")


## -----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

#lapply
##version 1
lm_result1 <- lapply(formulas,function(x)lm(x))
##version 2
lm_result1 <- lapply(formulas,lm,data=mtcars)

#for
lm_result2 <- vector("list",length(formulas))
for (i in seq_along(formulas))lm_result2[[i]] <- lm(formulas[[i]])

#coefficients
l_int <- round(sapply(lm_result1,function(x)x$coefficients[1]),3)
l_coef <- round(sapply(lm_result1,function(x)x$coefficients[2]),3)
f_int <- round(sapply(lm_result2,function(x)x$coefficients[1]),3)
f_coef <- round(sapply(lm_result2,function(x)x$coefficients[2]),3)


## -----------------------------------------------------------------------------
attach(mtcars)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

#for
results1 <- vector("list",length(bootstraps))
for(i in seq_along(bootstraps)){
  mt <- bootstraps[[i]]
  results1[[i]] <- lm(mt[,1]~mt[,3])
}

#lapply
results2 <- lapply(bootstraps,function(x)lm(x[,1]~x[,3]))

#coefficients
f_intercpt <- f_coef <- numeric(length(results1))
for(i in seq_along(bootstraps)){
  f_intercpt[i] <- results1[[i]]$coefficients[1]
  f_coef[i] <- results1[[i]]$coefficients[2]
}
round(f_intercpt,3)
round(f_coef,3)

## -----------------------------------------------------------------------------
#without anonymous function
results4 <- lapply(bootstraps,lm,formula=mpg~disp)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

#T3
t3 <- sapply(lm_result1,rsq)
round(t3,3)

## -----------------------------------------------------------------------------
#T4
t4 <- sapply(results1,rsq)
round(t4,3)

## -----------------------------------------------------------------------------
trials <- replicate(100,t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
p_result1 <- sapply(trials,function(x)x$p.value)
round(p_result1[1:10],3)

## -----------------------------------------------------------------------------
p_result2 <- sapply(trials,"[[","p.value")
round(p_result2[1:10],3)

## -----------------------------------------------------------------------------
library(parallel)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(getOption("cl.cores", cores))
parSapply(cl, 1:500, get("+"), 3)[1:10]

## -----------------------------------------------------------------------------
lapply2 <- function(x, f, ...) {
  out <- vector("list", length(x))
  for (i in seq_along(x))out[[i]] <- f(x[[i]], ...)
  out
}

## -----------------------------------------------------------------------------
lapply3 <- function(x, f, ...) {
  out <- vector("list", length(x))
  for (i in sample(seq_along(x)))out[[i]] <- f(x[[i]], ...)
  out
}

## -----------------------------------------------------------------------------
x <- 1:10
y <- rep("a",10)
df1 <- data.frame(x,y)
df2 <- data.frame(y,x)
vapply(df1,class,FUN.VALUE=character(1))
vapply(df2,class,FUN.VALUE=character(1))

## -----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
## r
library(GeneralizedHyperbolic)
rw.Metropolis <- function(x0,N,sigma){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N){
    y <- rnorm(1,x[i-1],sigma)
    if (u[i] <= dskewlap(y) / dskewlap(x[i-1]))
      x[i] <- y
    else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x=x,k=k))
}

## -----------------------------------------------------------------------------
## C++
library(Rcpp)

### generate random numbers
cppFunction('NumericVector rwMetropolisC(double x0, int N, double sigma){
  NumericVector xC(N);
  xC[0] = x0;
  NumericVector u = runif(N);
  int k = 0;
  for(int i = 1; i < N; i++) {
    double y = (rnorm(1,xC[i-1],sigma))[0];
    double res = exp(abs(xC[i-1]) - abs(y));
    if(u[i] < res){
      xC[i] = y;
    }
    else{
      xC[i] = xC[i-1];
      k ++;}
    };
  return xC;}')

### generate k
cppFunction('int rw_k_C(double x0, int N, double sigma){
  NumericVector xC(N);
  xC[0] = x0;
  NumericVector u = runif(N);
  int k = 0;
  for(int i = 1; i < N; i++) {
    double y = (rnorm(1,xC[i-1],sigma))[0];
    double res = exp(abs(xC[i-1]) - abs(y));
    if(u[i] < res){
      xC[i] = y;
    }
    else{
      xC[i] = xC[i-1];
      k ++;}
    };
  return k;}')

## -----------------------------------------------------------------------------
rw.Metropolis <- function(x0,N,sigma){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N){
    y <- rnorm(1,x[i-1],sigma)
    if (u[i] <= dskewlap(y) / dskewlap(x[i-1]))
      x[i] <- y
    else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x=x,k=k))
}

sigma <- c(.1,.5,2,5,10)
n <- length(sigma)
N <- 2000
x0 <- rnorm(1)
rwM <- matrix(nrow=N,ncol=n)
rwM_C <- matrix(nrow=N,ncol=n)
Num_of_Reject <- numeric(n)
Num_of_Reject_C <- numeric(n)
for (i in 1:n){
 rwM[,i] <- rw.Metropolis(x0,N,sigma[i])$x
 Num_of_Reject[i] <- rw.Metropolis(x0,N,sigma[i])$k
 rwM_C[,i] <- rwMetropolisC(x0,N,sigma[i])
 Num_of_Reject_C[i] <- rw_k_C(x0,N,sigma[i])
}
Accept_rate_R <- round(1-Num_of_Reject/N,3)
Accept_rate_Cpp <- round(1-Num_of_Reject_C/N,3)

#compare random numbers
qqplot(rwM, rwM_C, xlab = deparse(substitute(rwM)),ylab = deparse(substitute(rwM_C)),main="rwMetropolis with R and Cpp",col="blue")
x <- seq(-6,6,.05)
lines(x,x,col="red")

#compare accept rate
knitr::kable(rbind(Accept_rate_R,Accept_rate_Cpp),formate="html",col.names=c("0.1","0.5","2","5","10"))


## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- list(n)
for (i in 1:n){
 times <- microbenchmark(rwM=rw.Metropolis(x0,N,sigma[i])$x,
 rwM_C=rwMetropolisC(x0,N,sigma[i]))
 ts[[i]] <- summary(times)[,c(1,3,5,6)]
}
ts

