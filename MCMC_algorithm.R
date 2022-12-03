#Metropolis-Hastings MCMC algorithm implementation
#https://bigiskander.github.io/mcmc.html
#author: Sultanov Iskander
# BigIskander@gmail.com

#aprior likelyhood
prior <- function(b, sig_sq, b_TH, alpha, beta, stdp)
{
  dn <- dnorm(b, mean=b_TH, sd=stdp, log=T)
  dg <- dgamma(1/sig_sq, shape=alpha, scale=1/beta, log=T)
  return(sum(dn)+dg)
}

#likelyhood function
likelyhood <- function(b, Y, pY)
{
  st <- sqrt(var(pY))
  dn <- dnorm(Y, mean = pY, sd = st, log = T)
  return(sum(dn))
}

#aposterior likelyhood
aposterior <- function(b, sig_sq, b_TH, alpha, beta, stdp, Y, pY)
{
  return(prior(b, sig_sq, b_TH, alpha, beta, stdp)+likelyhood(b, Y, pY))
}

#generate new varibles
propose <- function(b, st)
{
  return(rnorm(length(b),mean=b,sd=st))
}

#mcmc simulation
#iter - number of iterations
#X, Y - curent data
#b_TH - vector of coeffients from previous data estimation
#alpha, beta - gamma distribution parameters from previous data estimation
#h, L - accuracy variables from previous data estimation
mcmc <- function(iter, X, Y, b_TH, alpha, beta, h, L)
{
  #preliminary estimation 1
  B <- h*L
  C <- matrix(0, length(b_TH), 1)
  for(i in 1:length(b_TH))
  {
    C[i] <- B[i,i]-B[i,-i]%*%solve(B[-i,-i],B[-i,i])
  }
  stdp <- 1/sqrt(C)
  
  #preliminary estimation 2
  iter <- iter*1.1
  b <- b_TH
  pY <- X%*%b
  resid <- Y-pY
  sig_sq <- as.numeric(var(resid))
  output=matrix(0,ncol=ncol(X)+1,nrow=iter)	
  output[1,]=c(sig_sq,b)
  
  prevp <- aposterior(b, sig_sq, b_TH, alpha, beta, stdp, Y, pY)
  
  #mcmc algorithm main cycle
  i <- 2
  while(i<=iter)
  {
    bn <- propose(b, stdp)
    pY <- X%*%bn
    residn <- Y-pY
    sig_sqn <- as.numeric(var(residn))
    ace <- aposterior(bn, sig_sqn, b_TH, alpha, beta, stdp, Y, pY) 
    if(log(runif(1))<(ace-prevp)) {
      b <- bn
      sig_sq <- sig_sqn
      prevp <- ace	
    } 
    output[i,]=c(sig_sq,b)
    i <- i+1
  }
  return(output[-c(1:(iter/11)),])
}

#dencity plot
dence <- function(data)
{
  cc <- ncol(data)
  if(cc>3)
  {
    mm <- round(cc/3,0)
    if(mm<cc/3) mm <- mm+1 
    par(mfrow=c(mm,3))
  } else {
    par(mfrow=c(1,cc))
  }
  #dencity
  for(i in 2:cc)
  {
    hist(data[,i],breaks=150,main=paste0("th",(i-2)),freq=FALSE)
    lines(density(data[,i],bw="nrd",
                  kernel="gaussian"),col="red",lwd=2)
  }
  i <- 1
  hist(data[,i],breaks=150,main="sig_sq",freq=FALSE)
  lines(density(data[,1],bw="nrd",kernel="gaussian"),col="red",lwd=2)
}

#example how you can calculate variables necessary for the algorithm
#and make calculations using this implementation of MCMC algorithm


#load the data
library(foreign)
#file with data saved in stata format (not provided)
mydata <- read.dta("D:/Учеба/ВШЭ Аспирантура/Science/Calc in Stata/470.dta") 
#econometric model for estimation (for example only)
f <- "YS~Leverage+LN_Maturity+LN_Issue" 

# 2014 excluded because of too smal number of observations
#estimate linear regressions
m1 <- lm(f, mydata, year==2007)
m2 <- lm(f, mydata, year==2008)
m3 <- lm(f, mydata, year==2009)
m4 <- lm(f, mydata, year==2010)
m5 <- lm(f, mydata, year==2011)
m6 <- lm(f, mydata, year==2012)
m7 <- lm(f, mydata, year==2013)
#m8 <- lm(f, mydata, year==2014)
  
#get coefficients (theta)
cf1 <- coef(m1)
cf2 <- coef(m2)
cf3 <- coef(m3)
cf4 <- coef(m4)
cf5 <- coef(m5)
cf6 <- coef(m6)
cf7 <- coef(m7)
#cf8 <- coef(m8)

#primary estimations  
hh <- 1/c(var(residuals(m1)),var(residuals(m2)),
          var(residuals(m3)),var(residuals(m4)),
          var(residuals(m5)),var(residuals(m6)),
          var(residuals(m7)))#,var(residuals(m8)))
h <- mean(hh)
st_h <- sqrt(var(hh))
  
#calc the distribution params
alpha <- (h^2)/(st_h^2)
beta <- h/(st_h^2)

ll <- length(cf1)
  
cfl <- matrix(0, ll, 7) #8
th <- matrix(0, ll, 1)
st_th <- matrix(0, ll, 1)
#L matrix calculation
L <- matrix(0, ll, ll)
for(i in 1:ll)
{
  cfl[i,] <- c(cf1[i],cf2[i],cf3[i],cf4[i],cf5[i],
               cf6[i],cf7[i])#,cf8[i])
  th[i] <- mean(cfl[i,])
  st_th[i] <- sqrt(var(cfl[i,]))
  L[i,i] <- (1/(st_th[i]^2))*(beta/(alpha-1))
}
  
#get X and Y
ndata <- subset(mydata, year==2015)
M <- as.matrix(model.frame(formula=f, data=ndata))
Y <- M[,1]
X <- cbind(1,M[,-1])
n <- length(Y)

#estimate using MCMC
iter <- 1000000
data <- mcmc(iter, X, Y, th, alpha, beta, h, L)

#get model parameters from MCMC estimation
th_sim <- matrix(0,ncol(data)-1,1)
for(i in 2:ncol(data)) th_sim[i-1] <- rbind(mean(data[,i]))
print(th_sim)

#plot dencity graphs
dence(data)
