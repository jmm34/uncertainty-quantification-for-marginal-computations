library(parallel)
library(doParallel)
library(doRNG)
library(mvtnorm)

ncores <- 50
cl <- makeCluster(ncores)
registerDoParallel(cl)
randomseedpara <- 1974

p <- 10
y <- rep(5,p)
target <- -sum(dnorm(y,0,sqrt(2),log=TRUE))
m <- rep(3,p)
sig <- diag(rep(0.8,p))


mqprime <- function(const,Nsim,lh)
{
  -Nsim+sum(exp(target+lh)/(exp(target+lh)+exp(target-const)))
}

mqprimew <- function(const,Nsim,lh,w)
{
  -sum(w[1:Nsim])+sum(w*exp(target+lh)/(exp(target+lh)+exp(target-const)))
}

Nsim <- 10^5
Nwlb <- 1000
Nrepcouv <- 1000
wlb <- rep(0,Nwlb)
# res <- rep(0,Nrepcouv)
# tp <- txtProgressBar(min = 1, max = Nrepcouv, style = 3, char = "*")
# for (j in 1:Nrepcouv)
res <- foreach(j=1:Nrepcouv, .combine="c",.options.RNG = randomseedpara, .packages="mvtnorm") %dorng% {
thetapost <- rmvnorm(Nsim,mean=y/2,sigma=diag(rep(1/2,p)))
thetag <- rmvnorm(Nsim,mean=m,sigma=sig)
zeta <- rbind(thetapost,thetag)
lh  <- dmvnorm(zeta,mean=y,log=TRUE)+dmvnorm(zeta,log=TRUE)-
  dmvnorm(zeta,mean=m,sigma=sig,log=TRUE)
for (i in 1:Nwlb)
{
    w1 <- rexp(Nsim)
    w1 <- w1/sum(w1)
    w2 <- rexp(Nsim)
    w2 <- w2/sum(w2)
    w <- c(w1,w2)
    wlb[i] <- uniroot(mqprimew,Nsim=Nsim,lh=lh,w=w,
                      c(target+10,target-10),tol=.Machine$double.eps^0.5)$root
  }
  quant <- quantile(wlb,c(0.025,0.975),names=FALSE)
  return(target >= quant[1] & target <=quant[2])
#  res[j] <- target >= quant[1] & target <=quant[2]
#  setTxtProgressBar(tp, j)
}
stopCluster(cl)
save(res,file="res10.RData")
