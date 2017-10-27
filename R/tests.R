library("mvtnorm")
set.seed(1453)
n=1e4
sigmax=rbind(c(0.6,0.3,0.4),c(0.3,0.5,0.2),c(0.4,0.2,0.3))
#sigmax=diag(c(0.6,0.5,0.3))
x=rmvnorm(n,sigma=sigmax)
sigmay=rbind(c(0.4,0.2,0.5),c(0.2,0.5,0.4),c(0.5,0.4,0.7))
#sigmay=diag(c(0.4,0.5,0.7))
y=x+rmvnorm(n,sigma=sigmay)
# now change x to random predictions
#x=rmvnorm(n,sigma=sigmax)

#x=rmvnorm(n,sigma=rbind(c(0.6,0,0),c(0,0.5,0),c(0,0,0.5)))
#y=x+rmvnorm(n,sigma=rbind(c(0.4,0,0),c(0,0.5,0),c(0,0,0.5)))

d=cbind(y[,1]>qnorm(0.8),y[,2]>qnorm(0.7),y[,3]>qnorm(0.3))
risk=pnorm(qnorm(0.8),mean=x[,1],sd=sqrt(0.4),lower=F)
risk=cbind(risk,pnorm(qnorm(0.7),mean=x[,2],sd=sqrt(0.5),lower=F))
risk=cbind(risk,pnorm(qnorm(0.3),mean=x[,3],sd=sqrt(0.7),lower=F))
w=which(risk[,1]>0.15 & risk[,2]>0.35 & risk[,3]>0.6)
