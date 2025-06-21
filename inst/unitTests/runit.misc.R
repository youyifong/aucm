library("RUnit")
library("aucm")

test.misc <- function() {

RNGkind("Mersenne-Twister", "Inversion")
#RNGkind("Marsaglia-Multicarry", "Kinderman-Ramage") 

tol=1e-6

dat = sim.dat.1(n=200,seed=1)
checkEqualsNumeric(mean(as.matrix(dat)), 0.1411382, tol=tol)


set.seed(1)
x0 <- matrix(rnorm(100,1))
y  <- as.numeric(runif(100)>0.5)        # numeric(runif(100)>0.5)
dat=data.frame(y=y, x=x0)

fit=rlogit(y~x, dat)
checkEqualsNumeric(coef(fit), c(-0.51703584,  0.02472164), tol=tol)

fit=rlogit(y~-1+x, dat)
checkEqualsNumeric(coef(fit), c(-0.51703584,  0.02472164), tol=tol)


}
