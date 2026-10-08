########################################
#### Author: Artur J. Lemonte 
#### Office: UFRN, Natal/RN, Brazil
#### Version: 3.14 (October 08, 2026)
#### Time: 10:02
########################################

rm(list=ls(all=T))

##############################
#### Loading packages
require(bbreg)
require(fBasics)

#### Auxiliar functions
D2Q_Obs_Fisher_bet = function(theta,z,x,v,link.mean,link.precision){
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  nlam = length(lam)
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap) # mean parameter.
  phi = link_precision$linkinv(v%*%lam) # phi precision parameter.

  dmudeta = link_mean$mu.eta(x%*%kap)
  dphideta = link_precision$mu.eta(v%*%lam)
  d2mu = d2mudeta2(link.mean,mu)
  d2phi = d2phideta2(link.precision,phi)

  auxKK1 = ( trigamma(mu*phi)+trigamma((1-mu)*phi) ) * phi^2 * dmudeta^2
  auxKK2 = phi*(log(z)-log(1-z)-digamma(mu*phi)+digamma((1-mu)*phi)) * d2mu
  KK = diag(c(auxKK1 - auxKK2))

  auxLL1 = ( (mu^2)*trigamma(mu*phi)+((1-mu)^2)*trigamma((1-mu)*phi)) * dphideta^2
  auxLL2 = ( mu*(log(z)-log(1-z))+digamma(phi)+log(1-z)-mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)) * d2phi
  LL = diag(c(auxLL1 + auxLL2))

  auxKL = (-log(z)+log(1-z)+digamma(mu*phi)-digamma((1-mu)*phi)+mu*phi*trigamma(mu*phi) -(1-mu)*phi*trigamma((1-mu)*phi))*dmudeta*dphideta
  KL = diag(c(auxKL))

  D2QKappa = (t(x)%*%KK)%*%x
  D2QKL = (t(x)%*%KL)%*%v
  D2QKappa = cbind(D2QKappa, D2QKL)
  D2QLambda = (t(v)%*%LL)%*%v
  D2QLambda = cbind(t(D2QKL), D2QLambda)
  D2Q = rbind(D2QKappa, D2QLambda)

  return(D2Q)
}
d2mudeta2 = function(link.mean, mu){
  d2mu = switch(link.mean,
                "logit" = {
                  mu*(1-mu)*(1-2*mu)
                },
                "probit" = {
                  (-mu/sqrt(2*pi)) * exp(-mu^2/2)
                },
                "cloglog" = {
                  -(1-mu)*log(1-mu)*(1+log(1-mu))
                },
                "cauchit" = {
                  (-2/pi) * tan(pi*(mu - 1/2))/ ( (1 + tan(pi*(mu - 1/2))^2 )^2 ) 
                }
  )
  return(d2mu)
}
d2phideta2 = function(link.precision, phi){
  d2phi = switch(link.precision,
                "identity" = {
                  0
                },
                "log" = {
                  phi
                },
                "sqrt" = {
                  2
                },
                "inverse" = {
                  2*phi^3
                },
                "1/precision^2" = {
                  3*phi^5/4
                }
  )
  return(d2phi)
}


##############################
dados <- read.table("fooddata.txt", head=T)
attach(dados)
head(dados)

vy <- dados[,3]/100 ##= response variable on (0,1)
income <- dados[,1]
persons <- dados[,2]


#########################################
#### Fitting the Beta regression model
#########################################

fit.model <- bbreg(vy ~ income + persons, link.precision = c("log"), model="beta")#, residual = "quantile", prob = 0.95, envelope = 1000)
summary(fit.model)
#plot(fit.model)


##############################
#### Global influence
##############################

n <- NROW(vy)
mX <- model.matrix(~ income + persons)
mV <- matrix(1,n,1)
p <- NCOL(mX)
q <- NCOL(mV) 
vkappa <- coef(fit.model)[1:p]
vlambda <- coef(fit.model)[(p+1):(p+q)]
veta <- mX%*%vkappa
vtau <- mV%*%vlambda
vmu <- exp(veta)/(1 + exp(veta)) 
vphi <- exp(vtau)
vdotg <- 1/(vmu*(1-vmu))
vdoth <- 1/vphi
mddQ <- -D2Q_Obs_Fisher_bet(coef(fit.model),vy,mX,mV,link.mean="logit",link.precision="log")
mddQkk <- mddQ[1:p,1:p]
mddQkl <- mddQ[1:p,(p+1):(p+q)]
mddQll <- mddQ[(p+1):(p+q),(p+1):(p+q)]
va <- vphi*( log(vy/(1-vy)) - digamma(vmu*vphi) + digamma((1-vmu)*vphi) )*(1/vdotg)
vb <- ( vmu*log(vy/(1-vy)) + digamma(vphi) + log(1-vy) - vmu*digamma(vmu*vphi) - 
              (1-vmu)*digamma((1-vmu)*vphi) )*(1/vdoth)
mDdiag <- matrix(0,n,p+q)
for (i in 1:n){
 mDdiag[i,] <- diag( rbind(cbind((va[i]^2)*mX[i,]%*%t(mX[i,]),va[i]*vb[i]*mX[i,]%*%t(mV[i,])),
                           cbind(va[i]*vb[i]*mV[i,]%*%t(mX[i,]),(vb[i]^2)*mV[i,]%*%t(mV[i,])) )%*%solve(-mddQ) )
}

## theta
GDi <- apply(mDdiag,1,sum)
plot(GDi, pch=20,ylab="global influence",xlab="index") 
identify(GDi)

## kappa
GDi.k <- apply(mDdiag[,1:p],1,sum)
plot(GDi.k, pch=20,ylab="global influence",xlab="index")
mtext("(a)", cex=1.5, line=0.5)
identify(GDi.k)

## lambda
GDi.l <- mDdiag[,(p+1):(p+q)]  
plot(GDi.l, pch=20,ylab="global influence",xlab="index")
mtext("(b)", cex=1.5, line=0.5)
identify(GDi.l)


##############################
#### Local influence
##############################

#### Case-weight perturbation
mA <- diag(as.vector(va))
mB <- diag(as.vector(vb))
mDeltak.cwp <- t(mX)%*%mA
mDeltal.cwp <- t(mV)%*%mB
mDelta.cwp <- rbind(mDeltak.cwp, mDeltal.cwp)
mddQinv <- solve(mddQ)
mddQinv.k <- rbind( cbind(solve(mddQkk),matrix(0,p,q)),cbind(matrix(0,q,p),matrix(0,q,q)) )
mddQinv.l <- rbind( cbind(matrix(0,p,p),matrix(0,p,q)),cbind(matrix(0,q,p),solve(mddQll)) )

## theta 
vBh.cwp <- numeric()
for (i in 1:n){
 vBh.cwp[i] <- 2*t(mDelta.cwp)[i,]%*%mddQinv%*%mDelta.cwp[,i]
}
Bi.cwp <- abs(vBh.cwp)/tr(-2*t(mDelta.cwp)%*%mddQinv%*%mDelta.cwp)
plot(Bi.cwp, pch=20,ylab="local influence",xlab="index")
identify(Bi.cwp)

## kappa
vBh.k.cwp <- numeric()
for (i in 1:n){
 vBh.k.cwp[i] <- 2*t(mDelta.cwp)[i,]%*%(mddQinv - mddQinv.l)%*%mDelta.cwp[,i]
}
Bi.k.cwp <- abs(vBh.k.cwp)/tr(-2*t(mDelta.cwp)%*%(mddQinv - mddQinv.l)%*%mDelta.cwp)
plot(Bi.k.cwp, pch=20,ylab="local influence",xlab="index") 
mtext("(a)", cex=1.5, line=0.5)
identify(Bi.k.cwp)

## lambda
vBh.l.cwp <- numeric()
for (i in 1:n){
 vBh.l.cwp[i] <- 2*t(mDelta.cwp)[i,]%*%(mddQinv - mddQinv.k)%*%mDelta.cwp[,i]
}
Bi.l.cwp <- abs(vBh.l.cwp)/tr(-2*t(mDelta.cwp)%*%(mddQinv - mddQinv.k)%*%mDelta.cwp)
plot(Bi.l.cwp, pch=20,ylab="local influence",xlab="index") 
mtext("(b)", cex=1.5, line=0.5)
identify(Bi.l.cwp)


#### Response perturbation
Sy <- sd(vy)
mT1 <- diag(as.vector(1/vdotg))
mT2 <- diag(as.vector(1/vdoth))
mPhi <- diag(as.vector(vphi))
mU <- diag(as.vector(vmu))
mMz <- diag(as.vector(1/(vy*(1-vy))))
mDeltak.resp <- Sy * t(mX)%*%mT1%*%mU%*%mMz
mDeltal.resp <- Sy * t(mV)%*%mT2%*%mPhi%*%mMz
mDelta.resp <- rbind(mDeltak.resp, mDeltal.resp)

## theta
vBh.resp <- numeric()
for (i in 1:n){
 vBh.resp[i] <- 2*t(mDelta.resp)[i,]%*%mddQinv%*%mDelta.resp[,i]
}
Bi.resp <- abs(vBh.resp)/tr(-2*t(mDelta.resp)%*%mddQinv%*%mDelta.resp)
plot(Bi.resp, pch=20,ylab="local influence",xlab="index") 
identify(Bi.resp)

## kappa
vBh.k.resp <- numeric()
for (i in 1:n){
 vBh.k.resp[i] <- 2*t(mDelta.resp)[i,]%*%(mddQinv - mddQinv.l)%*%mDelta.resp[,i]
}
Bi.k.resp <- abs(vBh.k.resp)/tr(-2*t(mDelta.resp)%*%(mddQinv - mddQinv.l)%*%mDelta.resp)
plot(Bi.k.resp, pch=20,ylab="local influence",xlab="index") 
mtext("(a)", cex=1.5, line=0.5)
identify(Bi.k.resp)

## lambda
vBh.l.resp <- numeric()
for (i in 1:n){
 vBh.l.resp[i] <- 2*t(mDelta.resp)[i,]%*%(mddQinv - mddQinv.k)%*%mDelta.resp[,i]
}
Bi.l.resp <- abs(vBh.l.resp)/tr(-2*t(mDelta.resp)%*%(mddQinv - mddQinv.k)%*%mDelta.resp)
plot(Bi.l.resp, pch=20,ylab="local influence",xlab="index") 
mtext("(b)", cex=1.5, line=0.5)
identify(Bi.l.resp)


#### Mean covariate perturbation: income covariate
Sj <- sd(income)
vdelta.cov <- c(0,1,0)
kappa.p <- vkappa[2]
vg1line <- 1/(vmu*(1-vmu))
vg1lineline <- - (1 - 2*vmu)/( (vmu*(1 - vmu))^2 )
vg2line <- 1/vphi                   
vg2lineline <- -1/(vphi^2)
vai <- trigamma(vmu*vphi) + trigamma((1-vmu)*vphi)
vbi <- (vmu^2)*trigamma(vmu*vphi) + ((1-vmu)^2)*trigamma((1-vmu)*vphi)
vci <- vmu*log(vy/(1-vy)) + digamma(vphi) + log(1-vy) - vmu*digamma(vmu*vphi) - 
       (1-vmu)*digamma((1-vmu)*vphi)
vGkk <- vphi*( vphi*vai*(1/vg1line^2) - (log(vy/(1-vy)) - digamma(vmu*vphi) + 
                              digamma((1-vmu)*vphi))*(-vg1lineline/(vg1line^3)) )
vGkl <- ( vphi*( vmu*vai - trigamma((1-vmu)*vphi) ) - log(vy/(1-vy)) + digamma(vmu*vphi) - 
                              digamma((1-vmu)*vphi) )*(1/vg1line)*(1/vg2line)
vGll <- vbi*(1/vg2line^2) - vci*(-vg2lineline/(vg2line^3))
mGkk <- diag( as.vector(vGkk) )
mGkl <- diag( as.vector(vGkl) )
mGll <- diag( as.vector(vGll) )
mDeltak.cov <- -Sj * kappa.p*t(mX)%*%mGkk
mDeltal.cov <- -Sj * kappa.p*t(mV)%*%t(mGkl)
mDelta.cov <- rbind(mDeltak.cov, mDeltal.cov)

## theta 
vBh.cov <- numeric()
for (i in 1:n){
 vBh.cov[i] <- 2*t(mDelta.cov)[i,]%*%mddQinv%*%mDelta.cov[,i]
}
Bi.cov <- abs(vBh.cov)/tr(-2*t(mDelta.cov)%*%mddQinv%*%mDelta.cov)
plot(Bi.cov, pch=20,ylab="local influence",xlab="index") 
identify(Bi.cov)

## kappa
vBh.k.cov <- numeric()
for (i in 1:n){
 vBh.k.cov[i] <- 2*t(mDelta.cov)[i,]%*%(mddQinv - mddQinv.l)%*%mDelta.cov[,i]
}
Bi.k.cov <- abs(vBh.k.cov)/tr(-2*t(mDelta.cov)%*%(mddQinv - mddQinv.l)%*%mDelta.cov)
plot(Bi.k.cov, pch=20,ylab="local influence",xlab="index") 
mtext("(a)", cex=1.5, line=0.5)
identify(Bi.k.cov)

## lambda
vBh.l.cov <- numeric()
for (i in 1:n){
 vBh.l.cov[i] <- 2*t(mDelta.cov)[i,]%*%(mddQinv - mddQinv.k)%*%mDelta.cov[,i]
}
Bi.l.cov <- abs(vBh.l.cov)/tr(-2*t(mDelta.cov)%*%(mddQinv - mddQinv.k)%*%mDelta.cov)
plot(Bi.l.cov, pch=20,ylab="local influence",xlab="index") 
mtext("(b)", cex=1.5, line=0.5)
identify(Bi.l.cov)



##############################
#### Confirmatory analysis
##############################

## obs #4
fit.model.s04 <- bbreg(vy[-c(4)] ~ income[-c(4)] + persons[-c(4)], link.precision = c("log"), model="beta")
summary(fit.model.s04)
round(((coef(fit.model)-coef(fit.model.s04))/coef(fit.model)) * 100, 2)

## obs #9
fit.model.s09 <- bbreg(vy[-c(9)] ~ income[-c(9)] + persons[-c(9)], link.precision = c("log"), model="beta")
summary(fit.model.s09)
round(((coef(fit.model)-coef(fit.model.s09))/coef(fit.model)) * 100, 2) 

## obs #10
fit.model.s10 <- bbreg(vy[-c(10)] ~ income[-c(10)] + persons[-c(10)], link.precision = c("log"), model="beta")
summary(fit.model.s10)
round(((coef(fit.model)-coef(fit.model.s10))/coef(fit.model)) * 100, 2)

## obs #11
fit.model.s11 <- bbreg(vy[-c(11)] ~ income[-c(11)] + persons[-c(11)], link.precision = c("log"), model="beta")
summary(fit.model.s11)
round(((coef(fit.model)-coef(fit.model.s11))/coef(fit.model)) * 100, 2) 

## obs #20
fit.model.s20 <- bbreg(vy[-c(20)] ~ income[-c(20)] + persons[-c(20)], link.precision = c("log"), model="beta")
summary(fit.model.s20)
round(((coef(fit.model)-coef(fit.model.s20))/coef(fit.model)) * 100, 2) 

## obs #25
fit.model.s25 <- bbreg(vy[-c(25)] ~ income[-c(25)] + persons[-c(25)], link.precision = c("log"), model="beta")
summary(fit.model.s25)
round(((coef(fit.model)-coef(fit.model.s25))/coef(fit.model)) * 100, 2) 


