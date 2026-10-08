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
  auxLL2 = ( mu*(log(z)-log(1-z))+digamma(phi)+log(1-z)-
			mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)) * d2phi
  LL = diag(c(auxLL1 + auxLL2))

  auxKL = (-log(z)+log(1-z)+digamma(mu*phi)-digamma((1-mu)*phi)+mu*phi*trigamma(mu*phi) -
			(1-mu)*phi*trigamma((1-mu)*phi))*dmudeta*dphideta
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
#### DataSet: Water Quality Index 
#### Bonat, W.H., Petterle, R.R., Hinde, J., Demétrio, C.G.B. (2019). 
#### 	 Flexible quasi-beta regression models for continuous bounded data. 
#### 	 Statistical Modelling 19, 617-633.

dataset <- read.table("iqadata.txt", header = TRUE)
attach(dataset)
head(dataset)
vy <- IQA ## response variable on (0,1)
TRIM <- factor(TRIM)
LOCAL <- factor(LOCAL, levels=c("MONT","RESER","JUSA"))
levels(LOCAL) <- c("Upstream","Reservoir","Downstream")

#########################################
#### Fitting the Beta regression model
#########################################

fit.model <- bbreg(vy ~ LOCAL + TRIM | LOCAL + TRIM, model="beta")#, residual = "quantile", prob = 0.95, envelope = 1000)
summary(fit.model)
#plot(fit.model)


##############################
#### Global influence
##############################

mX <- model.matrix(~ LOCAL + TRIM)
mV <- model.matrix(~ LOCAL + TRIM)
p <- NCOL(mX)
q <- NCOL(mV) 
n <- NROW(vy)
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
GDi.l <- apply(mDdiag[,(p+1):(p+q)],1,sum)
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
#abline(h=mean(Bi.l.cwp) + 3*sd(Bi.l.cwp),lty=2)
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
plot(Bi.resp, pch=20,ylab="local influence",xlab="index")# , las=1,ylim=c(0,1))
identify(Bi.resp)

## kappa
vBh.k.resp <- numeric()
for (i in 1:n){
 vBh.k.resp[i] <- 2*t(mDelta.resp)[i,]%*%(mddQinv - mddQinv.l)%*%mDelta.resp[,i]
}
Bi.k.resp <- abs(vBh.k.resp)/tr(-2*t(mDelta.resp)%*%(mddQinv - mddQinv.l)%*%mDelta.resp)
plot(Bi.k.resp, pch=20,ylab="local influence",xlab="index")# , las=1,ylim=c(0,1))
mtext("(a)", cex=1.5, line=0.5)
identify(Bi.k.resp)

## lambda
vBh.l.resp <- numeric()
for (i in 1:n){
 vBh.l.resp[i] <- 2*t(mDelta.resp)[i,]%*%(mddQinv - mddQinv.k)%*%mDelta.resp[,i]
}
Bi.l.resp <- abs(vBh.l.resp)/tr(-2*t(mDelta.resp)%*%(mddQinv - mddQinv.k)%*%mDelta.resp)
plot(Bi.l.resp, pch=20,ylab="local influence",xlab="index")# , las=1,ylim=c(0,1))
mtext("(b)", cex=1.5, line=0.5)
identify(Bi.l.resp)



##############################
#### Confirmatory analysis
##############################

## obs #4
fit.model.s04 <- bbreg(vy[-c(4)] ~ LOCAL[-c(4)] + TRIM[-c(4)] | LOCAL[-c(4)] + TRIM[-c(4)], model="beta")
summary(fit.model.s04)
round(((coef(fit.model)-coef(fit.model.s04))/coef(fit.model)) * 100, 2)

## obs #6
fit.model.s06 <- bbreg(vy[-c(6)] ~ LOCAL[-c(6)] + TRIM[-c(6)] | LOCAL[-c(6)] + TRIM[-c(6)], model="beta")
summary(fit.model.s06)
round(((coef(fit.model)-coef(fit.model.s06))/coef(fit.model)) * 100, 2)

## obs #18
fit.model.s18 <- bbreg(vy[-c(18)] ~ LOCAL[-c(18)] + TRIM[-c(18)] | LOCAL[-c(18)] + TRIM[-c(18)], model="beta")
summary(fit.model.s18)
round(((coef(fit.model)-coef(fit.model.s18))/coef(fit.model)) * 100, 2)

## obs #20
fit.model.s20 <- bbreg(vy[-c(20)] ~ LOCAL[-c(20)] + TRIM[-c(20)] | LOCAL[-c(20)] + TRIM[-c(20)], model="beta")
summary(fit.model.s20)
round(((coef(fit.model)-coef(fit.model.s20))/coef(fit.model)) * 100, 2)

## obs #39
fit.model.s39 <- bbreg(vy[-c(39)] ~ LOCAL[-c(39)] + TRIM[-c(39)] | LOCAL[-c(39)] + TRIM[-c(39)], model="beta")
summary(fit.model.s39)
round(((coef(fit.model)-coef(fit.model.s39))/coef(fit.model)) * 100, 2)

## obs #43
fit.model.s43 <- bbreg(vy[-c(43)] ~ LOCAL[-c(43)] + TRIM[-c(43)] | LOCAL[-c(43)] + TRIM[-c(43)], model="beta")
summary(fit.model.s43)
round(((coef(fit.model)-coef(fit.model.s43))/coef(fit.model)) * 100, 2)

## obs #55
fit.model.s55 <- bbreg(vy[-c(55)] ~ LOCAL[-c(55)] + TRIM[-c(55)] | LOCAL[-c(55)] + TRIM[-c(55)], model="beta")
summary(fit.model.s55)
round(((coef(fit.model)-coef(fit.model.s55))/coef(fit.model)) * 100, 2)

## obs #67
fit.model.s67 <- bbreg(vy[-c(67)] ~ LOCAL[-c(67)] + TRIM[-c(67)] | LOCAL[-c(67)] + TRIM[-c(67)], model="beta")
summary(fit.model.s67)
round(((coef(fit.model)-coef(fit.model.s67))/coef(fit.model)) * 100, 2)

## obs #68
fit.model.s68 <- bbreg(vy[-c(68)] ~ LOCAL[-c(68)] + TRIM[-c(68)] | LOCAL[-c(68)] + TRIM[-c(68)], model="beta")
summary(fit.model.s68)
round(((coef(fit.model)-coef(fit.model.s68))/coef(fit.model)) * 100, 2)

## obs #78
fit.model.s78 <- bbreg(vy[-c(78)] ~ LOCAL[-c(78)] + TRIM[-c(78)] | LOCAL[-c(78)] + TRIM[-c(78)], model="beta")
summary(fit.model.s78)
round(((coef(fit.model)-coef(fit.model.s78))/coef(fit.model)) * 100, 2)

## obs #91
fit.model.s91 <- bbreg(vy[-c(91)] ~ LOCAL[-c(91)] + TRIM[-c(91)] | LOCAL[-c(91)] + TRIM[-c(91)], model="beta")
summary(fit.model.s91)
round(((coef(fit.model)-coef(fit.model.s91))/coef(fit.model)) * 100, 2)

## obs #92
fit.model.s92 <- bbreg(vy[-c(92)] ~ LOCAL[-c(92)] + TRIM[-c(92)] | LOCAL[-c(92)] + TRIM[-c(92)], model="beta")
summary(fit.model.s92)
round(((coef(fit.model)-coef(fit.model.s92))/coef(fit.model)) * 100, 2)

## obs #108
fit.model.s108 <- bbreg(vy[-c(108)] ~ LOCAL[-c(108)] + TRIM[-c(108)] | LOCAL[-c(108)] + TRIM[-c(108)], model="beta")
summary(fit.model.s108)
round(((coef(fit.model)-coef(fit.model.s108))/coef(fit.model)) * 100, 2)

## obs #126
fit.model.s126 <- bbreg(vy[-c(126)] ~ LOCAL[-c(126)] + TRIM[-c(126)] | LOCAL[-c(126)] + TRIM[-c(126)], model="beta")
summary(fit.model.s126)
round(((coef(fit.model)-coef(fit.model.s126))/coef(fit.model)) * 100, 2)

## obs #138
fit.model.s138 <- bbreg(vy[-c(138)] ~ LOCAL[-c(138)] + TRIM[-c(138)] | LOCAL[-c(138)] + TRIM[-c(138)], model="beta")
summary(fit.model.s138)
round(((coef(fit.model)-coef(fit.model.s138))/coef(fit.model)) * 100, 2)

## obs #158
fit.model.s158 <- bbreg(vy[-c(158)] ~ LOCAL[-c(158)] + TRIM[-c(158)] | LOCAL[-c(158)] + TRIM[-c(158)], model="beta")
summary(fit.model.s158)
round(((coef(fit.model)-coef(fit.model.s158))/coef(fit.model)) * 100, 2)


