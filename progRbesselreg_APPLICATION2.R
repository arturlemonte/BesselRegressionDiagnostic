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

##############################
dados <- read.table("fooddata.txt", head=T)
attach(dados)
head(dados)

vy <- dados[,3]/100 ##= response variable on (0,1)
income <- dados[,1]
persons <- dados[,2]


#########################################
#### Fitting the Bessel regression model
#########################################

fit.model <- bbreg(vy ~ income + persons, link.precision = c("log"), model="bessel")#, residual = "quantile", prob = 0.95, envelope = 1000)
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
vddotg <- - (1 - 2*vmu)/( (vmu*(1 - vmu))^2 )
vdoth <- 1/vphi
vddoth <- -1/(vphi^2)
vxi <- sqrt( 1 + ((vy-vmu)^2)/(vy*(1-vy)) ) 
vpsi <- (besselK(vphi*vxi, -2)/besselK(vphi*vxi, -1))/(vphi*vxi)
vw1 <- - ( (1-2*vmu)/((vmu^2)*((1-vmu)^2)) + 2/((1-vmu)^2) + 
	     vpsi*(vphi^2)/(vy*(1-vy)) )*(1/vdotg) - ( (1-2*vmu)/(vmu*(1-vmu)) + 
            vpsi*(vphi^2)*(vy-vmu)/(vy*(1-vy)) )*(vddotg/vdotg)
vw2 <- 2*vpsi*vphi*(vy-vmu)/(vy*(1-vy))
vw3 <- - ( 2/(vphi^2) + vpsi*((vmu^2)/vy + ((1-vmu)^2)/(1-vy) ) )*(1/vdoth) - ( 2/vphi + 1 -
            vpsi*vphi*((vmu^2)/vy + ((1-vmu)^2)/(1-vy)) )*(vddoth/vdoth)
mG <- diag(as.vector(1/vdotg))
mH <- diag(as.vector(1/vdoth))
mW1 <- diag(as.vector(vw1))
mW2 <- diag(as.vector(vw2))
mW3 <- diag(as.vector(vw3))
mddQkk <- t(mX)%*%mW1%*%mG%*%mX
mddQkl <- t(mX)%*%mW2%*%mG%*%mH%*%mV
mddQll <- t(mV)%*%mW3%*%mH%*%mV
mddQ <- rbind(cbind(mddQkk,mddQkl),cbind(t(mddQkl),mddQll))
va <- ( (1-2*vmu)/(vmu*(1-vmu)) + vpsi*(vphi^2)*(vy-vmu)/(vy*(1-vy)) )*(1/vdotg)
vb <- ( 2/vphi + 1 - vpsi*vphi*((vmu^2)/vy + ((1-vmu)^2)/(1-vy)) )*(1/vdoth)
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
plot(Bi.cwp, pch=20,ylab="local influence",xlab="index")# , las=1,ylim=c(0,1))
identify(Bi.cwp)

## kappa
vBh.k.cwp <- numeric()
for (i in 1:n){
 vBh.k.cwp[i] <- 2*t(mDelta.cwp)[i,]%*%(mddQinv - mddQinv.l)%*%mDelta.cwp[,i]
}
Bi.k.cwp <- abs(vBh.k.cwp)/tr(-2*t(mDelta.cwp)%*%(mddQinv - mddQinv.l)%*%mDelta.cwp)
plot(Bi.k.cwp, pch=20,ylab="local influence",xlab="index")# , las=1,ylim=c(0,1))
mtext("(a)", cex=1.5, line=0.5)
identify(Bi.k.cwp)

## lambda
vBh.l.cwp <- numeric()
for (i in 1:n){
 vBh.l.cwp[i] <- 2*t(mDelta.cwp)[i,]%*%(mddQinv - mddQinv.k)%*%mDelta.cwp[,i]
}
Bi.l.cwp <- abs(vBh.l.cwp)/tr(-2*t(mDelta.cwp)%*%(mddQinv - mddQinv.k)%*%mDelta.cwp)
plot(Bi.l.cwp, pch=20,ylab="local influence",xlab="index")# , las=1,ylim=c(0,1))
mtext("(b)", cex=1.5, line=0.5)
identify(Bi.l.cwp)

#### Response perturbation
Sy <- sd(vy)
vPhi <- vpsi*(vphi^2)*( vmu/(vy^2) + (1-vmu)/((1-vy)^2) )/vdotg
vGamma <- -vpsi*vphi*( (vmu^2)/(vy^2) + ((1-vmu)^2)/((1-vy)^2) )/vdoth
mPhi <- diag(as.vector(vPhi))
mGamma <- diag(as.vector(vGamma))
mDeltak.resp <- Sy * t(mX)%*%mPhi
mDeltal.resp <- Sy * t(mV)%*%mGamma 
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


#### Mean covariate perturbation: income covariate
Sj <- sd(income)
v1.kappa.p <- c(0,1,0)
kappa.p <- vkappa[2]
vrho <- ( (1-2*vmu)/(vmu*(1-vmu)) + vpsi*(vphi^2)*(vy-vmu)/(vy*(1-vy)) )
mDeltak.cov <- Sj * ( kappa.p*t(mX)%*%mW1 + v1.kappa.p%*%t(vrho) ) %*% mG
mDeltal.cov <- Sj * kappa.p * t(mV)%*%mH%*%mG%*%mW2
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
#### Generalized leverage
##############################

vSigma.k <- vpsi*(vphi^2)*( vmu/(vy^2) + (1-vmu)/((1-vy)^2) )
vSigma.l <- -vpsi*vphi*( (vmu^2)/(vy^2) + ((1-vmu)^2)/((1-vy)^2) )
mSigma.k <- diag(as.vector(vSigma.k))
mSigma.l <- diag(as.vector(vSigma.l))
mQtheta.mu <- rbind(t(mX)%*%mG%*%mSigma.k, t(mV)%*%mH%*%mSigma.l)
mDmu <- cbind(mG%*%mX, matrix(0,n,q))
mL <- mDmu%*%solve(-mddQ)%*%mQtheta.mu
vleverage <- diag( mL )
plot(vleverage, pch=20,ylab="generalized leverage",xlab="index")
identify(vleverage)
mtext("(a)", cex=1.5, line=0.5) 

quantRESID <- fit.model$residuals
plot(quantRESID, vleverage, pch=20,ylab="generalized leverage",xlab="quantile residuals")
identify(quantRESID,vleverage)
mtext("(b)", cex=1.5, line=0.5)



##############################
#### Confirmatory analysis
##############################

## obs #4
fit.model.s04 <- bbreg(vy[-c(4)] ~ income[-c(4)] + persons[-c(4)], link.precision = c("log"), model="bessel")
summary(fit.model.s04)
round(((coef(fit.model)-coef(fit.model.s04))/coef(fit.model)) * 100, 2)

## obs #5
fit.model.s05 <- bbreg(vy[-c(5)] ~ income[-c(5)] + persons[-c(5)], link.precision = c("log"), model="bessel")
summary(fit.model.s05)
round(((coef(fit.model)-coef(fit.model.s05))/coef(fit.model)) * 100, 2)

## obs #8
fit.model.s08 <- bbreg(vy[-c(8)] ~ income[-c(8)] + persons[-c(8)], link.precision = c("log"), model="bessel")
summary(fit.model.s08)
round(((coef(fit.model)-coef(fit.model.s08))/coef(fit.model)) * 100, 2) 

## obs #9
fit.model.s09 <- bbreg(vy[-c(9)] ~ income[-c(9)] + persons[-c(9)], link.precision = c("log"), model="bessel")
summary(fit.model.s09)
round(((coef(fit.model)-coef(fit.model.s09))/coef(fit.model)) * 100, 2) 

## obs #11
fit.model.s11 <- bbreg(vy[-c(11)] ~ income[-c(11)] + persons[-c(11)], link.precision = c("log"), model="bessel")
summary(fit.model.s11)
round(((coef(fit.model)-coef(fit.model.s11))/coef(fit.model)) * 100, 2) 

## obs #20
fit.model.s20 <- bbreg(vy[-c(20)] ~ income[-c(20)] + persons[-c(20)], link.precision = c("log"), model="bessel")
summary(fit.model.s20)
round(((coef(fit.model)-coef(fit.model.s20))/coef(fit.model)) * 100, 2) 

## obs #38
fit.model.s38 <- bbreg(vy[-c(38)] ~ income[-c(38)] + persons[-c(38)], link.precision = c("log"), model="bessel")
summary(fit.model.s38)
round(((coef(fit.model)-coef(fit.model.s38))/coef(fit.model)) * 100, 2) 

