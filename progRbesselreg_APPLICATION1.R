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
#### DataSet: Water Quality Index 
#### Bonat, W.H., Petterle, R.R., Hinde, J., Dem�trio, C.G.B. (2019). 
#### 	 Flexible quasi-beta regression models for continuous bounded data. 
#### 	 Statistical Modelling, 19, 617-633.

dataset <- read.table("iqadata.txt", header = TRUE)
attach(dataset)
head(dataset)
vy <- IQA ##= response variable on (0,1)
TRIM <- factor(TRIM)
LOCAL <- factor(LOCAL, levels=c("MONT","RESER","JUSA"))
levels(LOCAL) <- c("Upstream","Reservoir","Downstream")

#########################################
#### Fitting the Bessel regression model
#########################################

fit.model <- bbreg(vy ~ LOCAL + TRIM | LOCAL + TRIM, model="bessel")#, residual = "quantile", prob = 0.95, envelope = 1000)
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
mtext("(a)", cex=1.5, line=0.5)
identify(GDi)

## kappa
GDi.k <- apply(mDdiag[,1:p],1,sum)
plot(GDi.k, pch=20,ylab="global influence",xlab="index")
mtext("(a)", cex=1.5, line=0.5)
#abline(h=mean(GDi.k) + 3*sd(GDi.k),lty=2)
#abline(h=median(GDi.k) + 6*mad(GDi.k),lty=2)
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
fit.model.s04 <- bbreg(vy[-c(4)] ~ LOCAL[-c(4)] + TRIM[-c(4)] | LOCAL[-c(4)] + TRIM[-c(4)], model="bessel")
summary(fit.model.s04)
round(((coef(fit.model)-coef(fit.model.s04))/coef(fit.model)) * 100, 2)

## obs #6
fit.model.s06 <- bbreg(vy[-c(6)] ~ LOCAL[-c(6)] + TRIM[-c(6)] | LOCAL[-c(6)] + TRIM[-c(6)], model="bessel")
summary(fit.model.s06)
round(((coef(fit.model)-coef(fit.model.s06))/coef(fit.model)) * 100, 2)

## obs #18
fit.model.s18 <- bbreg(vy[-c(18)] ~ LOCAL[-c(18)] + TRIM[-c(18)] | LOCAL[-c(18)] + TRIM[-c(18)], model="bessel")
summary(fit.model.s18)
round(((coef(fit.model)-coef(fit.model.s18))/coef(fit.model)) * 100, 2)

## obs #26
fit.model.s26 <- bbreg(vy[-c(26)] ~ LOCAL[-c(26)] + TRIM[-c(26)] | LOCAL[-c(26)] + TRIM[-c(26)], model="bessel")
summary(fit.model.s26)
round(((coef(fit.model)-coef(fit.model.s26))/coef(fit.model)) * 100, 2)

## obs #39
fit.model.s39 <- bbreg(vy[-c(39)] ~ LOCAL[-c(39)] + TRIM[-c(39)] | LOCAL[-c(39)] + TRIM[-c(39)], model="bessel")
summary(fit.model.s39)
round(((coef(fit.model)-coef(fit.model.s39))/coef(fit.model)) * 100, 2)

## obs #68
fit.model.s68 <- bbreg(vy[-c(68)] ~ LOCAL[-c(68)] + TRIM[-c(68)] | LOCAL[-c(68)] + TRIM[-c(68)], model="bessel")
summary(fit.model.s68)
round(((coef(fit.model)-coef(fit.model.s68))/coef(fit.model)) * 100, 2)

## obs #78
fit.model.s78 <- bbreg(vy[-c(78)] ~ LOCAL[-c(78)] + TRIM[-c(78)] | LOCAL[-c(78)] + TRIM[-c(78)], model="bessel")
summary(fit.model.s78)
round(((coef(fit.model)-coef(fit.model.s78))/coef(fit.model)) * 100, 2)

## obs #91
fit.model.s91 <- bbreg(vy[-c(91)] ~ LOCAL[-c(91)] + TRIM[-c(91)] | LOCAL[-c(91)] + TRIM[-c(91)], model="bessel")
summary(fit.model.s91)
round(((coef(fit.model)-coef(fit.model.s91))/coef(fit.model)) * 100, 2)

## obs #92
fit.model.s92 <- bbreg(vy[-c(92)] ~ LOCAL[-c(92)] + TRIM[-c(92)] | LOCAL[-c(92)] + TRIM[-c(92)], model="bessel")
summary(fit.model.s92)
round(((coef(fit.model)-coef(fit.model.s92))/coef(fit.model)) * 100, 2)

## obs #93
fit.model.s93 <- bbreg(vy[-c(93)] ~ LOCAL[-c(93)] + TRIM[-c(93)] | LOCAL[-c(93)] + TRIM[-c(93)], model="bessel")
summary(fit.model.s93)
round(((coef(fit.model)-coef(fit.model.s93))/coef(fit.model)) * 100, 2)

## obs #108
fit.model.s108 <- bbreg(vy[-c(108)] ~ LOCAL[-c(108)] + TRIM[-c(108)] | LOCAL[-c(108)] + TRIM[-c(108)], model="bessel")
summary(fit.model.s108)
round(((coef(fit.model)-coef(fit.model.s108))/coef(fit.model)) * 100, 2)

## obs #126
fit.model.s126 <- bbreg(vy[-c(126)] ~ LOCAL[-c(126)] + TRIM[-c(126)] | LOCAL[-c(126)] + TRIM[-c(126)], model="bessel")
summary(fit.model.s126)
round(((coef(fit.model)-coef(fit.model.s126))/coef(fit.model)) * 100, 2)

## obs #160
fit.model.s160 <- bbreg(vy[-c(160)] ~ LOCAL[-c(160)] + TRIM[-c(160)] | LOCAL[-c(160)] + TRIM[-c(160)], model="bessel")
summary(fit.model.s160)
round(((coef(fit.model)-coef(fit.model.s160))/coef(fit.model)) * 100, 2)

## obs #161
fit.model.s161 <- bbreg(vy[-c(161)] ~ LOCAL[-c(161)] + TRIM[-c(161)] | LOCAL[-c(161)] + TRIM[-c(161)], model="bessel")
summary(fit.model.s161)
round(((coef(fit.model)-coef(fit.model.s161))/coef(fit.model)) * 100, 2)

## obs #158
fit.model.s158 <- bbreg(vy[-c(158)] ~ LOCAL[-c(158)] + TRIM[-c(158)] | LOCAL[-c(158)] + TRIM[-c(158)], model="bessel")
summary(fit.model.s158)
round(((coef(fit.model)-coef(fit.model.s158))/coef(fit.model)) * 100, 2)



