#### Bias Adjustment in Recruitment
#### 1. Steepness in Ricker model (AMAK VS. BAM and SS) ####
bam_ss_h <- seq(from=0, to=2, by=0.01)
amak_h <- exp(bam_ss_h) / (4 + exp(bam_ss_h))
plot(x=bam_ss_h, y=amak_h,
     xlab="h from BAM, SS, and OM", ylab="h' from AMAK",
     xlim=range(amak_h, bam_ss_h), ylim=range(amak_h, bam_ss_h),
     pch=16, cex=0.5)

om_h <- 0.75
abline(v=om_h, col="gray20", lty=2)
abline(h=exp(om_h) / (4 + exp(om_h)), col="gray50", lty=3)
legend("topright",
       legend=c(paste("OM h = BAM h = SS h =", om_h),
                paste("AMAK h =", round(exp(om_h) / (4 + exp(om_h)), digits=2))),
       lty=c(2,3),
       col=c("gray20", "gray50"),
       bty="n")

#### 2. Compare BAM and SS estimates following Erik's method (Ricker model) ####
R0_med_true=1000     # true median R0 value
h_med_true=0.75   # vector of median steepness (h) values
phi0=10

sd.val=0.6
BC=(sd.val^2)/2

n=100000 # use a fairly large number because fits need be very tight

S0_med_true=R0_med_true*phi0

#S.vec=runif(n,min=0,max=S0)
S.vec=seq(0.001,S0_med_true, length=n)

R.vec.true=S.vec/phi0*exp(h_med_true*(1-S.vec/(R0_med_true*phi0)))
lnormdev.vec=rlnorm(n,0,sd.val)
R.vec.obs=R.vec.true*lnormdev.vec

# median S-R fit (produces median R0 and h, used in BAM)
SR.fit.med=function(par){
  R0.par=par[1]
  h.par=par[2]
  R.vec.pred=S.vec/phi0*exp(h.par*(1-S.vec/(R0.par*phi0)))
  resid=(log(R.vec.obs)-log(R.vec.pred))^2 # median likelihood
  resid.out=sum(resid)
  resid.out}

# SS S-R fit (produces mean R0 and h, used in SS)
SR.fit.SS=function(par){
  R0.par=par[1]
  h.par=par[2]
  R.vec.pred=S.vec/phi0*exp(h.par*(1-S.vec/(R0.par*phi0)))*exp(-BC)
  resid=(log(R.vec.obs)-log(R.vec.pred))^2 # median likelihood
  resid.out=sum(resid)
  resid.out}

#optimization for SR parameters using mean and median fit criteria
fit.med=optim(c(R0_med_true,h_med_true),SR.fit.med,lower=c(0,0.01),method="L-BFGS-B",upper=c(3*R0_med_true,10))
fit.SS=optim(c(R0_med_true,h_med_true),SR.fit.SS,lower=c(0,0.01),method="L-BFGS-B",upper=c(3*R0_med_true,10))

R0.BAM=fit.med$par[1]
h.BAM=fit.med$par[2]

R0.SS=fit.SS$par[1]
h.SS=fit.SS$par[2]
S0.SS=R0.SS*phi0

phiF=seq(from=5,to=10,by=0.1)

# SS method
Seq.SS=S0.SS*(1.+(log(R0.SS/S0.SS)+log(phiF))/h.SS);
Req.SS=R0.SS*Seq.SS/S0.SS * exp(h.SS*(1.-Seq.SS/S0.SS));
# BAM method
Req.BAM=R0.BAM/(phiF/phi0)*(1+log(exp(BC)*phiF/phi0)/h.BAM)
Seq.BAM=Req.BAM*phiF

cat(paste("True_Median_R0=", R0_med_true, sep=""),
    paste("True_Median_h=", h_med_true, sep=""),
    paste("BAM_Median_R0=", R0.BAM, sep=""),
    paste("BAM_Median_h=", h.BAM, sep=""),
    paste("SS_Mean_R0=", R0.SS, sep=""),
    paste("SS_Mean_h=", h.SS, sep=""),
    sep="\n")

plot(phiF, Req.BAM, col="red", pch=1, type="o", lty=1, xlab="phiF", ylab="Req")
lines(phiF, Req.SS, col="deepskyblue3", pch=3, type="o", lty=2)
legend("bottomright",
       legend=c("BAM_Req", "SS_Req"),
       col=c("red", "deepskyblue3"),
       lty=c(1, 2),
       pch=c(1, 3),
       bty="n")

#### 3. Conversion function following Chris's method (Ricker model) ####
convertSRparms <- function(R0, h, phi, sigmaR, mean2med, model){
  ## model=1: B-H model
  ## model=2: Ricker model
  BC <- ifelse(mean2med == TRUE, exp(-0.5 * sigmaR^2), exp(0.5 * sigmaR^2))
  if (model==1) {
    S0BC <- (BC * 0.8 * R0 * h * phi - 0.2 * phi * R0 * (1 - h)) / (h - 0.2)
    R0BC <-  S0BC / phi
    Rnew <- BC * 0.8 * R0 * h * 0.2 * S0BC / (0.2 * phi * R0 * (1 - h) + (h - 0.2) * 0.2 * S0BC)
    hBC <- Rnew / R0BC
  }

  if (model==2) {
    S0BC <- R0 * phi * (1 + log(BC) / h)
    R0BC <- S0BC / phi
    Rnew <- BC * 0.2 *S0BC / phi * exp(h * (1 - 0.2 * S0BC/ (R0 * phi)))
    hBC <- 1.25*log(5*Rnew /R0BC)
  }

  return(list(S0BC = S0BC, R0BC = R0BC, hBC = hBC))
}

R0 <- R0_med_true
h <- h_med_true
phi <- phi0
sigmaR <- sd.val
mean2med <- FALSE
model <- 2
p <- convertSRparms(R0, h, phi, sigmaR, mean2med, model)

cat(paste("Median R0=", R0, sep=""),
    paste("Median h=", h, sep=""),
    paste("phi=", phi, sep=""),
    paste("sigmaR=", sigmaR, sep=""),
    paste("Mean R0 (R0BC)=", p$R0BC, sep=""),
    paste("Mean h (hBC)=", p$hBC, sep=""),
    sep="\n")

BC <- ifelse(mean2med == TRUE, exp(-0.5 * sigmaR^2), exp(0.5 * sigmaR^2))
S.vec <- seq(0, R0 * phi, length.out = 1000)
R.vec <- S.vec / phi * exp(h * (1 - S.vec/ (R0 * phi)))
R.vecBC <- R.vec * BC
R.calc <- S.vec / phi * exp(p$hBC * (1 - S.vec / (p$R0BC * phi)))
diff <- R.calc - R.vecBC
range(diff)
plot(S.vec, R.vecBC, col="red", xlab="SSB", ylab="R", type="l", lty=1)
lines(S.vec, R.calc, col="deepskyblue3", lty=2)
legend("topleft",
       legend=c("BAM_R", "SS_R"),
       col=c("red", "deepskyblue3"),
       lty=c(1,2),
       bty="n")

#### 4. Conversion plot over a range of steepness and sigmaR (Ricker model) ####
R0=1000000
phi=0.01025625
sigmaR=seq(0.01, 2, by=0.01)
h=seq(0.21, 1, by=0.01)

R0diff <- hdiff <- matrix(NA, nrow=length(sigmaR), ncol=length(h))

for(i in 1:length(sigmaR)){
  for(j in 1:length(h)){
    R0diff[i,j] <- (convertSRparms(R0=R0, h=h[j], phi=phi, sigmaR=sigmaR[i], mean2med=FALSE, model=2)$R0BC-R0)/R0*100
    hdiff[i,j] <- (convertSRparms(R0=R0, h=h[j], phi=phi, sigmaR=sigmaR[i], mean2med=FALSE, model=2)$hBC-h[j])
  }
}

par(mfrow=c(1,2))
plot(NA,xlim=range(sigmaR),
     ylim=range(h),xlab="SigmaR",ylab="Median h",
     frame=FALSE,axes=F,xaxs="i",yaxs="i", main="Relative Diff in R0 (%)")
contour(sigmaR, h, R0diff,
        col = "black",
        add = TRUE,
        levels = c(10, 50, 100, 300, 500),
        method="edge")
axis(1); axis(2);
legend("topleft", "A)", bty="n")
box()

plot(NA,xlim=range(sigmaR),
     ylim=range(h),xlab="SigmaR",ylab="Median h",
     frame=FALSE,axes=F,xaxs="i",yaxs="i", main="Diff in h")
contour(sigmaR, h, hdiff,
        col = "black",
        add = TRUE,
        levels = c(0.1, 0.4, 0.8, 1.2, 1.6),
        method="edge")
axis(1); axis(2);
legend("topleft", "B)", bty="n")
box()

#### 5. Mean h VS. Median h (B-H model and Ricker model) ####
# B-H model
h_vec <- seq(0.21, 1, by=0.01)
logR_sd_vec <- seq(from=0.2, to=2, by=0.2)

medianh_vec <- h_vec
logR_sd_vec <- logR_sd_vec
forloop_id <- expand.grid(medianh_vec = medianh_vec, logR_sd_vec = logR_sd_vec)

meanh_vec <- c()
meanR0_vec <- c()

for (i in 1:nrow(forloop_id)){
  meanh_vec[i] <- convertSRparms(R0=1000000,
                                 h=forloop_id$medianh_vec[i],
                                 phi=0.01025625,
                                 sigmaR=forloop_id$logR_sd_vec[i],
                                 mean2med=FALSE,
                                 model=1)$hBC
  meanR0_vec[i] <- convertSRparms(R0=1000000,
                                  h=forloop_id$medianh_vec[i],
                                  phi=0.01025625,
                                  sigmaR=forloop_id$logR_sd_vec[i],
                                  mean2med=FALSE,
                                  model=1)$R0BC
}
forloop_id$medianR0_vec <- rep(1000000, nrow(forloop_id))
forloop_id$converted_meanR0 <- meanR0_vec
forloop_id$converted_meanh  <- meanh_vec

par(mfrow=c(1,1), mar=c(4,4,0.5, 0.5))
col=rainbow(n=length(logR_sd_vec))
plot(forloop_id$medianh_vec[which(forloop_id$logR_sd_vec==logR_sd_vec[1])], forloop_id$converted_meanh[which(forloop_id$logR_sd_vec==logR_sd_vec[1])], type="l", xlab="Median h", ylab="Mean h", xlim=c(0.2,1), ylim=c(0.2,1), lty=1, col=col[1])
for(i in 1:length(logR_sd_vec)){
  lines(forloop_id$medianh_vec[which(forloop_id$logR_sd_vec==logR_sd_vec[i])], forloop_id$converted_meanh[which(forloop_id$logR_sd_vec==logR_sd_vec[i])], lty=i, col=col[i])
}
legend("bottomright",
       paste("sigmaR =", logR_sd_vec),
       col=col,
       lty=1:length(logR_sd_vec),
       bty="n")

# Ricker model
h_vec <- seq(0.01, 3, by=0.01)
logR_sd_vec <- seq(from=0.2, to=2, by=0.2)

medianh_vec <- h_vec
logR_sd_vec <- logR_sd_vec
forloop_id <- expand.grid(medianh_vec = medianh_vec, logR_sd_vec = logR_sd_vec)

meanh_vec <- c()
meanR0_vec <- c()

for (i in 1:nrow(forloop_id)){
  meanh_vec[i] <- convertSRparms(R0=1000000,
                                 h=forloop_id$medianh_vec[i],
                                 phi=0.01025625,
                                 sigmaR=forloop_id$logR_sd_vec[i],
                                 mean2med=FALSE,
                                 model=2)$hBC
  meanR0_vec[i] <- convertSRparms(R0=1000000,
                                  h=forloop_id$medianh_vec[i],
                                  phi=0.01025625,
                                  sigmaR=forloop_id$logR_sd_vec[i],
                                  mean2med=FALSE,
                                  model=2)$R0BC
}
forloop_id$medianR0_vec <- rep(1000000, nrow(forloop_id))
forloop_id$converted_meanR0 <- meanR0_vec
forloop_id$converted_meanh  <- meanh_vec

par(mfrow=c(1,1), mar=c(4,4,0.5, 0.5))
col=rainbow(n=length(logR_sd_vec))
plot(forloop_id$medianh_vec[which(forloop_id$logR_sd_vec==logR_sd_vec[1])], forloop_id$converted_meanh[which(forloop_id$logR_sd_vec==logR_sd_vec[1])], type="l", xlab="Median h", ylab="Mean h", xlim=c(0.2,3), ylim=c(0.2,3), lty=1, col=col[1])
for(i in 1:length(logR_sd_vec)){
  lines(forloop_id$medianh_vec[which(forloop_id$logR_sd_vec==logR_sd_vec[i])], forloop_id$converted_meanh[which(forloop_id$logR_sd_vec==logR_sd_vec[i])], lty=i, col=col[i])
}
legend("bottomright",
       paste("sigmaR =", logR_sd_vec),
       col=col,
       lty=1:length(logR_sd_vec),
       bty="n")
