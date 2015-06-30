# Demo 2 code #
library(ENERGYRISK)
library(zoo)
data(data_NBP)

#Mean Reversion Process
#data_NBP <- as.zoo(data_NBP)
plot(data_NBP, type = "l", xlab = "Time (Month)", ylab = "Spot Price", main = "Spot Price Aug 2007-2008")
log_ret = diff(log(data_NBP$St))

lnSt <- log(data_NBP$St[-nrow(data_NBP)])
plot(lnSt, log_ret)
abline(0,0)
dt = 1/365
nObs = length(lnSt)
MR_reg = lm(log_ret~lnSt)
results = summary(MR_reg)
a0 = MR_reg$coefficients[1]; a1 = MR_reg$coefficients[2]
SE = results$coefficients[,"Std. Error"]
fit = fitted(MR_reg)
fit = as.numeric(a0 + a1*lnSt)
lines(lnSt, fit, col = "red", lwd = 2, lty = 3)

alpha = -a1/dt
halflife = log(2)/alpha

#Capture Jumps in Spot Price
std = sd(log_ret)
limit = 3*std
jumps = log_ret[which(abs(log_ret)>limit)]
nonjumps = log_ret[-which(abs(log_ret) > limit)]
jumpStorage = list("jumps" = jumps, "nonjumps" = nonjumps, "nJumps" = length(jumps))

phi = jumpStorage$nJumps/(nObs/365) #how many jumps to expect per year
sigma = sd(jumpStorage$nonjumps)*sqrt(365) #annualize standard deviation into years
kappa_bar = sum(jumpStorage$jumps)/jumpStorage$nJumps #usually taken to be zero...

