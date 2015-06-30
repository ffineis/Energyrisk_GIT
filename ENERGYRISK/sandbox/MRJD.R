library(ENERGYRISK)
### For Mean Reversion and Jump Diffusion
## Load Data from ICENBP.txt file in directory
ICENBP = read.delim("~/Google Drive/UWGARP_FRM/EnergyRisk/ENERGYRISK/data/csv/ICENBP.txt")
# If you want to convert a .txt file into .rda then save and reload data using data()
ICENBP$Date = as.Date(ICENBP$Date, "%m/%d/%y")
data_NBP = ICENBP
save(data_NBP,file="data/data_NBP.rda")
data(data_NBP)
# Illustarte first and last 6 obs to get a quick view of the data
head(data_NBP)
tail(data_NBP)
# Quick Plot of Loaded Spot Prices
plot(data_NBP,type="l", xlab="Time in Months", ylab="Spot Price")
title("Spot Price Aug 2007-2008")

## Mean Reversion
# Estimate log returns
# This can be done in two steps or just one.
# The idea is to 1) Estimate the log of the spot prices and then 
# 2) take the difference (log-return) 
# (Making sure to shift/index the time series by one to start)
nRow = nrow(data_NBP)
# Match prices to expected returns
lnSt = log(data_NBP$St[2:nRow])[-(nRow-1)]
head(lnSt)
returns = log(data_NBP$St[2:nRow]) - log(data_NBP$St[1:(nRow-1)])
returns = returns[-1]
head(returns)
plot(lnSt,returns,type="p")
abline(0,0)
# Initialze Starting Assumptions
# Time Step
dt = 1/365
nObs = length(returns)
# Estimate Parameters using Least Square to get MR_reg Object
MR_reg = lm(returns~lnSt)
# View Summary Results from regression
result = summary(MR_reg)
# Illustrating two methods for retrieve 1) coefs and/or 2) s.e.
# Retrieve intercept and slope from summary
alpha0 = coef(MR_reg)[1]
alpha1 = coef(MR_reg)[2]
c(alpha0,alpha1)
# Retrieve standard error
result$coef[,2]
# Fitting Graph: Either use alpha0 and alpha1 to estimate fit or fitted() 
fit = data.frame(alpha0 + alpha1*lnSt)
head(fit)
# OR
fit = fitted(MR_reg)
head(fit)
# Fit line on graph
lines(lnSt,fit,type="l", col="red", lwd=2, lty=3)

# Mean Reversion Parameter Specific Results 
# Underlying Mean and Alpha
UM = exp(-alpha0/alpha1)
alpha = -alpha1/dt
UMAlpha = data.frame(c(UM,alpha))
colnames(UMAlpha)=c("Results")
rownames(UMAlpha) = c("Underlying Mean","Alpha")
t(UMAlpha)
# Estimate Half-Life Exercise
t_Onehalf = data.frame((log(2)/alpha)*365)
colnames(t_Onehalf) = c("Half-Life")
rownames(t_Onehalf) = c("Results")
t_Onehalf

## Jump Diffusion
# Initilize starting variable to start recursively filtering
# Remember to set a limit of 3 times the standard deviation
std = sd(returns)
limit = 3*std
c(std,limit)
# Now Estimate the jump diffusion time-series
jumps = returns[which(abs(returns)>limit)]
nonJumps = returns[-which(abs(returns)>limit)]
nJumps = length(jumps)
# Organize data into a list to be drawn later
jumpStorage = list(jumps,nonJumps,nJumps)
names(jumpStorage) = c("jumps","nonJumps","nJumps")
jumpStorage$nJumps
jumpStorage$jumps
# Only look at the first 6 rows b/c the key still represents 365 obs post jump removal
head(jumpStorage$nonJumps)
# Estimate phi, sigma, \bar{kappa} & gamma
# Sigma is expressed in annual volatility terms
phi = jumpStorage$nJumps/(nObs/365)
sigma = sd(jumpStorage$nonJumps) * sqrt(365)
kappa_bar = sum(jumpStorage$jumps)/jumpStorage$nJumps
gamma = sd(jumpStorage$jumps)
# Organize data
JD = data.frame(c(phi,sigma,kappa_bar,gamma))
rownames(JD) = c("phi","sigma","kappa_bar","gamma")
colnames(JD) = c("Iteration 1: ")
t(JD)

# Deriving JDit R
tol0 = -Inf
tolVar = "phi"
tol = 1
done = 0
it = 1
nJumps = 0
limit = 3
R = returns
table = matrix(0,5,1)
if(tolVar == "phi"){index = 1}else if(tolVar=="sigma"){index=2}else if(toldVar=="kappa_bar"){index=3
}else if(toldVar == "gamma"){index = 4}else{stop("check tolVar specification")}
while(!done){
  std = sd(R)
  threshold = limit*std
  jumps = returns[which(abs(returns)>threshold)]
  nonJumps = R[-which(abs(R)>threshold)]
  nJumps = length(jumps)
  jumpStorage = list(jumps,nonJumps,nJumps)
  names(jumpStorage) = c("jumps","nonJumps","nJumps")
  
  phi = jumpStorage$nJumps/(nObs/365)
  sigma = sd(jumpStorage$nonJumps) * sqrt(365)
  kappa_bar = sum(jumpStorage$jumps)/jumpStorage$nJumps
  gamma = sd(jumpStorage$jumps)
  
  JD = data.frame(c(nJumps,phi,sigma,kappa_bar,gamma))
  rownames(JD) = c("nJumps","phi","sigma","kappa_bar","gamma")
  colnames(JD) = paste("Iteration",it,":")
  
  R = jumpStorage$nonJumps
  it = it + 1
  table = cbind(table,JD)
  if((JD[index,]-tol0)<=tol){done = 1}
  tol0 = JD[index,]
}
return(table[,2:ncol(table)])
#Generalize result in functional form
JPit(returns, 3, "phi", 1)
