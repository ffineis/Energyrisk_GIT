library(ENERGYRISK)
library(lubridate)
library(matrixStats)
library(zoo)
# If you want to convert a .txt file into .rda then save and reload data using data()
# Load Data
#S9_ForwardCurve$Date = as.Date(S9_ForwardCurve$Date,"%m/%d/%y") 
#save.image(file="data/data_S9.rda")
data(data_S9)
# Illustrate Forward Curves
S9_ForwardCurve
# Illustrate Model Parameters
S9_ModelParams
# Setup costs (strike), heat rate and interest rates
# Number of simulations
N = 200
# Set start and end dates for theforward curve
valDate = as.Date("01/01/09","%m/%d/%y")
valEnd = as.Date("12/31/09","%m/%d/%y")
c(valDate,valEnd)
# Estimate number of days for simulation: beg-end of forward strip date 
nbDays = as.numeric(valEnd-valDate)
# Don't forget to load library(zoo) to manage dates
storDates = as.Date(S9_ForwardCurve[1,1]:S9_ForwardCurve[length(S9_ForwardCurve[,1]),1])
# Costs (strike)
KCap=90;KFloor=70;
VOM=1.50;Emission=10;heatRate=7.25;K=VOM+Emission;r=0.05
scaleFactor=0.1598;
# Create a list of simulation variables
param = list(KCap,KFloor,K,heatRate,r,scaleFactor,N)
names(param) = c("KCap","KFloor","K","heatRate","r","scaleFactor","N")
t(param)
######################################################
## Simple Cap/Floor Example
# Define the Period/Power/PumpGen/PumpingCost/GenRevenue matrix
valueStor = matrix(0,nbDays,5)
colnames(valueStor) = c("Power","Gas","HO",
                        "Cap","Floor")

# Bucketing of cashflows
indexStartDates = c(1,which(as.numeric(format(storDates[-length(storDates)],"%m"))!=
                              as.numeric(format(storDates[-1],"%m")))+1)[-13]
bucketStartDates = storDates[indexStartDates]
indexEndDates = c(which(as.numeric(format(storDates[-length(storDates)],"%m"))!=
                          as.numeric(format(storDates[-1],"%m")))[-12],364)
bucketEndDates = storDates[indexEndDates]
head(cbind(indexStartDates,indexEndDates))
# Setup matrix to hold buckets
bucketStor = matrix(0,N,length(indexStartDates))

# Generate the correlation structure of the error terms
# Initialize a=0(no correlation Power/Gas);b=0 (no correlation HO/Power)
a=0;b=0;
if(S9_ModelParams[2,"rho_Power"]<1){a = (S9_ModelParams[2,"rho_HO"]-
                                           S9_ModelParams[3,"rho_Power"]*
                                           S9_ModelParams[2,"rho_Power"])/
                                      sqrt(1 -S9_ModelParams[2,"rho_Power"]^2 )}
if((S9_ModelParams[3,"rho_Power"]^2+a^2)<=1){b = sqrt(1-S9_ModelParams[3,"rho_Power"]^2-
                                                        a^2)}

# Time step dt
dt =1/365
for(j in 1:N){
  # Initialize starting values & time step dt
  valueStor[1,"Power"] = S9_ForwardCurve$Power[1]
  valueStor[1,"Gas"] = S9_ForwardCurve$Gas[1]
  valueStor[1,"HO"] = S9_ForwardCurve$Heating.Oil[1]
  # Initialize the contracts: Cap/Floor/Swap/FueldSwitching/Collar
  ## Cap/Floor/Swap/FuelSwitching/Collar Cash-Flow
  ## Cap/Floor on Power
  valueStor[1,"Cap"] = max(0,valueStor[1,"Power"]-param$KCap)*
    exp(-param$r*(as.numeric(storDates[1]-valDate))/365)
  valueStor[1,"Floor"] = -max(0,param$KFloor- valueStor[1,"Power"])*
    exp(-param$r*(as.numeric(storDates[1]-valDate))/365)
  for(i in 1:(nbDays-1)){
    # Keep track of the month for forward curve between month estimation
    trackMonth = which(format(S9_ForwardCurve$Date,"%m%y")==format(storDates[i+1],"%m%y"))
    F_t1 = c(S9_ForwardCurve[trackMonth,"Power"],S9_ForwardCurve[trackMonth,"Gas"],
             S9_ForwardCurve[trackMonth,"Heating.Oil"])
    F_t2 = c(S9_ForwardCurve[trackMonth+1,"Power"],S9_ForwardCurve[trackMonth+1,"Gas"],
             S9_ForwardCurve[trackMonth+1,"Heating.Oil"])
    # Step matrix for the whole period of the model
    t1 = 1/nbDays
    # Use qnorm: give it a probability and returns the number whose 
    # cumulative distribution matches the probability
    rand1 = qnorm(runif(1,min=0,max=1))
    rand2 = qnorm(runif(1,min=0,max=1))
    rand3 = qnorm(runif(1,min=0,max=1))
    # Stochastic component estimation
    dz1 = sqrt(dt)*rand1
    dz2 = S9_ModelParams[2,"rho_Power"]*sqrt(dt)*rand1+
      sqrt(1 -S9_ModelParams[3,"rho_Power"]^2)*sqrt(dt)*rand2
    dz3 = S9_ModelParams[3,"rho_Power"]*sqrt(dt)*rand1+
      a*sqrt(dt)*rand2+b*sqrt(dt)*rand3
    
    # Estimate mu from model 
    muPower = (1/S9_ModelParams[1,"alpha"])*
      (log(F_t2[1])-log(F_t1[1]))/(1/12) +
      log(F_t1[1]) + S9_ModelParams[1,"sigma"]^2 /
      (4*S9_ModelParams[1,"alpha"] ) *(1 - exp(-2*S9_ModelParams[1,"alpha"]*t1)) 
    muGas = (1/S9_ModelParams[2,"alpha"])*
      (log(F_t2[2])-log(F_t1[2]))/(1/12) +
      log(F_t1[2]) + S9_ModelParams[2,"sigma"]^2 /
      (4*S9_ModelParams[2,"alpha"] ) *(1 - exp(-2*S9_ModelParams[2,"alpha"]*t1)) 
    muHO = (1/S9_ModelParams[3,"alpha"])*
      (log(F_t2[3])-log(F_t1[3]))/(1/12) +
      log(F_t1[3]) + S9_ModelParams[3,"sigma"]^2 /
      (4*S9_ModelParams[3,"alpha"] ) *(1 - exp(-2*S9_ModelParams[3,"alpha"]*t1)) 
    
    tempPower = log(valueStor[i,"Power"]) + (S9_ModelParams[1,"alpha"]*
              (muPower-log(valueStor[i,"Power"])) -
              0.5 * S9_ModelParams[1,"sigma"]^2) * dt + S9_ModelParams[1,"sigma"]* dz1
    tempGas = log(valueStor[i,"Gas"]) + (S9_ModelParams[2,"alpha"]*
              (muGas-log(valueStor[i,"Gas"])) -
              0.5 * S9_ModelParams[2,"sigma"]^2) * dt + S9_ModelParams[2,"sigma"]* dz2
    tempHO = log(valueStor[i,"HO"]) + (S9_ModelParams[3,"alpha"]*
              (muHO-log(valueStor[i,"HO"])) -
              0.5 * S9_ModelParams[3,"sigma"]^2) * dt + S9_ModelParams[3,"sigma"]* dz3
    # Store values to be used on the next iteration
    valueStor[i+1,"Power"] = exp(tempPower)
    valueStor[i+1,"Gas"] = exp(tempGas)
    valueStor[i+1,"HO"] = exp(tempHO)
    ## Cap/Floor/Swap/FuelSwitching/Collar Cash-Flow
    ## Cap/Floor on Power
    valueStor[i+1,"Cap"] = max(0,valueStor[i+1,"Power"]-param$KCap)*
      exp(-param$r*(as.numeric(storDates[i+1]-valDate))/365)
    valueStor[i+1,"Floor"] = -max(0,param$KFloor- valueStor[i+1,"Power"])*
      exp(-param$r*(as.numeric(storDates[i+1]-valDate))/365)
  }
  # Sum all the bucket cash flows together
  for(k in 1:length(indexStartDates)){bucketStor[j,k] = sum(valueStor[indexStartDates[k]:
                                                                        indexEndDates[k],4:5])}
}
# Illustrate the payoff storage for all products
head(valueStor)
# Illustrate aggregation of the bucket cash-flows for all buckets
head(bucketStor)
# Estimate Average per Bucket
averageValue = colMeans(bucketStor)
# Standard Error per Bucket (library(matrixStats))
se = colSds(bucketStor)/sqrt(N)
Result = rbind(averageValue,se)
rownames(Result) = c("Mean_Value","se")
Result
# Estimate the Simulation Quartiles: 19 quartiles [0.05:0.05:0.95]
percentile = matrix(0,19,length(indexEndDates))
for(t in 1:length(indexEndDates)){percentile[,t] = quantile(bucketStor[,t],
                                                            c(seq(0.05,0.95,0.05)))}

# Plot quartiles EaR
plot(c(percentile[1,]),type="l", ylab="Earning at Risk", xlab="Months",
     ylim=c(min(percentile),max(percentile)))
for(v in 2:19){lines(percentile[v,],col=v,lwd=1, lty=v)}
legend("topright",legend=c("5%","10%","90%","95%"),
       col = c(1,2,18,19), lwd=c(1.5,1.5,1.5,1.5), lty=c(1,2,18,19), cex=0.8)