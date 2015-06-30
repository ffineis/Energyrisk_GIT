library(ENERGYRISK)
library(lubridate)
library(zoo)
# If you want to convert a .txt file into .rda then save and reload data using data()
#S6_ForwardCurve$Date = as.Date(S6_ForwardCurve$Date,"%m/%d/%y") 
#save.image(file="data/data_S6.rda")
# Load storage contract specs and forward curve
data(data_S6)
# Illustrate Forward Curves
S6_ForwardCurve
# Illustrate Model Parameters
S6_ModelParams
# Setup costs (strike), heat rate and interest rates
# Number of simulations
N = 250
# Estimate number of days for simulation: beg-end of forward strip date 
nbDays = as.numeric(S6_ForwardCurve[length(S6_ForwardCurve[,1]),1]-
                      S6_ForwardCurve[1,1])
# Don't forget to load library(zoo) to manage dates
storDates = as.Date(S6_ForwardCurve[1,1]:S6_ForwardCurve[length(S6_ForwardCurve[,1]),1])
# Set start and end dates for theforward curve
valDate = S6_ForwardCurve[1,1]
valEnd = S6_ForwardCurve[length(S6_ForwardCurve[,1]),1]
c(valDate,valEnd)
# Costs (strike)
VOM = 1.50; Emission = 10; heatRate=7.25; K = VOM + Emission;r=0.05
# Create a list of simulation variables
param = list(K,VOM,Emission,heatRate,r,N)
names(param) = c("K","VOM", "Emission","heatRate","r","N")
t(param)

# Define the Power/GasPrice/CT matrix
valueStor = matrix(0,nbDays,3)
colnames(valueStor) = c("Power","Gas","CT")
# Number of simulation matrix #sim/OptionValue
simMatrix = matrix(0,param$N,2)
# Initialize starting values & time step dt
valueStor[1,c("Power","Gas")] = c(S6_ForwardCurve$Power[1], S6_ForwardCurve$Gas[1])
# step matrix for the whole period of the model
dt =1/365
for(j in 1:N){
  ## Model price changes
  for(i in 1:(nbDays-1)){
    # Keep track of the month for forward curve between month estimation
    trackMonth = which(format(S6_ForwardCurve$Date,"%m%y")==format(storDates[i+1],"%m%y"))
    F_t1 = c(S6_ForwardCurve[trackMonth,"Power"],S6_ForwardCurve[trackMonth,"Gas"])
    F_t2 = c(S6_ForwardCurve[trackMonth+1,"Power"],S6_ForwardCurve[trackMonth+1,"Gas"])
    
    # step matrix for the whole period of the model
    t1 = i/nbDays
    # Set seed so that the random variable does not change
    # Use qnorm: give it a probability and returns the number whose 
    # cumulative distribution matches the probability
    rand1 = qnorm(runif(1,min=0,max=1))
    rand2 = qnorm(runif(1,min=0,max=1))
    # stochastic component estimation
    dz1 = sqrt(dt)*rand1
    # Capture the stochastic gas compnent which is correlated to power demand
    dz2 = S6_ModelParams$rho[1]*dz1 + sqrt(1-(S6_ModelParams$rho[1])^2)*sqrt(dt)*rand2
      
    # Estimate mu from model 
    muPower = (1/S6_ModelParams[1,"alpha"])*
      (log(F_t2[1])-log(F_t1[1]))/(1/12) +
      log(F_t1[1]) + S6_ModelParams[1,"sigma"]^2 /
      (4*S6_ModelParams[1,"alpha"] ) *(1 - exp(-2*S6_ModelParams[1,"alpha"]*t1)) 
    muGas = (1/S6_ModelParams[2,"alpha"])*
      (log(F_t2[2])-log(F_t1[2]))/(1/12) +
      log(F_t1[2]) + S6_ModelParams[2,"sigma"]^2 /
      (4*S6_ModelParams[2,"alpha"] ) *(1 - exp(-2*S6_ModelParams[2,"alpha"]*t1)) 
    
    tempPower = log(valueStor[i,"Power"]) + (S6_ModelParams[1,"alpha"]*
                                              (muPower-log(valueStor[i,"Power"])) -
      0.5 * S6_ModelParams[1,"sigma"]^2) * dt + S6_ModelParams[1,"sigma"]* dz1
    tempGas = log(valueStor[i,"Gas"]) + (S6_ModelParams[2,"alpha"]*
                                              (muGas-log(valueStor[i,"Gas"])) -
      0.5 * S6_ModelParams[2,"sigma"]^2) * dt + S6_ModelParams[2,"sigma"]* dz2
    # Store values to be used on the next iteration
    valueStor[i+1,"Power"] = exp(tempPower)
    valueStor[i+1,"Gas"] = exp(tempGas)
    # Option value payoff for plant
    valueStor[i+1,"CT"] = max(0,valueStor[i+1,"Power"]-param$heatRate*valueStor[i+1,"Gas"]-
                                param$K)*
                                exp(-param$r*(as.numeric(storDates[i]-valDate)/365))
  }
  # option value
  simMatrix[j,] = c(j,sum(valueStor[,"CT"]))
}
# Illustrate the result from all N simulations
head(simMatrix)

# Standard Error
se = sd(simMatrix[,2])/sqrt(N)
# 4) Call Value discounted back
Option_Value = mean(simMatrix[,2])
Result = data.frame(c(Option_Value,se))
colnames(Result) = c("Result")
rownames(Result) = c("Call_Value","se")
Result

# Estimate cumulative average payoff over simulated paths 
avgOptionValue = matrix(0,N,1)
for(k in 1:N){avgOptionValue[k] = mean(simMatrix[1:k,2])}
head(simMatrix)
# Plot Convergence path for MC N simulations
plot(c(avgOptionValue),type="l", ylab="Payoff", xlab="Number of Simulations")
legend("topright",legend=c("Simple Method"),col = c(1), lwd=c(1.5), lty=c(1), cex=0.8)
title("Monte Carlo Convergence of Option Value")
# Plot frequency distribution of the simulated payoffs: breaks = 20
hist(simMatrix[,2],xlab="Payoff Bin",main=paste("Frequency of",N,"Simulations"),breaks=20)