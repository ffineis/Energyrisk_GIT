# Demo 6 code. July 8, 2015
library(ENERGYRISK)
library(lubridate)
library(zoo)
set.seed(1)

data(data_S6) #S6_ForwardCurve has forward prices for Power, Gas, and Heating.Oil
              #S6_ModelParams have alpha, sigma, and rho for Power and Gas mean-reverting model
M = 250 #number of MC simulations
nbDays = as.numeric(S6_ForwardCurve$Date[nrow(S6_ForwardCurve)]-S6_ForwardCurve$Date[1])
storDates = as.Date(x = S6_ForwardCurve$Date[1]:S6_ForwardCurve$Date[nrow(S6_ForwardCurve)]) #every day btwn Jan 1 2k9 and Jan1 2k10
valDate = storDates[1]; valEnd = storDates[length(storDates)] #valuation dates
VOM = 1.50; Emission = 10; heatRate=7.25; K = VOM + Emission; r=0.05 #parameters
param = list("K" = K, "VOM" = 1.5, "heatRate" = heatRate, "Emission" = Emission, "r" = r, "M" = M)

#simulate fuel, power prices and option payoff:
valueStor = matrix(0, nbDays, 3)
colnames(valueStor) = c("Power", "Gas", "CT")
#storage for each simulation:
simStor = matrix(0, param$M, 2)
valueStor[1,c("Power","Gas")] = c(S6_ForwardCurve$Power[1], S6_ForwardCurve$Gas[1]) #today's spot prices = today's forward
dt = 1/(as.numeric(storDates[length(storDates)]-storDates[1])) #dt = 1/365

for (j in 1:param$M){ #250 simulations
  for(i in 1:(nbDays-1)){   #each simulation simulates 364 days
    trackMonth = which(format(S6_ForwardCurve$Date,"%m%y")==format(storDates[i+1],"%m%y")) #keep track of month
    #get current month and next month's forward prices
    F_t1 = c(S6_ForwardCurve[trackMonth,"Power"],S6_ForwardCurve[trackMonth,"Gas"])
    F_t2 = c(S6_ForwardCurve[trackMonth+1,"Power"],S6_ForwardCurve[trackMonth+1,"Gas"])
    t1 = i/nbDays
    
    #standard normal random variables... get number between 0 and 1 and find value giving that quantile from normal dist.
    rand1 = qnorm(runif(1,min=0,max=1)) #epsilon 1
    rand2 = qnorm(runif(1,min=0,max=1)) #epsilon 2
    
    dz1 = sqrt(dt)*rand1
    # stochastic gas component which is correlated to power demand
    dz2 = S6_ModelParams$rho[1]*dz1 + sqrt(1-(S6_ModelParams$rho[1])^2)*sqrt(dt)*rand2
    
    # Estimate mu from mean-reverting model and forward prices?
    muPower = (1/S6_ModelParams[1,"alpha"])*
      (log(F_t2[1])-log(F_t1[1]))/(1/12) +
      log(F_t1[1]) + S6_ModelParams[1,"sigma"]^2 /
      (4*S6_ModelParams[1,"alpha"] ) *(1 - exp(-2*S6_ModelParams[1,"alpha"]*t1)) 
    muGas = (1/S6_ModelParams[2,"alpha"])*
      (log(F_t2[2])-log(F_t1[2]))/(1/12) +
      log(F_t1[2]) + S6_ModelParams[2,"sigma"]^2 /
      (4*S6_ModelParams[2,"alpha"] ) *(1 - exp(-2*S6_ModelParams[2,"alpha"]*t1)) 
    
    # delta_x = (alpha(mu-x) - .5*sigma^2)*dt - sigma*sqrt(t)*epsilon
    tempPower = log(valueStor[i,"Power"]) + (S6_ModelParams[1,"alpha"]*
                                               (muPower-log(valueStor[i,"Power"])) -
                                               0.5 * S6_ModelParams[1,"sigma"]^2) * dt + S6_ModelParams[1,"sigma"]* dz1
    tempGas = log(valueStor[i,"Gas"]) + (S6_ModelParams[2,"alpha"]*
                                           (muGas-log(valueStor[i,"Gas"])) -
                                           0.5 * S6_ModelParams[2,"sigma"]^2) * dt + S6_ModelParams[2,"sigma"]* dz2
    #transform back into price
    valueStor[i+1,"Power"] = exp(tempPower)
    valueStor[i+1,"Gas"] = exp(tempGas)
    
    # Option value payoff for plant: df*max(0, P-HR*G-K)
    valueStor[i+1,"CT"] = max(0,valueStor[i+1,"Power"]-param$heatRate*valueStor[i+1,"Gas"]-
                                param$K)*
      exp(-param$r*(as.numeric(storDates[i]-valDate)/365))#discount_factor
  }
  simStor[j,] = c(j,sum(valueStor[,"CT"])) #value = sum of daily payoffs. Store MC valuations.
}
#Formatting results:
colnames(simStor) = c("MC_iteration", "Thermal_valuation")
se = sd(simStor[,"Thermal_valuation"])/sqrt(param$M) #standard error in MC valuation
thermal_val = mean(simStor[,"Thermal_valuation"])
result = data.frame(t(c(thermal_val, se)))
colnames(result) = c("Call_value", "std_error")
rownames(result) = "Result"

#Plotting valuation progression:
avgOptionValue = matrix(0, param$M, 1) #store cumulative mean payoff
for (kk in simStor[,"MC_iteration"]){avgOptionValue[kk] = mean(simStor[1:kk,"Thermal_valuation"])}
plot(c(avgOptionValue),type="l", ylab="Payoff", xlab="Number of Simulations", main = "Monte Carlo Convergence of Option Value", cex.main = .9)
legend("topright",legend=c("Simple Method"),col = c(1), lwd=c(1.5), lty=c(1), cex=0.8) #we haven't implemented antithetic valuation method.
abline(h = result$Call_value, col = "red")

#Distribution of MC payoffs
hist(simStor[,"Thermal_valuation"],xlab="Payoff Bin",main=paste("Frequency of",param$M,"Simulations"), breaks = 20)



