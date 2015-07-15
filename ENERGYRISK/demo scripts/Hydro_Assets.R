library(ENERGYRISK)
library(lubridate)
library(zoo)
# If you want to convert a .txt file into .rda then save and reload data using data()
#S6_ForwardCurve$Date = as.Date(S6_ForwardCurve$Date,"%m/%d/%y") 
#save.image(file="data/data_S6.rda")
# Load storage contract specs and forward curve
data(data_S7)
# Illustrate Model Parameters
S7_ModelParams
# Setup costs (strike), heat rate and interest rates
# Number of simulations
N = 100
# Set start and end dates for theforward curve
valDate = as.Date("01/01/09","%m/%d/%y")
valEnd = as.Date("03/31/09","%m/%d/%y")
c(valDate,valEnd)
# Store the number of days the study will cover
storDates =  as.Date(valDate:valEnd)
# Store the day for the week to be used for off-Peak/Peak
storWkday = as.POSIXlt(storDates,format="%Y-%m-%d")
storWkday = storWkday$wday

# Estimate number of days for simulation: beg-end of forward strip date 
nbDays = as.numeric(valEnd-valDate)
# 24 hours in a day 48 (2*24) 30 minute intervals
nb30mn = nbDays *24*2

# Costs (strike)
Effeciency = 0.7222; F_Peak = 55; F_OPeak=25;r=0.05
# Create a list of simulation variables
param = list(Effeciency,F_Peak,F_OPeak,r,N)
names(param) = c("Effeciency","F_Peak", "F_OPeak","r","N")
t(param)

# Define the Period/Power/PumpGen/PumpingCost/GenRevenue matrix
valueStor = matrix(0,nb30mn,7)
colnames(valueStor) = c("Period","Power","PumpGen","PumpingCost",
                        "GenRevenue","volPumped","CT")
# Setup Period 1:48 to cover all 30 minute intervals in a day
valueStor[,"Period"] = rep(1:48,nb30mn/48)
# Number of simulation matrix #sim/OptionValue
simMatrix = matrix(0,param$N,2)
# step matrix for the whole period of the model
dt =1/(365*2*24)
for(j in 1:N){
  # Define the Period/Power/PumpGen/PumpingCost/GenRevenue matrix
  valueStor = matrix(0,nb30mn,7)
  colnames(valueStor) = c("Period","Power","PumpGen","PumpingCost",
                          "GenRevenue","volPumped","CT")
  # Setup Period 1:48 to cover all 30 minute intervals in a day
  valueStor[,"Period"] = rep(1:48,nb30mn/48)
  set.seed(j)
  # Initialize starting values & time step dt
  valueStor[1,"Power"] = param$F_OPeak
  F_t1 = param$F_OPeak
  valueStor[1,"PumpGen"] = -1
  offsetDate = 1
  # Indexing for a new day: used for 3) & 4) etc.
  offset = -48
  trackDate = 1
  ## Model price changes
  for(i in 1:(nb30mn-1)){
    # Indexing for a new day
    if(valueStor[i,"Period"] ==1){offset =offset+48}
    
    # Capture Holiday/Weekend & time of day o.w. peak price is charged
    if(storDates[offsetDate] =="2009-01-01"){muPower=param$F_OPeak
    }else if(storWkday[trackDate]==6||storWkday[trackDate]==0){muPower = param$F_OPeak
    }else if(valueStor[i,"Period"]>44 || valueStor[i,"Period"]<15){muPower = param$F_OPeak
    }else{muPower = param$F_Peak}
    
    # Keep track of estimated F_Price
    F_t2 = muPower
    
    # step matrix for the whole period of the model
    t1 = i/(365*2*24)
    # Set seed so that the random variable does not change
    # Use qnorm: give it a probability and returns the number whose 
    # cumulative distribution matches the probability
    rand1 = randPolarRejc()
    rand2 = randPolarRejc()
    # stochastic component estimation
    dz1 = sqrt(dt)*rand1
    # Capture the stochastic gas compnent which is correlated to power demand
    dz2 = sqrt(dt)*rand2
    
    # Test whether the random draw allow for a jump
    jFrequency = as.numeric(S7_ModelParams[1,"jFrequency"]*dt <qnorm(runif(1,min=0,max=1)))
    # Estimate mu from model 
    muPower = (1/S7_ModelParams[1,"alpha"])*
      (log(F_t2)-log(F_t1))/(dt) +
      log(F_t1) + S7_ModelParams[1,"sigma"]^2 /
      (4*S7_ModelParams[1,"alpha"] ) *(1 - exp(-2*S7_ModelParams[1,"alpha"]*t1)) 
    
    tempPower = log(valueStor[i,"Power"]) + (S7_ModelParams[1,"alpha"]*
    (muPower-log(valueStor[i,"Power"])) - 0.5 * S7_ModelParams[1,"sigma"]^2)*
      dt + S7_ModelParams[1,"sigma"]* dz1 +
      jFrequency*S7_ModelParams[1,"volJump"]* dz2
    
    # 1) Store values to be used on the next iteration
    valueStor[i+1,"Power"] = exp(tempPower)
    # 2) Pumping/Generating: determines when the unit should pump or generate 
    # given the prices and assuming a daily cycle
    # Utilize the off-peak critirion
    if(i==(48+offset-1)){
    valueStor[(1:14)+offset,"PumpGen"] = -1; valueStor[(45:48)+offset,"PumpGen"]=-1;
    # Use 30 b/c that is the amount we want to index by in the peak 
    # period minus current period
    # Capture the 13th largest number to match the effeciency rate of
    # the plant 13 (peak) /18(oPeak)=72.22%. Decide when in the peak period to generate power.
    valueStor[(15:44)+offset,"PumpGen"] = as.numeric(valueStor[(15:44)+offset,"Power"]>=
    valueStor[order(valueStor[(15:44)+offset,"Power"], decreasing=TRUE)[13]+offset+14,"Power"])
    }

    # 3) Pumping Cost: pumping cost is calculated as the average overnight price 
    # Use 14 capture the overnight rate of 22:00-6:30: everytime pumping is set to -1
    if(i==(48-1)){
      valueStor[(1:48)+offset,"PumpingCost"] = sum(valueStor[which(valueStor[(1:14)+offset,
      "PumpGen"]==-1),"Power"])/length(which(valueStor[(1:14)+offset,"PumpGen"]==-1))
      }else if(i==(48+offset-1)){
        valueStor[(1:48)+offset,"PumpingCost"] = sum(valueStor[c((-3:14)+offset),"Power"])/
              length(c((-3:14)+offset))}

    # 4) Generating Revenue: everytime generation is set to 1
    if(i==(48+offset-1)){
        diviser = length(which(valueStor[(1:48)+offset,"PumpGen"]==1))
        # Catch the case when there is no generating occuring in a day
        if(length(which(valueStor[(1:48)+offset,"PumpGen"]==1))==0){diviser=1}
        valueStor[(1:48)+offset,"GenRevenue"] = sum(valueStor[which(valueStor[(1:48)+
        offset,"PumpGen"]==1)+offset,"Power"])/diviser
        # 5) Volume Pumped/Generated
        if(all((valueStor[(1:48)+offset,"PumpingCost"]/valueStor[(1:48)+offset,"GenRevenue"])<param$Effeciency)){
          valueStor[(1:48)+offset,"volPumped"] = valueStor[(1:48)+offset,"PumpGen"]*100
        }
        
        # 6) Option Value
        valueStor[(1:48+offset),"CT"] = valueStor[(1:48+offset),"Power"]*valueStor[(1:48)+offset,"volPumped"]*
          exp(-param$r*(as.numeric(as.Date(storDates[trackDate]) -as.Date(valDate))+
                          (valueStor[(1:48+offset),"Period"]-1)/48)/365)
        trackDate = trackDate +1
    }

    
    # Used to capture the Jan holiday structure 01/01/2009
    if(i==48){offsetDate = offsetDate +1}
    # Reset the initial starting futures price
    F_t1 = F_t2

    # A way to check your result is to estimate: Sum realeased/Sume Pumped
    # if effeciencyCheck is higher than 72.22% then there is a mistake
    sum(valueStor[which(valueStor[,"volPumped"]>0),"volPumped"])
    sum(valueStor[which(valueStor[,"volPumped"]<0),"volPumped"])
    effeciencyCheck = abs(sum(valueStor[which(valueStor[,"volPumped"]>0),"volPumped"])/ 
                            sum(valueStor[which(valueStor[,"volPumped"]<0),"volPumped"]))
  }
  # Option value
  simMatrix[j,] = c(j,sum(valueStor[,"CT"]))

}
# Illustrate how valueStor looks
head(valueStor)

#Illustrate the result from all N simulations
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
# Plot electricity prices for a simulated path: MRJD
plot(valueStor[,"Power"],type="l", col=2,ylab="Power Prices", xlab="Number Time Step")

# Plot CT and volume pumped/generated
# Optimal position
par(mar=c(4.5,3.5,3,6)+0.1)
plot(valueStor[,"volPumped"],type="l",xlab="Time",ylab="",col=2)
title("Volume Pumped/Generated & CT")
mtext("Volume",side=2,line=2,cex=0.75)
par(new=TRUE)
plot(valueStor[,"CT"],type="l",xlab="Time",ylab="",col=1,lwd=1, lty=1,yaxt="n")
axis(4, ylim=c(min(valueStor[,"CT"]),max(valueStor[,"CT"])),line=0,cex.axis=0.75)
mtext("CT Value",side=4,line=1.75,cex=0.75)
legend("bottomleft",legend=c("CT","Volume PumpGen"),
       col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)