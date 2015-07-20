# Demo 7 code: Valuing a pump hydro power plant

library(ENERGYRISK)
library(lubridate)
library(zoo)
source("./ENERGYRISK/R/ES1_F.R")
data(data_S7)
#get model parameters to project Power and rainfall evolution
head(S7_ModelParams)

N = 100 # number of simulations
valDate = as.Date("01/01/09", format = "%m/%d/%y")
valEnd = as.Date("03/31/09", format = "%m/%d/%y")
storDates = as.Date(valDate:valEnd) #dates included in valuation

storWkday = as.POSIXlt(storDates, format = "%Y-%m-%d")
storWkday = storWkday$wday
nbDays = as.numeric(valEnd-valDate) #89 days to determine prices
nb30min = nbDays*24*2 #number of 30-min intervals in day

Efficiency = 0.7222 # 72.22% of energy used to pump water up will be regained
F_Peak = 55; F_OPeak = 25; r = 0.05
param = list(Efficiency = Efficiency, F_Peak = F_Peak, F_OPeak = 25, r = r, N = N)

# store period, power, pump generation, pumping cost, and revenues:
valueStor = matrix(0, nb30min, 7) #one observation for every 30min interval in 3 month span.
#pump gen is -1, 0, or 1
#volume pumped is function of having a 100 MW unit
colnames(valueStor) = c("Period", "Power", "PumpGen", "PumpingCost", "GenRevenue", "volPumped", "CT")
simMatrix = matrix(0, param$N, 2)
dt = 1/(365*2*24)

for (j in 1:param$N){
  set.seed(j)
  valueStor[,"Period"] = rep(1:48, nb30min/48)
  valueStor[1,"Power"] = param$F_OPeak
  F_t1 = param$F_OPeak ##### IMPORTANT ######
  valueStor[1, "PumpGen"] = -1 #offpeak
  offsetDate = 1
  offset = -48 #index for new day
  trackDate = 1
  
  #simulate price changes
  for (i in 1:(nb30min-1)){
    if(valueStor[i, "Period"] == 1){offset = offset+48} #indexing can't start at -48
    if(storDates[offsetDate]=="2009-01-01"){muPower=param$F_OPeak}
    else if(storWkday[trackDate]==6||storWkday[trackDate]==0){muPower=param$F_OPeak}
    else if(valueStor[i,"Period"]>44||valueStor[i,"Period"]<15){muPower=param$F_OPeak}
    else{muPower=param$F_Peak}
    
    F_t2 = muPower
    
    t1=i/(365*2*24)
    rand1=randPolarRejc()
    rand2=randPolarRejc()
    
    dz1 = sqrt(dt)*rand1
    dz2 = sqrt(dt)*rand2
    
    #TRUE = 1, FALSE = 0
    jFrequency = as.numeric(S7_ModelParams$jFrequency[1]*dt < qnorm(runif(1,min = 0 ,max = 1)))
    #get mu from model:
    muPower = (1/S7_ModelParams$alpha[1])*(log(F_t2)-log(F_t1))/dt + log(F_t1) + S7_ModelParams$sigma[1]^2/(4*S7_ModelParams$alpha[1])*(1-exp(-2*S7_ModelParams$alpha[1]*t1))
    tempPower = log(valueStor[i,"Power"])+(S7_ModelParams[1,"alpha"]*(muPower-log(valueStor[i,"Power"]))-0.5*S7_ModelParams[1,"sigma"]^2)*dt + S7_ModelParams[1,"sigma"]* dz1 +
      jFrequency*S7_ModelParams[1,"volJump"]* dz2
    
    valueStor[i+1, "Power"] = exp(tempPower)
    if(i == (48+offset-1)){
      #offpeak (nighttime)
      valueStor[(1:14)+offset,"PumpGen"]=-1;valueStor[(45:48)+offset,"PumpGen"]=-1;
      #onpeak (daytime)
      #Because you can only generate for 13 segments in a day, choose top 13 highest prices to produce on.
      valueStor[(15:44)+offset,"PumpGen"]=as.numeric(valueStor[(15:44)+offset,"Power"]>=
        valueStor[order(valueStor[(15:44)+offset,"Power"],decreasing=TRUE)[13]+offset+14,"Power"])
    }
    #Pumping cost: average of overnight price (8:00 PM to 6:30 AM)
    if(i ==(47)){
      valueStor[(1:48)+offset, "PumpingCost"] = sum(valueStor[which(valueStor[(1:14)+offset, "PumpGen"]==-1), "Power"])/
        length(which(valueStor[(1:14)+offset, "PumpGen"]==-1))
    }
    else if(i == (48+offset-1)){
      valueStor[(1:48)+offset, "Pumping Cost"] = sum(valueStor[c((-3:14)+offset),"Power"])/length(c((-3:14)+offset))
    }

    #Calculate generating Revenue: every time generation is 1, add revenue
  if(i==(48+offset-1)){
    diviser = length(which(valueStor[(1:48)+offset,"PumpGen"]==1))
    #Catch the case when there is no generating occuring in a day
    if(diviser == 0){diviser <- 1}
    valueStor[(1:48)+offset,"GenRevenue"] = sum(valueStor[which(valueStor[(1:48)+offset,"PumpGen"]==1)+offset,"Power"])/diviser
    #5) Volume Pumped/Generated
    if(all((valueStor[(1:48)+offset,"PumpingCost"]/valueStor[(1:48)+offset,"GenRevenue"])< param$Effeciency)){
      valueStor[(1:48)+offset,"volPumped"] = valueStor[(1:48)+offset,"PumpGen"]*100
    }
    #6) Option Value
    valueStor[(1:48+offset),"CT"] = valueStor[(1:48+offset),"Power"]*valueStor[(1:48)+offset,"volPumped"]*
      exp(-param$r*(as.numeric(as.Date(storDates[trackDate]) - as.Date(valDate))+(valueStor[(1:48+offset),"Period"]-1)/48)/365)
    trackDate = trackDate +1
  }
   
  if(i == 48){offsetDate = offsetDate + 1}
  F_t1 = F_t2 #reset initial starting futures price once we've moved up a day
  #Check to see that we're meeting the efficiency cutoff of 72.22% (sum(released energy)/sum(pumped energy))
  Echeck = abs(sum(valueStor[which(valueStor[,"volPumped"]>0),"volPumped"]))/abs(sum(valueStor[which(valueStor[,"volPumped"]<0),"volPumped"]))
  #if (Echeck > param$Efficiency){print("WARNING: Efficiency is too high")}
  }
  
  #End up with N option values:
  simMatrix[j,]=c(j,sum(valueStor[,"CT"]))
}

