### CFaR demo 9 code
library(ENERGYRISK)
library(zoo)

S9_ForwardCurve = read.csv("/Users/fineiskid/Desktop/AMATH_Summer_UW/CFRM_520/EnergyRisk/Energyrisk_GIT/ENERGYRISK/data/csv/S9_ForwardCurve.csv")
S9_Params = read.csv("/Users/fineiskid/Desktop/AMATH_Summer_UW/CFRM_520/EnergyRisk/Energyrisk_GIT/ENERGYRISK/data/csv/S9_ModelParams.csv")
S9_ForwardCurve$Date <- as.Date(S9_ForwardCurve$Date, "%m/%d/%y")
set.seed(1)

N = 200
valDate = as.Date("01/01/09","%m/%d/%y")
valEnd = as.Date("12/31/09","%m/%d/%y")

# 364 days for simulation
nbDays = as.numeric(valEnd-valDate)
#we'll need daily Power, Gass, and Heating Oil prices
storDates = as.Date(S9_ForwardCurve$Date[1]:S9_ForwardCurve$Date[nrow(S9_ForwardCurve)])

#More parameters:
KCap = 90; KFloor = 70; KSwap = 80;
VOM = 1.50; Emission = 10; heatRate = 7.25; r = 0.05
K = VOM+Emission
scaleFactor = 0.1538 #convert dollars/barrel to $/MMBtu


############################### Simulations for Problems 1, 2, 3, 4 #################################

#store simulated prices, payoffs
valueStor = matrix(0, nbDays, 8)
colnames(valueStor) = c("Power", "Gas", "HO", "Cap", "Floor", "Swap", "Fuel_Switching", "Collar")

#beginning days of months: for month bucketz
indexStartDates = c(1, which(as.numeric(format(storDates[-length(storDates)],"%m"))!=as.numeric(format(storDates[-1],"%m")))+1)[-13]
bucketStartDates = storDates[indexStartDates]

#ending days of months: for month bucketzzzz
indexEndDates = c(which(as.numeric(format(storDates[-length(storDates)],"%m"))!= as.numeric(format(storDates[-1],"%m")))[-12],364)
bucketEndDates = storDates[indexEndDates]
head(cbind(indexStartDates,indexEndDates))


#store the buckets, 1 bucket for each month for each simulation, each problem gets N = 200 rows:
bucketStor = matrix(0,4*N,length(indexStartDates)) #simulation (rows) by month bucket (columns)
a = 0; b = 0;
if(S9_Params[2,"rho_Power"]<1){
  a = (S9_Params[2,"rho_HO"]-S9_Params[3,"rho_Power"]*S9_Params[2,"rho_Power"])/sqrt(1 -S9_Params[2,"rho_Power"]^2 )
}
if((S9_Params[3,"rho_Power"]^2+a^2)<=1){b = sqrt(1-S9_Params[3,"rho_Power"]^2-a^2)}

dt = 1/365

for(j in 1:N){
  #beginning prices: S0 = first forward price
  valueStor[1,"Power"] = S9_ForwardCurve$Power[1]
  valueStor[1,"Gas"] = S9_ForwardCurve$Gas[1]
  valueStor[1,"HO"] = S9_ForwardCurve$Heating.Oil[1]
  valueStor[1,"Cap"] = max(0, valueStor[1,"Power"] - KCap) #<- no need to present-value the present
  valueStor[1,"Floor"] = -max(0, KFloor- valueStor[1,"Power"]) #payoff on floor when price below Floor
  valueStor[1,"Swap"] = (valueStor[1,"Power"] -KSwap) #futures payoff: make money when spot > strike, losing o/w
  valueStor[1,"Fuel_Switching"] = max(0, valueStor[1,"Power"] - (heatRate*min(valueStor[1,"Gas"], scaleFactor*valueStor[1,"HO"]))-K)
  #collar: make money if S>90, S<70
  if(valueStor[1,"Power"] < KFloor){
      valueStor[1, "Collar"] = -max(0,KFloor- valueStor[1,"Power"])*
      exp(-r*(as.numeric(storDates[i+1]-valDate))/365)}
  if(valueStor[1, "Power"] > KCap){
      valueStor[1, "Collar"] = max(0,valueStor[1,"Power"]-KCap)*
      exp(-r*(as.numeric(storDates[i+1]-valDate))/365)}
  
  for(i in 1:(nbDays-1)){
    trackMonth = which(format(S9_ForwardCurve$Date,"%m%y")==format(storDates[i+1],"%m%y"))
    #current month's forward price
    F_t1 = c(S9_ForwardCurve[trackMonth,"Power"],S9_ForwardCurve[trackMonth,"Gas"],S9_ForwardCurve[trackMonth,"Heating.Oil"])
    #next month's forward price:
    F_t2 = c(S9_ForwardCurve[trackMonth+1,"Power"],S9_ForwardCurve[trackMonth+1,"Gas"],S9_ForwardCurve[trackMonth+1,"Heating.Oil"])
    t1 = i/nbDays
    
    #standard normal RV's:
    rand1 = qnorm(runif(1,min=0,max=1))
    rand2 = qnorm(runif(1,min=0,max=1))
    rand3 = qnorm(runif(1,min=0,max=1))
    dz1 = sqrt(dt)*rand1
    dz2 = S9_Params[2, "rho_Power"]*sqrt(dt)*rand1+
      sqrt(1-S9_Params[2, "rho_Power"]^2)*sqrt(dt)*rand2
    dz3 = S9_Params[3, "rho_Power"]*sqrt(dt)*rand1 + a*sqrt(dt)*rand2+b*sqrt(dt)*rand3
    
    #get mu:
    muPower = (1/S9_Params[1,"alpha"])*
      (log(F_t2[1])-log(F_t1[1]))/(1/12) +
      log(F_t1[1]) + S9_Params[1,"sigma"]^2 /
      (4*S9_Params[1,"alpha"] ) *(1 - exp(-2*S9_Params[1,"alpha"]*t1))
    muGas = (log(F_t2[2])-log(F_t1[2]))/(1/12) +
      log(F_t1[2]) + S9_Params[2,"sigma"]^2 /
      (4*S9_Params[2,"alpha"] ) *(1 - exp(-2*S9_Params[2,"alpha"]*t1))
    muHO = 
      (log(F_t2[3])-log(F_t1[3]))/(1/12) +
      log(F_t1[3]) + S9_Params[3,"sigma"]^2 /
      (4*S9_Params[3,"alpha"] ) *(1 - exp(-2*S9_Params[3,"alpha"]*t1))
    
    tempPower = log(valueStor[i,"Power"]) + (S9_Params[1,"alpha"]*
                                               (muPower-log(valueStor[i,"Power"])) -
                                               0.5 * S9_Params[1,"sigma"]^2) * dt + S9_Params[1,"sigma"]* dz1
    tempGas = log(valueStor[i,"Gas"]) + (S9_Params[2,"alpha"]*
                                           (muGas-log(valueStor[i,"Gas"])) -
                                           0.5 * S9_Params[2,"sigma"]^2) * dt + S9_Params[2,"sigma"]* dz2
    tempHO = log(valueStor[i,"HO"]) + (S9_Params[3,"alpha"]*
                                         (muHO-log(valueStor[i,"HO"])) -
                                         0.5 * S9_Params[3,"sigma"]^2) * dt + S9_Params[3,"sigma"]* dz3
    
    # Store values to be used on the next iteration
    valueStor[i+1,"Power"] = exp(tempPower)
    valueStor[i+1,"Gas"] = exp(tempGas)
    valueStor[i+1,"HO"] = exp(tempHO)
    valueStor[i+1,"Cap"] = max(0,valueStor[i+1,"Power"]-KCap)*
        exp(-r*(as.numeric(storDates[i+1]-valDate))/365)
    valueStor[i+1,"Floor"] = -max(0,KFloor- valueStor[i+1,"Power"])*
        exp(-r*(as.numeric(storDates[i+1]-valDate))/365)
    valueStor[i+1,"Swap"] = (valueStor[i+1,"Power"] -KSwap)*
      exp(-r*(as.numeric(storDates[i+1]-valDate))/365)
    valueStor[i+1,"Fuel_Switching"] = max(0, valueStor[i+1,"Power"] - (heatRate*min(valueStor[i+1,"Gas"], scaleFactor*valueStor[i+1,"HO"]))-K)*
      exp(-r*(as.numeric(storDates[i+1]-valDate))/365)
      
    #collar is long a cap (max(0, S-K)) and short a floor (-max(0, K-S))
    if(valueStor[i+1,"Power"] < KFloor){
        valueStor[i+1, "Collar"] = -max(0,KFloor- valueStor[i+1,"Power"])*
        exp(-r*(as.numeric(storDates[i+1]-valDate))/365)}
    if(valueStor[i+1, "Power"] > KCap){
        valueStor[i+1, "Collar"] = max(0,valueStor[i+1,"Power"]-KCap)*
        exp(-r*(as.numeric(storDates[i+1]-valDate))/365)}
    
  }
  # Sum all the bucket cash flows for (1) Swap, (2) Fuel_Switching, (3) Collar, (4) All assets combined
  for(k in 1:length(indexStartDates)){
      bucketStor[j,k] = sum(valueStor[indexStartDates[k]:indexEndDates[k],6])
      bucketStor[j+N, k] = sum(valueStor[indexStartDates[k]:indexEndDates[k],7])
      bucketStor[j+(2*N), k] = sum(valueStor[indexStartDates[k]:indexEndDates[k],8])
      bucketStor[j+(3*N), k] = sum(valueStor[indexStartDates[k]:indexEndDates[k],4:8])
  }
}

# check out prices for power, gas, heating oil, and call/floor options on power for 1 simulation:
head(valueStor)

#checkout simulated cashflows based on months (there are 12 months)
head(bucketStor)

############################################ PAYOFFS ###############################################

#considering that each column of bucketStor has 200 simulated monthly payoffs,
#then column averages represent expected monthly payoff:
Prob1 = t(as.data.frame(colMeans(bucketStor[1:N,])))
Prob2 = t(as.data.frame(colMeans(bucketStor[(N+1):(2*N),])))
Prob3 = t(as.data.frame(colMeans(bucketStor[((2*N)+1):(3*N),])))
Prob4 = t(as.data.frame(colMeans(bucketStor[((3*N)+1):(4*N),])))
results = list(Prob1 = Prob1, Prob2 = Prob2, Prob3 = Prob3, Prob4 = Prob4)
for (k in 1:4){rownames(results[[k]]) <- "Cashflow"; colnames(results[[k]]) <- unique(months(storDates))}

# Standard Error per Bucket in portfolio:
se = apply(bucketStor[bucketStor((3*N)+1):(4*N), MARGIN = 2, FUN = sd)/sqrt(N)

#Plotting EaR:
percentile = matrix(0,10,length(indexEndDates)) #storage matrix
for(t in 1:length(indexEndDates)){
  percentile[,t] = quantile(bucketStor[((3*N)+1):(4*N),t],c(seq(0.05,0.95,0.10)))
}

#plot 5% empirical quantile:
plot(c(percentile[1,]),type="l", ylab="Earning at Risk", xlab="Months",ylim=c(min(percentile),max(percentile)),
     main = "Portfolio Earnings at Risk\nJan 1, 2009 to Jan 31, 2009")
#plot all other empirical quantiles:
for(v in 2:10){lines(percentile[v,],col=v,lwd=1, lty=v)}
legend("topright",legend=c("5%","15%","85%","95%"),
       col = c(1,2,9,10), lwd=c(1.5,1.5,1.5,1.5), lty=c(1,2,9,10), cex=0.8)

##################################################################################################





