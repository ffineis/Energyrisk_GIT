library(ENERGYRISK)
### The spot energy price and the process parameters are chosen:
param = list(21.05,21.05,1.196,3.053,0.529,0.5,1,0.1,10)
names(param) = c("S","K","alpha","mu","sigma","T","s","r","N")
t(param)

# Estimate variable from parameters
lnS = log(param$S)
# Estimate discount factor for a 6 month option
P = exp(-param$r*param$T)
# The forward price (F(T,s)) from equation (3)
F = FT(param$S, param$alpha, param$mu, param$sigma, param$s)
dt = param$T/param$N
sigrdt = param$sigma*sqrt(dt)
var = list(lnS,P,F,dt,sigrdt)
names(var) = c("lnS","P","F","dt","sigrdt")
t(var)

# Generate spot price
# 1) Generate error term using the polar rejection method
x2_available <<- 0
nSims = 100
# Set seed so that part 1, 2 and 3 share the same randomized numbers
set.seed(1)
time = as.matrix(seq(0,0.5,0.05))
nObs = length(seq(0,0.5,0.05)) # or param$N
lnS1 = matrix(0,nSims,nObs)
colnames(lnS1) = t(time)

lnS1[,1]= lnS
# 2) Loop for number of simulations
for(i in 1:nSims){
  # Simulate ln(spot price) path
  for(j in 1:param$N){
    setRand = randPolarRejc()
    # Use a quick lnS path function to create the price path
    lnS1[i,j+1] = lnS(lnS1[i,j], param, sigrdt, dt, setRand)    
  }
}
# Spot price path
S1 = exp(lnS1)
head(S1)
# Plot the first three simulated price paths
plot(time,S1[1,],type="l",xlab="Time",ylab="S1",ylim=c(min(S1[1:4,]),max(S1[1:4,])))
lines(time,S1[2,],type="l", col="red", lwd=2, lty=3)
lines(time,S1[3,],type="l", col="blue", lwd=3, lty=4)
lines(time,S1[4,],type="l", col="green", lwd=4, lty=5)
title("Simulated Price Paths")

# Initialize the max b/w zero and Terminal (S-K)
zero = matrix(0,nSims,1)
# 3) Estimate Payoff for every simulation
payoff = (apply(cbind(zero,S1[,11]-param$K), 1, max))
head(payoff)
# Standard Error
se = sd(payoff)/sqrt(nSims)
# 4) Call Value discounted back
Call_Value = mean(payoff)*P
Result = data.frame(c(Call_Value,se))
colnames(Result) = c("Result")
rownames(Result) = c("Call_Value","se")
Result

# Estimate cumulative average payoff over simulated paths 
storage = matrix(0,nSims,1)
for(k in 1:nSims){storage[k] = sum(payoff[1:k])/k}
plot(storage,type="l", ylab="Payoff", xlab="Number of Simulations")
title("Monte Carlo European Call Value")


### ANTITHETIC: Antithetic Approach
# Use Result from part 1
# Set seed so that part 1, 2 and 3 share the same randomized numbers
set.seed(1)
lnS2 = matrix(0,nSims,nObs)
colnames(lnS2) = t(time)
lnS2[,1]= lnS
# 2) Loop for number of simulations
for(i in 1:nSims){
  # Simulate ln(spot price) path
  for(j in 1:param$N){
    setRand = randPolarRejc()
    # Use a quick lnS path function to create the price path
    lnS2[i,j+1] = lnS(lnS2[i,j], param,-sigrdt, dt, setRand)
  }
}
# Spot price path
S2 = exp(lnS2)
head(S2)
# Initialize the max b/w zero and Terminal (S-K)
zero = matrix(0,nSims,1)
# Estimate Payoff for every simulation & antithetic: use average
payoff_S1 = apply(cbind(zero,S1[,11]-param$K), 1, max)
payoff_S2 = apply(cbind(zero,S2[,11]-param$K), 1, max)
avg_payoff = apply(cbind(payoff_S1,payoff_S2),1,mean)
head(avg_payoff)

# Standard Error
se = sd(avg_payoff)/sqrt(nSims)
# 4) Call Value discounted back
Call_Value_A = mean(avg_payoff)*P
Result = data.frame(c(Call_Value_A,se))
colnames(Result) = c("Result")
rownames(Result) = c("Call_Value_A","se")
Result
# Estimate cumulative average payoff over simulated paths 
avg_storage = matrix(0,nSims,1)
for(k in 1:nSims){avg_storage[k] = mean(avg_payoff[1:k])}
plot(storage,type="l", ylab="Payoff", xlab="Number of Simulations",
     ylim=c(min(avg_storage),max(storage)))
lines(avg_storage,type="l", col="red", lwd=2, lty=3)
legend("topright",legend=c("Simple","Antithetic"),
       col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)
title("Monte Carlo European Call Value")


### GREEKs: Estimate Greeks from part 1
# Add change of delta to parameter variable
# Set seed so that part 1, 2 and 3 share the same randomized numbers
set.seed(1)
param = c(param,deltaS = 0.01) #shouldn't this be 1% of S?
lnS_G = matrix(0,nSims,nObs);lnS_G_PDS = matrix(0,nSims,nObs);lnS_G_NDS = matrix(0,nSims,nObs)
# 1) Estimate S then +/- delta, the convert back to ln(S)
lnS_G[,1]= lnS;lnS_G_PDS[,1] = log(exp(lnS) + param$deltaS)
lnS_G_NDS[,1] = log(exp(lnS) - param$deltaS)
# Loop for number of simulations
for(i in 1:nSims){
  # Simulate ln(spot price) path
  for(j in 1:param$N){
    setRand = randPolarRejc()
    # Use a quick lnS path function to create the price path
    lnS_G[i,j+1] = lnS(lnS_G[i,j], param, sigrdt, dt, setRand)
    lnS_G_PDS[i,j+1] = lnS(lnS_G_PDS[i,j], param, sigrdt, dt, setRand)
    lnS_G_NDS[i,j+1] = lnS(lnS_G_NDS[i,j], param, sigrdt, dt, setRand) 
  }
}
# Initialize the max b/w zero and Terminal (S-K)
zero = matrix(0,nSims,1)
# Estimate Payoff for every simulation & antithetic: use average
payoff_G = ((exp(lnS_G[,11])-param$K))
payoff_G_PDS = ((exp(lnS_G_PDS[,11])-param$K))
payoff_G_NDS = ((exp(lnS_G_NDS[,11])-param$K))
# Estimate Call_Value, Delta and Gamma
# 2) estimate Payoff then Payoff w/ +/- delta
payoff_G = apply(cbind(zero,payoff_G), 1, max)
payoff_G_PDS = apply(cbind(zero,payoff_G_PDS), 1, max)
payoff_G_NDS = apply(cbind(zero,payoff_G_NDS), 1, max)
head(cbind(payoff_G,payoff_G_PDS,payoff_G_NDS))

# Standard Error
se = sd(payoff_G)/sqrt(nSims)
# 4) Call Value then  w/ +/- delta
Call_Value = mean(payoff_G)*P
Call_Value_PDS = mean(payoff_G_PDS)*P
Call_Value_NDS = mean(payoff_G_NDS)*P

Result = data.frame(c(Call_Value,Call_Value_NDS,Call_Value_PDS,se))
colnames(Result) = c("Result")
rownames(Result) = c("Call_Value","Call_Value_NDS","Call_Value_PDS","se")
Result
# Estimate cumulative average payoff over simulated paths 
greek_storage = matrix(0,nSims,1)
for(k in 1:nSims){greek_storage[k] = sum(payoff_G[1:k])/k}
plot(greek_storage,type="l", xlab="Payoff", ylab="Number of Simulations")
title("Monte Carlo European Call Value")

# 5) Estimate delta & gamma
delta = (Call_Value_PDS-Call_Value_NDS)/(2*param$deltaS)
gamma = (Call_Value_PDS - 2*Call_Value + Call_Value_NDS)/(param$deltaS^2)
Greeks = data.frame(c(delta,gamma))
colnames(Greeks) = c("Result")
rownames(Greeks) = c("Delta","Gamma")

# 6) Monte Carlo result for delta and gamma
delta_storage = matrix(0,nSims,1)
gamma_storage = matrix(0,nSims,1)
for(k in 1:nSims){
  delta_storage[k] = (mean(payoff_G_PDS[1:k])-mean(payoff_G_NDS[1:k]))/(2*param$deltaS)
  gamma_storage[k] = (mean(payoff_G_PDS[1:k]) 
                      - 2*mean(payoff_G[1:k]) + mean(payoff_G_NDS[1:k]))/(param$deltaS^2)
}
# Show the convergent 6 last observations for delta and gamma
Greek_ill = cbind(delta_storage,gamma_storage)
colnames(Greek_ill) = c("Delta","Gamma")
tail(Greek_ill)