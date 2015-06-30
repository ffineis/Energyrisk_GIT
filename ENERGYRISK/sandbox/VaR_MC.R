library(ENERGYRISK)
library(lubridate)
library(moments)

# Load Data
#save.image(file="data/S8_ForwardCurve.rda")
data(S8_ForwardCurve)
S8_ForwardCurve

# Specify model params for MC
set.seed(1)

# Create a list of simulation variables
# 1) Set option on future paramters for oil
n1=1;K1=20.11;T1=1;s1=2;F1=20.11;
C1_param = list(n1,K1,T1,s1,F1)
names(C1_param) = c("n1","K1", "T1","s1","F1")
t(C1_param)
# 2) Set option on future paramters for natural gas
n2=1;K2=2.35;T2=0.25;s2=0.5;F2=2.35;
C2_param = list(n2,K2,T2,s2,F2)
names(C2_param) = c("n2","K2", "T2","s2","F2")
t(C2_param)
# 3) Set model Parameters
nbMaturities =6;nbFactors=3; r=0.05; delta_t=1/252; Prob=0.95;N=500;
modelParams = list(nbMaturities,nbFactors,r,delta_t,Prob,N)

# Time step delta_t is daily, which allows for a day-ahead VaR estimation
names(modelParams) = c("nbMaturities","nbFactors", "r","delta_t","Prob","N")
t(modelParams)

# Number of standard deviation
NSD = qnorm(modelParams$Prob,0,1)

# Estimate future and option initial values
# i) Estimate 2 y maturity volatilty OIL
sigOil2y = sqrt(sum(S8_ForwardCurve[6,4:6]^2))
# ii) Estimate 6 m maturity volatilty GAS
sigGas6m = sqrt(sum(S8_ForwardCurve[10,4:6]^2))
# iii) Estimate correlation between these two contracts
rho = as.matrix(S8_ForwardCurve[6,4:6])%*%
  as.matrix(t(S8_ForwardCurve[10,4:6]))/(sigOil2y*sigGas6m)
c(sigOil2y,sigGas6m,rho)

# 1) Estimate the error terms for oil an gas
errorOilGas = matrix(0,N,2)
for(i in 1:N){for(j in 1:2){errorOilGas[i,j] = randPolarRejc()}}
# 2) Estimate futures paths (& returns)
F_sim = matrix(0,N,2) 
F_sim[,1] = F1*exp(-0.5*sigOil2y^2*delta_t+sigOil2y*sqrt(delta_t)*(errorOilGas[,1]))
F_sim[,2] = F2*exp(-0.5*sigGas6m^2*delta_t+sigGas6m*sqrt(delta_t)*(errorOilGas[,2]))
colnames(F_sim) = c("F1","F2")
head(F_sim)

# capture oil/gas risk factors
sigv1 = S8_ForwardCurve[1:6,4:6]
sigv2 = S8_ForwardCurve[7:12,4:6]
# capture maturities
matv1 = S8_ForwardCurve[1:6,2]
matv2 = S8_ForwardCurve[7:12,2]

# 2) Option Value
# For comparison purposes estimate the current option/portfolio values
# Substract delta_t because we are conduct a day-ahead VaR
# o is the type of option (call or put)
# g can either return the option payoff or the discounted total value that 
# is attributable to the states in which the forward has a value greater than 
# or equal to the exercise price
OptionValue_0 = cbind(E_Option(o=1,g=1,K1,T1-delta_t,s1-
                 delta_t,F1,sigv1,matv1,nbMaturities,nbFactors,r),
                     E_Option(o=1,g=1,K2,T2-delta_t,s2-
                 delta_t,F2,sigv2,matv2,nbMaturities,nbFactors,r))
colnames(OptionValue_0) = c("C1","C2")
OptionValue_0
PortValue_0 = sum(OptionValue_0)
PortValue_0

# 3) Simulate the 1-Days ahead horizon changes: T1/s1 - t is the time to maturity
OptionSim =  matrix(0,N,2) 
for(k in 1:N){OptionSim[k,1:2] = 
                cbind(E_Option(1,1,K1,T1-delta_t,s1-
                           delta_t,F_sim[k,1],sigv1,matv1,nbMaturities,nbFactors,r),
                  E_Option(1,1,K2,T2-delta_t,s2-
                             delta_t,F_sim[k,2],sigv2,matv2,nbMaturities,nbFactors,r))}
colnames(OptionSim)=c("C1","C2")
head(OptionSim)

# 4) Estimate portfolio value
port = matrix(0,N,2)
port[,1] = rowSums(OptionSim)
# Portfolio returns
port[,2] = log(port[,1])-log(PortValue_0)
colnames(port) = c("Port","Port Change")
head(port)

# 5) Estimate VaR @ 95% CI
VaR = port[order(port[,"Port Change"], decreasing=FALSE)[N*(1-Prob)],"Port Change"]
# To convert daily volatility into annual multiple by sqrt(252) or divid by sqrt(delta_t)
# Load library(moments) for skewness and kurtosis estimation
result = c(mean(port[,2]),var(port[,2]),sd(port[,2]),skewness(port[,2]),
           kurtosis(port[,2]), sd(port[,2])/sqrt(delta_t),abs(VaR))
names(result)=c("mean","var","stdev","skewness","kurtosis","Annual Vol","VaR")
result

# Plot frequency distribution of the simulated payoffs: breaks = 20
hist(port[,2],xlab="Port Returns",main=paste("Probability Based on",
                                             N,"Simulations"),breaks=20)

