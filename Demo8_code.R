#Code for Demo 8: Value at Risk
load("../data/S8_ForwardCurve.rda")
S8_ForwardCurve
library(moments)
library(ENERGYRISK)

set.seed(1)
#contract parameters
C1_param = list(n1 = 1, K1 = S8_ForwardCurve[6, "F_0T"], T1 = 1, s1 = 2, F1 = S8_ForwardCurve[6, "F_0T"])
C2_param = list(n2 = 1, K2 = S8_ForwardCurve[10, "F_0T"], T2 = .25, s2 = .5, F2 = S8_ForwardCurve[10, "F_0T"])

nbMaturities = 6; nbFactors = 3; r = 0.05; dt = 1/252; Prob = 0.95; N = 500;
modelParams = list(nbMaturities = nbMaturities, nbFactors =nbFactors, r = r, dt = dt, Prob = Prob, N = N)

NSD = qnorm(modelParams$Prob) #95% quantile of SN distribution
# 2yr maturity volatility of oil
sigOil2y = sqrt(sum(S8_ForwardCurve[6,4:6]^2))
# 6 month maturity volatility of gas
sigGas6mo = sqrt(sum(S8_ForwardCurve[10, 4:6]^2))
#dot product scaled by volatilities:
rho = as.matrix(S8_ForwardCurve[6,4:6])%*%as.matrix(t(S8_ForwardCurve[10,4:6]))/(sigOil2y*sigGas6mo)
c(sigOil2y, sigGas6mo, rho)

#Modeling Oil and Gas prices:
errorOilGas = matrix(0,N,2)
for (i in 1:N){
  for(j in 1:2){
    errorOilGas[i,j] = randPolarRejc()
  }
}

#simulate future forward prices:
F_sim = matrix(0, N, 2) #column1 = Oil column2 = Gas
F_sim[,1] = C1_param$F1*exp(-0.5*(sigOil2y^2)*dt+sigOil2y*sqrt(dt)*(errorOilGas[,1]))
F_sim[,2] = C2_param$F2*exp(-0.5*(sigGas6mo^2)*dt+sigGas6mo*sqrt(dt)*((errorOilGas[,2])*rho + sqrt(1-rho^2)*errorOilGas[,1]))

colnames(F_sim) = c("F1", "F2")

# capture oil/gas risk factors
sigv1 = S8_ForwardCurve[1:6,4:6]
sigv2 = S8_ForwardCurve[7:12,4:6]
matv1 = S8_ForwardCurve[1:6,2]
matv2 = S8_ForwardCurve[7:12, 2]

#take forward prices and value European (call) option prices on the forwards
OptionValue_O = cbind(E_Option(o=1,g=1, C1_param$K1,C1_param$T1-dt,C1_param$s1-dt,C1_param$F1,sigv1,matv1,modelParams$nbMaturities,modelParams$nbFactors,modelParams$r),
                      E_Option(o=1,g=1,C2_param$K2,C2_param$T2-dt,C2_param$s2-dt,C2_param$F2,sigv2,matv2,modelParams$nbMaturities,modelParams$nbFactors,modelParams$r))
colnames(OptionValue_O) = c("C1","C2")
head(OptionValue_O)

OptionSim = matrix(0,N,2)
for(k in 1:N){OptionSim[k,1:2] = cbind(E_Option(o=1,g=1, C1_param$K1,C1_param$T1-dt,C1_param$s1-dt,F_sim[k,1],sigv1,matv1,modelParams$nbMaturities,modelParams$nbFactors,modelParams$r),
                      E_Option(o=1,g=1,C2_param$K2,C2_param$T2-dt,C2_param$s2-dt,F_sim[k,2],sigv2,matv2,modelParams$nbMaturities,modelParams$nbFactors,modelParams$r))}

colnames(OptionSim)=c("C1","C2")

#store portfolio data:
PortValue_O = sum(OptionValue_O) #original portfolio value
port = matrix(0,N,2)
port[,1] = apply(X = OptionSim, MARGIN = 1, FUN = sum)
port[,2] = log(port[,1]) - log(PortValue_O)

colnames(port) = c("P_Value", "P_return")


# GET VaR:
#get the 25th (500*0.05) smallest portfolio return -
VaR = port[order(port[,"P_return"], decreasing = FALSE)[N*(1-modelParams$Prob)], "P_return"]
result = c("mean" = mean(port[,2]), "var" = var(port[,2]), "sd" = sd(port[,2]), "skew" = skewness(port[,2]),
           "kurt" = kurtosis(port[,2]), "SE" = sd(port[,2])/sqrt(dt), "VaR" = abs(VaR))
hist(port[,2],xlab="Port Returns",main=paste("Probability Based on", N, "Simulations"),breaks=20)
abline(v = VaR, col = 'red')

