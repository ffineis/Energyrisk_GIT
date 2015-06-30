## Structure of an Option
# Set initial parameters
alpha = 0.338; sigma = 0.305404;
# N is the number of periods in a year
N = 12
# In months we have
delta_t = 1/12
# The relationship between the space step and the time step is given by:
delta_x = sigma * sqrt(delta_t*3)
param = list(alpha,sigma,delta_t,delta_x)
names(param) = c("alpha","sigma","delta_t","delta_x")
t(param)
# Based on the delta_t from above load a data.frame() w/ parameters for estimation
# Load data from S4_Params.txt file in directory
# If save .txt file into .rda and reload it using folder path data()
data(S4_Params)
# Illustrate first and last 6 obs to get a quick view of the data
S4_Params

# 1) Setup and option object that will pass all the option information 
# Put into a list and return
S0 = S4_Params["F_Price",]$X0; K = 21; r = S4_Params["R_ts",]$X0;
vol = param$sigma; r = S4_Params["R_ts",]$X0;  ttm = S4_Params["Maturity",]$X12;
# Create optionSpec
optionSpec = list(S0, K, ttm, r, vol)
names(optionSpec) = c("S0", "K", "maturity", "r","volatility")
t(optionSpec)

# 2) Simplified tree
# Size of up move is $delta_x
# Setup Risk neutral probability for j = 0, i = 0
# (i represents the time step, j the spot price level)
j = 0 ; i = 0 ;
prob = prob(j, delta_t, alpha, delta_x, vol)
prob

# Discount factor
df = exp(-r * delta_t)
# At the terminal node, there are N/2 asset values (6-month of the year)
# Level_x represents the max/min price changes that can occur for the terminal node
level_x = -c(-delta_x*((N/2):1),delta_x*(0:(N/2)))
level_x
j.index = seq(from=0, to=N/2, by=1)
nbNodes = seq(from=1,to=length(level_x),by=2)
# a) Estimate initialize state price accumulation (t = 0)
# b) Preliminary steps to estimating a_i which are chosen to ensure that 
# the tree correctly returns the observed forward price curve:
a_i = a_i(prob, delta_x, j.index, df, nbNodes)
a_i

## Option Valuation: EUROPEAN
# d) Estimate Price Paths
# Set Call/Put multiplier +/-
mult = 1
# i) Initialize the backwardation apporach: start w/ terminal node
F = exp(as.numeric(a_i[length(j.index)]) + level_x)
V = pmax(0, mult * (F- K))
# The following print out is useful for a trinomial function
# It will print out the results as the function is running to keep track
cat("Time step: ", N/2, "\n", sep="")
cat("Prices:\n")
print(F)
cat("Option Values:\n")
print(V)

i.index = seq(from=N/2-1, to=0, by=-1)
offset = 1
# ii) Utilize backwardation technique to go from t=T to t=0
for (i in i.index) {
  level_xt = -c(-delta_x*((i):1),delta_x*(0:i))
  # Load up probabilities to estimate expected value
  j  = level_xt
  prob = prob(j, delta_t, alpha, delta_x, vol)
  
  # Sub bind the expectation values together
  E_V = cbind(V[1:(length(V)-2)],V[2:(length(V)-1)],V[3:length(V)])
  # F is the vector of prices at each time step and node
  F = exp(as.numeric(a_i[length(j.index)-offset]) + level_xt)
  # Update the V vector of option values at each time step and node
  V = pmax(0, df * diag(E_V%*%prob))
    
  if (i==0){F = F[1]}
  # The following print out is useful for a trinomial function
  # It will print out the results as the function is running to keep track
    cat("Time step: ", i, "\n", sep="")
    cat("Prices:\n")
    print(F)
    cat("Option Values:\n")
    print(V)
  
  offset = offset +1 
}

## Option Valuation: AMERICAN
# e) Estimate Price Paths
# Set Call/Put multiplier +/-
mult = -1
# i) Initialize the backwardation apporach: start terminal node
F = exp(as.numeric(a_i[length(j.index)]) + level_x)
V = pmax(0, mult * (F- K))
# The following print out is useful for a trinomial function
# It will print out the results as the function is running to keep track
cat("Time step: ", N/2, "\n", sep="")
cat("Prices:\n")
print(F)
cat("Option Values:\n")
print(V)

i.index = seq(from=N/2-1, to=0, by=-1)
offset = 1
# ii) Utilize backwardation technique to go from t=T to t=0
for (i in i.index) {
  level_xt = -c(-delta_x*((i):1),delta_x*(0:i))
  # Load up probabilities to estimate expected value
  j  = level_xt
  prob = prob(j, delta_t, alpha, delta_x, vol)
  
  # Sub bind the expectation values together
  E_V = cbind(V[1:(length(V)-2)],V[2:(length(V)-1)],V[3:length(V)])
  # F is the vector of prices at each time step and node
  F = exp(as.numeric(a_i[length(j.index)-offset]) + level_xt)
  
  # Primary difference between EUROPEAN and AMERICAN Options
  # Update the V vector of option values at each time step and node
  V = pmax(mult*(F-K), df * diag(E_V%*%prob))
  
  if (i==0){F = F[1]; V= V[1]}
  # The following print out is useful for a trinomial function
  # It will print out the results as the function is running to keep track
  cat("Time step: ", i, "\n", sep="")
  cat("Prices:\n")
  print(F)
  cat("Option Values:\n")
  print(V)
  
  offset = offset +1 
}