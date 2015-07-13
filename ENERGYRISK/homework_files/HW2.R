### CFRM 520 HW2.
library(ENERGYRISK)
source("./R/Lattices.R")
data(S4_Params)

tree_Params <- list("alpha" = 0.338, "delta_t" = 0.08333333, "delta_x" = 0.152702,
                    "S0" = 21.05, "K" = 21, "maturity" = 1, "r" = 0.06, "volatility" = 0.305404)
a_is <- t(data.frame(c(3.046901, 3.041576, 3.034532, 3.027756, 3.021138, 3.014664, 3.008326)))
colnames(a_is) <- c("X0", "X1", "X2", "X3", "X4", "X5", "X6")
rownames(a_is) <- "a_i" 

#set up tree characteristics for the preliminary tree.
discount <- exp(-tree_Params$r*tree_Params$delta_t)
level_x = -c(-tree_Params$delta_x*((6):1),tree_Params$delta_x*(0:(6))) #trinomial tree prices at time N/2. Symmetric about 0.
j.index = seq(from=0, to=6, by=1) #month number, 0 to 6 by 1's
nbNodes = seq(from=1,to=length(level_x),by=2) #number of nodes at each time point


##### Problem 1: Value EUROPEAN PUT Option ####
mult = -1 #put option...
# i) Fitted forward prices, spot price data fitted correctly
Forw = exp(as.numeric(a_is[length(j.index)]) + level_x)
V = pmax(0, mult * (Forw- tree_Params$K))

# forward prices and option values at time = 0.5 yr
cat("Time step: ", 6, "\n", sep="")
cat("Prices:\n")
print(Forw)
cat("European Put Values:\n")
print(V)

i.index = seq(from=6-1, to=0, by=-1)
offset = 1

# backward induction steps:
for (i in i.index) {
  # moving backwards, get spot prices on preliminary tree at previous time
  level_xt = -c(-tree_Params$delta_x*((i):1),tree_Params$delta_x*(0:i))
  # Get transition probabilities to estimate expected value
  j  = level_xt
  #print(sprintf("level_xt: %s", level_xt))
  prob = prob(j, tree_Params$delta_t, tree_Params$alpha, tree_Params$delta_x, tree_Params$volatility)
  
  # Sub bind the expectation values together
  E_V = cbind(V[1:(length(V)-2)],V[2:(length(V)-1)],V[3:length(V)])
  # F is the vector of prices at each time step and node
  Forw = exp(as.numeric(a_is[length(j.index)-offset]) + level_xt)
  
  # Primary difference between EUROPEAN and AMERICAN Options
  # Update the V vector of option values at each time step and node
  V = pmax(0, discount * diag(E_V%*%prob))
  
  if (i==0){
    Forw = Forw[1]
    V= V[1]
    result1 = as.data.frame(V); rownames(result1) = "Result 1:"; colnames(result1) = "Option Value"
  }
  # The following print out is useful for a trinomial function
  # It will print out the results as the function is running to keep track
  cat("Time step: ", i, "\n", sep="")
  cat("Prices:\n")
  print(Forw)
  cat("European Put Values:\n")
  print(V)
  
  offset = offset +1 
}


##### Problem 2: Value AMERICAN CALL Option ####
mult = 1 #call, not put. Put is valued at K-S.
# i) Fitted forward prices, spot price data fitted correctly
Forw = exp(as.numeric(a_is[length(j.index)]) + level_x)
V = pmax(0, mult * (Forw- tree_Params$K))

# forward prices and option values at time = 0.5 yr
cat("Time step: ", 6, "\n", sep="")
cat("Prices:\n")
print(Forw)
cat("American Call Values:\n")
print(V)

i.index = seq(from=6-1, to=0, by=-1)
offset = 1

# backward induction steps:
for (i in i.index) {
  # moving backwards, get spot prices on preliminary tree at previous time
  level_xt = -c(-tree_Params$delta_x*((i):1),tree_Params$delta_x*(0:i))
  # Get transition probabilities to estimate expected value
  j  = level_xt
  #print(sprintf("level_xt: %s", level_xt))
  prob = prob(j, tree_Params$delta_t, tree_Params$alpha, tree_Params$delta_x, tree_Params$volatility)
  
  # Sub bind the expectation values together
  E_V = cbind(V[1:(length(V)-2)],V[2:(length(V)-1)],V[3:length(V)])
  # F is the vector of prices at each time step and node
  Forw = exp(as.numeric(a_is[length(j.index)-offset]) + level_xt)
  
  # Now value is either exercising option or EV of not exercising
  V = pmax(mult*(Forw-tree_Params$K), discount * diag(E_V%*%prob))
  
  if (i==0){
    Forw = Forw[1]
    V= V[1]
    result2 = as.data.frame(V); rownames(result2) = "Result 2:"; colnames(result2) = "Option Value"
  }
  # The following print out is useful for a trinomial function
  # It will print out the results as the function is running to keep track
  cat("Time step: ", i, "\n", sep="")
  cat("Prices:\n")
  print(Forw)
  cat("American Call Values:\n")
  print(V)
  
  offset = offset +1 
}


##### Problem 3: Value AMERICAN PUT Double Barrier Knock-out Option ####
mult = -1 #put, not call. Put is valued at K-S.
barrier = c(16, 30)
# i) Fitted forward prices, spot price data fitted correctly
Forw = exp(as.numeric(a_is[length(j.index)]) + level_x)
KOed = which(Forw<barrier[1] | Forw>barrier[2])
V = pmax(0, mult * (Forw- tree_Params$K))
V[KOed] <- 0

# forward prices and option values at time = 0.5 yr
cat("Time step: ", 6, "\n", sep="")
cat("Prices:\n")
print(Forw)
cat("American Call Values:\n")
print(V)

i.index = seq(from=6-1, to=0, by=-1)
offset = 1

# backward induction steps:
for (i in i.index) {
  # moving backwards, get spot prices on preliminary tree at previous time
  level_xt = -c(-tree_Params$delta_x*((i):1),tree_Params$delta_x*(0:i))
  # Get transition probabilities to estimate expected value
  j  = level_xt
  #print(sprintf("level_xt: %s", level_xt))
  prob = prob(j, tree_Params$delta_t, tree_Params$alpha, tree_Params$delta_x, tree_Params$volatility)
  
  # Sub bind the expectation values together
  E_V = cbind(V[1:(length(V)-2)],V[2:(length(V)-1)],V[3:length(V)])
  # F is the vector of prices at each time step and node
  Forw = exp(as.numeric(a_is[length(j.index)-offset]) + level_xt)
  KOed = which(Forw<barrier[1] | Forw>barrier[2])
  
  # Now value is either exercising option or EV of not exercising
  V = pmax(mult*(Forw-tree_Params$K), discount * diag(E_V%*%prob))
  V[KOed] <- 0
  
  if (i==0){
    Forw = Forw[1]
    V= V[1]
    result3 = as.data.frame(V); rownames(result3) = "Result 3:"; colnames(result3) = "Option Value"
  }
  # Print out the results as the function is running to keep track
  cat("Time step: ", i, "\n", sep="")
  cat("Prices:\n")
  print(Forw)
  cat("American Put/Barrier Values:\n")
  print(V)
  
  offset = offset +1 
}


##### Problem 4: 3-level (indexing at 0) trinomial tree solution #####
#What do we need to fill out the 3rd level of the tree?
# simplified tree price levels
# Option values (i.e. discounted expected value of option values level i = 4)
# transition probabilities to go from level 3 to level i = 4
# State prices at level 3.

params = S4_Params[,c(1:4)]
alpha = 0.338; sigma = 0.305404
dt = 0.08333333; dx = 0.152702; r = params["R_ts",1]

#simplified tree prices at level 2:
level_x2 = dx*seq(from = -2, to = 2, by = 1)
#simplified tree prices at level 3:
level_x3 = dx*seq(from = -3, to = 3, by = 1)

#number of nodes at 3rd level
nbNodes = 7

#transition probabilities function:
prob <- function(x, dt, alpha, dx, vol){ 
  pu =  1/2*((((vol^2)*dt + (alpha^2)*(x^2)*dt^2)/(dx^2)) -
          ((alpha*x*dt)/dx))
  pd = 1/2*((((vol^2)*dt + (alpha^2)*(x^2)*dt^2)/(dx^2)) +
          ((alpha*x*dt)/dx))
  pm = 1 - pu - pd
  prob = list(pu = pu, pm = pm, pd = pd)
  return(prob)
}

# transition probs, from level 2 to level 3
P_2 = c(.1401, .6635, .1964, .1530, .6659, .1811,
        .1667, .6667, .1667, .1811, .6659, .1530,
        .1964, .1530, .6635, .1401)

# Calculate a_3 so we can get option values at j = 3
df = exp(-S4_Params["R_ts",1]*dt)

#state prices at level 2:
Q_2 = c(0.0252, 0.2199, 0.4998, 0.2199, 0.0252)

#Want Q_3: state prices at level 3 (will have 7 nodes)
Q_3 = numeric(length = nbNodes)

#according to "State prices continued" slide
Q_3[1] = Q_2[1]*P_2[1]*df
Q_3[nbNodes] = Q_3[1]
Q_3[2] = Q_2[1]*P_2[2]*df + Q_2[2]*P_2[4]*df
Q_3[nbNodes-1] = Q_3[2]
Q_3[3] = Q_2[1]*P_2[3]*df + Q_2[2]*P_2[5]*df + Q_2[3]*P_2[7]*df
Q_3[nbNodes-2] = Q_3[3]
Q_3[4] = Q_2[2]*P_2[6]*df + Q_2[3]*P_2[8]*df + Q_2[4]*P_2[10]*df

#E[x_(i+1),j] = \bar{x}_{i,j} - alpha*\bar{x}_{i,j}*dt
Exij = level_x3-(alpha*level_x3*dt)

# Transition probabilities from i = 3 to i = 4:
P_3 <- numeric(length = nbNodes*3)
offset = 1
for (j in 1:nbNodes){
  trans = prob(level_x3[j], dt, alpha, dx, sigma)
  P_3[offset:(offset+2)] = c(trans$pu, trans$pm, trans$pd)
  offset = offset + 3
}

#Assemble 3rd level of simplified tree:
L3 = as.data.frame(matrix(0, nrow = 6, ncol = 7))
colnames(L3) = c("-3", "-2", "-1", "0", "1", "2", "3")
rownames(L3) = c("xij", "E[xij]", "Qij", "pu", "pm", "pd")
offset = 1
for (k in 1:ncol(L3)){
  L3[1,k] = level_x3[k]
  L3[2,k] = Exij[k]
  L3[3,k] = Q_3[k]
  L3[4:6,k] = P_3[offset:(offset+2)]; offset = offset + 3
}

print(L3)

