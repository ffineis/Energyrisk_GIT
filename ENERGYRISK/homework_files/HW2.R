### CFRM 520 HW2. Frank Fineis.
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


##### Value EUROPEAN PUT Option ####
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
  
  if (i==0){Forw = Forw[1]; V= V[1]}
  # The following print out is useful for a trinomial function
  # It will print out the results as the function is running to keep track
  cat("Time step: ", i, "\n", sep="")
  cat("Prices:\n")
  print(Forw)
  cat("European Put Values:\n")
  print(V)
  
  offset = offset +1 
}


##### Value AMERICAN CALL Option ####
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
  
  if (i==0){Forw = Forw[1]; V= V[1]}
  # The following print out is useful for a trinomial function
  # It will print out the results as the function is running to keep track
  cat("Time step: ", i, "\n", sep="")
  cat("Prices:\n")
  print(Forw)
  cat("American Call Values:\n")
  print(V)
  
  offset = offset +1 
}


##### Value AMERICAN PUT Double Barrier Knock-out Option ####
mult = -1 #call, not put. Put is valued at K-S.
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
  
  if (i==0){Forw = Forw[1]; V= V[1]}
  # Print out the results as the function is running to keep track
  cat("Time step: ", i, "\n", sep="")
  cat("Prices:\n")
  print(Forw)
  cat("American Put/Barrier Values:\n")
  print(V)
  
  offset = offset +1 
}





