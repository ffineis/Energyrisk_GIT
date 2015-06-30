#Code from CFRM 520 Demo 4
# Trinomial tree model for fitting energy spot prices to forward curves

library(ENERGYRISK)

#load tree parameters: forward curve data, discount rate info
load("./data/S4_Params.rda")

#parameters
alpha = 0.338; sigma = 0.305404; dt = 1/12; N = 12
dx = sigma*sqrt(3*dt)
params = list("alpha" = alpha, "sigma" = sigma, "dt" = dt, "dx" = dx)
option_params = list("S0" = S4_Params["F_Price",1], "K" = 21, "ttm" = 1, "r" = S4_Params["R_ts",1], "vol" = params$alpha)

#transition probabilities between trinomial branches at starting node. From Thomas Fillebeen
prob <- function(j, delta_t, alpha, delta_x, vol){  
  i=0
  pu =  c(1/2*((vol^2*delta_t + alpha^2* j^2*delta_t^2)/
                 delta_x^2 + (i-i)^2 - (alpha* j*delta_t)/delta_x*
                 (1-2*(i-i)) - (i-i)))
  pd = c(1/2*((vol^2*delta_t + alpha^2* j^2*delta_t^2)/
                delta_x^2 + (i-i)^2 + (alpha* j*delta_t)/delta_x*
                (1-2*(i-i)) - (i-i)))
  pm = 1 - pu - pd
  prob = t(cbind(pu,pm,pd))
  return(prob)
}

#Estimate mean reversion parameters to fit spot price to forward curve. From Thomas Fillebeen.
a_i <- function(prob, delta_x, j.index, df, nbNodes){
  # a) Estimate initialize state price accumulation (t = 0)
  Q = 1; sum_Q = exp(j)*Q;
  # b) Preliminary steps to estimating a_i which are chosen to ensure that 
  # the tree correctly returns the observed forward price curve:
  # -1 b/c we initalized in a) already
  for (k in 1:(length(j.index)-1)){ 
    overLap = 0
    offset = 0
    # Count the number of overlapping possibilities
    nb3Legs = nbNodes[k+1] - 4
    if(nb3Legs<0){nb3Legs =1}
    level_xt = -c(-delta_x*((k):1),delta_x*(0:(k)))  
    for(z in 1:nb3Legs){
      if(length(prob)==3){str_Q = df*Q*prob[1:3]
      }else{
        overLap = c(overLap,sum(prob[c(3,5,7)+offset]*Q[c(1:3)+(z-1)]))
        # Make sure to capture all overlapping possibilities
        if(z==nb3Legs){
          str_Q = df*c(Q[1]*prob[1], sum(Q[1:2]*prob[c(2,4)]),overLap[2:length(overLap)],
                       sum(Q[length(Q):(length(Q)-1)]*prob[c(length(prob)-1,length(prob)-3)]),
                       Q[length(Q)]*prob[length(prob)])
        }
        offset = offset + 3
      }
    }
    Q = str_Q #value at t = 0 of a security that pays $1 if on node (i,j), $0.00 otherwise
    sum_Q = c(sum_Q,(sum(Q*exp(level_xt))))
    # Re-estimate the probabilities for a given j
    j  = level_xt
    prob = prob(j, delta_t, alpha, delta_x, vol)
  }
  # c) Estimate a_i
  a_i = log(S4_Params["P_ts",]* S4_Params["F_Price",]/sum_Q)
  rownames(a_i) = "a_i"
  # Since we are only interested in the first 6-months (index starts 0)
  temp = a_i[j.index+1]
  return(temp)
}

df = exp(-option_params$r * params$dt)
level_x = -c(-params$dx*((N/2):1),params$dx*(0:(N/2))) #trinomial tree prices at time N/2. Symmetric about 0.
j.index = seq(from=0, to=N/2, by=1) 
nbNodes = seq(from=1,to=length(level_x),by=2) #number of nodes at each time point

