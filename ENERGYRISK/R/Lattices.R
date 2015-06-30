#' Probbility estimation for Simplified Trees: Mean Reversion
#' 
#' @param j is the spot price level
#' @param delta_t is the mean reversion rate, 
#' @param alpha is the mean reversion rate
#' @param delta_x price steps
#' @param vol is the spot price volatility
#' @author Thomas Fillebeen
#' @export
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
#' Estimating a_i which are chosen to ensure that 
#' the tree correctly returns the observed forward price curve
#' 
#' @param prob 
#' @param delta_x
#' @param j.index  
#' @param df 
#' @param nbNodes 
#' @author Thomas Fillebeen
#' @export
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
    Q = str_Q
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