#' Schwartz one factor model forward curve function
#' 
#' @param S is the spot price
#' @param alpha is the mean reversion rate, 
#' @param mu is the long term level to wich S reverts
#' @param sig is the spot price volatility
#' @param s is the original maturity of 1 year
#' @author Thomas Fillebeen
#' @export
FT <- function(S, alpha, mu, sig, s){  
  eaT = exp(-alpha * s)
  sig2a = sig * sig / (2 * alpha)
  FT = exp(log(S) * eaT + (mu - sig2a) * (1 - eaT) + sig2a * (1 - eaT * eaT) / 2)
  return(FT)
}
#' Create Mean Reverting Price Path
#' 
#' @param logPrice is the log of the spot price
#' @param param is the simulation process parameters
#' @param sigrdt ia the spot price volatility in years 
#' @param dt is the time step
#' @param setRand is the polar method random number generator
#' @author Thomas Fillebeen
#' @export
lnS <- function(logPrice, param, sigrdt, dt, setRand){  
  lnS = logPrice + (param$alpha*(param$mu-logPrice) 
                             - 0.5*param$sigma^2)*dt + sigrdt*setRand
  return(lnS)
}
#' Rand_polar_rejc
#' 
#' The polar rejection method is stated algorithmically
#' 
#' @author Thomas Fillebeen
#' @export
randPolarRejc <- function(){  
  x2_available = FALSE
  x1=0;x2=0;u1=0;u2=0;w=0;c=0
  if(x2_available){
    x2_available = FALSE
    randPolarRejc = x2
  }else{
    repeat{
    u1 = runif(1) * 2 - 1
    u2 = runif(1) * 2 - 1
    w = u1 * u1 + u2 * u2
    if(w <= 1){
      c = sqrt(-2 * log(w) / w)
      x1 = c * u1
      x2 = c * u2
      x2_available <<- TRUE
      randPolarRejc = x1
      break
    }
  }
}
  return(randPolarRejc)
}