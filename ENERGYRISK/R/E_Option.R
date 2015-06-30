#' Standard European Option on a Forward
#' 
#' @param o set to 1 for call and 2 for put
#' @param g set to 1 
#' @param K strike
#' @param T option maturity
#' @param s future maturity
#' @param F future price
#' @param sigv risk factor volatility (volatility function)
#' @param matv maturity volatility
#' @param nd number of maturities
#' @param nf number of risk factors
#' @param r interest rate
#' @author Thomas Fillebeen
#' @export
E_Option <- function(o,g,K,T,s,F,sigv,matv,nd,nf,r){
  if(K <= 1e-06){K = 1e-06}
  if(T <= 1e-06){T = 1e-06}
  if(s < T){s = T + 1e-06}
  if(F <= 1e-06){Fs = 1e-06}
  PT = exp(-r * T)
  
  var = cov_pl(T, s, s, sigv, matv, nd, nf)
  svar = sqrt(var)
  
  h1 = (log(F / K) + 0.5 * var) / svar
  h2 = h1 - svar
  Nh1 = Cumn(h1)
  Nh2 = Cumn(h2)
  
  if(o == 1){
    if(g == 1){
    val = PT * (F * Nh1 - K * Nh2)
    }else if(g == 2){
    val = PT * Nh1}
  }else if(o == 2){
    if(g == 1){
    val = PT * (K * (1 - Nh2) - F * (1 - Nh1))
    }else if(g == 2){
    val = -PT * (1 - Nh1)
    }
  }
  
  return(val)
}
#' Covariance function of ln forwards (peicewise linear) 
#' approximation of the definite integral of a function,
#' usually stated as a weighted sum of function values at 
#' specified points within the domain of integration
#' 
#' @param T option maturity
#' @param s1 future maturity oil
#' @param s2 future maturity gas
#' @param sigv risk factor volatility (volatility function)
#' @param matv maturity volatility
#' @param nd number of maturities
#' @param nf number of risk factors
#' @author Thomas Fillebeen
cov_pl <- function(T,s1, s2,sigv,matv,nd,nf){
  
  # 10 point Gaussian Quadrature abscissa and weights
  np = 5
  # Abscissa
  x = c(0.1488743389,0.4333953941,0.679409568,0.8650633666,0.9739065285)
  # Weights
  w = c(0.2955242247,0.2692667193,0.2190863625,0.1494513491,0.0666713443)
  sum = 0
  
  # For nf factors
  for(i in 1:nf){ 
  # Integrate using Gaussian Quadrature
  a = 0; b = T
  xm = 0.5 * (b + a)
  xr = 0.5 * (b - a)
  sumi = 0
    for(j in 1:np){
    dx = xr * x[j]
    sumi = sumi + w[j] *
      ((MF_vpl(i, s1 - (xm + dx), sigv, matv, nd)) *
         (MF_vpl(i, s2 - (xm + dx), sigv, matv, nd)) +
         (MF_vpl(i, s1 - (xm - dx), sigv, matv, nd)) *
         (MF_vpl(i, s2 - (xm - dx), sigv, matv, nd)))
    }
  
  sumi = sumi * xr
  sum = sum + sumi
  
  }
   cov_pl = sum
  return(cov_pl)
}
#' Multi-factor piecewise linear volatility function
#'  
#' @param i option maturity
#' @param T future maturity oil/gas
#' @param sigv risk factor volatility (volatility function)
#' @param matv maturity volatility
#' @param nd number of maturities
#' @author Thomas Fillebeen
MF_vpl <- function(i,T,sigv, matv,nd){
for(j in 1:nd){if (T <= matv[j]){break}}
if(j == 1) {v = sigv[1, i]
}else{
v = sigv[j - 1, i] + (sigv[j, i] - sigv[j - 1, i]) /
  (matv[j] - matv[j - 1]) * (T - matv[j - 1])}

MF_vpl = v
return(MF_vpl)
}

#' Cumulative normal distribution function
#' Used in the BS model N(d1)/N(d2) is the cumulative 
#' distribution function of the standard normal distribution
#' 
#' @param dx is the stock state variable
#' @author Thomas Fillebeen
Cumn <- function (dx){
  dg = (2.71828182845905 ^ -(dx * dx / 2)) / 2.50662827463
  # five terms
  dT = 1 / (1 + 0.2316419 * abs(dx))
  dx1 = dg * dT * (0.31938153 - 0.356563782 * dT + 1.781477937 * dT * dT -
                     1.821255978 * dT * dT * dT + 1.330274429 * dT * dT * dT * dT)
  
   if(dx < 0){Cumn = dx1}else{Cumn = 1 - dx1}
  return(Cumn)
}