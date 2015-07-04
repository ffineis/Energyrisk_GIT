# HW 1 file: CFRM 520
# Frank Fineis
library(ENERGYRISK)

#### Problem 1: Recursive Filter, estimating Jump diffusion parameters ###

JPit <- function(R = returns, limit = 3, tolVar = "phi", tol = 1){
  tolVarnames = c("phi", "gamma", "kappa_bar", "nJumps")
  if(!(tolVar %in% tolVarnames)){stop("tolVar must be phi, gamma, or kappa_bar")}
  if(tolVar == "nJumps" & !(class(tol)=="integer" | class(tol)=="numeric")){stop("tolVar = nJumps must be an integer")}
  
  #analyze data for first iteration:
  cutoff = sd(returns)*limit
  jumps = returns[which(abs(returns)>cutoff)]
  nojumps = returns[which(abs(returns)<= cutoff)]
  nJ = length(jumps)
  jumpStore = list("jumps" = jumps, "nojumps" = nojumps, "nJ" = nJ)
  nObs = length(returns)
  phi = jumpStore$nJ/(nObs/365)
  sigma = sd(jumpStore$nojumps)*(sqrt(365))
  kappa_bar = mean(jumpStore$jumps)
  gamma = sd(jumpStore$jumps)
  JDiters = data.frame(c(jumpStore$nJ, phi, sigma, kappa_bar, gamma))
  rownames(JDiters) <- c("nJumps", "phi", "sigma", "kappa_bar", "gamma")
  colnames(JDiters) <- "Iteration 1"
  conv = F
  ctr = 1
  
  while(!conv){
    ctr = ctr+1
    cutoff = sd(jumpStore$nojumps)*limit 
    jumpStore$jumps <- returns[which(abs(returns)>cutoff)]
    jumpStore$nojumps <- returns[which(abs(returns)<=cutoff)]
    jumpStore$nJ <- length(jumpStore$jumps)
    # calculate new jump diffusion parameters
    phi = jumpStore$nJ/(nObs/365)
    sigma = sd(jumpStore$nojumps)*(sqrt(365))
    kappa_bar = mean(jumpStore$jumps)
    gamma = sd(jumpStore$jumps)
    # store these new parameters
    temp = data.frame(c(jumpStore$nJ, phi, sigma, kappa_bar, gamma))
    rownames(temp) <- c("nJumps", "phi", "sigma", "kappa_bar", "gamma")
    colnames(temp) <- paste("Iteration ", ctr, sep = "")
    JDiters <- cbind(JDiters, temp)
    # identify which variable is the tolerance parameter,
    # check to see if change in tolerance parameter is below convergence threshold (tol) 
    tvOld <- JDiters[which(rownames(JDiters)==tolVar),ctr-1]
    tvNew <- JDiters[which(rownames(JDiters)==tolVar),ctr]
    if(abs(tvOld-tvNew)<tol){conv=T}
  }
  return(JDiters)
}

data <- read.table('./data/csv/ICENBP.txt', sep = "\t", header = T)
returns <- diff(log(data$St))[-1]
results1 <- JPit(R = returns, limit = 3, tolVar = "phi", tol = 1)
print(results1)


################## Problem 2: value Asian option ###################
set.seed(1)
#Polar Rejection Random Number Generator:
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
#Spot price iteration function
lnS <- function(logPrice, param, sigrdt, dt, setRand){  
  lnS = logPrice + (param$alpha*(param$mu-logPrice) 
                    - 0.5*param$sigma^2)*dt + sigrdt*setRand
  return(lnS)
}
params = list("S" = 21.05, "alpha" = 1.196, "K" = 21.05, "T" = 0.5, "mu" = 3.053,
              "sigma" = 0.529, "s" = 1, "N" = 10, "M" = 100)
discrate = exp(-0.1*params$T)
dt = params$T/params$N
time <- seq(0, params$T, dt)
lnSmat1 = matrix(0, nrow = params$M, ncol = length(time))
lnSmat1[,1] <- log(params$S)
lnSmat2 = lnSmat1
sigrdt = params$sigma*sqrt(dt)
meanStorage1 <- exp(lnSmat1); meanStorage2 <- exp(lnSmat2)
for (ii in 1:params$M){
  for (jj in 2:ncol(lnSmat1)){
    rnd = randPolarRejc()
    lnSmat1[ii,jj] = lnS(lnSmat1[ii, (jj-1)], params, sigrdt, dt, rnd)
    lnSmat2[ii,jj] = lnS(lnSmat2[ii, (jj-1)], params, sigrdt, dt, -rnd)
  }
}
# transform from log prices to regular spot prices.
S1 <- as.data.frame(exp(lnSmat1))
S2 <- as.data.frame(exp(lnSmat2))
colnames(S1) <- time; colnames(S2) <- time; 
# get row means of simulated spot prices for valuing Azn opts
AznS1 <- apply(S1[,2:11], 1, mean); AznS2 <- apply(S2[,2:11], 1, mean)
zero = rep(0, params$M)
# payoff for regular and antithetic Azn opts
payoff_S1 = apply(cbind(zero,AznS1-params$K), 1, max)
payoff_S2 = apply(cbind(zero,AznS2-params$K), 1, max)
# average the two payoffs 
payoff = apply(cbind(payoff_S1, payoff_S2), MARGIN = 1, FUN = mean)
call_value2 = discrate*mean(payoff)
se2 = sd(payoff)/sqrt(params$M)
results2 = list("call_value" = call_value2, "se"  = se2)
print(results2)

#plot evolution of cumulative average payoff.
avg_storage = matrix(0,params$M,1)
for(k in 1:params$M){avg_storage[k] = discrate*mean(payoff[1:k])}
plot(avg_storage, type = "l", lwd = 3, xlab = "Number of Simulations",
     ylab = "Avg-ed, Discounted Payoff", main = "MC Asian Call Option, Antithetics", cex.lab = .85)
abline(h = results2$call_value)


###################### Problem 3: value knock-out barrier option ######################
barrier <- c(0.9*params$S,  1.1*params$S)
time = as.matrix(seq(0,0.5,0.05))
S1_3mo <- S1[,6:11]; S2_3mo <- S2[,6:11]
# store survival of individual price paths in tempB1, tempB2
tempB1 <- numeric(params$M); tempB2 <- tempB1
for (ii in 1:nrow(S1_3mo)){
  if(any(S1_3mo[ii,] <= barrier[1] | S1_3mo[ii,] >= barrier[2])){
    tempB1[ii] <- 0
  }
  else{
    tempB1[ii] <- S1_3mo[ii,ncol(S1_3mo)]
  }
  if(any(S2_3mo[ii,] <= barrier[1] | S2_3mo[ii,] >= barrier[2])){
    tempB2[ii] <- 0
  }
  else{
    tempB2[ii] <- S2_3mo[ii,ncol(S2_3mo)]
  }
}
#calculate payoffs: max(0, survival_price-K)
payoff_B1 = apply(cbind(zero,tempB1-params$K), 1, max)
payoff_B2 = apply(cbind(zero,tempB2-params$K), 1, max)
payoff_barrier <- apply(cbind(payoff_B1, payoff_B2), 1, mean)
call_value3 <- discrate*mean(payoff_barrier)
se3 <- sd(payoff_barrier)/sqrt(params$M)
results3 <- list("call_value" = call_value3, "se" = se3)
#Because of thin barrier knock-out range, value of option is diminished. 
print(results3)


####### Problem 4: value 6 mo option on on 1yr & 6 mo forward price spreads ########

#From Thomas Fillebeen: Schwartz one factor model forward curve function
FT <- function(S, alpha, mu, sig, s){  #s = time until maturity, alpha = mrr, S = spot price, sig = sd(Spots)
  eaT = exp(-alpha * s)
  sig2a = sig * sig / (2 * alpha)
  FT = exp(log(S) * eaT + (mu - sig2a) * (1 - eaT) + sig2a * (1 - eaT * eaT) / 2)
  return(FT)
}

#price the option on the spread itself (so K is $0.00)
K <- 0
#apply FT function to spot prices at maturity date to get forward prices in future
ft05 <- sapply(X = S1[,ncol(S1)], FUN = FT, params$alpha, params$mu, params$sigma, .5)
ft1 <- sapply(X = S1[,ncol(S1)], FUN = FT, params$alpha, params$mu, params$sigma, 1)
spread = ft1-ft05
payoffSpread <- apply(cbind(zero,spread-K), MARGIN = 1, max)
call_value4 <- discrate*mean(payoffSpread)
se4 <- sd(payoffSpread)/sqrt(params$M)
results4 <- list("call value" = call_value4, "se" = se4)
print(results4)
  
