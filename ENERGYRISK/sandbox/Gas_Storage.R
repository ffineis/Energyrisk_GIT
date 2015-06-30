library(ENERGYRISK)
library(lubridate)
library(optimx)
library(Rglpk)

# If you want to convert a .txt file into .rda then save and reload data using data()
#S5_ForwardCurve$Date = as.Date(S5_ForwardCurve$Date,"%m/%d/%y") 
#save.image(file="data/data_S5.rda")
# Load storage contract specs and forward curve
data(data_S5)
# Illustarte both data sets upon loading
S5_ForwardCurve
# Forward (Mid; pence/therm)
S5_StorageContract
#Set Valuation date
valDate = as.Date("03/15/07","%m/%d/%y")
# For rounding precision use the following parameter specification
S5_StorageContract$withRate =  16393.4426229508
  
# Compute the forward spreads for injection and withrawal
striplength = dim(S5_ForwardCurve)[1]
# B/c April Injection/April Withdrawal nothing happens we skip it
injectionL = seq(2,striplength,1)
fSpread = matrix(0, striplength-1,striplength-1)
# Create a topright triangular matrix
for (j in injectionL){
  fSpread[(j-1),(j-1):(striplength-1)] = 
    ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
      -S5_StorageContract$withCost)*
        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
        (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
}
colnames(fSpread)= paste(S5_ForwardCurve[2:striplength,1])
rownames(fSpread) = paste(S5_ForwardCurve[1:(striplength-1),1])
# Illustrate the spread matrix: Forward Spread (pence/therm)
# it will be vectorized to fit the Rglpk package
# Note: injection (rows) and withdrawal (columns)
fSpread
# Monthly 2x Injection/Withdrawal Constraint (MMBtu)
injWith = matrix(0,2,striplength)
# Estimate the number of days between months
ForwardCurve_AddMonth = S5_ForwardCurve[,1]
month(ForwardCurve_AddMonth) = month(S5_ForwardCurve[,1]) + 1
daysInMonth = ForwardCurve_AddMonth-S5_ForwardCurve[,1]
# inj/with result: Injection/Withdrawal Constraint (MMBtu)
injWith = 
  round(rbind(pmin(as.numeric(daysInMonth*S5_StorageContract$injRate) ,S5_StorageContract$capacity),
        pmin(as.numeric(daysInMonth*S5_StorageContract$withRate) ,S5_StorageContract$capacity)))
colnames(injWith) = paste(S5_ForwardCurve[1:striplength,1])
colnames(injWith) = paste(S5_ForwardCurve[1:striplength,1])
rownames(injWith) = c("Injection","Withdrawal")
injWith

# Generate solution array (use array to ) and then a guess (upper right triangle)
position = array(0, c(striplength-1,striplength-1))
# Imagine that the postion (guess) vector looks like
guess = position[upper.tri(position,diag=TRUE)]
guess

# Setup Constraint to be concatenated
# a) Storage Constraint
matrixIndex = seq(1,striplength-1,1)
# Setup numbered matrix to track combinations through time
offset = 0
#for(k in 1:11){matrixNb[k,] = seq(1,11,1) +offset; offset =offset +11}
matrixNb = matrix(0, striplength-1,striplength-1) 
for(v in 1:11){matrixNb[v,v:11] = c(1:(12-v)) + offset; offset = offset + 12- length(1:v)}
# Constraint Matrix
strMatrix = matrix(0,striplength-1,length(guess))
for(k in 1:11){strMatrix[k,] = seq(1,length(guess),1)}
withrawalRemoval=0
storedAmount = 0
for (i in matrixIndex){
  tempIndex = c(storedAmount,matrixNb[i,])
  tempIndex = tempIndex[-which(tempIndex==0)]
  if (i>1){tempIndex = tempIndex[-which(storedAmount %in% withrawalRemoval)]}
  # Remove all injections so that can set everything else to zero
  setOV = strMatrix[i,-match(tempIndex,strMatrix[i,])]
  strMatrix[i,match(tempIndex,strMatrix[i,])] =1
  strMatrix[i,setOV] =0
  #strMatrix[tempIndex] = 1
  withrawalRemoval = matrixNb[,i]
  # Remove zeros so indexing stays correct
  withrawalRemoval = withrawalRemoval[-which(withrawalRemoval==0)]
  storedAmount = c(storedAmount,matrixNb[i,])
  # Remove zeros so indexing stays correct
  storedAmount = storedAmount[-which(storedAmount==0)]
  # Update accumulation
  storedAmount = tempIndex
}  
head(strMatrix)

# b) Injection Constraint
strInjectC = matrix(0,striplength-1,length(guess))
offset = 0
for (i in matrixIndex){strInjectC[i,(c(1:(12-i)) + offset)] = 1; offset = offset + 12- length(1:i)}
# c) Withrawal Constraint
strWithC = matrix(0,striplength-1,length(guess))
offset = 0
for (i in matrixIndex){strWithC[i,matrixNb[,i]] = 1}

# Now we have a) strMatrix (11x66), b) strInjectC (11x66) & c) strWithC (11x66)
# Stacking a), b) & c) leads to a 33x66 A matrix
dim(strMatrix);dim(strInjectC);dim(strWithC);
A = rbind(strMatrix,strInjectC,strWithC)
x = guess

# e) Constraints 1MMBtu, Inject & Withrawal
b = matrix(0,3*dim(strMatrix)[1],1)
b[1:11]=S5_StorageContract$capacity;b[12:22]=t(injWith["Injection",1:11]);
b[23:33]=t(injWith["Withdrawal",2:12])
# e) Negativity constraint
negConstraint = rep(0,33)

# VfSpread = matrix(0,66,1);offset=0;a=1;rowUp=11;
# # Vectorize fSpreads
# for (j in 2:12){
#   VfSpread[a:(rowUp-(j-2))] = 
#     ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
#       -S5_StorageContract$withCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
#        (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
#   a = a+11-offset; rowUp=rowUp+11-offset;offset = offset +1
# }
# 
# # Objective Coefficients
# monthRevenue = pmax(VfSpread,0)*S5_StorageContract$convRate
# # Specify gt, gte, lt, lte, eq for each constraint
# dir = c(rep("<=",33))
# sol = Rglpk_solve_LP(obj=c(monthRevenue),mat=A,dir=dir,rhs=c(b),types = rep("I"),max=TRUE)
# 
# # Position to be taken in the spreads
# sol$solution
# # Maximized revenue result
# sol$optimum
# # Check Constraints
# A%*%sol$solution
# 
# # Optimal position
# optPosition = sol$solution
# storageLevel = rbind(0,strMatrix%*%sol$solution)
# # Plot the first three simulated price paths
# par(mar=c(4.5,3.5,3,6)+0.1)
# plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
# title("Storage & Forward Strip")
# mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
# par(new=TRUE)
# plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
#      type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
# axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
# mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
# legend("topleft",legend=c("Storage","Foward Strip"),
#        col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)


# ## 1) Fuel cost of 1% to both inject and withdraw
# # Reset injection and withdrawal costs:
# S5_StorageContract$injFuelCost = 0.01
# S5_StorageContract$withFuelCost = 0.01
# VfSpread = matrix(0,66,1);offset=0;a=1;rowUp=11;
# # Vectorize fSpreads
# for (j in 2:12){
#   VfSpread[a:(rowUp-(j-2))] = 
#     ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
#       -S5_StorageContract$withCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
#        (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
#   a = a+11-offset; rowUp=rowUp+11-offset;offset = offset +1
# }
# # Objective Coefficients
# monthRevenue = pmax(VfSpread,0)*S5_StorageContract$convRate
# # Specify gt, gte, lt, lte, eq for each constraint
# dir = c(rep("<=",33))
# sol = Rglpk_solve_LP(obj=c(monthRevenue),mat=A,dir=dir,rhs=c(b),types = rep("I"),max=TRUE)
# # Position to be taken in the spreads
# sol$solution
# # Maximized revenue result
# sol$optimum
# # Check Constraints
# A%*%sol$solution

# # Optimal position
# optPosition = sol$solution
# storageLevel = rbind(0,strMatrix%*%sol$solution)
# # Plot the first three simulated price paths
# par(mar=c(4.5,3.5,3,6)+0.1)
# plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
# title("Storage & Forward Strip")
# mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
# par(new=TRUE)
# plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
#      type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
# axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
# mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
# legend("topleft",legend=c("Storage","Foward Strip"),
#        col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)


## 2) Interest rate curve of 5%
# Reset injection and withdrawal costs:
S5_StorageContract$injFuelCost = 0.01
S5_StorageContract$withFuelCost = 0.01
S5_StorageContract$r = 0.05
VfSpread = matrix(0,66,1);offset=0;a=1;rowUp=11;
# Vectorize fSpreads
for (j in 2:12){
  VfSpread[a:(rowUp-(j-2))] = 
    ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
      -S5_StorageContract$withCost)*
       exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
       (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
       exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
  a = a+11-offset; rowUp=rowUp+11-offset;offset = offset +1
}
# Objective Coefficients
monthRevenue = pmax(VfSpread,0)*S5_StorageContract$convRate
# Specify gt, gte, lt, lte, eq for each constraint
dir = c(rep("<=",33))
sol = Rglpk_solve_LP(obj=c(monthRevenue),mat=A,dir=dir,rhs=c(b),types = rep("I"),max=TRUE)
# Position to be taken in the spreads
sol$solution
# Maximized revenue result
sol$optimum
# Check Constraints
A%*%sol$solution

# Optimal position
optPosition = sol$solution
storageLevel = rbind(0,strMatrix%*%sol$solution)
# Plot the first three simulated price paths
par(mar=c(4.5,3.5,3,6)+0.1)
plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
title("Storage & Forward Strip")
mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
par(new=TRUE)
plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
     type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
legend("topleft",legend=c("Storage","Foward Strip"),
       col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)


## 3) The withdrawal rate is halved
# Reset injection and withdrawal costs:
S5_StorageContract$injFuelCost = 0.01
S5_StorageContract$withFuelCost = 0.01
S5_StorageContract$r = 0.05
S5_StorageContract$withRate = 8197
# Monthly 2x Injection/Withdrawal Constraint (MMBtu)
injWith = matrix(0,2,striplength)
# Estimate the number of days between months
ForwardCurve_AddMonth = S5_ForwardCurve[,1]
month(ForwardCurve_AddMonth) = month(S5_ForwardCurve[,1]) + 1
daysInMonth = ForwardCurve_AddMonth-S5_ForwardCurve[,1]
# inj/with result: Injection/Withdrawal Constraint (MMBtu)
injWith = 
  round(rbind(pmin(as.numeric(daysInMonth*S5_StorageContract$injRate) ,S5_StorageContract$capacity),
              pmin(as.numeric(daysInMonth*S5_StorageContract$withRate) ,S5_StorageContract$capacity)))
colnames(injWith) = paste(S5_ForwardCurve[1:striplength,1])
colnames(injWith) = paste(S5_ForwardCurve[1:striplength,1])
rownames(injWith) = c("Injection","Withdrawal")
injWith



### Re-estimate the monthly injection and withrawal constraint
# e) Constraints 1MMBtu, Inject & Withrawal
b = matrix(0,3*dim(strMatrix)[1],1)
b[1:11]=S5_StorageContract$capacity;b[12:22]=t(injWith["Injection",1:11]);
b[23:33]=t(injWith["Withdrawal",2:12])

# Initialize the spread vector
VfSpread = matrix(0,66,1);offset=0;a=1;rowUp=11;
# Vectorize fSpreads
for (j in 2:12){
  VfSpread[a:(rowUp-(j-2))] = 
    ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
      -S5_StorageContract$withCost)*
       exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
       (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
       exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
  a = a+11-offset; rowUp=rowUp+11-offset;offset = offset +1
}
# Objective Coefficients
monthRevenue = pmax(VfSpread,0)*S5_StorageContract$convRate
# Specify gt, gte, lt, lte, eq for each constraint
dir = c(rep("<=",33))
sol = Rglpk_solve_LP(obj=c(monthRevenue),mat=A,dir=dir,rhs=c(b),types = rep("I"),max=TRUE)
# Position to be taken in the spreads
sol$solution
# Maximized revenue result
sol$optimum
# Check Constraints
A%*%sol$solution

# Optimal position
optPosition = sol$solution
storageLevel = rbind(0,strMatrix%*%sol$solution)
# Plot the first three simulated price paths
par(mar=c(4.5,3.5,3,6)+0.1)
plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
title("Storage & Forward Strip")
mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
par(new=TRUE)
plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
     type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
legend("topleft",legend=c("Storage","Foward Strip"),
       col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)




#########################################################################
# position = c(rep(0,121))
# position[2]=34;position[8]=16393;position[10]=229483;position[13]=254107;
# position[42]=8180;position[43]=245927;position[53]=254107;position[64]=245910
# position[1,2] = 34;position[2,2] = 254100;position[1,8] = 16393;position[4,9] = 8180
# position[5,9] = 254107; position[6,9] = 245900;
# position[1,10] = 229483;position[4,10] = 245910
# # Create a function to call and optimize over
# maxRev = function(position,fSpread,S5_StorageContract,striplength,returnSt = FALSE){
#   S = diag(11);  S[upper.tri(S,diag=TRUE)] = position
#   injWith = 
#   round(rbind(pmin(as.numeric(daysInMonth*S5_StorageContract$injRate) ,S5_StorageContract$capacity),
#   pmin(as.numeric(daysInMonth*S5_StorageContract$withRate) ,S5_StorageContract$capacity)))
#   rownames(injWith) = c("Injection","Withdrawal")
#   
#   if(any(S<0)){expectedRevenue = -Inf; return(expectedRevenue)}
#   # Constraints
#   # 1) Withdrawal needs to be constraint by: [0,injWith["Withdrawal",]]
#   m_position = matrix(S,11,11)
#   posWith = colSums(m_position)
#   if(any(posWith>injWith["Withdrawal",2:12]) ||
#        any(posWith<0)){expectedRevenue = -Inf; return(expectedRevenue)}
#   # 2) Withdrawal needs to be constraint by: [0,injWith["Injection",]] 
#   posInj = rowSums(m_position)
#   if(any(posInj>injWith["Injection",1:11]) || 
#        any(posInj<0)){expectedRevenue = -Inf; return(expectedRevenue)}
#   
#   # Monthly Storage Level
#   storageLevel = matrix(0,1,striplength)
#   storageLevel[1] = S5_StorageContract$currStrgL
#   posWith =c(0,posWith); posInj =c(posInj,0)
#   for(j in 1:11){storageLevel[j+1] = storageLevel[j] + (posInj[j]-posWith[j])}
#   # 3) Need a recusive storage constraint: 0 =< storageLevel <= S5_StorageContract$capacity
#   if(any(storageLevel>S5_StorageContract$capacity) || 
#        any(storageLevel<0)){expectedRevenue = -Inf; return(expectedRevenue)}
#   # Monthly Revenue & objective function result
#   # convRate (multiple to convert pence/therm into pound/MMBtu)
#   monthRevenue = fSpread*S5_StorageContract$convRate*S
#   expectedRevenue = sum(monthRevenue)
#   if(returnSt){return(storageLevel)
#   }else{return(expectedRevenue)}
# }


# # Develop a function object to be called when optimizing
# fn = function(x) maxRev(x,fSpread,S5_StorageContract,striplength)
# # Use optim(.) function from library 'optimx' to optimize with the initial 
# # position guess.
# guess = position[upper.tri(position,diag=TRUE)]
# 
# result = optim(par=guess,fn=fn, method="Nelder-Mead",
#                control = list(fnscale = -1,maxit = 10000, reltol = 1e-6))
# # Control/decision variable 
# S = diag(11);  S[upper.tri(S,diag=TRUE)] = result$par
# S
# # Objective function 
# result$value
# 
# # Optimal position
# optPosition = result$par
# storageLevel = maxRev(optPosition,fSpread,S5_StorageContract,striplength, returnSt =TRUE)
# # Plot the first three simulated price paths
# par(mar=c(4.5,3.5,3,6)+0.1)
# plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
# title("Storage & Forward Strip")
# mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
# par(new=TRUE)
# plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
#      type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
# axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
# mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
# legend("topleft",legend=c("Storage","Foward Strip"),
#        col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)




## 1) Fuel cost of 1% to both inject and withdraw
# Reset injection and withdrawal costs:
# S5_StorageContract$injFuelCost = 0.01
# S5_StorageContract$withFuelCost = 0.01
# # Create a topright triangular matrix
# for (j in injectionL){
#   fSpread[(j-1),(j-1):(striplength-1)] = 
#     ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
#       -S5_StorageContract$withCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
#        (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
# }
# colnames(fSpread)= paste(S5_ForwardCurve[2:striplength,1])
# rownames(fSpread) = paste(S5_ForwardCurve[1:(striplength-1),1])
# 
# # Call function again
# fn = function(x) maxRev(x,fSpread,S5_StorageContract,striplength)
# # Use optim(.) function from library 'optimx' to optimize with the initial 
# # position guess.
# guess = position[upper.tri(position,diag=TRUE)]
# 
# result = optim(par=guess,fn=fn, method="Nelder-Mead",
#                control = list(fnscale = -1,maxit = 10000, reltol = 1e-6))
# # Control/decision variable 
# S = diag(11);  S[upper.tri(S,diag=TRUE)] = result$par
# S
# # Objective function 
# result$value
# optPosition = result$par
# storageLevel = maxRev(optPosition,fSpread,S5_StorageContract,striplength, returnSt =TRUE)  
# # Plot the first three simulated price paths
# par(mar=c(4.5,3.5,3,6)+0.1)
# plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
# title("Storage & Forward Strip")
# mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
# par(new=TRUE)
# plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
#      type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
# axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
# mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
# legend("topleft",legend=c("Storage","Foward Strip"),
#        col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)
# 
# 
# ## 2) Interest rate curve of 5%
# # Reset injection and withdrawal costs:
# S5_StorageContract$injFuelCost = 0.01
# S5_StorageContract$withFuelCost = 0.01
# S5_StorageContract$r = 0.05
# # Create a topright triangular matrix
# for (j in injectionL){
#   fSpread[(j-1),(j-1):(striplength-1)] = 
#     ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
#       -S5_StorageContract$withCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
#        (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
# }
# colnames(fSpread)= paste(S5_ForwardCurve[2:striplength,1])
# rownames(fSpread) = paste(S5_ForwardCurve[1:(striplength-1),1])
# # Use optim(.) function from library 'optimx' to optimize with the initial 
# # position guess.
# result = optim(par=test,fn=fn, method="Nelder-Mead",
#                control = list(fnscale = -1,maxit = 10000, reltol = 1e-6))
# # Control/decision variable 
# result$par
# # Objective function 
# result$value
# optPosition = result$par
# storageLevel = maxRev(optPosition,fSpread,S5_StorageContract,striplength, returnSt =TRUE)  
# # Plot the first three simulated price paths
# par(mar=c(4.5,3.5,3,6)+0.1)
# plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
# title("Storage & Forward Strip")
# mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
# par(new=TRUE)
# plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
#      type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
# axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
# mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
# legend("topleft",legend=c("Storage","Foward Strip"),
#        col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)


# ## 3) The withdrawal rate is halved
# # Reset injection and withdrawal costs:
# S5_StorageContract$injFuelCost = 0.01
# S5_StorageContract$withFuelCost = 0.01
# S5_StorageContract$r = 0.05
# S5_StorageContract$withRate = 8197
# # Create a topright triangular matrix
# for (j in injectionL){
#   fSpread[(j-1),(j-1):(striplength-1)] = 
#     ((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)
#       -S5_StorageContract$withCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
#        (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
#        exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
# }
# colnames(fSpread)= paste(S5_ForwardCurve[2:striplength,1])
# rownames(fSpread) = paste(S5_ForwardCurve[1:(striplength-1),1])
# 
# # Use optim(.) function from library 'optimx' to optimize with the initial 
# # position guess.
# result = optim(par=test,fn=fn, method="Nelder-Mead",
#                control = list(fnscale = -1,maxit = 10000, reltol = 1e-6))
# # Control/decision variable 
# result$par
# # Objective function 
# result$value
# optPosition = result$par
# storageLevel = maxRev(optPosition,fSpread,S5_StorageContract,striplength, returnSt =TRUE)  
# # Plot the first three simulated price paths
# par(mar=c(4.5,3.5,3,6)+0.1)
# plot(S5_ForwardCurve[,1],t(storageLevel),type="l",xlab="Time",ylab="")
# title("Storage & Forward Strip")
# mtext("Volume (MMBtu)",side=2,line=2,cex=0.75)
# par(new=TRUE)
# plot(S5_ForwardCurve[,1],S5_ForwardCurve[,2],
#      type="l",xlab="Time",ylab="",col=2,lwd=2, lty=2,yaxt="n")
# axis(4, ylim=c(min(S5_ForwardCurve[,2]),max(S5_ForwardCurve[,2])),line=0,cex.axis=0.75)
# mtext("Price (pence/therm)",side=4,line=1.75,cex=0.75)
# legend("topleft",legend=c("Storage","Foward Strip"),
#        col = c(1,2), lwd=c(1.5,1.5), lty=c(1,2), cex=0.8)


#expectedRevenue = maxRev(position,fSpread,S5_StorageContract,striplength)$expectedRevenue
#monthRevenue = maxRev(position,fSpread,S5_StorageContract,striplength)$monthRevenue
# col = injWith["Withdrawal",2:12]
# row = injWith["Injection",1:11]
# add = rep(S5_StorageContract$capacity,11)
# row.signs <- rep ("<", 11)
# col.signs <- rep ("<", 11)
# lp.transport(monthRevenue,"max",row.signs,row,col.signs,col)
#optim(position,fn,method = "SANN", control = list(fnscale = -1,maxit = 1000000, reltol = 1e-6))
#optim(position,fn,method = "SANN", control = list(fnscale = -1,maxit = 100, reltol = 1e-2))