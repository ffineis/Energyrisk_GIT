# Code for CFRM 520 Demo 5:
library(ENERGYRISK)
library(lubridate)
library(Rglpk)
data(data_S5)

forwards = S5_ForwardCurve
valDate = as.Date("03/15/07","%m/%d/%y") #valuation date
S5_StorageContract$withRate = 16393.4426229508

# intrinsic valuation of storage facility:
# determined by cashflow derived from daily trading decision, and the impact trading decision on 
# future expected value of the gas in storage. THINGS TO NOTE: you can't long and short a month, so same-month
# injections and withdrawals are prohibited
striplength = dim(forwards)[1]
injectionlength = seq(2,striplength,1) #cut off April withdraw/inject, can't happen
fSpread = matrix(0, striplength-1, striplength-1)
colnames(fSpread) = paste(forwards$Date[2:striplength]) #columns names are 1 month ahead of rownames (withdraw in future)
rownames(fSpread) = paste(forwards$Date[1:striplength-1])
for (j in injectionlength){
  fSpread[(j-1), (j-1):(striplength-1)] = ((S5_ForwardCurve[j:striplength, 2]*(1-S5_StorageContract$withFuelCost)- #make money: forward_price*(unit_price-withdrawal fuel cost)*discount_factor
    S5_StorageContract$withCost)*exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)) -
    ((S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)* #fill each row beginning at diagonal entry
    exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
}

###### CONSTRAINTS ########
#monthly injection/withdrawal constrains
injWith = matrix(0,2,striplength)
ForwardCurve_AddMonth = S5_ForwardCurve[,1] #dates (will add 1 month)
month(ForwardCurve_AddMonth) = month(S5_ForwardCurve[,1]) +1 #dates + 1month
daysInMonth = ForwardCurve_AddMonth-S5_ForwardCurve$Date #difftime object
#fill max injection volume per month and max widthrawal volume per month matrix
injWith =
  round(rbind(pmin(as.numeric(daysInMonth*S5_StorageContract$injRate),S5_StorageContract$capacity),
              pmin(as.numeric(daysInMonth*S5_StorageContract$withRate),S5_StorageContract$capacity)))
colnames(injWith) = paste(S5_ForwardCurve$Date)
rownames(injWith) = c("Max_injection", "Max_withdrawal")

#solution array storage. Guess is an "upper triangular" matrix.
position = array(0, c(striplength-1,striplength-1))
guess = position[upper.tri(position,diag=TRUE)] #need to fill 66 positions
#monthly storage constraint vector:
matrixIndex = seq(1,striplength-1,1)
offset = 0
matrixNb = matrix(0, striplength-1,striplength-1)
for(v in 1:11){matrixNb[v,v:11] = c(1:(12-v)) + offset; offset = offset + 12- length(1:v)}
#constraint matrix
strMatrix = matrix(0,striplength-1,length(guess))
for(k in 1:11){strMatrix[k,] = seq(1,length(guess),1)}
withrawalRemoval=0
storedAmount = 0
for (i in matrixIndex){
  tempIndex = c(storedAmount, matrixNb[i,])
  tempIndex = tempIndex[-which(tempIndex==0)]
  if (i>1){tempIndex = tempIndex[-which(storedAmount %in% withrawalRemoval)]}
  # Remove all injections so that can set everything else to zero
  setOV = strMatrix[i,-match(tempIndex,strMatrix[i,])]
  strMatrix[i,match(tempIndex,strMatrix[i,])] =1
  strMatrix[i,setOV] =0
  withrawalRemoval = matrixNb[,i]
  #Remove zeros so indexing stays correct
  withrawalRemoval = withrawalRemoval[-which(withrawalRemoval==0)]
  storedAmount = c(storedAmount,matrixNb[i,])
  # Remove zeros so indexing stays correct
  storedAmount = storedAmount[-which(storedAmount==0)]
  # Update accumulation
  storedAmount = tempIndex
}

# injection constraint
strInjectC = matrix(0,striplength-1, length(guess))
offset = 0
for (i in matrixIndex){strInjectC[i,(c(1:(12-i)) + offset)] = 1;
                       offset = offset + 12- length(1:i)}
strWithC = matrix(0,striplength-1,length(guess))
offset = 0
for (i in matrixIndex){strWithC[i,matrixNb[,i]] = 1}
A = rbind(strMatrix,strInjectC,strWithC) #storage constraints are stored on top of injection constraints are stored on top of withdrawal constraints
x = guess
b = matrix(0,3*dim(strMatrix)[1],1)
b[1:11]=S5_StorageContract$capacity;b[12:22]=t(injWith["Max_injection",1:11])
b[23:33] = t(injWith["Max_withdrawal", 1:11])
negConstraint = rep(0,33)
#Vectorize spreads in order specified in matrixNb. There has to be a better way of doing this than method below.
VfSpread = matrix(0,66,1);offset=0;a=1;rowUp=11;
for (j in 2:12){
  VfSpread[a:(rowUp-(j-2))] =((S5_ForwardCurve[j:striplength,2]*(1-S5_StorageContract$withFuelCost)-S5_StorageContract$withCost)*
                                exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j:striplength,1]-valDate)/365)-
                                (S5_ForwardCurve[j-1,2]*(1+S5_StorageContract$injFuelCost)+S5_StorageContract$injCost)*
                                exp(-S5_StorageContract$r*as.numeric(S5_ForwardCurve[j-1,1]-valDate)/365))
 a = a+11-offset; rowUp=rowUp+11-offset;offset = offset +1} #VfSpread has same values as fSpread

####### SOLVING
monthRevenue = pmax(VfSpread,0)*S5_StorageContract$convRate #c vector

