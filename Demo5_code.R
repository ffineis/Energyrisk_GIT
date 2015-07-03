# Code for CFRM 520 Demo 5:
library(ENERGYRISK)
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