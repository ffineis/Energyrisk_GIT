# CFRM 520, Homework 3. Frank Fineis
library(ENERGYRISK); library(zoo)
S9_Params <- read.csv("./data/csv/S9_ModelParams.csv", header = T)
S9_ForwardCurve <- read.csv("./data/csv/S9_ForwardCurve.csv", header = T)
S9_ForwardCurve$Date <- as.Date(S9_ForwardCurve$Date, format = "%m/%d/%y")


N = 200
valDate = S9_ForwardCurve$Date[1]
valEnd = as.Date("12/31/09","%m/%d/%y")
# 364 days in valuation:
nbDays = as.numeric(valEnd-valDate)
storeDates = as.Date(S9_ForwardCurve[1,1]:S9_ForwardCurve[nrow(S9_ForwardCurve),1])

KCap = 90; KFloor = 70; KSwap = 80;
VOM = 1.50; Emission = 10; heatRate = 7.25; r = 0.05
K = VOM+Emission
scaleFactor = 0.1538