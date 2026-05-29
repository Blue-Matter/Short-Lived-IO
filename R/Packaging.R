# Packaging slMSE

library(slMSE)


nYear = 25  # No. historical years
pYear = 15  # No. projection years
Seasons = 4 # time steps, subyears/seasons per year (two-monthly)
nAreas = 3
nAges = 8
nSim = 24
CurrentYear = 2026


source("C:/GitHub/Short-Lived-IO/R/Custom_stocks.R") # Three short-lived morphs
source("C:/GitHub/Short-Lived-IO/R/Custom_fleets.R") # Two fleets (identical fishing dynamics among morphs)

save(stock, file = "C:/GitHub/slMSE/data/stock.RData")
save(fleet, file = "C:/GitHub/slMSE/data/fleet.RData")
