library(slMSE)
packageVersion('slMSE')

library(MSEtool)
packageVersion('MSEtool')

#miceadds::source.all("C:/GitHub/slMSE/R")

# Single substock
om = slOM()
hist = Simulate(om)
slplot(hist)


# Multi stock

nYear = 15  # No. historical years
pYear = 15  # No. projection years
Seasons = 4 # time steps, subyears/seasons per year (two-monthly)
nAreas = 3
nAges = 8
nSim = 4
CurrentYear = 2026

source("C:/GitHub/Short-Lived-IO/R/Custom_stocks.R") # Three shortlived stocks

fleet = slFleet(nYear = nYear, pYear = pYear, Seasons = Seasons, nSim = nSim, nAges = nAges,
                sel50 = 2, selSLP = 4)

om = slOM(stock = list(small, medium, large), fleet = fleet,
          nSim = nSim, nYear = nYear, pYear = pYear, Seasons = Seasons, CurrentYear = CurrentYear,
          Interval = 6, Seed = 1)

hist = Simulate(om)


om = slOM(stock = small, fleet = fleet,
          nSim = nSim, nYear = nYear, pYear = pYear, Seasons = Seasons, CurrentYear = CurrentYear,
          Interval = 6, Seed = 1)









#mse = Project(hist, 'ECur')
#slPlot(mse)



# Multi stock with time-varying parameters

# L50

# M

# K




load("C:/GitHub/mahiMSE/data/smallOM.rda")
om = smallOM
