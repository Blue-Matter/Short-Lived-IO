
Long_highF = slFleet(Name = "Long_highF", sel50 = 4, selSLP = 5, relE_loc = 0.1, relE_freq = 1, maxF = 0.2, smu = 0.7, scv = 0.15, nYear = nYear, pYear = pYear, Seasons = Seasons, nSim = nSim, nAges = nAges, plot=T)
Short_lowF = slFleet(Name = "Short_lowF", sel50 = 3, selSLP = 4, relE_min = 0.5, relE_max = 0.75, maxF = 0.1, nYear = nYear, pYear = pYear, Seasons = Seasons, nSim = nSim, nAges = nAges, plot=T)

fleet = slCombineFleets(fleetlist = list(Long_highF = Long_highF, Short_lowF = Short_lowF), stock)

cat("Returing a StockFleetList object 'fleet' containing two identical fleets \n")



