# Test script

library(openMSE)
packageVersion('MSEtool')
source("C:/GitHub/Short-Lived-IO/R/SLSim.R")

stock = make_SL_stock()
fleet = make_SL_fleet()
om = make_SL_om(stock=stock, fleet=fleet)

hist = Simulate(om)


saveRDS(om, "C:/temp/om_agevecs.rds")

