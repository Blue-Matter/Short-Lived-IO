library(slMSE)
packageVersion('MSEtool')
packageVersion('slMSE')

stock = slStock()
fleet = slFleet()
om = slOM(stock = stock, fleet = fleet)

hist = Simulate(om)

load("C:/GitHub/mahiMSE/data/smallOM.rda")
om = smallOM
