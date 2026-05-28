library(MSEtool)

om1 = readRDS("C:/temp/OM1.rds")
hist = Simulate(om1)

om2 = readRDS("C:/temp/OM2.rds")
om2@Stock |> length()
om2@Complexes |> length()
om2@Obs <- list(om2@Obs[[1]])
hist = Simulate(om2)
