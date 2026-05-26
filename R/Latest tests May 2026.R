# ====================================================================================
# ===== Demonstration of Sim - Sam ===================================================
# ====================================================================================


library(slMSE)
packageVersion('slMSE')

#  miceadds::source.all("C:/GitHub/slMSE/R")

library(MSEtool)
packageVersion('MSEtool')



# === Single substock ==========================================================

om = slOM()
hist = Simulate(om)
slplot(hist)


# === Multi stock ==============================================================

# --- Historical Simulation --------------------------------------------

nYear = 15  # No. historical years
pYear = 15  # No. projection years
Seasons = 4 # time steps, subyears/seasons per year (two-monthly)
nAreas = 3
nAges = 8
nSim = 4
CurrentYear = 2026

source("C:/GitHub/Short-Lived-IO/R/Custom_stocks.R") # Three short-lived morphs
source("C:/GitHub/Short-Lived-IO/R/Custom_fleets.R") # Two fleets (identical fishing dynamics among morphs)

om = slOM(stock = stock, fleet = fleet, ComplexName = "Pacific JFS", 
          nSim = nSim, nYear = nYear, pYear = pYear, Seasons = Seasons, 
          CurrentYear = CurrentYear, Interval = 6, Seed = 1)

hist = Simulate(om)
slplot(hist)


# --- Sim Testing Stock Assessment Methods ------------------------------

simdata = slSimData(hist)





# ====================================================================================
# ====== END =========================================================================
# ====================================================================================






