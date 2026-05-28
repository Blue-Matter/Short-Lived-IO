# ====================================================================================
# ===== Demonstration of Sim - Sam ===================================================
# ====================================================================================

# remotes::install_github("DTUAqua/spict/spict")
# remotes::install_github("Blue-Matter/slMSE")


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
slplot(simdata)


# --- SPiCT ----------

library(spict)
sim = 1

do_spict(1,simdata)

#  r.pr = c(0.5,0.2,1); bk.pr = c(0.5,0.3,1);shape.pr = c(2, 0.001, 1);oe = c(0.2, 0.5, 1);pe = c(0.2,0.5,1);fdevs = c(4, 0.5, 1);ce = c(0.05, 0.001, 1);q.pr = NULL;timing = 0.625;dteuler = 0.25

do_spict=function(sim, simdata,
                  r.pr = c(0.5,0.2,1),
                  bk.pr = c(0.5,0.3,1),
                  shape.pr = c(2, 0.001, 1),
                  oe = c(0.2, 0.5, 1),
                  pe = c(0.2,0.5,1),
                  fdevs = c(4, 0.5, 1),
                  ce = c(0.05, 0.001, 1),
                  q.pr = NULL,
                  timing = 0.625,
                  dteuler = 0.25){

  Sdata = SPiCT_data(sim, simdata)
  Sinput = SPiCT_config(Sdata, r.pr, bk.pr, shape.pr, oe, pe, fdevs, ce,  q.pr, timing, dteuler)
  fit = fit.spict(Sinput)
  SPiCT_output(fit)

}

cbind(simdata$B[sim, ],Sout$B, Sout$B/simdata$B[sim, ])



# --- JABBA ----------


# --- RCM - ASPM -----


# --- RCM - SCAL -----



# ====================================================================================
# ====== END =========================================================================
# ====================================================================================






