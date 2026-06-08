# Packaging slMSE

library(slMSE)

# ============================ Make and example datafile ================================

nYear   = 25        # No. historical years
pYear   = 15        # No. projection years
Seasons = 4         # time steps, subyears/seasons per year (quarterly)
nAreas  = 3         # 3 areas used to simulate ontogeny and fleet distribution
nAges   = 8         # calculations run to 2 years
nSim    = 24        # a small number of simulations for demo purposes
CurrentYear = 2026  # 'Today'

# These functions invent 3-stock (phenotype), 2-fleet components
stock = demo_stocks(nYear, pYear, Seasons, nSim, nAges, CurrentYear) # invented
fleet = demo_fleets(nYear, pYear, Seasons, nSim, nAges, stock)       # invented

om = slOM(stock = stock, fleet = fleet,                  # combine in OM
          ComplexName = "Pacific JFS",
          nSim = nSim, nYear = nYear, pYear = pYear,
          Seasons = Seasons, CurrentYear = CurrentYear,
          Interval = 6, Seed = 1)

hist = Simulate(om, doMSYRefs=T)                 # Historical reconstruction
slplot(hist)                                     # plot simulated dynamics



simdata = slSimData(hist)                        # Extract simulated data
slplot(simdata)                                  # plot simulated data


OM = RCM_data(1, simdata, "Fit to sim data in lieu of real data",
              M = 0.5,                # Quarterly natural mortality rate
              R0init = 1E7,           # Initial value for unfished recruitment
              Steepness = 0.7,        # Assumed steepness of Ricker SRR
              i_oe  = 0.1,            # 10% CV on abundance indices
              ESS = 25,               # Effective Sample Size of 25 pa for Lcomp
              nSim = 16)              # Make OM 16 simulations big


fit = RCM(OM,                         # OM is an operating model
          OM@cpars$Data,              # That has the real data appended
          s_selectivity=c("B",1,2),   # Biomass survey, fleet 1 and 2 CPUE
          max_F = 3.0, mean_fit = T,  # report on a single fit to mean params
          condition = "catch",        # remove catch exactly
          cores = 8,                  # run in parallel
          resample = T,               # resample nsim from the var-covar matrix
          drop_nonconv = T)           # drop any non converged fits

myOM = Seasonality(fit@OM, Seasons = 4) # Add seasonality in dynamics

hist = Simulate(myOM)                 # Historical reconstruction


repGIR = GIR
formals(repGIR)$debugfile = "C:/GitHub/slMSE/data/Example_Data.rda"
class(repGIR) = 'mp'

Project(hist, "repGIR")


# ===================================================================================
