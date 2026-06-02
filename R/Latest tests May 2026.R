# ====================================================================================
# ===== Demonstration of slMSE: Sim > Sam > MSE ======================================
# ====================================================================================

# Tom Carruthers
# 29 May 2026

# A demo script showing:
#  (A) 'made up' multistock, multifleet, seasonal, spatial simulations;
#  (B) assessment data simulation
#  (C) annual and seasonal surplus production assessment testing (SPiCT)
#  (D) seasonal age-structured production model and seasonal SCAL (RCM)
#  (E) operating model specification (RCM)
#  (F) MP testing


# ==== Installation ============================================================

remotes::install_github("Blue-Matter/MSEtool", ref= "prelease")
remotes::install_github("Blue-Matter/slMSE")
remotes::install_github("DTUAqua/spict/spict")


# ==== Packages ================================================================

library(slMSE)
library(MSEtool)
library(spict)


# ==== A ==== Multi stock Simulation ===========================================

# ---- Historical Simulation --------------------------------------------

nYear = 25          # No. historical years
pYear = 15          # No. projection years
Seasons = 4         # time steps, subyears/seasons per year (two-monthly)
nAreas = 3          # 3 areas used to simulate ontogeny and fleet distribution
nAges = 8           # calculations run to 2 years
nSim = 24           # a small number of simulations for demo purposes
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



# Note that this operating model can be used directly to do MSE and MP testing:

myMSE = Project(hist, MPs = "CurrentEffort")     # Projection - current effort
B = Biomass(myMSE)                               # Extract Biomass

Bplot = do.call(data.frame,                      # Obtain quantiles
                aggregate(Value~Year,data = B, quantile,
                          probs=c(0.05,0.5,0.95)))

ggplot(Bplot) +                                  # Plot biomass
  geom_ribbon(aes(x=Year,ymin = Value.5.,ymax=Value.95.),fill = "steelblue2") +
  geom_line(aes(y=Value.50.,x=Year)) +
  geom_vline(xintercept = myMSE@OM@CurrentYear)



# ==== B ==== Simulated Data ===================================================

simdata = slSimData(hist)                        # Generate simulated data
slplot(simdata)                                  # plot simulated data


# ==== C ==== SPiCT Sim-Sam ====================================================

# ---- Annually ---------------------------------------------------------

# One simulation
Sout = do_spict(sim = 1, simdata)     # Fit spict (default settings) for sim 1
plot(Sout$fit)                        # Generic SPiCT fitting report

# Sim-Sam all simulations
SS_spict = SimSam_spict(simdata, timestep = "year",   # annual sim sam
                        parallel = T,                 # use parallel processing
                        r.pr = c(0.8,0.2,1),          # r prior
                        bk.pr = c(0.5,0.3,1),         # initial depletion prior
                        shape.pr = c(2, 0.001, 1),    # shape prior (about 0.4)
                        oe = c(0.2, 0.5, 1),          # Index observation error
                        pe = c(0.5, 0.5, 1),          # process error
                        fdevs = c(4, 0.5, 1),         # F deviation penalty
                        ce = c(0.05, 0.001, 1),       # Catch obs error
                        q.pr = NULL,                  # index q prior
                        timing = 0.01,                # index timing (within ts)
                        dteuler = 0.25)               # resolution of cont func.

slplot(SS_spict)                                      # plot sim sam results


# ---- Quarterly --------------------------------------------------------

# One simulation
Sout_q = do_spict(sim = 1, simdata, timestep = "quarter", dteuler = 0.05)
plot(Sout_q$fit)                                  # generic SPiCT fitting report

# Sim-Sam all simulations
SS_spict_q = SimSam_spict(simdata, timestep = "quarter",
                          parallel = T,               # use parallel processing
                          r.pr = c(0.8,0.2,1),        # r prior
                          bk.pr = c(0.5,0.3,1),       # initial depletion prior
                          shape.pr = c(2, 0.001, 1),  # shape prior (about 0.4)
                          oe = c(0.2, 0.5, 1),        # Index observation error
                          pe = c(0.5, 0.5, 1),        # process error
                          fdevs = c(4, 0.5, 1),       # F deviation penalty
                          ce = c(0.05, 0.001, 1),     # Catch obs error
                          q.pr = NULL,                # index q prior
                          timing = 0.01,              # index timing (within ts)
                          dteuler = 0.25)             # resolution of cont func.

slplot(SS_spict_q)                                    # plot sim sam results



# ==== D ==== RCM Sim-Sam ======================================================

# ---- Age-Structured Production Model (ASPM) -----------------------------

# One simulation
Rout = do_RCM(1, simdata)   # fit RCM default args in ASPM model for sim 1
plot(Rout$fit)

# Sim-Sam all simulations
SS_RCM_ASPM = SimSam_RCM(simdata, Name = "JFS demo",
                    mode = "ASPM",
                    c_oe = 0.05,
                    i_oe = 0.2,
                    ESS = 50,
                    C_eq_fac = 1,
                    C_eq_nyrs = 5,
                    nsubyr = 4,
                    R0init = 1E7,
                    M = 0.5,
                    Len_age = NA,
                    Wt_age = NA,
                    Mat_age = NA,
                    Sel_age = NA,
                    Steepness = 0.8,
                    SRrel = 2,
                    pe = 5.0,
                    max_F = 3.0)

slplot(SS_RCM_ASPM)



# ---- Statistical Catch-at-Length (SCAL) ---------------------------------

Rout = do_RCM(1, simdata, mode = "SCAL") # fit RCM SCAL model for sim 1
plot(Rout$fit)

# Sim-Sam all simulations
SS_RCM_SCAL = SimSam_RCM(simdata, Name = "JFS demo",
                    mode = "SCAL",
                    c_oe = 0.05,
                    i_oe = 0.15,
                    ESS = 5,
                    C_eq_fac = 1,
                    C_eq_nyrs = 5,
                    nsubyr = 4,
                    R0init = 1E7,
                    M = 0.5,
                    Len_age = NA,
                    Wt_age = NA,
                    Mat_age = NA,
                    Sel_age = NA,
                    Steepness = 0.9,
                    SRrel = 2,
                    pe = 5.0,
                    max_F = 3.0)

slplot(SS_RCM_SCAL)



# ===== Operating Model Creation ===============================================

# Normally we would fit the appropriate sim tested RCM to real data
# Here we just take the fit from simulation 1

fit = do_RCM(1, simdata, mode = "ASPM")$fit

myOM = ConvertOM(fit@OM)
myOM@Stock@Depletion = NULL # needed
Project(myOM,"")



# ===================================================================================
# ====== END ========================================================================
# ===================================================================================









# Appendix
# === Single substock ==========================================================

# om = slOM()
# hist = Simulate(om)
# slplot(hist)


