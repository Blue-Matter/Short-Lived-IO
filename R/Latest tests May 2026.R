#library(slMSE)
library(MSEtool)
miceadds::source.all("C:/GitHub/slMSE/R")
#packageVersion('MSEtool')
#packageVersion('slMSE')


# Single substock


om = slOM()
hist = Simulate(om)



mse = Project(hist, 'ECur')
slPlot(mse)


# Multi stock


nYear = 15  # No. historical years
pYear = 15  # No. projection years
Seasons = 6 # time steps, subyears/seasons per year (two-monthly)
nAreas = 3

small = slStock(Name = "Small Pacific JFS", Species = "Dosidicus gigas", CommonName = "Jumbo Flying Squid",
                nYear = nYear, pYear = pYear, Seasons = Seasons, # Historical years, projection years, subyears
                CurrentYear = 2026, nSim = 4,                    # Last historical year, number of simulations
                rec_age = 1, nages = 12, PlusGroup = F,          # Season / subyear to recruitment, max seasons, senescence after nages
                spawndist = c(0, 0, 0.4, 0.5, 0.1, 0),           # Spawning distribution among sub-years
                Len_age =  c(22.0, 28.0, 31.9, 34.6, 36.4, 37.6, # Length at age (cm)
                             38.4, 38.9, 39.3, 39.5, 39.7, 39.8),               # round(40 * (1-exp(-0.4*2:13)),1)
                Len_CV = 0.2,
                Wt_age = c(1.1, 2.2, 3.3, 4.1, 4.8, 5.3,         # Weight at age (kg)
                           5.6, 5.9, 6.1, 6.2, 6.2, 6.3),                       # round(1E-4* (40 * (1-exp(-0.4*2:13)))^3,1)
                M = 0.4,                                         # Instantaneous natural mortality rate (per season)
                Mat_age = c(0.03, 0.14, 0.50, 0.86, 0.97, 1.00,  # Spawning fraction at age ('maturity at age')
                            1.00, 1.00, 1.00, 1.00, 1.00, 1.00),                # avec = ((1 : 12)-3) * 1.8; round(exp(avec)/(1+exp(avec)),2)
                SR_type = "Ricker", h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                nareas = nAreas,
                Frac_area = matrix(c(0.95, 0.80, 0.75, 0.60, 0.45, 0.30,           # Area 1
                                     0.05, 0.15, 0.20, 0.25, 0.30, 0.40,           # Area 2
                                     0.00, 0.05, 0.05, 0.15, 0.25, 0.30),          # Area 3
                                   byrow = T, nrow=nAreas),
                prob_stay = 0.8)                                                               # Tendency to remain in area

# !!! need to extend Frac_area to age classes

medium = slStock(Name = "Medium Pacific JFS", Species = "Dosidicus gigas", CommonName = "Jumbo Flying Squid",
                nYear = nYear, pYear = pYear, Seasons = Seasons, # Historical years, projection years, subyears
                CurrentYear = 2026, nSim = 4,                    # Last historical year, number of simulations
                rec_age = 1, nages = 12, PlusGroup = F,          # Season / subyear to recruitment, max seasons, senescence after nages
                spawndist = c(0, 0.1, 0.4, 0.4, 0.1, 0),           # Spawning distribution among sub-years
                Len_age =  c(37.9, 46.6, 51.9, 55.1, 57.0, 58.2, # Length at age (cm)
                             58.9, 59.3, 59.6, 59.8, 59.9),                     # round(60 * (1-exp(-0.5*2:13)),1)
                Len_CV = 0.2,
                Wt_age = c(5.5,  10.1, 14.0, 16.7, 18.5, 19.7,   # Weight at age (kg)
                           20.4, 20.9, 21.2, 21.3, 21.4, 21.5),  # round(1E-4* (60 * (1-exp(-0.5*2:13)))^3,1)
                M = 0.35,                                         # Instantaneous natural mortality rate (per season)
                Mat_age = c(0.03, 0.14, 0.50, 0.86, 0.97, 1.00,  # Spawning fraction at age ('maturity at age')
                            1.00, 1.00, 1.00, 1.00, 1.00, 1.00),                # avec = ((1 : 12)-3) * 1.8; round(exp(avec)/(1+exp(avec)),2)
                SR_type = "Ricker", h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                nareas = nAreas, Frac_area = matrix(c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25,           # Area 1
                                                 0.50, 0.50, 0.50, 0.50, 0.50, 0.50,           # Area 2
                                                 0.25, 0.25, 0.25, 0.25, 0.25, 0.25), nrow=nAreas), # Area 3
                prob_stay = 0.8)                                                               # Tendency to remain in area



large = slStock(Name = "Large Pacific JFS", Species = "Dosidicus gigas", CommonName = "Jumbo Flying Squid",
                 nYear = nYear, pYear = pYear, Seasons = Seasons, # Historical years, projection years, subyears
                 CurrentYear = 2026, nSim = 4,                    # Last historical year, number of simulations
                 rec_age = 1, nages = 12, PlusGroup = F,          # Season / subyear to recruitment, max seasons, senescence after nages
                 spawndist = c(0, 0.2, 0.4, 0.3, 0.1, 0),         # Spawning distribution among sub-years
                 Len_age =  c(38.5, 48.9, 55.9, 60.5, 63.6, 65.7, # Length at age (cm)
                              67.1, 68.1, 68.7, 69.1, 69.4, 69.6),              # round(70 * (1-exp(-0.4*2:13)),1)
                 Len_CV = 0.2,
                 Wt_age = c(5.7,  11.7, 17.4, 22.2, 25.8, 28.4,   # Weight at age (kg)
                            30.3, 31.6, 32.4, 33.1, 33.5, 33.7),  # round(1E-4* (70 * (1-exp(-0.4*2:13)))^3,1)
                 M = 0.3,                                         # Instantaneous natural mortality rate (per season)
                 Mat_age = c(0.03, 0.14, 0.50, 0.86, 0.97, 1.00,  # Spawning fraction at age ('maturity at age')
                             1.00, 1.00, 1.00, 1.00, 1.00, 1.00),                # avec = ((1 : 12)-3) * 1.8; round(exp(avec)/(1+exp(avec)),2)
                 SR_type = "Ricker", h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                 nareas = nAreas, Frac_area = matrix(c(0.95, 0.80, 0.75, 0.60, 0.45, 0.30,           # Area 1
                                                  0.05, 0.15, 0.20, 0.25, 0.30, 0.40,           # Area 2
                                                  0.00, 0.05, 0.05, 0.15, 0.25, 0.30), nrow=nAreas), # Area 3
                 prob_stay = 0.8)                                                               # Tendency to remain in area


fleet = slFleet(nYear = nyear, pYear = pYear, Seasons = Seasons)

om = slOM(stock = list(small, medium, large), fleet = fleet)
hist = Simulate(om)












#mse = Project(hist, 'ECur')
#slPlot(mse)



# Multi stock with time-varying parameters

# L50

# M

# K




load("C:/GitHub/mahiMSE/data/smallOM.rda")
om = smallOM
