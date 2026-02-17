# Short-lived Sim

# library(MSEtool)
# https://ices-library.figshare.com/WGCEPH
# https://ices-library.figshare.com/articles/report/WGCEPH_2025_Data_call_for_2024_landings_discards_biological_effort_and_survey_data/28495505
# https://sprfmo.int/assets/Meetings/02-SC/13th-SC-2025/Working-Papers/SC13-WP17-SQUIDSIM-Papers-presented-by-Chile.pdf

# Name = "Short-lived simulation"; Agency = "A fishery agency"; Author = "A fishery analyst"; Species = "Shortus liveus"; CommonName = "Short-lived creature"
# nYear = 10; pYear = 10; Seasons = 12; CurrentYear = 2026; nSim = 4
# rec_age = 4; nages = 24; PlusGroup = F
# Linf = 1; K = 0.2; t0 = 0; Len_CV = 0.2; a = 1; b = 3
# M = 0.2; amat50 = 6; amatSLP = 2;
# h = 0.9; sigmaR = 1.0; trunc_sigmaR = 2.0; R_AC = 0.5
# nareas = 2;  prob_stay = 0.9
# Frac_area = matrix(c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.5,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05),nrow=1)
# prop_stay = 0.9
# sel50 = 6; selSLP = 2

SLSim = function(Name = "Short-lived simulation", Agency = "A fishery agency", Author = "A fishery analyst",
                 Species = "Shortus liveus", CommonName = "Short-lived creature",
                 nYear = 10, pYear = 10, Seasons = 12, CurrentYear = 2026, nSim = 4,
                 rec_age = 4, nages = 24, PlusGroup = F,
                 spawndist = c(0, 0, 0.1, 0.5, 0.3, 0.2, 0, 0, 0, 0, 0, 0),
                 Linf = 1, K = 0.2, t0 = 0, Len_CV = 0.2, a = 1, b = 3,
                 M = 0.2, amat50 = 6, amatSLP = 2,
                 h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5,
                 nareas = 2,
                 Frac_area = matrix(c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.5,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05),nrow=1),
                 prob_stay = 0.9,
                 Effort = array(NA, c(1,1,1)), sel50 = 6, selSLP = 2
                 ){


  om = new('om')     # Operating model object
  om@Name = Name
  om@Agency = Agency
  om@Author = Author
  om@Email = "An email address"
  om@Region = "A region"
  om@Latitude = 0
  om@Longitude = 0
  om@Sponsor = "A sponsor"
  om@nSim = nSim
  om@nYear = nYear
  om@pYear = pYear
  om@CurrentYear = CurrentYear
  om@Seasons = Seasons
  om@Interval = 1
  om@Seed = 1


  # --- Stock object ------------------------------------------------------------------

  na = nages - rec_age + 1
  stock = new('stock')
  stock@Name = Name
  stock@Species = Species
  stock@CommonName = CommonName
  stock@nYear = nYear
  stock@pYear = pYear
  stock@CurrentYear = CurrentYear
  stock@Years = CurrentYear+ (-(nYear-1):pYear)
  stock@Seasons = 12
  stock@nSim = nSim

  # Ages
  stock@Ages = Ages(MaxAge = nages, MinAge = rec_age, Units = "Month", PlusGroup = F)

  # Length
  Length = new('length')
  Length@MeanAtAge = Linf * 1-exp(-K * ((rec_age:nages)-t0))
  Length@CVatAge = rep(Len_CV, na)
  stock@Length = Length

  # Weight
  Weight = new('weight')
  Weight@MeanAtAge = a * Length@MeanAtAge ^ b
  Weight@Units = "Kg"
  stock@Weight = Weight

  # Natural Mortality
  stock@NaturalMortality = NaturalMortality(MeanAtAge = rep(M,na))

  # Maturity
  avec = ((rec_age : nages)-amat50)*amatSLP
  stock@Maturity = Maturity(MeanAtAge = exp(avec)/(1+exp(avec)))

  # Stock-Recruitment
  stock@SRR = SRR(Pars = list(h=0.9, R0 = 1e6), SD = sigmaR, AC = R_AC, TruncSD = trunc_sigmaR)

  # Spatial
  # nSim, nArea, nArea, nAge, and nTS,
  Movement = array(0,c(1,nareas,nareas,na))
  dist1 = c(Frac_area[,1],1-sum(Frac_area[,1]))
  Movement[1,,,aa] = matrix(dist,ncol=nareas,nrow=nareas,byrow=T) # Initial movement is fully mixed
  for(aa in 2:na){
    fracsin = c(Frac_area[,aa-1],1-sum(Frac_area[,aa-1]))
    fracsout = c(Frac_area[,aa],1-sum(Frac_area[,aa]))
    movout = get_mov_D(fracsin, fracsout, prob = rep(prob_stay,nareas))
    Movement[1,,,aa] = movout$mov
  }
  stock@Spatial = Spatial(Movement=Movement)

  # Optional slots
  # Fecundity: assumes proportional to SBiomass if left empty
  # Spatial: leave empty if no spatial structure
  # Depletion: leave empty if you don't need to optimize q for specific depeltion
  # Fleet - Retention: don't need it if all retained

  om@Stock = stock

  # --- Fleet object ------------------------------------------------------------------

  fleet = Fleet()
  fleet@Name = "Fleet1"
  fleet@nYear = nYear
  fleet@pYear = pYear
  fleet@CurrentYear = CurrentYear
  fleet@Years = CurrentYear+ (-(nYear-1):pYear)
  fleet@Seasons = 12
  fleet@nSim = nSim
  fleet@Effort = Effort_sim(nSim, nYear, Seasons, nArea, ymin = 0.25, yfac = 0.5, ECV = 0.15, maxF = 0.1, plot=T)
  fleet@Catchability = 1
  svec = ((rec_age : nages)-sel50)*selSLP
  fleet@Selectivity = Selectivity(MeanAtAge = exp(svec)/(1+exp(svec)))

  # Optional slots
  # Distribution: only needed if you don't want openMSE to solve the distribution according to vulnerable biomass

  om@Fleet = fleet

  # Observation model --------------------------------------------------------------

  obs = Obs()
  obs@Landings@CV = 0.025
  obs@Survey@CV = 0.1
  obs@CAL@ESS = 200

  # Implementation model -----------------------------------------------------------

  imp = Imp()
  imp@TAC@Mean=1
  imp@TAC@SD = 0.001

  hist = Simulate(om)

}


# Invent Effort

Effort_sim = function(nSim, nYear, Seasons, nArea, ymin = 0.25, yfac = 0.5, ECV = 0.15, maxF = 0.1, plot=F){
  nt =  nYear*Seasons
  Effort = array(0,c(nSim,nt))
  ye = 1 + ymin + (sin(((9+1:(nYear))/3))*0.5)
  se = dnorm(1:Seasons,Seasons/2,Seasons/5)
  seind = t(array(1:Seasons,c(nt,nSim)))
  yind = t(array(rep(1:nYear,each=Seasons),c(nt,nSim)))
  Effort[] = ye[yind] * se[seind] * trlnorm(nt*nSim,1,ECV)
  Effort = Effort / apply(Effort,1,max) * maxF
  if(plot)matplot(t(Effort),type="l",lty=1)
  Effort
}


dograv2 = function(log_visc,log_grav,fracsin){
  log_grav_1 = c(0,log_grav)
  nareas = length(fracsin)
  lmov = matrix(log_grav_1,byrow=T,nrow=nareas,ncol=nareas)
  diag(lmov) = diag(lmov) + log_visc
  emov = exp(lmov)
  mov = emov/apply(emov,1,sum)
  grav_out=list()
  grav_out$predfracsout = fracsin %*% mov
  grav_out$mov = mov
  grav_out$psum = mean(diag(mov))
  grav_out
}

opt_mov_D <- function(x, fracsin, fracsout, prob, probCV = 0.35, distCV = 0.01) {

  grav_out <- dograv2(log_visc = x[1:length(prob)], log_grav = x[(length(prob)+1):length(x)], fracsin=fracsin)

  nll_prior = dnorm(x,0,5,TRUE)
  nll_dist <- dnorm(log(grav_out$predfracsout), log(fracsout), distCV, TRUE)
  if(length(prob) == 1) {
    nll_stay <- dnorm(log(grav_out$psum), log(prob), probCV, TRUE)
  } else {
    nll_stay <- dnorm(log(diag(grav_out$mov)), log(prob), probCV, TRUE)
  }
  nll <- sum(nll_dist, nll_stay, nll_prior)
  return(-nll)
}

# fracsin = c(0.1, 0.2, 0.3, 0.4); fracsout = c(0.4,0.3,0.2,0.1); prob = c(0.5, 0.8, 0.9, 0.95)
get_mov_D <- function(fracsin = c(0.1, 0.2, 0.3, 0.4), fracsout = c(0.4,0.3,0.2,0.1), prob = c(0.5, 0.8, 0.9, 0.95)) {

  nareas <- length(fracsin)
  nprob <- length(prob)

  opt <- stats::nlminb(rep(0, nprob + nareas - 1), opt_mov_D,
                       fracsin = fracsin, fracsout=fracsout, prob = prob,
                       control = list(iter.max = 5e3, eval.max = 5e3))

  hess = numDeriv::hessian(opt_mov_D,opt$par, prob=prob, fracsin=fracsin,fracsout=fracsout)
  vcv = solve(hess)

  mov <- dograv2(log_visc = opt$par[1:length(prob)],
                 log_grav = opt$par[(length(prob)+1):length(opt$par)],
                 fracsin = fracsin)$mov

  out=list()
  out$opt = opt
  out$vcv = vcv
  out$par = opt$par
  out$mov = mov
  out$predfracsout = fracsin %*% mov
  out$predprobs = diag(mov)
  out
}


