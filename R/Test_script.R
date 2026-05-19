# Test script

library(openMSE)
packageVersion('MSEtool')
source("C:/GitHub/Short-Lived-IO/R/SLSim.R")

stock = make_SL_stock()
#all(apply(stock@Spatial@Movement[1,,,,1],c(1,3),sum)==1)

fleet = make_SL_fleet()
om = make_SL_om(stock=stock, fleet=fleet)
all(apply(om@Stock@Spatial@Movement[1,,,,1],c(1,3),sum)==1)

hist = Simulate(om)






la() # load all MSEtool
OM = om
parallel = FALSE
silent = FALSE
nSim = NULL
DoDynamicUnfished = FALSE
DoRefMSY = TRUE
DoRefLandings = TRUE
DoRefRemovals = FALSE
DoConditionObs = TRUE
DoGenerateData = TRUE
Reduce = TRUE


calc_seas_spat = function(InitAgeClassRecDevs, UnfishedDist, InitMovMat, R0, nSeasons){

  nSim = dim(InitAgeClassRecDevs)[1]
  nAge = dim(InitAgeClassRecDevs)[2]
  nArea = dim(InitMovMat)[2]

  # fill simulation dimension
  if(dim(InitMovMat)[1]==1){
    InitMovMat_s = array(rep(InitMovMat, each=nSim), c(nSim, nArea, nArea, nAge)) # automatically maps to nAge if dimlength=1
  } else{
    InitMovMat_s = array(InitMovMat, c(nSim, nArea, nArea, nAge)) # automatically maps to age if dimlength = 1
  }

  if(dim(UnfishedDist)[1]==1){
    UnfishedDist_s = array(rep(UnfishedDist, each=nSim),c(nSim, nAge, nArea))
  } else{
    UnfishedDist_s = UnfishedDist
  }

  nRemFirstSeasons = (ceiling(nAge/nSeasons)*nSeasons)-nAge # the number of initial seasons to remove to get nAge R0s of the correct indexing
  R0ind = nRemFirstSeasons+1:nAge # these are the nAge seasons prior to season 1 in year 1

  if(dim(R0)[1]==1){
    R0_s = array(rep(R0,each = nSim),c(nSim,dim(R0)[2]))[,R0ind]
  } else{
    R0_s = R0[,R0ind]
  }



  Ntemp = array(NA,c(nSim, nAge, nAge, nArea))

  # The initial distribution (age class, x axis)
  for(cc in 1:nAge){ # age class 1
    reclookup = nAge - cc + 1              # reverse Rec devs so first is the last being mapped (most recent)
    Ntemp[,1,cc,] = UnfishedDist_s[,1,] * InitAgeClassRecDevs[,reclookup] * R0_s[,cc] # Cohorts start at Age 1 distribution
  }

  for(aa in 2:nAge){  # loop through age classes (x axis)
    for(cc in aa:nAge){  # loop through cohorts (y axis) calculate triangle from previous age class / chort
      moved = array(Ntemp[,aa-1,cc-1,],c(nSim,nArea,nArea)) * InitMovMat_s[,,,aa]
      Ntemp[,aa,cc,] = apply(moved,c(1,3),sum) # sum over moved-to areas
    }
  }

  Ntemp[,,nAge,] # c(nSim, nAge, nArea)

}




matplot(t(hist@SBiomass[,1,]),type="l",lty=2)

dim(hist@Number$`A short-lived creature`)


dim(om@Stock@Spatial@Movement@)

