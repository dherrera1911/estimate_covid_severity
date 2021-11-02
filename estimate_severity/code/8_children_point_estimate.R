#############################################
#############################################
# 
# This script generates the fits for the main results of the
# paper. It fits the Bayesian model described in the manuscript
# and defined in 6_estimate_serology_outcome_rates.stan to the
# data collected in script 2_countries_values.R, and corrected
# in 4_countries_values_correction.R
# 
# This script also generates all the models reported in the
# supplementary, by taking subsets of the main dataset.
# 
# 
#############################################
#############################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
library(rriskDistributions)
source("./functions_auxiliary.R")
source("./stan_utility.R")
set.seed(2691)

##############
##############
#  Uncomment data to be used
##############
##############
for (i in 2) {
  if (i==1) {
#########
# Corrected full data
#########
savingName <- "../data/processed_data/8_serology_fits_corrected.RDS"
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100)
englandRepDeathsInd <- which(with(countryData, Location=="England" &
                         !(is.na(Severe) & is.na(Critical))))
countryData$Deaths[englandRepDeathsInd] <- NA
  }

  if (i==2) {
#########
# Uncorrected data
#########
savingName <- "../data/processed_data/8_serology_fits.RDS"
countryData <- read.csv("../data/collected_data/locations_serology_data.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100)
  }


  if (i==3) {
#########
# Uncorrected, with older ages filtered out
#########
savingName <- "../data/processed_data/6_serology_fits_u40.RDS"
countryData <- read.csv("../data/collected_data/locations_serology_data.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100) %>%
  dplyr::filter(., meanAge<=40)
  }

  if (i==4) {
#########
# Uncorrected, with older ages filtered out 2
#########
savingName <- "../data/processed_data/6_serology_fits_u50.RDS"
countryData <- read.csv("../data/collected_data/locations_serology_data.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100) %>%
  dplyr::filter(., meanAge<=50)
  }


  if (i==5) {
#########
# Corrected data, without data from comprehensive testing
#########
savingName <- "../data/processed_data/6_serology_fits_corrected_noTest_corrected.RDS"
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100) %>%
  #dplyr::filter(., Type != "Testing") 
  dplyr::filter(., Type != "Testing") 
englandRepDeathsInd <- which(with(countryData, Location=="England" &
                         !(is.na(Severe) & is.na(Critical))))
countryData$Deaths[englandRepDeathsInd] <- NA
  }


#########
# Corrected data, without data from convenience seroprevalence
#########
  if (i==6) {
savingName <- "../data/processed_data/6_serology_fits_corrected_noConvenience_corrected.RDS"
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100) %>%
  #dplyr::filter(., Type != "Testing") 
  dplyr::filter(., Type != "Seroprevalence_convenience") 
englandRepDeathsInd <- which(with(countryData, Location=="England" &
                         !(is.na(Severe) & is.na(Critical))))
countryData$Deaths[englandRepDeathsInd] <- NA
  }


#########
# Corrected data without most changing locations
#########
  if (i==7) {
savingName <- "../data/processed_data/6_serology_fits_corrected_lessChanging_corrected.RDS"
fastestCountries <- c("Atlanta", "Belgium", "NYC", "Sweden", "Iceland")
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100) %>%
  dplyr::filter(., !(Location %in% fastestCountries))
englandRepDeathsInd <- which(with(countryData, Location=="England" &
                         !(is.na(Severe) & is.na(Critical))))
countryData$Deaths[englandRepDeathsInd] <- NA
  }


#########
# Corrected data without locations where the epidemic peak was
# most in the past
#########
  if (i==8) {
savingName <- "../data/processed_data/6_serology_fits_corrected_lessPast_corrected.RDS"
mostRecentCountries <- c("Iceland", "Belgium", "Geneva", "NYC", "Atlanta",
                         "Netherlands", "Sweden", "South_Korea",
                         "Spain")
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100) %>%
  dplyr::filter(., (Location %in% mostRecentCountries))
  }

##############
# Fit beta distribution to seroprevalences
##############
fittedBetas <- fit_beta_ci(meanEstimate=countryData$Prevalence,
                            lower=countryData$PrevalenceL,
                            upper=countryData$PrevalenceH)
countryData$seroprevShape1 <- fittedBetas$shape1
countryData$seroprevShape2 <- fittedBetas$shape2


childrenData <- dplyr::filter(countryData, Age=="0-9")

# compile model
outcome_reg <- rstan::stan_model("./8_children_point_estimate.stan")

####################
# define function to return initial values list
####################
initial_values <- function(nChains, paramListName, lowerUni, upperUni,
                           prevalenceVals) {
  initList <- list()
  for (ch in c(1:nChains)) {
    initList[[ch]] <- list()
    for (p in c(1:length(paramListName))) {
      initList[[ch]][[paramListName[p]]] <- runif(1, min=lowerUni[p],
                                                  max=upperUni[p])
    }
    randomVec <- runif(length(prevalenceVals), min=1, max=1.1)
    initList[[ch]][["prevalence_raw"]] <-prevalenceVals*randomVec
  }
  return(initList)
}

nChains <- 4
paramListName <- c("intercept", "interceptSigma")
lowerUni <- c(-10, 0.1)
upperUni <- c(-8, 1)

########
# define some lists to fill
########
outcome <- c("Severe", "Critical", "Deaths")
outcomeData <- list()
priorReg <- list()
outcomeDataList <- list()
model <- list()
locationKey <- list()

########
# Fit the models
########
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  # Select countries with outcome data
  selectInd <-  !is.na(childrenData[[oStr]])
  outcomeData <- childrenData[which(selectInd),] %>%
    dplyr::mutate(., locationNum=as.integer(factor(Location)))

  # Put data input in list
  outcomeDataList <- list(N=nrow(outcomeData),
                    location=outcomeData$locationNum,
                    population=outcomeData$Population,
                    seroprevShape1=outcomeData$seroprevShape1,
                    seroprevShape2=outcomeData$seroprevShape2,
                    outcomes=outcomeData[[oStr]])

  # Sample initial values
  initList <- initial_values(nChains, paramListName, lowerUni, upperUni,
    outcomeData$Prevalence)

  # Fit model
  print(paste("Fitting:", oStr))
  model[[oStr]] <- rstan::sampling(outcome_reg, data=outcomeDataList,
                             chains=nChains, iter=5000, refresh=0,
                             verbose=TRUE, cores=1, init=initList,
                             control=list(adapt_delta=0.99,
                                          max_treedepth=15))

  # save the code of location number and location name
  locationKey[[oStr]] <- unique(dplyr::select(outcomeData, Location, locationNum))
}

modelList <- list(model=model, locationKey=locationKey)

saveRDS(modelList, savingName)
}


posteriorTraces <- NULL
for (oStr in outcome) {
  posteriorTemp <- tidybayes::gather_draws(model[[oStr]],
                                           intercept, interceptSigma) %>%
    dplyr::mutate(., type=oStr)
  posteriorTraces <- rbind(posteriorTraces, posteriorTemp)
}

serologyTrace <- dplyr::mutate(posteriorTraces, .chain=factor(.chain))  %>%
  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
  geom_line() +
  facet_grid(type~.variable, scales="free_y") +
  theme_bw()

#pairs(model[["Critical"]], pars=c("ageSlope", "ageSlopeSigma", "intercept",
#                                "interceptSigma"))
#pairs(model[["Deaths"]], pars=c("locationSlope"))
#
