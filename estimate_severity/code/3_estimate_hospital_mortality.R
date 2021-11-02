###################
###################
#
# This script fits the bayesian model described
# in the methods, and defined in 3_estimate_hospital_mortality.stan
# to the collected data of hospital and ICU mortality
#
###################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
source("./functions_auxiliary.R")
source("./stan_utility.R")
set.seed(2691)

mortalityData <- read.csv("../data/collected_data/hospitalized_patient_studies.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.)

outcome_reg <- rstan::stan_model("./3_estimate_hospital_mortality.stan")


####################
# define function to return initial values list
####################
initial_values <- function(nChains, paramListName, lowerUni, upperUni) {
  initList <- list()
  for (ch in c(1:nChains)) {
    initList[[ch]] <- list()
    for (p in c(1:length(paramListName))) {
      initList[[ch]][[paramListName[p]]] <- runif(1, min=lowerUni[p],
                                                  max=upperUni[p])
    }
  }
  return(initList)
}
nChains <- 4
nCores <- 2
paramListName <- c("ageSlope", "intercept", "interceptSigma", "ageSlopeSigma")
lowerUni <- c(1, -4, 0, 0)
upperUni <- c(2, -2, 0.5, 0.5)


#############################
# Fit hospital mortality by age
#############################

# Select studies of hospitalization mortality
hospMort <- dplyr::filter(mortalityData, Type=="Hospitalized") %>%
  dplyr::mutate(., locationNum=as.integer(factor(Location)))

meanHospAge <- 0 # mean(hospMort$meanAge)
sdHospAge <- sd(hospMort$meanAge)
hospMort$stdAge <- (hospMort$meanAge-meanHospAge)/sdHospAge

# Sample initial values
initList <- initial_values(nChains, paramListName, lowerUni, upperUni)

# Put data input in list
hospMortList <- list(N=nrow(hospMort),
                  K=length(unique(hospMort$Location)),
                  location=hospMort$locationNum,
                  ageVec=hospMort$stdAge,
                  cases=hospMort$Patients,
                  outcomes=hospMort$Deaths)

# Compile and run model
hospMortFit <- rstan::sampling(outcome_reg, data=hospMortList,
                               chains=4, iter=5000, refresh=0,
                               cores=nCores, init=initList,
                               control=list(adapt_delta=0.999,
                                            max_treedepth=15))

#############################
# Fit ICU mortality by age
#############################

# Select studies of hospitalization mortality
icuMort <- dplyr::filter(mortalityData, Type=="ICU") %>%
  dplyr::mutate(., locationNum=as.integer(factor(Study)))

meanICUAge <- 0 # mean(icuMort$meanAge)
sdICUAge <- sd(icuMort$meanAge)
icuMort$stdAge <- (icuMort$meanAge-meanICUAge)/sdICUAge

# Sample initial values
initList <- initial_values(nChains, paramListName, lowerUni, upperUni)

# Put data input in list
icuMortList <- list(N=nrow(icuMort),
                  K=length(unique(icuMort$Location)),
                  location=icuMort$locationNum,
                  ageVec=icuMort$stdAge,
                  cases=icuMort$Patients,
                  outcomes=icuMort$Deaths)

# Compile and run model
icuMortFit <- rstan::sampling(outcome_reg, data=icuMortList,
                              chains=4, iter=5000, refresh=0,
                              cores=4, init=initList,
                              control=list(adapt_delta=0.999,
                                           max_treedepth=15))

#pairs(icuMortFit, pars=c("ageSlope", "ageSlopeSigma", "intercept", "interceptSigma"))

##################################
# Put everything in a list to export
##################################
hospLocationDf <- unique(dplyr::select(hospMort, locationNum, Location))
icuLocationDf <- unique(dplyr::select(icuMort, locationNum, Location))

models <- list(Hospitalized=hospMortFit,
               ICU=icuMortFit)
locationKey <- list(Hospitalized=hospLocationDf,
                    ICU=icuLocationDf)
meanAge <- list(Hospitalized=meanHospAge,
                ICU=meanICUAge)
sdAge <- list(Hospitalized=sdHospAge,
                ICU=sdICUAge)
modelList <- list(model=models, locationKey=locationKey,
                  meanAge=meanAge, sdAge=sdAge)

saveRDS(modelList, "../data/processed_data/3_hospital_mortality_fit.RDS")

