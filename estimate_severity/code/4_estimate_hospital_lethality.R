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

lethalityData <- read.csv("../data/collected_data/hospitalized_patient_studies.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.)

outcome_reg <- rstan::stan_model("./4_estimate_hospital_lethality.stan")

#############################
# Fit hospital lethality by age
#############################

# Select studies of hospitalization lethality
hospLeth <- dplyr::filter(lethalityData, Type=="Hospitalized") %>%
  dplyr::mutate(., locationNum=as.integer(factor(Location)))

meanHospAge <- mean(hospLeth$meanAge)
sdHospAge <- sd(hospLeth$meanAge)
hospLeth$stdAge <- (hospLeth$meanAge-meanHospAge)/sdHospAge

# Put data input in list
hospLethList <- list(N=nrow(hospLeth),
                  K=length(unique(hospLeth$Location)),
                  location=hospLeth$locationNum,
                  ageVec=hospLeth$stdAge,
                  cases=hospLeth$Patients,
                 outcomes=hospLeth$Deaths)

# Compile and run model
hospLethFit <- rstan::sampling(outcome_reg, data=hospLethList,
                               chains=4, iter=5000, refresh=0)

#############################
# Fit ICU lethality by age
#############################

# Select studies of hospitalization lethality
icuLeth <- dplyr::filter(lethalityData, Type=="ICU") %>%
  #dplyr::filter(., Type != "Testing") %>%
  dplyr::mutate(., locationNum=as.integer(factor(Location)))

meanICUAge <- mean(icuLeth$meanAge)
sdICUAge <- sd(icuLeth$meanAge)
icuLeth$stdAge <- (icuLeth$meanAge-meanICUAge)/sdICUAge

# Put data input in list
icuLethList <- list(N=nrow(icuLeth),
                  K=length(unique(icuLeth$Location)),
                  location=icuLeth$locationNum,
                  ageVec=icuLeth$stdAge,
                  cases=icuLeth$Patients,
                  outcomes=icuLeth$Deaths)

# Compile and run model
icuLethFit <- rstan::sampling(outcome_reg, data=icuLethList,
                               chains=4, iter=5000, refresh=0)

##################################
# Put everything in a list to export
##################################
hospLocationDf <- unique(dplyr::select(hospLeth, locationNum, Location))
icuLocationDf <- unique(dplyr::select(icuLeth, locationNum, Location))

models <- list(Hospitalized=hospLethFit,
               ICU=icuLethFit)
locationKey <- list(Hospitalized=hospLocationDf,
                    ICU=icuLocationDf)
meanAge <- list(Hospitalized=meanHospAge,
                ICU=meanICUAge)
sdAge <- list(Hospitalized=sdHospAge,
                ICU=sdICUAge)
modelList <- list(model=models, locationKey=locationKey,
                  meanAge=meanAge, sdAge=sdAge)

saveRDS(modelList, "../data/processed_data/4_hospital_lethality_fit.RDS")

