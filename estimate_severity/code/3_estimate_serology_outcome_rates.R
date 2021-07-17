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
set.seed(2961)

##############
# Uncomment data to be used
##############

#savingName <- "../data/processed_data/3_serology_fits.RDS"
#countryData <- read.csv("../data/collected_data/locations_serology_data.csv",
#                        stringsAsFactors=FALSE) %>%
#  as_tibble(.) %>%
#  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
#                PrevalenceH=PrevalenceH/100)
#

##savingName <- "../data/processed_data/3_serology_fits_corrected_Phi.RDS"
#countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
#                        stringsAsFactors=FALSE) %>%
#  as_tibble(.) %>%
#  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
#                PrevalenceH=PrevalenceH/100)
#englandRepDeathsInd <- which(with(countryData, Location=="England" &
#                         !(is.na(Hospitalized) & is.na(ICU))))
#countryData$Deaths[englandRepDeathsInd] <- NA

#savingName <- "../data/processed_data/3_serology_fits_u40.RDS"
#countryData <- read.csv("../data/collected_data/locations_serology_data.csv",
#                        stringsAsFactors=FALSE) %>%
#  as_tibble(.) %>%
#  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
#                PrevalenceH=PrevalenceH/100) %>%
#  dplyr::filter(., meanAge<=40)

savingName <- "../data/processed_data/3_serology_fits_corrected_noTest_corrected.RDS"
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence=Prevalence/100, PrevalenceL=PrevalenceL/100,
                PrevalenceH=PrevalenceH/100) %>%
  dplyr::filter(., Type != "Testing") 
  #dplyr::filter(., Type != "Seroprevalence_convenience") 


##############
# Fit gamma distribution to seroprevalences
##############

outcome_reg <- rstan::stan_model("./3_estimate_serology_outcome_rates.stan")
outcome <- c("Hospitalized", "ICU", "Deaths")

outcomeData <- list()
meanAge <- list()
sdAge <- list()
priorReg <- list()
outcomeDataList <- list()
model <- list()
locationKey <- list()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  # Select countries with outcome data
  selectInd <-  !is.na(countryData[[oStr]])
  outcomeData[[oStr]] <- countryData[which(selectInd),] %>%
    dplyr::mutate(., locationNum=as.integer(factor(Location)))
  meanAge[[oStr]] <- mean(outcomeData[[oStr]]$meanAge)
  sdAge[[oStr]] <- sd(outcomeData[[oStr]]$meanAge)
  outcomeData[[oStr]]$stdAge <- (outcomeData[[oStr]]$meanAge-meanAge[[oStr]])/sdAge[[oStr]]

  lowerBoundPrev <- outcomeData[[oStr]][[oStr]]/outcomeData[[oStr]]$Population
  # Put data input in list
  outcomeDataList[[oStr]] <- list(N=nrow(outcomeData[[oStr]]),
                    K=length(unique(outcomeData[[oStr]]$Location)),
                    location=outcomeData[[oStr]]$locationNum,
                    ageVec=outcomeData[[oStr]]$stdAge,
                    population=outcomeData[[oStr]]$Population,
                    seroprevShape=outcomeData[[oStr]]$gammaShape,
                    seroprevRate=outcomeData[[oStr]]$gammaRate,
                    outcomes=outcomeData[[oStr]][[oStr]])

  # Fit model
  print(paste("Fitting:", oStr))
  model[[oStr]] <- rstan::sampling(outcome_reg, data=outcomeDataList[[oStr]],
                             chains=4, iter=5000, refresh=0,
                             verbose=TRUE, cores=4)

  locationKey[[oStr]] <- unique(dplyr::select(outcomeData[[oStr]],
                                              Location, locationNum))
}

modelList <- list(model=model, locationKey=locationKey,
                  meanAge=meanAge, sdAge=sdAge)

saveRDS(modelList, savingName)


#posteriorTemp <- tidybayes::gather_draws(model$Hospitalized, ageSlope,
#                                         ageSlopeSigma, intercept, interceptSigma) 
#
#serologyTrace <- dplyr::mutate(posteriorTemp, .chain=factor(.chain))  %>%
#  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
#  geom_line() +
#  facet_grid(.~.variable, scales="free_y") +
#  theme_bw()
#

