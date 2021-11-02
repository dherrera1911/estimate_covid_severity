###########################################
###########################################
#
#
# This script uses the IFR estimates collected in
# script 1_literature_values_covid.R, together with
# the models fitted in script 3_estimate_hospital_mortality.R
# to perform the ratio-of-ratios method described in
# the manuscript. This generates estimates of ISR and ICR
# for each IFR study
#
#
###########################################
###########################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
library(matrixStats)
source("./functions_auxiliary.R")
source("./stan_utility.R")
set.seed(2691)

# load literature IFR estimate
literatureIFR <- read.csv("../data/collected_data/literature_rates_estimations.csv",
                           stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::filter(., Type=="IFR") %>%
  dplyr::mutate(., Type="Deaths",
  Proportion=Proportion/100, Proportion_L=Proportion_L/100,
  Proportion_H=Proportion_H/100)

# so that the betas can fit, take the value with lower bound = mean
# and set the lower bound a bit smaller
equalProps <- with(literatureIFR, which(Proportion==Proportion_L))
literatureIFR$Proportion_L[equalProps] <- literatureIFR$Proportion_L[equalProps]*0.9 

# load fitted mortality data
mortalityFit <- readRDS("../data/processed_data/3_hospital_mortality_fit.RDS")

##############
# Fit beta distribution IFR estimates
##############
fittedBetas <- fit_beta_ci(meanEstimate=literatureIFR$Proportion,
                            lower=literatureIFR$Proportion_L,
                            upper=literatureIFR$Proportion_H)

#plot(log(fittedBetas$meanEstimate), log(fittedBetas$priorMean)); abline(0,1)
literatureIFR$betaShape1 <- fittedBetas$shape1
literatureIFR$betaShape2 <- fittedBetas$shape2

brazeauInd <- which(literatureIFR$Study=="Brazeau")
literatureIFR$betaShape1[brazeauInd] <- NA
literatureIFR$betaShape1[brazeauInd] <- NA

ifrSamples <- 500
outcomes1 <- c("Hospitalized", "ICU")
outcomes2 <- c("Severe", "Critical")
ageVec <- literatureIFR$meanAge
mortalityPosterior <- list()
outputDf <- dplyr::select(literatureIFR, Age, Proportion, Proportion_L,
                          Proportion_H, Study, Type, meanAge)

for (no in c(1:length(outcomes1))) {
  print(paste("Outcome:", no))
  oStr <- outcomes1[no]
  lowValsList <- list()
  highValsList <- list()
  outcomeFit <- mortalityFit$model[[oStr]]
  meanAge <- mortalityFit$meanAge[[oStr]]
  sdAge <- mortalityFit$sdAge[[oStr]]
  for (r in c(1:nrow(literatureIFR))) {
    print(paste("Row:", r))
    rowAge <- literatureIFR$meanAge[r]
    mortalityPosterior <- proportion_samples(model=outcomeFit,
                                                   ageVec=rowAge,
                                                   meanAge=meanAge,
                                                   sdAge=sdAge)
    # sample the IFR, to do the ratio of ratios
    if (!is.na(literatureIFR$betaShape1[r])){
      sampleIFR <- rbeta(1, shape1=literatureIFR$betaShape1[r],
                            shape2=literatureIFR$betaShape2[r],
                            n=length(mortalityPosterior$samples$sample))
      ratioOfRatios <- sampleIFR/mortalityPosterior$samples$proportion
      meanProp <- mean(ratioOfRatios)
      ciProp <- quantile(ratioOfRatios, probs=c(0.025, 0.975))
    } else {
      sampleIFR <- rep(literatureIFR$Proportion[r],
                       length(mortalityPosterior$samples$sample))
      ratioOfRatios <- sampleIFR/mortalityPosterior$samples$proportion
      meanProp <- mean(ratioOfRatios)
      ciProp <- rep(NA, 2)
    }
    tempDf <- data.frame(Age=literatureIFR$Age[r],
                        Proportion=meanProp,
                        Proportion_L=ciProp[1],
                        Proportion_H=ciProp[2],
                        Study=literatureIFR$Study[r],
                        Type=outcomes2[no],
                        meanAge=literatureIFR$meanAge[r])
    outputDf <- rbind(outputDf, tempDf)
  }
}

outputDf <- dplyr::mutate(outputDf, Proportion=Proportion*100,
                          Proportion_L=Proportion_L*100,
                          Proportion_H=Proportion_H*100)

fileName <- "../data/processed_data/5_literature_outcome_estimates.csv"
write.csv(outputDf, fileName, row.names=FALSE)

