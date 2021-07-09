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
  dplyr::mutate(., Type="Deaths")

# load fitted lethality data
lethalityFit <- readRDS("../data/processed_data/4_hospital_lethality_fit.RDS")

#outcomes <- c("Hospitalized", "ICU")
#ageVec <- literatureIFR$meanAge
#lethalityPosterior <- list()
#outputDf <- literatureIFR
#for (no in c(1:length(outcomes))) {
#  oStr <- outcomes[no]
#  outcomeFit <- lethalityFit$model[[oStr]]
#  meanAge <- lethalityFit$meanAge[[oStr]]
#  sdAge <- lethalityFit$sdAge[[oStr]]
#  lethalityPosterior[[oStr]] <- proportion_samples(model=outcomeFit,
#                                                   ageVec=ageVec,
#                                                   meanAge=meanAge,
#                                                   sdAge=sdAge)
#  ifrVec <- rep(literatureIFR$Proportion, max(lethalityPosterior[[oStr]]$sample))
#  outcomeProp <- ifrVec/lethalityPosterior[[oStr]]$samples$proportion
#  # Find mean and quantiles of samples
#  outcomePropMat <- matrix(outcomeProp, nrow=nrow(literatureIFR),
#            ncol=max(lethalityPosterior[[oStr]]$sample))
#  prop_mean <- rowMeans(outcomePropMat)
#  ciProp <- matrixStats::rowQuantiles(outcomePropMat, probs=c(0.025, 0.975))
#  # Put into data frame
#  tempDf <- data.frame(Age=literatureIFR$Age,
#                      Proportion=prop_mean,
#                      Proportion_L=ciProp[,1],
#                      Proportion_H=ciProp[,2],
#                      Study=literatureIFR$Study,
#                      Type=oStr,
#                      meanAge=ageVec)
#  outputDf <- rbind(outputDf, tempDf)
#}
#
#fileName <- "../data/processed_data/5_literature_outcome_estimates.csv"
#write.csv(outputDf, fileName, row.names=FALSE)
#

# Do the same but sampling from the IFR distributios
outcomes <- c("Hospitalized", "ICU")
ageVec <- literatureIFR$meanAge
lethalityPosterior <- list()


##############
# Fit gamma distribution IFR estimates
##############
library(rriskDistributions)
# number of cases the values of a gamma distribution over the cases
gammaPars <- function(casesQuantiles, quants=c(0.025, 0.5, 0.975)) { #
  #pars <- get.gamma.par(p=quants, q=casesQuantiles,
  #                      plot=FALSE, verbose=FALSE, show.output=FALSE)
  pars <- get.lnorm.par(p=quants, q=casesQuantiles,
                        plot=FALSE, verbose=FALSE, show.output=FALSE)
}

meanlog <- NULL
sdlog <- NULL
for (r in c(1:nrow(literatureIFR))) {
  row <- literatureIFR[r,]
  pars <- NA
  if (row$Proportion_L!=row$Proportion & row$Proportion_L!=0) {
    pars <- with(row, gammaPars(c(Proportion_L, Proportion, Proportion_H)))
  }
  if (is.na(pars) & row$Proportion_L==0) {
    pars <- with(row, gammaPars(c(Proportion_L+10^(-6), Proportion, Proportion_H),
                                quants=c(0.025, 0.5, 0.975)))
  }
  if (is.na(pars) & row$Proportion_L!=0) {
    pars <- with(row, gammaPars(c(Proportion_L, Proportion, Proportion_H*1.01),
                                quants=c(0.025, 0.5, 0.975)))
  }
  if (is.na(pars) & (row$Proportion==row$Proportion_H) & row$Proportion_L==0) {
    pars <- with(row, gammaPars(c(Proportion_L+10^(-6), Proportion, Proportion_H*1.05),
                                quants=c(0.025, 0.5, 0.975)))
  }
  meanlog[r] <- pars[1]
  sdlog[r] <- pars[2]
}

# Just check that the gammas and the quantiles match
#meanPrev <- NULL
#lPrev <- NULL
#hPrev <- NULL
#for (r in c(1:length(meanlog))) {
##  meanPrev[r] <- qgamma(p=c(0.5), rate=gammaRate[r], shape=gammaShape[r])
##  lPrev[r] <- qgamma(p=c(0.025), rate=gammaRate[r], shape=gammaShape[r])
##  hPrev[r] <- qgamma(p=c(0.975), rate=gammaRate[r], shape=gammaShape[r])
#  meanPrev[r] <- qlnorm(p=c(0.5), sdlog=sdlog[r], meanlog=meanlog[r])
#  lPrev[r] <- qlnorm(p=c(0.025), sdlog=sdlog[r], meanlog=meanlog[r])
#  hPrev[r] <- qlnorm(p=c(0.975), sdlog=sdlog[r], meanlog=meanlog[r])
#}
#plot(log(literatureIFR$Proportion), log(meanPrev))
#abline(0,1)
#plot(log(literatureIFR$Proportion_L), log(lPrev))
#abline(0,1)
#plot(log(literatureIFR$Proportion_H), log(hPrev))
#abline(0,1)
#
#plot(log(b1$Proportion), log(meanB))
#abline(0,1)
#plot(log(b1$Proportion_L), log(lowerB))
#abline(0,1)
#plot(log(b1$Proportion_H), log(upperB))
#abline(0,1)
#lele <- which(literatureIFR$Proportion_L/lPrev>1.5)


ifrSamples <- 500
outcomes <- c("Hospitalized", "ICU")
ageVec <- literatureIFR$meanAge
lethalityPosterior <- list()
outputDf <- literatureIFR

for (no in c(1:length(outcomes))) {
  print(paste("Outcome:", no))
  oStr <- outcomes[no]
  lowValsList <- list()
  highValsList <- list()
  outcomeFit <- lethalityFit$model[[oStr]]
  meanAge <- lethalityFit$meanAge[[oStr]]
  sdAge <- lethalityFit$sdAge[[oStr]]
  lethalityPosterior[[oStr]] <- proportion_samples(model=outcomeFit,
                                                   ageVec=ageVec,
                                                   meanAge=meanAge,
                                                   sdAge=sdAge)
  # sample the IFR vec many times, to do the ratio of ratios
  outcomeProp <- NULL
  propMeanMat <- NULL
  for (s in c(1:ifrSamples)) {
    print(paste("IFR sample:", s))
    sampleIFR <- NULL
    for (r in c(1:length(meanlog))) {
      sampleIFR[r] <- rlnorm(1, meanlog=meanlog[r], sdlog=sdlog[r])
    }
    ifrVec <- rep(sampleIFR, max(lethalityPosterior[[oStr]]$sample))
    tempOutcomeProp <- ifrVec/lethalityPosterior[[oStr]]$samples$proportion
    tempOutcomeProp <- pmin(tempOutcomeProp, 100)
    # Find mean and quantiles of samples
    outcomePropMat <- matrix(tempOutcomeProp, nrow=nrow(literatureIFR),
              ncol=max(lethalityPosterior[[oStr]]$sample))
    ind1 <- round(dim(outcomePropMat)[2]*0.025)
    ind2 <- round(dim(outcomePropMat)[2]*0.975)
    # collect the extreme values to compute CI. due to memory insufficiency
    for (r in c(1:length(meanlog))) {
      first025 <- sort(outcomePropMat[r,])[1:ind1]
      last025 <- sort(outcomePropMat[r,])[ind2:ncol(outcomePropMat)]
      if (s==1) {
        lowValsList[[r]] <- first025
        highValsList[[r]] <- last025
      } else {
        lowValsList[[r]] <- c(lowValsList[[r]], first025)
        highValsList[[r]] <- c(highValsList[[r]], last025)
      }
    }
    propMeanMat <- cbind(propMeanMat, rowMeans(outcomePropMat))
  }
  meanProp <- rowMeans(propMeanMat)
  lowerProp <- NULL
  upperProp <- NULL
  for (r in c(1:nrow(literatureIFR))) {
    lowerProp[r] <- max(lowValsList[[r]])
    upperProp[r] <- min(highValsList[[r]])
  }
  # Put into data frame
  tempDf <- data.frame(Age=literatureIFR$Age,
                      Proportion=meanProp,
                      Proportion_L=lowerProp,
                      Proportion_H=upperProp,
                      Study=literatureIFR$Study,
                      Type=oStr,
                      meanAge=ageVec)
  outputDf <- rbind(outputDf, tempDf)
}

fileName <- "../data/processed_data/5_literature_outcome_estimates.csv"
write.csv(outputDf, fileName, row.names=FALSE)

