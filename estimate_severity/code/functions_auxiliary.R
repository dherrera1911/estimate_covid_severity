library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidybayes)
library(rstan)
library(matrixStats)

# Return a vector with the mean age in each bin
mid_bin_age <- function(binsVec) {
  defaultW <- getOption("warn")
  options(warn = -1)
  ageBins <- strsplit(binsVec, "-") %>%
    lapply(., as.numeric) %>%
    lapply(., mean) %>%
    unlist(.)
  for (naInd in which(is.na(ageBins))) {
    naVal <- as.numeric(substr(binsVec[naInd], start=1, stop=2))
    ageBins[naInd] <- mean(c(naVal, 90))
  }
  # extract last value
  options(warn = defaultW)
  return(ageBins)
}

# Take a demography dataframe and a vector indicating a
# new subdivission, and return demography for the new division
change_demography_bins <- function(demographyDf, newBins) {
  midAges <- mid_bin_age(as.character(demographyDf$age))
  newDemRow <- bin_ages(midAges, newBins)
  newDemography <- dplyr::mutate(demographyDf, newBin = newDemRow) %>%
    group_by(newBin) %>%
    summarize(., proportion = sum(proportion)) %>%
    ungroup(.) %>%
    dplyr::mutate(age = newBins) %>%
    dplyr::select(-newBin)
  return(newDemography)
}


# take a vector with numbers and bin them as given in
# binsVec. binsVec can be a character vector indicating
# the ranges as "Xi-Xf" or a numeric vector indicating
# the superior limits of the bins
bin_ages <- function(agesVec, binsVec) {
  if (is.character(binsVec)) {
    ageBins <- strsplit(binsVec, "-") %>%
      lapply(., as.numeric) %>%
      unlist(.)
    supInd <-seq(2, length(ageBins), by = 2)
    ageBins <- ageBins[supInd]
    if (any(is.na(ageBins))) {
      ageBins <- ageBins[-which(is.na(ageBins))]
    }
    ageBins <- c(-1, ageBins, 200)
  } else {
    ageBins <- binsVec
  }
  indVec <- .bincode(agesVec, ageBins)
  return(indVec)
}



get_bins_limits <- function(ageBins) {
  ageList <- strsplit(ageBins, "-") %>%
    lapply(., as.numeric)
  ageLow <- NULL
  ageHigh <- NULL
  for (l in c(1:length(ageList))) {
    if (is.na(ageList[[l]])) {
      naVal <- as.numeric(substr(ageBins[l], start=1, stop=2))
      ageList[[l]] <- c(naVal, 300)
    }
    ageLow[l] <- ageList[[l]][1]
    ageHigh[l] <- ageList[[l]][2]
  }
  return(list(lower=ageLow, upper=ageHigh))
}

binomial_confint <- function(countTotal, occurrences, input="count"){
  lower <- NULL
  upper <- NULL
  if (input=="count") {
    countOccurrences <- occurrences
  } else {
    countOccurrences <- round(occurrences * countTotal)
  }
  for (i in c(1:length(countTotal))) {
    confint <- binom.test(countOccurrences[i], countTotal[i])$conf.int
    lower[i] <- confint[1]
    upper[i] <- confint[2]
  }
  confintList <- list(lower=lower, upper=upper)
}


extract_country_population <- function(popM, popF, countryName, ageBins) {
  countryPopM <- dplyr::filter(popM, name==countryName) %>%
    dplyr::select(., age, "2020")
  countryPopF <- dplyr::filter(popF, name==countryName) %>%
    dplyr::select(age, "2020")
  countryPop <- data.frame(age=countryPopM$age,
                         pop=(countryPopM[["2020"]]+countryPopF[["2020"]])*1000)
  countryPop$proportion <- countryPop$pop # rename pop to proportion so function works 
  countryPop <- change_demography_bins(countryPop, ageBins) %>%
    rename(., pop=proportion)
  return(round(countryPop$pop))
}



# function used for plotting log axis
scaleFun <- function(x) sprintf("%1g", x)


# Get posterior samples of bayesian fit
proportion_samples <- function(model, ageVec,
                               slopeName="ageSlope",
                               interceptName="intercept",
                               meanAge=0,
                               sdAge=1) {
  stdAgeVec <- (ageVec-meanAge)/sdAge
  posterior <- rstan::extract(model)
  fitSampleMat <- matrix(nrow=length(ageVec), ncol=0)
  for (n in c(1:length(posterior[[1]]))) {
    intSample <- posterior[[interceptName]][n]
    if (is.na(slopeName)) {
      slopeSample <- 0
    } else {
      slopeSample <- posterior[[slopeName]][n]
    }
    lin <- intSample + slopeSample * stdAgeVec
    fitProp <- exp(lin)/(1+exp(lin))
    fitSampleMat <- cbind(fitSampleMat, as.matrix(fitProp))
  }
  meanProp <- rowMeans(fitSampleMat) 
  ciProp <- matrixStats::rowQuantiles(fitSampleMat, probs=c(0.025, 0.975))
  if (!is.matrix(ciProp)) { # if only one param, this is needed
    ciProp <- t(as.matrix(ciProp))
  }
  sampleVec <- sort(rep(c(1:ncol(fitSampleMat)), nrow(fitSampleMat)))
  ageVecLong <- rep(ageVec, ncol(fitSampleMat))
  ageIndVec <- rep(c(1:length(ageVec)), ncol(fitSampleMat))

  samplesDf <- data.frame(sample=sampleVec,
                          proportion=as.vector(fitSampleMat),
                          age=ageVecLong, ageInd=ageIndVec)
  sampleList <- list(samples=samplesDf, prop_mean=meanProp,
                     prop_L=ciProp[,1], prop_H=ciProp[,2])
  return(sampleList)
}
  


## Get posterior samples of bayesian fit and put into matrix
#proportion_samples <- function(model, ageVec,
#                               slopeName="ageSlope",
#                               interceptName="intercept",
#                               meanAge=0,
#                               sdAge=1,
#                               link="logit") {
#  stdAgeVec <- (ageVec-meanAge)/sdAge
#  posterior <- rstan::extract(model)
#  fitSampleMat <- matrix(nrow=length(ageVec), ncol=0)
#  for (n in c(1:length(posterior[[1]]))) {
#    intSample <- posterior[[interceptName]][n]
#    slopeSample <- posterior[[slopeName]][n]
#    lin <- intSample + slopeSample * stdAgeVec
#    if (link=="logit") {
#      fitProp <- exp(lin)/(1+exp(lin))
#    } else {
#      fitProp <- VGAM::probitlink(theta=lin, inverse=TRUE)
#    }
#
#    fitSampleMat <- cbind(fitSampleMat, as.matrix(fitProp))
#  }
#  meanProp <- rowMeans(fitSampleMat) 
#  ciProp <- matrixStats::rowQuantiles(fitSampleMat, probs=c(0.025, 0.975))
#  sampleVec <- sort(rep(c(1:ncol(fitSampleMat)), nrow(fitSampleMat)))
#  ageVecLong <- rep(ageVec, ncol(fitSampleMat))
#  ageIndVec <- rep(c(1:length(ageVec)), ncol(fitSampleMat))
#
#  samplesDf <- data.frame(sample=sampleVec,
#                          proportion=as.vector(fitSampleMat),
#                          age=ageVecLong, ageInd=ageIndVec)
#  sampleList <- list(samples=samplesDf, prop_mean=meanProp,
#                     prop_L=ciProp[,1], prop_H=ciProp[,2])
#  return(sampleList)
#}


# estimate the number of out of hospital (or ICU) deaths
# from the hospital mortality fit, hospitalizations (or ICU) and deaths
ooh_deaths_estimation <- function(mortalitySamples,
                                  hospitalized,
                                  deaths) {
  sampleMat <- matrix(nrow=length(hospitalized), ncol=0)
  for (s in unique(mortalitySamples$samples$sample)) {
    sampleMort <- dplyr::filter(mortalitySamples$samples, sample==s)
    oohDeaths <- deaths - hospitalized*sampleMort$proportion
    sampleMat <- cbind(sampleMat, oohDeaths)
  }
  oohDeaths <- pmax(0, round(rowMeans(sampleMat)))
  oohDeaths_ci <- round(matrixStats::rowQuantiles(sampleMat, probs=c(0.025, 0.975)))
  oohDeaths_ci[,1] <- pmax(0, oohDeaths_ci[,1])
  oohDeaths_ci[,2] <- pmax(0, oohDeaths_ci[,2])
  oohDeathsDf <- data.frame(mean=oohDeaths, lower=oohDeaths_ci[,1],
                          upper=oohDeaths_ci[,2])
  return(oohDeathsDf)
}



# beta fitting by moment matching
fit_beta_ci <- function(meanEstimate, lower, upper) {
  var <- ((upper-lower)/4)^2
  commonFactor <- meanEstimate*(1-meanEstimate)/var - 1
  shape1 <- meanEstimate*commonFactor
  shape2 <- (1-meanEstimate)*commonFactor
  fittingDf <- data.frame(meanEstimate=meanEstimate, lower=lower, upper=upper,
                    shape1=NA, # stores estimates of shape1
                    shape2=NA, # same for shape 2
                    priorMean=NA, # the mean of the fitted prior
                    priorLB=NA,  # the 0.025 quantile of the fitted prior
                    priorUB=NA) # the 0.975 quantile of the fitted prior
  fittingDf$shape1 <- shape1
  fittingDf$shape2 <- shape2
  fittingDf$priorMean <- shape1/(shape1+shape2)
  for (i in c(1:nrow(fittingDf))) {
    shape1i <- fittingDf$shape1[i]
    shape2i <- fittingDf$shape2[i]
    fittingDf$priorLB[i] <- qbeta(p=.025, shape1=shape1i, shape2=shape2i)
    fittingDf$priorUB[i] <- qbeta(p=.975, shape1=shape1i, shape2=shape2i)
  }
  return(fittingDf)
}


# Extract stan params means and CrI undoing predictor normalization
extract_model_params_norm <- function(model, xCenter, xSd) {
  posterior <- rstan::extract(model)
  slopeSamples <- posterior$ageSlope
  interceptSamples <- posterior$intercept
  slopeSamplesOrig <- slopeSamples/xSd 
  interceptSamplesOrig <- interceptSamples - slopeSamples*xCenter/xSd
  paramNames <- c("ageSlope", "intercept")
  paramMeans <- c(mean(slopeSamplesOrig), mean(interceptSamplesOrig))
  ciL <- c(quantile(slopeSamplesOrig, probs=0.025),
              quantile(interceptSamplesOrig, probs=0.025))
  ciH <- c(quantile(slopeSamplesOrig, probs=0.975),
              quantile(interceptSamplesOrig, probs=0.975))
  paramDf <- data.frame(param=paramNames, meanVal=paramMeans,
                        lower=ciL, higher=ciH)
  return(paramDf)
}


