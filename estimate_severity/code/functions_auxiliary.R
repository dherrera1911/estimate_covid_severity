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
  ageBins <- strsplit(binsVec, "-") %>%
    lapply(., as.numeric) %>%
    lapply(., mean) %>%
    unlist(.)
  for (naInd in which(is.na(ageBins))) {
    naVal <- as.numeric(substr(binsVec[naInd], start=1, stop=2))
    ageBins[naInd] <- mean(c(naVal, 90))
  }
  # extract last value
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
    slopeSample <- posterior[[slopeName]][n]
    lin <- intSample + slopeSample * stdAgeVec
    fitProp <- exp(lin)/(1+exp(lin))
    fitSampleMat <- cbind(fitSampleMat, as.matrix(fitProp))
  }
  meanProp <- rowMeans(fitSampleMat) 
  ciProp <- matrixStats::rowQuantiles(fitSampleMat, probs=c(0.025, 0.975))
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
  


# Get posterior samples of bayesian fit and put into matrix
proportion_samples <- function(model, ageVec,
                               slopeName="ageSlope",
                               interceptName="intercept",
                               meanAge=0,
                               sdAge=1,
                               link="logit") {
  stdAgeVec <- (ageVec-meanAge)/sdAge
  posterior <- rstan::extract(model)
  fitSampleMat <- matrix(nrow=length(ageVec), ncol=0)
  for (n in c(1:length(posterior[[1]]))) {
    intSample <- posterior[[interceptName]][n]
    slopeSample <- posterior[[slopeName]][n]
    lin <- intSample + slopeSample * stdAgeVec
    if (link=="logit") {
      fitProp <- exp(lin)/(1+exp(lin))
    } else {
      fitProp <- VGAM::probitlink(theta=lin, inverse=TRUE)
    }

    fitSampleMat <- cbind(fitSampleMat, as.matrix(fitProp))
  }
  meanProp <- rowMeans(fitSampleMat) 
  ciProp <- matrixStats::rowQuantiles(fitSampleMat, probs=c(0.025, 0.975))
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


# Fit gamma distribution to CIs
gammaPars <- function(casesQuantiles, quants=c(0.025, 0.5, 0.975)) { #
  pars <- get.gamma.par(p=quants, q=casesQuantiles,
                        plot=FALSE, verbose=FALSE, show.output=FALSE)
}

# Fit gamma distribution to CIs
fit_gamma_ci <- function(meanEstimate, lower, upper) {
  gammaShape <- NULL
  gammaRate <- NULL
  for (r in c(1:length(meanEstimate))) {
    pars <- NA
    if (lower[r]!=meanEstimate[r] & upper[r]!=meanEstimate[r]) {
      if (meanEstimate[r]<0.01) {
        pars <- gammaPars(c(lower[r], meanEstimate[r], upper[r])*100)
        pars[2] <- pars[2]*100
      } else if (meanEstimate[r]>10) {
        pars <- gammaPars(c(lower[r], meanEstimate[r], upper[r])/10)
        pars[2] <- pars[2]/10
      } else {
        pars <- gammaPars(c(lower[r], meanEstimate[r], upper[r]))
      }
    }
    if (lower[r]==meanEstimate[r] | upper[r]==meanEstimate[r]) {
      pars <- gammaPars(c(lower[r], upper[r])*100, quants=c(0.025, 0.975))
      pars[2] <- pars[2]*100
    }
    gammaShape[r] <- pars[1]
    gammaRate[r] <- pars[2]
  }
  return(list(gammaShape=gammaShape, gammaRate=gammaRate))
}


## Fit beta distribution to CIs
betaPars <- function(casesQuantiles, quants=c(0.025, 0.5, 0.975)) { #
  pars <- get.beta.par(p=quants, q=casesQuantiles,
                       plot=FALSE, verbose=FALSE, show.output=FALSE,
                       tol=0.01)
}

fit_beta_ci <- function(meanEstimate, lower, upper) {
  shape1 <- NULL
  shape2 <- NULL
  for (r in c(1:length(meanEstimate))) {
    pars <- NA
    if (lower[r]!=meanEstimate[r] & upper[r]!=meanEstimate[r]) {
      pars <- betaPars(c(lower[r], meanEstimate[r], upper[r]))
    }
    if (lower[r]==meanEstimate[r] | is.na(pars[1])) {
      pars <- betaPars(c(lower[r], upper[r]), quants=c(0.025, 0.975))
    }
    shape1[r] <- pars[1]
    shape2[r] <- pars[2]
  }
  return(list(shape1=shape1, shape2=shape2))
}


########
# Home made function for getting beta parameters
# https://stats.stackexchange.com/questions/112614/determining-beta-distribution-parameters-alpha-and-beta-from-two-arbitrary
#######
#
# Logistic transformation of the Beta CDF.
#
f.beta <- function(alpha, beta, x, lower=0, upper=1) {
  p <- pbeta((x-lower)/(upper-lower), alpha, beta)
  log(p/(1-p))
}

# Sums of squares.

delta <- function(fit, actual) sum((fit-actual)^2)

# The objective function handles the transformed parameters `theta` and
# uses `f.beta` and `delta` to fit the values and measure their discrepancies.

objective <- function(theta, x, prob, ...) {
  ab <- exp(theta) # Parameters are the *logs* of alpha and beta
  fit <- f.beta(ab[1], ab[2], x, ...)
  return (delta(fit, prob))
}



#
# Solve two problems.
#
par(mfrow=c(1,2))
alpha <- 15; beta <- 22 # The true parameters
for (x in list(c(1e-3, 2e-3), c(1/3, 2/3))) {
  x.p <- f.beta(alpha, beta, x)        # The correct values of the p_i
  start <- log(c(1e1, 1e1))            # A good guess is useful here
  sol <- nlm(objective, start, x=x, prob=x.p, lower=0, upper=1,
             typsize=c(1,1), fscale=1e-12, gradtol=1e-12)
  parms <- exp(sol$estimate)           # Estimates of alpha and beta
  #
  # Display the actual and estimated values.
  #
  print(rbind(Actual=c(alpha=alpha, beta=beta), Fit=parms))
  #
  # Plot the true and estimated CDFs.
  #      
  curve(pbeta(x, alpha, beta), 0, 1, n=1001, lwd=2)
  curve(pbeta(x, parms[1], parms[2]), n=1001, add=TRUE, col="Red")
  points(x, pbeta(x, alpha, beta))
}



