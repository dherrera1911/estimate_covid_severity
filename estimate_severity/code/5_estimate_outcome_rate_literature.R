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

outcomes <- c("Hospitalized", "ICU")
ageVec <- literatureIFR$meanAge
lethalityPosterior <- list()
outputDf <- literatureIFR
for (no in c(1:length(outcomes))) {
  oStr <- outcomes[no]
  outcomeFit <- lethalityFit$model[[oStr]]
  meanAge <- lethalityFit$meanAge[[oStr]]
  sdAge <- lethalityFit$sdAge[[oStr]]
  lethalityPosterior[[oStr]] <- proportion_samples(model=outcomeFit,
                                                   ageVec=ageVec,
                                                   meanAge=meanAge,
                                                   sdAge=sdAge)
  ifrVec <- rep(literatureIFR$Proportion, max(lethalityPosterior[[oStr]]$sample))
  outcomeProp <- ifrVec/lethalityPosterior[[oStr]]$samples$proportion
  # Find mean and quantiles of samples
  outcomePropMat <- matrix(outcomeProp, nrow=nrow(literatureIFR),
            ncol=max(lethalityPosterior[[oStr]]$sample))
  prop_mean <- rowMeans(outcomePropMat)
  ciProp <- matrixStats::rowQuantiles(outcomePropMat, probs=c(0.025, 0.975))
  # Put into data frame
  tempDf <- data.frame(Age=literatureIFR$Age,
                      Proportion=prop_mean,
                      Proportion_L=ciProp[,1],
                      Proportion_H=ciProp[,2],
                      Study=literatureIFR$Study,
                      Type=oStr,
                      meanAge=ageVec)
  outputDf <- rbind(outputDf, tempDf)
}

fileName <- "../data/processed_data/5_literature_outcome_estimates.csv"
write.csv(outputDf, fileName, row.names=FALSE)

