library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
source("./functions_auxiliary.R")

regLineSize=1.2
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1

############################
############################
###
### Plot serology data fit
###
############################
############################

# load data and models
#countryData <- read.csv("../data/collected_data/locations_serology_data.csv",
#                        stringsAsFactors=FALSE) %>%
#  as_tibble(.)
#serologyModels <- readRDS("../data/processed_data/3_serology_fits.RDS")
#serologyPlotName <- "../data/plots/3_serology_regression.png"
#serologySamplesPlotName <- "../data/plots/3_serology_samples.png"
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.)
serologyModels <- readRDS("../data/processed_data/3_serology_fits_corrected.RDS")
serologyPlotName <- "../data/plots/3_serology_regression_corrected.png"
serologySamplesPlotName <- "../data/plots/3_serology_samples_corrected.png"
serologyCsvName <- "../data/processed_data/3_serology_fits_corrected.csv"
#
# predefine some variables
outcome <- c("Hospitalized", "ICU", "Deaths")

ageVec <- seq(2.5, 90, 5)
serologyPosterior <- list()
outcomeFitDf <- NULL
serologySamplesDf <- tibble()
# extract posteriors and put into data frame
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  stdAgeVec <- (ageVec-serologyModels$meanAge[[oStr]])/serologyModels$sdAge[[oStr]]
  serologyPosterior[[oStr]] <- proportion_samples(model=serologyModels$model[[oStr]],
                                                  ageVec=stdAgeVec)
  tempFitDf <- data.frame(meanAge=ageVec,
                          outcomeProp=serologyPosterior[[oStr]]$prop_mean,
                          outcome_L=serologyPosterior[[oStr]]$prop_L,
                          outcome_H=serologyPosterior[[oStr]]$prop_H,
                          Outcome_type=oStr)
  outcomeFitDf <- rbind(outcomeFitDf, tempFitDf)
  tempSamplesDf <- serologyPosterior[[oStr]]$samples
  tempSamplesDf$meanAge <- rep(ageVec, max(tempSamplesDf$sample))
  tempSamplesDf$Outcome_type <- oStr
  serologySamplesDf <- as_tibble(rbind(tempSamplesDf, serologySamplesDf))
}

outcome2 <- c("Severe", "Critical", "Fatal")
longCountryData <- tidyr::pivot_longer(data=countryData, 
                                       cols=all_of(outcome),
                                       names_to="Outcome_type",
                                       values_to="Outcomes") %>%
  dplyr::filter(., !is.na(Outcomes)) %>%
  dplyr::mutate(., Location=factor(Location))
longCountryData$Outcome_type <- factor(longCountryData$Outcome_type,
                                       levels=outcome, labels=outcome2)

outcomeFitDf$Outcome_type <- factor(outcomeFitDf$Outcome_type,
                                       levels=outcome, labels=outcome2)


# https://www.sciencemag.org/news/2021/06/israel-reports-link-between-rare-cases-heart-inflammation-and-covid-19-vaccination
myocarditisProp <- c(1/6000, 1/3000)
ageRangeMyo <- c(16, 24)
ageMeanMyo <- mean(ageRangeMyo)
myoTextX <- 55
myoTextY <- 0.005

serologyPlot <- longCountryData %>%
  ggplot(., aes(x=meanAge, y=Outcomes*100/(Population*Prevalence/100), color=Location,
                linetype=Type, facet=Outcome_type)) +
  geom_point(aes(shape=Type), size=locPointSize) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  facet_grid(.~Outcome_type) +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  geom_line(data=outcomeFitDf, aes(y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize) +
  geom_ribbon(data=outcomeFitDf,
              aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  geom_segment(aes(x=ageMeanMyo, xend=ageMeanMyo, y=myocarditisProp[1]*100,
                   yend=myocarditisProp[2]*100), color="black", size=1.5,
               inherit.aes=FALSE) +
  geom_segment(aes(x=ageMeanMyo+2, xend=myoTextX-10, y=100/4000, yend=myoTextY*1.5),
               arrow=arrow(length=unit(0.2, "cm")), color="black",
               inherit.aes=FALSE, size=0.3) +
  geom_text(aes(x=myoTextX, y=myoTextY), label="Myocarditis", color="black",
            inherit.aes=FALSE, size=3.1) +
  theme_bw() +
  xlab("Age") +
  ylab("% outcome")

ggsave(serologyPlotName, serologyPlot, width=25, height=10, units="cm")


serologySamplesDf$Outcome_type <- factor(serologySamplesDf$Outcome_type,
                                       levels=outcome)
serologySamplesPlot <- serologySamplesDf %>%
  ggplot(., aes(x=meanAge, y=proportion*100, group=sample, facet=Outcome_type)) +
  geom_line(alpha=sampleLineAlpha, size=sampleLineSize, color="#33ADFF") +
  facet_grid(.~Outcome_type) +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  geom_line(data=outcomeFitDf, aes(x=meanAge, y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize, inherit.aes=FALSE) +
#  geom_ribbon(data=outcomeFitDf,
#              aes(x=meanAge,ymin=outcome_L*100, ymax=outcome_H*100),
#              alpha=0.2, colour=NA, show.legend=FALSE,
#              inherit.aes=FALSE) +
  theme_bw() +
  xlab("Age") +
  ylab("% outcome")

ggsave(serologySamplesPlotName, serologySamplesPlot,
       width=29, height=10, units="cm")

exportFit <- dplyr::mutate(outcomeFitDf, Percentage=signif(outcomeProp*100, digits=3),
                           Percentage_L=signif(outcome_L*100, digits=3),
                           Percentage_H=signif(outcome_H*100, digits=3)) %>%
  dplyr::select(., -outcomeProp, -outcome_L, -outcome_H)
write.csv(exportFit, file=serologyCsvName, row.names=FALSE)

############################
############################
###
### Plot hospital lethality data fit
###
############################
############################

lethalityData <- read.csv("../data/collected_data/hospitalized_patient_studies.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.)
lethalityModels <- readRDS("../data/processed_data/4_hospital_lethality_fit.RDS")

# get fitted line
ageVec <- seq(2.5, 90, 5)
letType <- c("Hospitalized", "ICU")
lethalityPosterior <- list()
lethalityFitDf <- NULL
lethalitySamplesDf <- tibble()
for (no in c(1:length(letType))) {
  oStr <- letType[no]
  # extract model fit results
  stdAgeVec <- (ageVec-lethalityModels$meanAge[[oStr]])/lethalityModels$sdAge[[oStr]]
  lethalityPosterior[[oStr]] <- proportion_samples(model=lethalityModels$model[[oStr]],
                                                  ageVec=stdAgeVec)
  tempFitDf <- data.frame(meanAge=ageVec,
                          outcomeProp=lethalityPosterior[[oStr]]$prop_mean,
                          outcome_L=lethalityPosterior[[oStr]]$prop_L,
                          outcome_H=lethalityPosterior[[oStr]]$prop_H,
                          Type=oStr)
  lethalityFitDf <- rbind(lethalityFitDf, tempFitDf)
  tempSamplesDf <- lethalityPosterior[[oStr]]$samples
  tempSamplesDf$meanAge <- rep(ageVec, max(tempSamplesDf$sample))
  tempSamplesDf$Type <- oStr
  lethalitySamplesDf <- as_tibble(rbind(tempSamplesDf, lethalitySamplesDf))
}

letType2 <- c("Hospital", "ICU")
lethalityData$Type <- factor(lethalityData$Type, levels=letType, labels=letType2)
lethalityData$Location <- factor(lethalityData$Location)
lethalityFitDf$Type <- factor(lethalityFitDf$Type, levels=letType, labels=letType2)

lethalityPlot <- lethalityData %>%
  ggplot(., aes(x=meanAge, y=Deaths/Patients*100, color=Location, facet=Type)) +
  geom_point(size=2) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  facet_grid(Type~.) +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  geom_line(data=lethalityFitDf, aes(y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize) +
  geom_ribbon(data=lethalityFitDf,
              aes(x=meanAge,ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  xlab("Age") +
  ylab("Mortality (%)")

ggsave("../data/plots/4_hospital_lethality_regression.png", lethalityPlot,
       width=30, height=11, units="cm")


lethalitySamplesDf$Type <- factor(lethalitySamplesDf$Type,
                                       levels=letType)
lethalitySamplesPlot <- lethalitySamplesDf %>%
  ggplot(., aes(x=meanAge, y=proportion*100, group=sample, facet=Type)) +
  geom_line(alpha=sampleLineAlpha, size=sampleLineSize, color="#33ADFF") +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  geom_line(data=lethalityFitDf, aes(x=meanAge, y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize, inherit.aes=FALSE) +
  facet_grid(.~Type) +
#  geom_ribbon(data=outcomeFitDf,
#              aes(x=meanAge,ymin=outcome_L*100, ymax=outcome_H*100),
#              alpha=0.2, colour=NA, show.legend=FALSE,
#              inherit.aes=FALSE) +
  theme_bw() +
  xlab("Age") +
  ylab("% outcome")

ggsave("../data/plots/4_lethality_samples.png", lethalitySamplesPlot,
       width=29, height=10, units="cm")


############################
############################
###
### Plot outcome proportions estimated from literature
###
############################
############################

outcomePropLit <- read.csv("../data/processed_data/5_literature_outcome_estimates.csv",
                           stringsAsFactors=FALSE) %>%
  as_tibble(.)

outcomePropLit$Type <- factor(outcomePropLit$Type, levels=outcome,
  labels=outcome2)
outcomePropLit$Outcome_type <- outcomePropLit$Type

literaturePropPlot <- outcomePropLit %>%
  ggplot(., aes(x=meanAge, y=Proportion, color=Study)) +
  geom_point(size=locPointSize) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  facet_grid(.~Outcome_type) +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  geom_line(data=outcomeFitDf, aes(y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize) +
  geom_ribbon(data=outcomeFitDf,
              aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  theme_bw() +
  geom_segment(aes(x=ageMeanMyo, xend=ageMeanMyo, y=myocarditisProp[1]*100,
                   yend=myocarditisProp[2]*100), color="black", size=2) +
  geom_segment(aes(x=ageMeanMyo+2, xend=myoTextX-10, y=100/4000, yend=myoTextY*1.5),
               arrow=arrow(length=unit(0.2, "cm")), color="black",
               inherit.aes=FALSE, size=0.3) +
  geom_text(aes(x=myoTextX, y=myoTextY), label="Myocarditis", color="black",
            inherit.aes=FALSE, size=3.1) +
  xlab("Age") +
  ylab("% Outcome")


ggsave("../data/plots/5_literature_outcome_estimates.png", literaturePropPlot,
       width=17, height=10, units="cm")


figure2 <- ggarrange(plotlist=list(lethalityPlot, literaturePropPlot),
          ncol=2, widths=c(0.33, 0.66), labels=c("(A)", "(B)"))

ggsave("../data/plots/figure2.png", figure2, width=25, height=11, units="cm")

############################
############################
###
### Plot models params
###
############################
############################

# prior functions
gaussian <- function(x, mu, sigma){
  f <- (1/sigma*sqrt(2*pi))*exp(-0.5*((x-mu)/sigma)^2)
}
exponential <- function(x, lambda){lambda*exp(-lambda*x)}


mainParams <- c("ageSlope", "ageSlopeSigma", "intercept", "interceptSigma")

############
# plot serology Priors and trace
############
outcome <- c("Hospitalized", "ICU", "Deaths")
x1 <- seq(-1, 5, 0.02)
prior_ageSlope <- data.frame(value=x1, dens=gaussian(x1, mu=2, sigma=1),
                             .variable="ageSlope")
x2 <- seq(0, 4.5, 0.02)
prior_ageSlopeSigma <- data.frame(value=x2, dens=exponential(x2, lambda=0.5),
                                  .variable="ageSlopeSigma")
x3 <- seq(-8, 0, 0.02)
prior_intercept <- data.frame(value=x3, dens=gaussian(x3, mu=-6, sigma=2),
                              .variable="intercept")
x4 <- seq(0, 4.5, 0.02)
prior_interceptSigma <- data.frame(value=x4, dens=exponential(x4, lambda=0.5),
                                   .variable="interceptSigma")
priorDf <- rbind(prior_ageSlope, prior_ageSlopeSigma, prior_intercept,
                 prior_interceptSigma)

posteriorSerology <- data.frame()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  posteriorTemp <- tidybayes::gather_draws(serologyModels$model[[oStr]], ageSlope,
                                           ageSlopeSigma, intercept, interceptSigma) %>%
    dplyr::mutate(., Outcome_type=oStr)
  posteriorSerology <- rbind(posteriorSerology, posteriorTemp)
}
posteriorSerology$Outcome_type <- factor(posteriorSerology$Outcome_type,
                                         levels=outcome)

posteriorSerologyPlot <- posteriorSerology %>%
  ggplot(., aes(x=.value, color=Outcome_type, fill=Outcome_type,
                facet=.variable)) +
  #geom_density(alpha=0.6) +
  stat_density(geom="area", position="identity", alpha=0.6) +
  stat_density(geom="line", position="identity") +
  geom_line(data=priorDf, aes(x=value, y=dens), size=1.5, inherit.aes=FALSE) +
  facet_grid(.~.variable, scales="free") +
  theme_bw()

serologyTrace <- dplyr::mutate(posteriorSerology, .chain=factor(.chain))  %>%
  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
  geom_line() +
  facet_grid(.variable~Outcome_type, scales="free_y") +
  theme_bw()

ggsave("../data/plots/3_serology_param_posterior.png", posteriorSerologyPlot, 
       width=30, height=12, units="cm", device="png")

ggsave("../data/plots/3_serology_chains.png", serologyTrace, 
       width=30, height=20, units="cm")



############
# plot hospital lethality Priors and trace
############
posteriorLethality <- data.frame()
lethOutcome <- c("Hospitalized", "ICU")
for (no in c(1:length(lethOutcome))) {
  oStr <- lethOutcome[no]
  posteriorTemp <- tidybayes::gather_draws(lethalityModels$model[[oStr]], ageSlope,
                                           ageSlopeSigma, intercept, interceptSigma) %>%
    dplyr::mutate(., Outcome_type=oStr)
  posteriorLethality <- rbind(posteriorLethality, posteriorTemp)
}
posteriorLethality$Outcome_type <- factor(posteriorLethality$Outcome_type,
                                         levels=outcome)

posteriorLethalityPlot <- posteriorLethality %>%
  ggplot(., aes(x=.value, color=Outcome_type, fill=Outcome_type,
                facet=.variable)) +
  stat_density(geom="area", position="identity", alpha=0.6) +
  stat_density(geom="line", position="identity") +
  geom_line(data=priorDf, aes(x=value, y=dens), size=1.5, inherit.aes=FALSE) +
  facet_grid(.~.variable, scales="free") +
  theme_bw()

lethalityTrace <- dplyr::mutate(posteriorLethality, .chain=factor(.chain))  %>%
  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
  geom_line() +
  facet_grid(.variable~Outcome_type, scales="free_y") +
  theme_bw()

ggsave("../data/plots/4_lethality_param_posterior.png", posteriorLethalityPlot, 
       width=30, height=12, units="cm", device="png")

ggsave("../data/plots/4_lethality_chains.png", lethalityTrace, 
       width=30, height=20, units="cm")

