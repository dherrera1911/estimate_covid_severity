library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(lemon)
source("./functions_auxiliary.R")

regLineSize=0.8
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5

############################
############################
###
### Plot corrected vs uncorrected data
###
############################
############################

# load data and models
# Uncorrected data
countryDataUncorrected <- read.csv("../data/collected_data/locations_serology_data.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., processing="Uncorrected") %>%
  dplyr::arrange(., Location, Age, Hospitalized, ICU)

serologyModelsUncorrected <- readRDS("../data/processed_data/3_serology_fits.RDS")

# Corrected data
countryDataCorrected <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., processing="Corrected") %>%
  dplyr::arrange(., Location, Age, Hospitalized, ICU)

serologyModelsCorrected <- readRDS("../data/processed_data/3_serology_fits_corrected.RDS")

countryDataUncorrected$oohDeaths_source <- countryDataCorrected$oohDeaths_source
countryDataUncorrected$ooiDeaths_source <- countryDataCorrected$ooiDeaths_source


countryData <- dplyr::select(countryDataCorrected,
                             all_of(names(countryDataUncorrected))) %>%
  rbind(., countryDataUncorrected)

serologyModels <- list(Corrected=serologyModelsCorrected,
                    Uncorrected=serologyModelsUncorrected)


# predefine some variables
outcome <- c("Hospitalized", "ICU", "Deaths")
outcome2 <- c("Severe disease", "Critical disease", "Fatal disease")
labelsType <- c("Representative seroprevalence", "Convenience seroprevalence",
                "Comprehensive testing")

# enlongate data
longCountryData <- tidyr::pivot_longer(data=countryData, 
                                       cols=all_of(outcome),
                                       names_to="Outcome_type",
                                       values_to="Outcomes") %>%
  dplyr::filter(., !is.na(Outcomes)) %>%
  dplyr::mutate(., Location=factor(Location))

longCountryData$Outcome_type <- factor(longCountryData$Outcome_type,
                                       levels=outcome, labels=outcome2)
levsType <- levels(factor(longCountryData$Type))
longCountryData$Type <- factor(longCountryData$Type,
                                       levels=levsType,
                                       labels=labelsType)

england <- filter(longCountryData, Location=="England" & Outcome_type=="Critical disease")

longCountryData$correction_source <- NA
criticalRows <- which(longCountryData$Outcome_type == "Critical disease")
longCountryData$correction_source[criticalRows] <-
  longCountryData$ooiDeaths_source[criticalRows]
severeRows <- which(longCountryData$Outcome_type == "Severe disease")
longCountryData$correction_source[severeRows] <-
  longCountryData$oohDeaths_source[severeRows]

dataCorrectedInd <- with(longCountryData, which(!is.na(correction_source) &
                                                correction_source=="Data"))
longCountryData$availableData_oohD <- "No"
longCountryData$availableData_oohD[dataCorrectedInd] <- "Yes"

# pivot longer the serology data
outcomeFitDf <- NULL
for (n in names(serologyModels)) {
  # extract model posteriors and put into data frame
  ageVec <- seq(2.5, 90, 5)
  serologyPosterior <- list()
  serologySamplesDf <- tibble()
  for (no in c(1:length(outcome))) {
    oStr <- outcome[no]
    stdAgeVec <- (ageVec-serologyModels[[n]]$meanAge[[oStr]])/serologyModels[[n]]$sdAge[[oStr]]
    serologyPosterior[[oStr]] <- proportion_samples(model=serologyModels[[n]]$model[[oStr]],
                                                    ageVec=stdAgeVec)
    tempFitDf <- data.frame(meanAge=ageVec,
                            outcomeProp=serologyPosterior[[oStr]]$prop_mean,
                            outcome_L=serologyPosterior[[oStr]]$prop_L,
                            outcome_H=serologyPosterior[[oStr]]$prop_H,
                            Outcome_type=oStr,
                            processing=n)
    outcomeFitDf <- rbind(outcomeFitDf, tempFitDf)
  }
}

levsNames <- c("Hospitalized", "ICU", "Deaths")
labelsNames <- c("Severe disease", "Critical disease", "Deaths")
outcomeFitDf$Outcome_type <- factor(outcomeFitDf$Outcome_type,
                                       levels=levsNames, labels=outcome2)

corrected_uncorrectedPlot <- longCountryData %>%
  dplyr::filter(., Outcome_type!="Fatal disease") %>%
  droplevels(.) %>%
  ggplot(., aes(x=meanAge, y=Outcomes*100/(Population*Prevalence/100),
                color=Type, group=as.character(Location), facet=Outcome_type,
                linetype=availableData_oohD, shape=availableData_oohD)) +
  geom_point(size=locPointSize, alpha=dataAlpha) +
  geom_line(alpha=locLineAlpha, size=locLineSize, alpha=dataAlpha) +
  facet_rep_grid(Outcome_type~processing, repeat.tick.labels="left") +
  geom_line(data=droplevels(dplyr::filter(outcomeFitDf, Outcome_type!="Fatal disease")),
            aes(x=meanAge, y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize,
            inherit.aes=FALSE) +
  geom_ribbon(data=droplevels(dplyr::filter(outcomeFitDf, Outcome_type!="Fatal disease")),
              aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  scale_y_continuous(trans='log10', labels=scaleFun,
                     breaks=10^c(-3, -2, -1, 0, 1, 2)) +
  theme(strip.background=element_rect(fill="white", color="white"),
        strip.text=element_text(face="bold"),
        legend.position="top",
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour = "black"),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid"),
        legend.title=element_blank()) +
  guides(shape=FALSE, linetype=FALSE) +
  xlab("Age") +
  ylab("% Infected with outcome")

ggsave("../data/plots/figureS1.png", corrected_uncorrectedPlot,
       width=18, height=16, units="cm")


############################
############################
###
### Plot the proportion change by the correction method
###
############################
############################

countryCorrLong <- dplyr::mutate(countryDataCorrected,
                                      propCorrICU=ooiDeaths/ICU*100,
                                      propCorrHosp=oohDeaths/Hospitalized*100) %>%
  dplyr::select(., -Hospitalized, -ICU, -Deaths) %>%
  tidyr::pivot_longer(data=.,
                      cols=c("propCorrICU", "propCorrHosp"),
                      names_to="Outcome_type",
                      values_to="propCorrected")

propLabels <- c("Severe", "Critical")
propLevels <- levels(factor(countryCorrLong$Outcome_type))
countryCorrLong$Outcome_type <- factor(countryCorrLong$Outcome_type,
                                       levels=propLevels, labels=propLabels)

levsType <- levels(factor(countryCorrLong$Type))
countryCorrLong$Type <- factor(countryCorrLong$Type,
                                       levels=levsType,
                                       labels=labelsType)

countryCorrLong$correction_source <- NA
criticalRows <- which(countryCorrLong$Outcome_type=="Critical")
countryCorrLong$correction_source[criticalRows] <-
  countryCorrLong$ooiDeaths_source[criticalRows]
severeRows <- which(countryCorrLong$Outcome_type=="Severe")
countryCorrLong$correction_source[severeRows] <-
  countryCorrLong$oohDeaths_source[severeRows]

dataCorrectedInd <- with(countryCorrLong, which(!is.na(correction_source) &
                                                correction_source=="Data"))
countryCorrLong$availableData_oohD <- "Estimated"
countryCorrLong$availableData_oohD[dataCorrectedInd] <- "Data"


correctionPropPlot <- list()
for (no in c(1:length(propLabels))) {
  outcomeName <- propLabels[no]
  correctionPropPlot[[outcomeName]] <- countryCorrLong %>%
    dplyr::filter(., Outcome_type==outcomeName & !is.na(propCorrected) &
                  Location != "South Korea") %>%
    droplevels(.) %>%
    ggplot(., aes(x=meanAge, y=propCorrected,
                  color=factor(availableData_oohD), group=as.character(Location),
                  facet=Outcome_type)) +
    geom_point(size=locPointSize, alpha=dataAlpha) +
    geom_line(alpha=locLineAlpha, size=locLineSize, alpha=dataAlpha) +
  #  facet_rep_grid(.~Outcome_type, repeat.tick.labels="left")
    #scale_color_brewer(palette="Dark2") +
    theme_bw() +
#    scale_y_continuous(trans='log10', labels=scaleFun,
#                       limits=c(0,100)) +
    scale_color_brewer(palette="Set1") +
    theme(strip.background=element_rect(fill="white", color="white"),
          strip.text=element_text(face="bold"),
          legend.position="top",
          panel.border=element_blank(),
          #panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line=element_line(colour = "black"),
          axis.line.x=element_line(size=0.5, linetype="solid"),
          axis.line.y=element_line(size=0.5, linetype="solid"),
          legend.title=element_blank()) +
    guides(shape=FALSE, linetype=FALSE) +
      xlab("Age") +
      ylim(0, 100) +
      NULL
}

correctionPropPlot[[1]] <- correctionPropPlot[[1]] +
  ylab(expression(over("Out of hospital deaths", "Severe cases")%*%100))
correctionPropPlot[[2]] <- correctionPropPlot[[2]] +
  ylab(expression(over("Out of ICU deaths", "Critical cases")%*%100))

correctionPropPlot <- ggpubr::ggarrange(plotlist=correctionPropPlot, ncol=2,
                          labels=c("A)", "B)"), common.legend=TRUE)

ggsave("../data/plots/figureS2.png", correctionPropPlot,
       width=20, height=10, units="cm")



############################
############################
###
### Plot fits to younger ages
###
############################
############################

serologyModelsU50 <- readRDS("../data/processed_data/3_serology_fits_u50.RDS")
serologyModelsU40 <- readRDS("../data/processed_data/3_serology_fits_u40.RDS")

serologyModels2 <- list("Uncorrected_U50"=serologyModelsU50,
                        "Uncorrected_U40"=serologyModelsU40)

# pivot longer the serology data
outcomeFit2Df <- NULL
for (n in names(serologyModels2)) {
  # extract model posteriors and put into data frame
  ageVec <- seq(2.5, 90, 5)
  serologyPosterior <- list()
  serologySamplesDf <- tibble()
  for (no in c(1:length(outcome))) {
    oStr <- outcome[no]
    stdAgeVec <- (ageVec-serologyModels2[[n]]$meanAge[[oStr]])/serologyModels2[[n]]$sdAge[[oStr]]
    serologyPosterior[[oStr]] <- proportion_samples(model=serologyModels2[[n]]$model[[oStr]],
                                                    ageVec=stdAgeVec)
    tempFitDf <- data.frame(meanAge=ageVec,
                            outcomeProp=serologyPosterior[[oStr]]$prop_mean,
                            outcome_L=serologyPosterior[[oStr]]$prop_L,
                            outcome_H=serologyPosterior[[oStr]]$prop_H,
                            Outcome_type=oStr,
                            processing=n)
    outcomeFit2Df <- rbind(outcomeFit2Df, tempFitDf)
  }
}

levsNames <- c("Hospitalized", "ICU", "Deaths")
labelsNames <- c("Severe disease", "Critical disease", "Deaths")
outcomeFit2Df$Outcome_type <- factor(outcomeFit2Df$Outcome_type,
                                       levels=levsNames, labels=outcome2)

allModels <- rbind(outcomeFitDf, outcomeFit2Df)

colorMal <- c("black", brewer.pal(3, "Set1"))

differentModelsPlot <- allModels %>%
  dplyr::filter(., Outcome_type!="Fatal disease") %>%
  droplevels(.) %>%
  ggplot(., aes(x=meanAge, y=outcomeProp*100,
                color=processing, fill=processing, facet=Outcome_type)) +
  facet_rep_grid(.~Outcome_type, repeat.tick.labels="left") +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha*0.7, show.legend=FALSE, color=NA) +
  #scale_color_brewer(palette="Dark2") +
  scale_color_manual(values=colorMal) +
  scale_fill_manual(values=colorMal) +
  theme_bw() +
  scale_y_continuous(trans='log10', labels=scaleFun,
                     breaks=10^c(-3, -2, -1, 0, 1, 2)) +
  theme(strip.background=element_rect(fill="white", color="white"),
        strip.text=element_text(face="bold"),
        legend.position="top",
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour = "black"),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid"),
        legend.title=element_blank()) +
  guides(shape=FALSE, linetype=FALSE) +
  xlab("Age") +
  ylab("% Infected with outcome")


ggsave("../data/plots/figureS3.png", differentModelsPlot,
       width=16, height=10, units="cm")


############################
############################
###
### Plot fits without Testing or without convenience serosampling
###
############################
############################

serologyModels_noTest <- readRDS("../data/processed_data/3_serology_fits_corrected_noTest.RDS")
serologyModels_noConv <- readRDS("../data/processed_data/3_serology_fits_corrected_noConvenience.RDS")

serologyModels3 <- list(Corrected=serologyModelsCorrected,
                        No_test=serologyModels_noTest,
                        No_convenience=serologyModels_noConv)

# pivot longer the serology data
outcomeFit3Df <- NULL
for (n in names(serologyModels3)) {
  # extract model posteriors and put into data frame
  ageVec <- seq(2.5, 90, 5)
  serologyPosterior <- list()
  serologySamplesDf <- tibble()
  for (no in c(1:length(outcome))) {
    oStr <- outcome[no]
    stdAgeVec <- (ageVec-serologyModels3[[n]]$meanAge[[oStr]])/serologyModels3[[n]]$sdAge[[oStr]]
    serologyPosterior[[oStr]] <- proportion_samples(model=serologyModels3[[n]]$model[[oStr]],
                                                    ageVec=stdAgeVec)
    tempFitDf <- data.frame(meanAge=ageVec,
                            outcomeProp=serologyPosterior[[oStr]]$prop_mean,
                            outcome_L=serologyPosterior[[oStr]]$prop_L,
                            outcome_H=serologyPosterior[[oStr]]$prop_H,
                            Outcome_type=oStr,
                            processing=n)
    outcomeFit3Df <- rbind(outcomeFit3Df, tempFitDf)
  }
}

outcomeFit3Df <- rbind(outcomeFit3Df, outcomeFitDf


levsNames <- c("Hospitalized", "ICU", "Deaths")
labelsNames <- c("Severe disease", "Critical disease", "Deaths")
outcomeFit2Df$Outcome_type <- factor(outcomeFit2Df$Outcome_type,
                                       levels=levsNames, labels=outcome2)

allModels <- rbind(outcomeFitDf, outcomeFit2Df)

colorMal <- c("black", brewer.pal(3, "Set1"))

differentModelsPlot <- allModels %>%
  dplyr::filter(., Outcome_type!="Fatal disease") %>%
  droplevels(.) %>%
  ggplot(., aes(x=meanAge, y=outcomeProp*100,
                color=processing, fill=processing, facet=Outcome_type)) +
  facet_rep_grid(.~Outcome_type, repeat.tick.labels="left") +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha*0.7, show.legend=FALSE, color=NA) +
  #scale_color_brewer(palette="Dark2") +
  scale_color_manual(values=colorMal) +
  scale_fill_manual(values=colorMal) +
  theme_bw() +
  scale_y_continuous(trans='log10', labels=scaleFun,
                     breaks=10^c(-3, -2, -1, 0, 1, 2)) +
  theme(strip.background=element_rect(fill="white", color="white"),
        strip.text=element_text(face="bold"),
        legend.position="top",
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour = "black"),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid"),
        legend.title=element_blank()) +
  guides(shape=FALSE, linetype=FALSE) +
  xlab("Age") +
  ylab("% Infected with outcome")


ggsave("../data/plots/figureS3.png", differentModelsPlot,
       width=16, height=10, units="cm")


############################
############################
###
### Plot death dynamics around outcome collection point
###
############################
############################

deathDynamics <- read.csv("../data/collected_data/death_dynamics_countries.csv",
  stringsAsFactors=FALSE) %>%
  dplyr::mutate(deathProp=NA)

deathDynamics2 <- NULL
for (locN in c(1:length(unique(deathDynamics$Location)))) {
  loc <- unique(deathDynamics$Location)[locN]
  locDf <- dplyr::filter(deathDynamics, Location==loc)
  deathsRef <- locDf$deaths[locDf$daysSinceData==0]
  locDf$deathProp <- 100*locDf$deaths/deathsRef
  head(locDf)
  deathDynamics2 <- rbind(deathDynamics2, locDf)
}
deathDynamics2$Location <- factor(deathDynamics2$Location)

deathsDynPlot <- deathDynamics2 %>%
  ggplot(., aes(x=daysSinceData, y=deathProp, color=Location)) +
  geom_line(size=locLineSize) +
  theme_bw() +
  geom_segment(aes(x=-30, y=120, xend=15, yend=120), linetype=2,
               color="black") +
  theme(strip.background=element_rect(fill="white", color="white"),
        strip.text=element_text(face="bold"),
        legend.position="top",
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour = "black"),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid"),
        legend.title=element_blank()) +
  xlab("Days since outcome data date") +
  ylab("% cumulative deaths (reference to data date)") +
  ylim(0, 150)

ggsave("../data/plots/figureS4.png", deathsDynPlot,
       width=18, height=15, units="cm")

lastDayDyn <- dplyr::filter(deathDynamics2, daysSinceData==15) %>%
  dplyr::arrange(., deathProp)
nLoc <- nrow(lastDayDyn)
fastestCountries <- droplevels(lastDayDyn[(nLoc-2):nLoc,"Location"])











#########################
## Plot serology model params
#########################
#mainParams <- c("ageSlope", "ageSlopeSigma", "intercept", "interceptSigma")
#
## Get param summary statistics
#modelSummary <- list()
#for (no in c(1:length(outcome))) {
#  oStr <- outcome[no]
#  modelSummary[[oStr]] <- summary(serologyModels$model[[oStr]],
#                          pars=mainParams)
#}
#
## Compare ages
#ratiosAges <- data.frame()
#for (no in c(1:length(outcome))) {
#  oStr <- outcome[no]
#  agesComparison <- proportion_samples(model=serologyModels$model[[oStr]],
#           ageVec=c(22.5, 72.5),
#           meanAge=serologyModels$meanAge[[oStr]],
#           sdAge=serologyModels$sdAge[[oStr]])
#  ratiosTemp <- agesComparison$samples %>%
#    dplyr::select(., -age) %>%
#    tidyr::pivot_wider(., names_from=ageInd, values_from=proportion,
#                       names_prefix="age") %>%
#    dplyr::mutate(., ratio=age2/age1)
#  tempDf <- data.frame(meanRatio=mean(ratiosTemp$ratio),
#                       lowerRatio=quantile(ratiosTemp$ratio, p=0.025),
#                       upperRatio=quantile(ratiosTemp$ratio, p=0.975))
#  ratiosAges <- rbind(ratiosAges, tempDf)
#}
#
#
## Compare to vaccine
#vaccineRatio <- data.frame()
#for (no in c(1:length(outcome))) {
#  oStr <- outcome[no]
#  agesComparison <- proportion_samples(model=serologyModels$model[[oStr]],
#           ageVec=c(20, 20),
#           meanAge=serologyModels$meanAge[[oStr]],
#           sdAge=serologyModels$sdAge[[oStr]])
#  ratiosTemp <- agesComparison$samples %>%
#    dplyr::select(., -age) %>%
#    dplyr::filter(., ageInd==1) %>%
#    dplyr::mutate(., propRatio=proportion/(1/4000))
#  tempDf <- data.frame(meanRatio=mean(ratiosTemp$propRatio),
#                       lowerRatio=quantile(ratiosTemp$propRatio, p=0.025),
#                       upperRatio=quantile(ratiosTemp$propRatio, p=0.975))
#  vaccineRatio <- rbind(vaccineRatio, tempDf)
#}
#
#
## prior functions
#gaussian <- function(x, mu, sigma){
#  f <- (1/sigma*sqrt(2*pi))*exp(-0.5*((x-mu)/sigma)^2)
#}
#exponential <- function(x, lambda){lambda*exp(-lambda*x)}
#outcome <- c("Hospitalized", "ICU", "Deaths")
#x1 <- seq(-1, 5, 0.02)
#prior_ageSlope <- data.frame(value=x1, dens=gaussian(x1, mu=2, sigma=1),
#                             .variable="ageSlope")
#x2 <- seq(0, 4.5, 0.02)
#prior_ageSlopeSigma <- data.frame(value=x2, dens=exponential(x2, lambda=0.5),
#                                  .variable="ageSlopeSigma")
#x3 <- seq(-8, 0, 0.02)
#prior_intercept <- data.frame(value=x3, dens=gaussian(x3, mu=-6, sigma=2),
#                              .variable="intercept")
#x4 <- seq(0, 4.5, 0.02)
#prior_interceptSigma <- data.frame(value=x4, dens=exponential(x4, lambda=0.5),
#                                   .variable="interceptSigma")
#priorDf <- rbind(prior_ageSlope, prior_ageSlopeSigma, prior_intercept,
#                 prior_interceptSigma)
#
#
#posteriorSerology <- data.frame()
#for (no in c(1:length(outcome))) {
#  oStr <- outcome[no]
#  posteriorTemp <- tidybayes::gather_draws(serologyModels$model[[oStr]], ageSlope,
#                                           ageSlopeSigma, intercept, interceptSigma) %>%
#    dplyr::mutate(., Outcome_type=oStr)
#  posteriorSerology <- rbind(posteriorSerology, posteriorTemp)
#}
#posteriorSerology$Outcome_type <- factor(posteriorSerology$Outcome_type,
#                                         levels=outcome)
#
#
#posteriorSerologyPlot <- posteriorSerology %>%
#  ggplot(., aes(x=.value, color=Outcome_type, fill=Outcome_type,
#                facet=.variable)) +
#  #geom_density(alpha=0.6) +
#  stat_density(geom="area", position="identity", alpha=0.6) +
#  stat_density(geom="line", position="identity") +
#  geom_line(data=priorDf, aes(x=value, y=dens), size=1.5, inherit.aes=FALSE) +
#  facet_grid(.~.variable, scales="free") +
#  theme_bw()
#
#serologyTrace <- dplyr::mutate(posteriorSerology, .chain=factor(.chain))  %>%
#  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
#  geom_line() +
#  facet_grid(.variable~Outcome_type, scales="free_y") +
#  theme_bw()
#
#ggsave("../data/plots/3_serology_param_posterior.png", posteriorSerologyPlot, 
#       width=30, height=12, units="cm", device="png")
#
#ggsave("../data/plots/3_serology_chains.png", serologyTrace, 
#       width=30, height=20, units="cm")
#
#
#
##
##serologySamplesDf$Outcome_type <- factor(serologySamplesDf$Outcome_type,
##                                       levels=outcome)
##serologySamplesPlot <- serologySamplesDf %>%
##  ggplot(., aes(x=meanAge, y=proportion*100, group=sample, facet=Outcome_type)) +
##  geom_line(alpha=sampleLineAlpha, size=sampleLineSize, color="#33ADFF") +
##  facet_grid(.~Outcome_type) +
##  scale_y_continuous(trans='log10', labels=scaleFun) +
##  geom_line(data=outcomeFitDf, aes(x=meanAge, y=outcomeProp*100),
##            color="black", linetype="solid", size=regLineSize, inherit.aes=FALSE) +
###  geom_ribbon(data=outcomeFitDf,
###              aes(x=meanAge,ymin=outcome_L*100, ymax=outcome_H*100),
###              alpha=0.2, colour=NA, show.legend=FALSE,
###              inherit.aes=FALSE) +
##  theme_bw() +
##  xlab("Age") +
##  ylab("% outcome")
##
##ggsave(serologySamplesPlotName, serologySamplesPlot,
##       width=29, height=10, units="cm")
##
#
#exportFit <- dplyr::mutate(outcomeFitDf, Percentage=signif(outcomeProp*100, digits=3),
#                           Percentage_L=signif(outcome_L*100, digits=3),
#                           Percentage_H=signif(outcome_H*100, digits=3)) %>%
#  dplyr::select(., -outcomeProp, -outcome_L, -outcome_H)
#write.csv(exportFit, file=serologyCsvName, row.names=FALSE)
#
#############################
#############################
####
#### Plot hospital lethality data fit
####
#############################
#############################
#
#lethalityData <- read.csv("../data/collected_data/hospitalized_patient_studies.csv",
#                        stringsAsFactors=FALSE) %>%
#  as_tibble(.)
#lethalityModels <- readRDS("../data/processed_data/4_hospital_lethality_fit.RDS")
#
## get fitted line
#ageVec <- seq(2.5, 90, 5)
#letType <- c("Hospitalized", "ICU")
#lethalityPosterior <- list()
#lethalityFitDf <- NULL
#lethalitySamplesDf <- tibble()
#for (no in c(1:length(letType))) {
#  oStr <- letType[no]
#  # extract model fit results
#  stdAgeVec <- (ageVec-lethalityModels$meanAge[[oStr]])/lethalityModels$sdAge[[oStr]]
#  lethalityPosterior[[oStr]] <- proportion_samples(model=lethalityModels$model[[oStr]],
#                                                  ageVec=stdAgeVec)
#  tempFitDf <- data.frame(meanAge=ageVec,
#                          outcomeProp=lethalityPosterior[[oStr]]$prop_mean,
#                          outcome_L=lethalityPosterior[[oStr]]$prop_L,
#                          outcome_H=lethalityPosterior[[oStr]]$prop_H,
#                          Type=oStr)
#  lethalityFitDf <- rbind(lethalityFitDf, tempFitDf)
#  tempSamplesDf <- lethalityPosterior[[oStr]]$samples
#  tempSamplesDf$meanAge <- rep(ageVec, max(tempSamplesDf$sample))
#  tempSamplesDf$Type <- oStr
#  lethalitySamplesDf <- as_tibble(rbind(tempSamplesDf, lethalitySamplesDf))
#}
#
#letType2 <- c("Hospital", "ICU")
#lethalityData$Type <- factor(lethalityData$Type, levels=letType, labels=letType2)
#lethalityData$Location <- factor(lethalityData$Location)
#lethalityFitDf$Type <- factor(lethalityFitDf$Type, levels=letType, labels=letType2)
#
#lethalityPlot <- lethalityData %>%
#  ggplot(., aes(x=meanAge, y=Deaths/Patients*100, color=Location, facet=Type)) +
#  geom_point(size=locPointSize) +
#  geom_line(alpha=locLineAlpha, size=locLineSize) +
#  facet_rep_wrap(Type~., nrow=2, strip.position="top") +
#  geom_line(data=lethalityFitDf, aes(y=outcomeProp*100),
#            color="black", linetype="solid", size=regLineSize) +
#  geom_ribbon(data=lethalityFitDf,
#              aes(x=meanAge,ymin=outcome_L*100, ymax=outcome_H*100),
#              alpha=ribbonAlpha*0.8, colour=NA, show.legend=FALSE,
#              inherit.aes=FALSE) +
#  scale_y_continuous(trans='log10', labels=scaleFun) +
#  theme_bw() +
#  theme(strip.background =element_rect(fill="white", color="white"),
#        strip.text = element_text(face="bold"),
#        legend.position = "top",
#        panel.border = element_blank(),
#        #panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        strip.text.y = element_text(angle=0),
#        axis.line = element_line(colour = "black"),
#        axis.line.x = element_line(size=0.5, linetype="solid"),
#        axis.line.y = element_line(size=0.5, linetype="solid"),
#        legend.title=element_blank()) +
#  theme(legend.title = element_blank()) +
#  theme(legend.position = "none") +
#  xlab("Age") +
#  ylab("Mortality (%)")
#
#ggsave("../data/plots/4_hospital_lethality_regression.png", lethalityPlot,
#       width=10, height=14, units="cm")
#
#
##lethalitySamplesDf$Type <- factor(lethalitySamplesDf$Type,
##                                       levels=letType)
##lethalitySamplesPlot <- lethalitySamplesDf %>%
##  ggplot(., aes(x=meanAge, y=proportion*100, group=sample, facet=Type)) +
##  geom_line(alpha=sampleLineAlpha, size=sampleLineSize, color="#33ADFF") +
##  scale_y_continuous(trans='log10', labels=scaleFun) +
##  geom_line(data=lethalityFitDf, aes(x=meanAge, y=outcomeProp*100),
##            color="black", linetype="solid", size=regLineSize, inherit.aes=FALSE) +
##  facet_grid(.~Type) +
##  geom_ribbon(data=outcomeFitDf,
##              aes(x=meanAge,ymin=outcome_L*100, ymax=outcome_H*100),
##              alpha=0.2, colour=NA, show.legend=FALSE,
##              inherit.aes=FALSE) +
##  theme_bw() +
##  xlab("Age") +
##  ylab("% outcome")
##
##ggsave("../data/plots/4_lethality_samples.png", lethalitySamplesPlot,
##       width=29, height=10, units="cm")
#
#
#############################
#############################
####
#### Plot outcome proportions estimated from literature
####
#############################
#############################
#
#outcomePropLit <- read.csv("../data/processed_data/5_literature_outcome_estimates.csv",
#                           stringsAsFactors=FALSE) %>%
#  as_tibble(.)
#
#outcomePropLit$Type <- factor(outcomePropLit$Type, levels=outcome,
#  labels=outcome2)
#outcomePropLit$Outcome_type <- outcomePropLit$Type
#
#literaturePropPlot <- outcomePropLit %>%
#  ggplot(., aes(x=meanAge, y=Proportion, color=Study)) +
#  geom_point(size=locPointSize) +
#  geom_errorbar(aes(x=meanAge, ymin=Proportion_L, ymax=Proportion_H)) +
#  geom_line(alpha=locLineAlpha, size=locLineSize) +
#  facet_rep_grid(.~Outcome_type, repeat.tick.labels="left") +
#  geom_line(data=outcomeFitDf, aes(y=outcomeProp*100),
#            color="black", linetype="solid", size=regLineSize) +
#  geom_ribbon(data=outcomeFitDf,
#              aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
#              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
#              inherit.aes=FALSE) +
#  geom_segment(aes(x=ageMeanMyo, xend=ageMeanMyo, y=myocarditisProp[1]*100,
#                   yend=myocarditisProp[2]*100), color="black", size=1) +
#  geom_segment(aes(x=ageMeanMyo+2, xend=myoTextX-20, y=100/4000, yend=myoTextY*1.5),
#               arrow=arrow(length=unit(0.2, "cm")), color="black",
#               inherit.aes=FALSE, size=0.3) +
#  geom_text(aes(x=myoTextX-5, y=myoTextY*0.5), label="Israel vaccine\nmyocarditis rate",
#            color="black", inherit.aes=FALSE, size=3.1) +
#  scale_color_brewer(palette="Set1") +
#  scale_y_continuous(trans='log10', labels=scaleFun,
#                     breaks=10^c(-3, -2, -1, 0, 1, 2),
#                     limits=10^c(-3.5, 2)) +
#  theme_bw() +
#  theme(strip.background =element_rect(fill="white", color="white"),
#        strip.text=element_text(face="bold"),
#        legend.position="top",
#        panel.border=element_blank(),
#        #panel.grid.major=element_blank(),
#        panel.grid.minor=element_blank(),
#        axis.line=element_line(colour="black"),
#        axis.line.x=element_line(size=0.5, linetype="solid"),
#        axis.line.y=element_line(size=0.5, linetype="solid"),
#        legend.title=element_blank()) +
#  xlab("Age") +
#  ylab("% Infected with outcome")
#
#
#ggsave("../data/plots/5_literature_outcome_estimates.png", literaturePropPlot,
#       width=18, height=9, units="cm")
#
#
#figure2 <- ggarrange(plotlist=list(lethalityPlot, literaturePropPlot),
#          ncol=2, widths=c(0.25, 0.74), labels=c("(A)", "(B)"))
#
#ggsave("../data/plots/figure2.png", figure2, width=20, height=10, units="cm")
#
#############################
#############################
####
#### Plot models params
####
#############################
#############################
#
#
#############
## plot hospital lethality Priors and trace
#############
#posteriorLethality <- data.frame()
#lethOutcome <- c("Hospitalized", "ICU")
#for (no in c(1:length(lethOutcome))) {
#  oStr <- lethOutcome[no]
#  posteriorTemp <- tidybayes::gather_draws(lethalityModels$model[[oStr]], ageSlope,
#                                           ageSlopeSigma, intercept, interceptSigma) %>%
#    dplyr::mutate(., Outcome_type=oStr)
#  posteriorLethality <- rbind(posteriorLethality, posteriorTemp)
#}
#posteriorLethality$Outcome_type <- factor(posteriorLethality$Outcome_type,
#                                         levels=outcome)
#
#posteriorLethalityPlot <- posteriorLethality %>%
#  ggplot(., aes(x=.value, color=Outcome_type, fill=Outcome_type,
#                facet=.variable)) +
#  stat_density(geom="area", position="identity", alpha=0.6) +
#  stat_density(geom="line", position="identity") +
#  geom_line(data=priorDf, aes(x=value, y=dens), size=1.5, inherit.aes=FALSE) +
#  facet_grid(.~.variable, scales="free") +
#  theme_bw()
#
#lethalityTrace <- dplyr::mutate(posteriorLethality, .chain=factor(.chain))  %>%
#  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
#  geom_line() +
#  facet_grid(.variable~Outcome_type, scales="free_y") +
#  theme_bw()
#
#ggsave("../data/plots/4_lethality_param_posterior.png", posteriorLethalityPlot, 
#       width=30, height=12, units="cm", device="png")
#
#ggsave("../data/plots/4_lethality_chains.png", lethalityTrace, 
#       width=30, height=20, units="cm")
#
#
