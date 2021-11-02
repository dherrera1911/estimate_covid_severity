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

# predefine some variables
outcome <- c("Severe", "Critical", "Deaths")
outcome2 <- c("Severe disease", "Critical disease", "Fatal disease")
labelsType <- c("Representative seroprevalence", "Convenience seroprevalence",
                "Comprehensive testing")


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
  dplyr::arrange(., Location, Age, Severe, Critical)

serologyModelsUncorrected <- readRDS("../data/processed_data/6_serology_fits.RDS")

# Corrected data
countryDataCorrected <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., processing="Corrected") %>%
  dplyr::arrange(., Location, Age, Severe, Critical)


serologyModelsCorrected <- readRDS("../data/processed_data/6_serology_fits_corrected.RDS")

countryDataUncorrected$oohDeaths_source <- countryDataCorrected$oohDeaths_source
countryDataUncorrected$ooiDeaths_source <- countryDataCorrected$ooiDeaths_source

countryData <- dplyr::select(countryDataCorrected,
                             all_of(names(countryDataUncorrected))) %>%
  rbind(., countryDataUncorrected)

serologyModels <- list(Corrected=serologyModelsCorrected,
                    Uncorrected=serologyModelsUncorrected)


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
    stdAgeVec <- (ageVec-serologyModels[[n]]$meanAge)/serologyModels[[n]]$sdAge
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

levsNames <- c("Severe", "Critical", "Deaths")
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
### Plot the proportion of severe/critical that are ooh/ooi deaths
###
############################
############################

critCorrPlot <- dplyr::mutate(countryDataCorrected,
                            propCorrCrit=ooiDeaths/Critical*100) %>%
  dplyr::filter(., !is.na(propCorrCrit)) %>%
  droplevels(.) %>%
  ggplot(., aes(x=meanAge, y=propCorrCrit,
                color=factor(ooiDeaths_source), group=as.character(Location))) +
  geom_point(size=locPointSize, alpha=dataAlpha) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
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
  ylab(expression(over("Out of ICU deaths", "Critical cases")%*%100)) +
    NULL

severeCorrPlot <- dplyr::mutate(countryDataCorrected,
                            propCorrSevere=oohDeaths/Severe*100) %>%
  dplyr::filter(., !is.na(propCorrSevere)) %>%
  droplevels(.) %>%
  ggplot(., aes(x=meanAge, y=propCorrSevere,
                color=factor(oohDeaths_source), group=as.character(Location))) +
  geom_point(size=locPointSize, alpha=dataAlpha) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
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
  ylab(expression(over("Out of hospital deaths", "Severe cases")%*%100)) +
    NULL

correctionPropPlot <- ggpubr::ggarrange(plotlist=list(severeCorrPlot, critCorrPlot),
                                        ncol=2, labels=c("A)", "B)"), common.legend=TRUE)

ggsave("../data/plots/figureS2.png", correctionPropPlot,
       width=20, height=10, units="cm")


############################
############################
###
### Plot fits to younger ages
###
############################
############################

serologyModelsU50 <- readRDS("../data/processed_data/6_serology_fits_u50.RDS")
serologyModelsU40 <- readRDS("../data/processed_data/6_serology_fits_u40.RDS")

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
    stdAgeVec <- (ageVec-serologyModels2[[n]]$meanAge)/serologyModels2[[n]]$sdAge
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

levsNames <- c("Severe", "Critical", "Deaths")
labelsNames <- c("Severe disease", "Critical disease", "Deaths")
outcomeFit2Df$Outcome_type <- factor(outcomeFit2Df$Outcome_type,
                                       levels=levsNames, labels=outcome2)

allModels <- rbind(outcomeFitDf, outcomeFit2Df)
u40ind <- allModels$processing == "Uncorrected_U40"
u50ind <- allModels$processing == "Uncorrected_U50"
allModels$processing[u40ind] <- "Under 40 (uncorrected)"
allModels$processing[u50ind] <- "Under 50 (uncorrected)"

colorMal <- c("black", brewer.pal(3, "Set1"))

differentModelsPlot <- allModels %>%
  dplyr::filter(., processing!="Uncorrected") %>%
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

rm(serologyModelsU50, serologyModelsU40, serologyModels2)

############################
############################
###
### Plot fits without Testing or without convenience serosampling
###
############################
############################

serologyModels_noTest <- readRDS("../data/processed_data/6_serology_fits_corrected_noTest_corrected.RDS")
serologyModels_noConv <- readRDS("../data/processed_data/6_serology_fits_corrected_noConvenience_corrected.RDS")

serologyModels3 <- list(No_test=serologyModels_noTest,
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
    stdAgeVec <- (ageVec-serologyModels3[[n]]$meanAge)/serologyModels3[[n]]$sdAge
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

levsNames <- c("Severe", "Critical", "Deaths")
labelsNames <- c("Severe disease", "Critical disease", "Deaths")
outcomeFit3Df$Outcome_type <- factor(outcomeFit3Df$Outcome_type,
                                       levels=levsNames, labels=outcome2)
outcomeFit3Df$processing[outcomeFit3Df$processing=="No_convenience"] <- "No convenience samples"
outcomeFit3Df$processing[outcomeFit3Df$processing=="No_test"] <- "No comprehensive testing"

allModels <- dplyr::filter(outcomeFitDf, processing=="Corrected") %>%
  dplyr::mutate(., processing="Full dataset") %>%
  rbind(., outcomeFit3Df)

colorMal <- c("black", brewer.pal(3, "Set1"))

differentModelsPlot <- allModels %>%
  #dplyr::filter(., Outcome_type!="Fatal disease") %>%
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


ggsave("../data/plots/figureS4.png", differentModelsPlot,
       width=16, height=10, units="cm")


rm(serologyModels_noConv, serologyModels_noTest, serologyModels3)

############################
############################
###
### Plot fits without fastest changing locations
###
############################
############################
serologyModels_slow <- readRDS("../data/processed_data/6_serology_fits_corrected_lessChanging_corrected.RDS")

# extract model posteriors and put into data frame
ageVec <- seq(2.5, 90, 5)
serologyPosterior <- list()
serologySamplesDf <- tibble()
outcomeFit4Df <- NULL
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  stdAgeVec <- (ageVec-serologyModels_slow$meanAge)/serologyModels_slow$sdAge
  serologyPosterior[[oStr]] <- proportion_samples(model=serologyModels_slow$model[[oStr]],
                                                  ageVec=stdAgeVec)
  tempFitDf <- data.frame(meanAge=ageVec,
                          outcomeProp=serologyPosterior[[oStr]]$prop_mean,
                          outcome_L=serologyPosterior[[oStr]]$prop_L,
                          outcome_H=serologyPosterior[[oStr]]$prop_H,
                          Outcome_type=oStr,
                          processing="Without fastest growing")
  outcomeFit4Df <- rbind(outcomeFit4Df, tempFitDf)
}
outcomeFit4Df$Outcome_type <- factor(outcomeFit4Df$Outcome_type,
                                       levels=levsNames, labels=outcome2)

allModels <- dplyr::filter(outcomeFitDf, processing=="Corrected") %>%
  dplyr::mutate(., processing="Full dataset") %>%
  rbind(., outcomeFit4Df)

colorMal <- c("black", brewer.pal(3, "Set1"))

differentModelsPlot <- allModels %>%
  #dplyr::filter(., Outcome_type!="Fatal disease") %>%
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

ggsave("../data/plots/figureS6.png", differentModelsPlot,
       width=16, height=10, units="cm")


############################
############################
###
### Plot fits without most seroconverted changing locations
###
############################
############################
serologyModels_sero <- readRDS("../data/processed_data/6_serology_fits_corrected_lessPast_corrected.RDS")

# extract model posteriors and put into data frame
ageVec <- seq(2.5, 90, 5)
serologyPosterior <- list()
serologySamplesDf <- tibble()
outcomeFit5Df <- NULL
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  stdAgeVec <- (ageVec-serologyModels_sero$meanAge)/serologyModels_sero$sdAge
  serologyPosterior[[oStr]] <- proportion_samples(model=serologyModels_sero$model[[oStr]],
                                                  ageVec=stdAgeVec)
  tempFitDf <- data.frame(meanAge=ageVec,
                          outcomeProp=serologyPosterior[[oStr]]$prop_mean,
                          outcome_L=serologyPosterior[[oStr]]$prop_L,
                          outcome_H=serologyPosterior[[oStr]]$prop_H,
                          Outcome_type=oStr,
                          processing="Without most seroreverted")
  outcomeFit5Df <- rbind(outcomeFit5Df, tempFitDf)
}
outcomeFit5Df$Outcome_type <- factor(outcomeFit5Df$Outcome_type,
                                       levels=levsNames, labels=outcome2)

allModels <- dplyr::filter(outcomeFitDf, processing=="Corrected") %>%
  dplyr::mutate(., processing="Full dataset") %>%
  rbind(., outcomeFit5Df)

colorMal <- c("black", brewer.pal(3, "Set1"))

differentModelsPlot <- allModels %>%
  #dplyr::filter(., Outcome_type!="Fatal disease") %>%
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

ggsave("../data/plots/figureS7.png", differentModelsPlot,
       width=16, height=10, units="cm")

#seroCont <- dplyr::filter(allModels, processing=="Without most seroconverted" &
#                          Outcome_type=="Fatal disease")
#original <- dplyr::filter(allModels, processing=="Full dataset",
#                          Outcome_type=="Fatal disease")
#plot(original$outcomeProp/seroCont$outcomeProp)

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
  dplyr::filter(., daysSinceData>=-15 & daysSinceData<=15) %>%
  ggplot(., aes(x=daysSinceData, y=deathProp, color=Location)) +
  geom_line(size=locLineSize) +
  theme_bw() +
  geom_segment(aes(x=-15, y=110, xend=15, yend=110), linetype=2,
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

ggsave("../data/plots/figureS5.png", deathsDynPlot,
       width=18, height=15, units="cm")

lastDayDyn <- dplyr::filter(deathDynamics2, daysSinceData==15) %>%
  dplyr::arrange(., deathProp)
nLoc <- nrow(lastDayDyn)

# see proportion of deaths acquired in las 45 days
previousDay50 <- dplyr::mutate(deathDynamics2, newerDeaths=100-deathProp) %>%
  dplyr::filter(., daysSinceData==-45 & newerDeaths>50)

recentDeathsDf <- NULL
for (l in unique(previousDay50$Location)) {
  locDf <- dplyr::filter(previousDay50, Location==l)
  maxDate <- which(locDf$date == max(locDf$date))
  recentDeathsDf <- rbind(recentDeathsDf, locDf[maxDate,])
}

# Literature plots
literatureRates <- read.csv("../data/collected_data/literature_rates_estimations.csv",
                            stringsAsFactors=FALSE)

ihrLit <- dplyr::filter(literatureRates, Type=="IHR")

literatureIHR <- ihrLit %>%
  dplyr::filter(., meanAge < 60) %>%
  ggplot(., aes(x=meanAge, y=Proportion, color=Study)) +
  geom_point(size=locPointSize, alpha=dataAlpha)


##########
# Point estimates of children rates
##########
childrenFits <- readRDS("../data/processed_data/8_serology_fits.RDS")
childrenData <- dplyr::filter(countryDataUncorrected, Age=="0-9")

childrenPosterior <- list()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  childrenPosterior[[oStr]] <- proportion_samples(model=childrenFits$model[[oStr]],
                                                  ageVec=0, slopeName=NA)
}

