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

regLineSize=1
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
### Plot serology data fit
###
############################
############################

# Corrected data
countryData <- read.csv("../data/collected_data/locations_serology_data_corrected.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.)
serologyModels <- readRDS("../data/processed_data/3_serology_fits_corrected.RDS")
#serologyModels <- readRDS("../data/processed_data/3_serology_fits_corrected_Phi.RDS")
serologyPlotName <- "../data/plots/figure1.png"
serologySamplesPlotName <- "../data/plots/3_serology_samples_corrected.png"
serologyCsvName <- "../data/processed_data/3_serology_fits_corrected.csv"

# predefine some variables
outcome <- c("Hospitalized", "ICU", "Deaths")
outcome2 <- c("Severe disease", "Critical disease", "Fatal disease")

# pivot longer the serology data
longCountryData <- tidyr::pivot_longer(data=countryData, 
                                       cols=all_of(outcome),
                                       names_to="Outcome_type",
                                       values_to="Outcomes") %>%
  dplyr::filter(., !is.na(Outcomes)) %>%
  dplyr::mutate(., Location=factor(Location))
longCountryData$Outcome_type <- factor(longCountryData$Outcome_type,
                                       levels=outcome, labels=outcome2)

levsType <- levels(factor(longCountryData$Type))
labelsType <- c("Representative seroprevalence", "Convenience seroprevalence",
                "Comprehensive testing")
longCountryData$Type <- factor(longCountryData$Type,
                                       levels=levsType,
                                       labels=labelsType)

# extract model posteriors and put into data frame
ageVec <- seq(2.5, 90, 5)
serologyPosterior <- list()
outcomeFitDf <- NULL
serologySamplesDf <- tibble()
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
outcomeFitDf$Outcome_type <- factor(outcomeFitDf$Outcome_type,
                                       levels=outcome, labels=outcome2)


# https://www.sciencemag.org/news/2021/06/israel-reports-link-between-rare-cases-heart-inflammation-and-covid-19-vaccination
myocarditisProp <- c(1/6000, 1/3000)
ageRangeMyo <- c(16, 24)
ageMeanMyo <- mean(ageRangeMyo)
myoTextX <- 68
myoTextY <- 0.005

serologyPlot <- longCountryData %>%
  ggplot(., aes(x=meanAge, y=Outcomes*100/(Population*Prevalence/100),
                color=Type, group=as.character(Location), facet=Outcome_type)) +
  geom_point(size=locPointSize, alpha=dataAlpha) +
  geom_line(alpha=locLineAlpha, size=locLineSize, alpha=dataAlpha) +
  facet_rep_grid(.~Outcome_type, repeat.tick.labels="left") +
  geom_line(data=outcomeFitDf, aes(x=meanAge, y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize,
            inherit.aes=FALSE) +
  geom_ribbon(data=outcomeFitDf,
              aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  geom_segment(aes(x=ageMeanMyo, xend=ageMeanMyo, y=myocarditisProp[1]*100,
                   yend=myocarditisProp[2]*100), color="black", size=1,
               inherit.aes=FALSE) +
  geom_segment(aes(x=ageMeanMyo+2, xend=myoTextX-20, y=100/4000, yend=myoTextY*1.5),
               arrow=arrow(length=unit(0.2, "cm")), color="black",
               inherit.aes=FALSE, size=0.3) +
  geom_text(aes(x=myoTextX, y=myoTextY), label="Israel vaccine\nmyocarditis rate", color="black",
            inherit.aes=FALSE, size=3.1) +
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
  xlab("Age") +
  ylab("% Infected with outcome")

ggsave(serologyPlotName, serologyPlot, width=18, height=9.5, units="cm")


##############
# Extract serology data summary statistics
##############

# Get param summary statistics
mainParams <- c("ageSlope", "ageSlopeSigma", "intercept", "interceptSigma")
modelSummary <- list()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  modelSummary[[oStr]] <- summary(serologyModels$model[[oStr]],
                          pars=mainParams)
}

# Compare ages
ratiosAges <- data.frame()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  agesComparison <- proportion_samples(model=serologyModels$model[[oStr]],
           ageVec=c(22.5, 72.5),
           meanAge=serologyModels$meanAge[[oStr]],
           sdAge=serologyModels$sdAge[[oStr]])
  ratiosTemp <- agesComparison$samples %>%
    dplyr::select(., -age) %>%
    tidyr::pivot_wider(., names_from=ageInd, values_from=proportion,
                       names_prefix="age") %>%
    dplyr::mutate(., ratio=age2/age1)
  tempDf <- data.frame(meanRatio=mean(ratiosTemp$ratio),
                       lowerRatio=quantile(ratiosTemp$ratio, p=0.025),
                       upperRatio=quantile(ratiosTemp$ratio, p=0.975))
  ratiosAges <- rbind(ratiosAges, tempDf)
}


# Compare to vaccine
vaccineRatio <- data.frame()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  agesComparison <- proportion_samples(model=serologyModels$model[[oStr]],
           ageVec=c(20, 20),
           meanAge=serologyModels$meanAge[[oStr]],
           sdAge=serologyModels$sdAge[[oStr]])
  ratiosTemp <- agesComparison$samples %>%
    dplyr::select(., -age) %>%
    dplyr::filter(., ageInd==1) %>%
    dplyr::mutate(., propRatio=proportion/(1/4000))
  tempDf <- data.frame(meanRatio=mean(ratiosTemp$propRatio),
                       lowerRatio=quantile(ratiosTemp$propRatio, p=0.025),
                       upperRatio=quantile(ratiosTemp$propRatio, p=0.975))
  vaccineRatio <- rbind(vaccineRatio, tempDf)
}


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
  geom_point(size=locPointSize) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  facet_rep_wrap(Type~., nrow=2, strip.position="top") +
  geom_line(data=lethalityFitDf, aes(y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize) +
  geom_ribbon(data=lethalityFitDf,
              aes(x=meanAge,ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha*0.8, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white", color="white"),
        strip.text = element_text(face="bold"),
        legend.position = "top",
        panel.border = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle=0),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid"),
        legend.title=element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  xlab("Age") +
  ylab("Mortality (%)")

ggsave("../data/plots/4_hospital_lethality_regression.png", lethalityPlot,
       width=10, height=14, units="cm")


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

brazeauIFRind <- with(outcomePropLit, which(Study=="Brazeau" &
                                            Outcome_type=="Fatal disease"))
outcomePropLit$Proportion_L[brazeauIFRind] <- outcomePropLit$Proportion[brazeauIFRind]
outcomePropLit$Proportion_H[brazeauIFRind] <- outcomePropLit$Proportion[brazeauIFRind]


literaturePropPlot <- outcomePropLit %>%
  ggplot(., aes(x=meanAge, y=Proportion, color=Study)) +
  geom_point(size=locPointSize) +
  geom_errorbar(aes(x=meanAge, ymin=Proportion_L, ymax=Proportion_H), size=0.5,
                alpha=locLineAlpha*0.7) +
  #geom_line(alpha=locLineAlpha, size=locLineSize) +
  facet_rep_grid(.~Outcome_type, repeat.tick.labels="left") +
  geom_line(data=outcomeFitDf, aes(y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize) +
  geom_ribbon(data=outcomeFitDf,
              aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  geom_segment(aes(x=ageMeanMyo, xend=ageMeanMyo, y=myocarditisProp[1]*100,
                   yend=myocarditisProp[2]*100), color="black", size=1) +
  geom_segment(aes(x=ageMeanMyo+2, xend=myoTextX-20, y=100/4000, yend=myoTextY*1.5),
               arrow=arrow(length=unit(0.2, "cm")), color="black",
               inherit.aes=FALSE, size=0.3) +
  geom_text(aes(x=myoTextX-5, y=myoTextY*0.5), label="Israel vaccine\nmyocarditis rate",
            color="black", inherit.aes=FALSE, size=3.1) +
  scale_color_brewer(palette="Set1") +
  scale_y_continuous(trans='log10', labels=scaleFun,
                     breaks=10^c(-3, -2, -1, 0, 1, 2),
                     limits=10^c(-3.5, 2)) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white", color="white"),
        strip.text=element_text(face="bold"),
        legend.position="top",
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid"),
        legend.title=element_blank()) +
  xlab("Age") +
  ylab("% Infected with outcome")


#ggsave("../data/plots/5_literature_outcome_estimates.png", literaturePropPlot,
#       width=18, height=9, units="cm")


figure2 <- ggarrange(plotlist=list(lethalityPlot, literaturePropPlot),
          ncol=2, widths=c(0.25, 0.74), labels=c("(A)", "(B)"))

ggsave("../data/plots/figure2.png", figure2, width=20, height=10, units="cm")


