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
library(cowplot)
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
serologyModels <- readRDS("../data/processed_data/6_serology_fits_corrected.RDS")
serologyPlotName <- "../data/plots/figure1.png"
serologySamplesPlotName <- "../data/plots/3_serology_samples_corrected.png"
serologyCsvName <- "../data/processed_data/3_serology_fits_corrected.csv"

# predefine some variables
outcome <- c("Severe", "Critical", "Deaths")
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
ageVec <- seq(5, 90, 10)
serologyPosterior <- list()
outcomeFitDf <- NULL
serologySamplesDf <- tibble()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  stdAgeVec <- (ageVec-serologyModels$meanAge)/serologyModels$sdAge
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


serologyPlot <- longCountryData %>%
  ggplot(., aes(x=meanAge, y=Outcomes*100/(Population*Prevalence/100),
                color=Type, group=as.character(Location), facet=Outcome_type)) +
  geom_point(size=locPointSize, alpha=dataAlpha) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  facet_rep_grid(.~Outcome_type, repeat.tick.labels="left") +
  geom_line(data=outcomeFitDf, aes(x=meanAge, y=outcomeProp*100),
            color="black", linetype="solid", size=regLineSize,
            inherit.aes=FALSE) +
  geom_ribbon(data=outcomeFitDf,
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
  xlab("Age") +
  ylab("% Infected with outcome")

ggsave(serologyPlotName, serologyPlot, width=18, height=9.5, units="cm")


fitPlotLin <- outcomeFitDf %>%
#  dplyr::mutate(., old = meanAge>=40) %>%
  ggplot(., aes(x=meanAge, y=outcomeProp*100,
                color=Outcome_type, fill=Outcome_type)) +
  geom_line(size=regLineSize*0.8) +
  geom_ribbon(aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, color=NA, show.legend=FALSE) +
#  facet_wrap(.~old, scales="free") +
  theme_bw() +
  scale_x_continuous(expand=c(0,0)) +
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

fitPlotExp <- outcomeFitDf %>%
  ggplot(., aes(x=meanAge, y=outcomeProp*100,
                color=Outcome_type, fill=Outcome_type)) +
  geom_line(size=regLineSize*0.8) +
  geom_ribbon(aes(x=meanAge, ymin=outcome_L*100, ymax=outcome_H*100),
              alpha=ribbonAlpha, color=NA, show.legend=FALSE) +
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


fitsPlot <- ggpubr::ggarrange(plotlist=list(fitPlotExp, fitPlotLinSplit),
                  ncol=2, common.legend=TRUE)

figure1 <- ggpubr::ggarrange(plotlist=list(serologyPlot, fitsPlot),
                           nrow=2, labels=c("(A)", "(B)"),
                           heights=c(1.1, 0.9))

ggsave("../data/plots/figure1_new.png", figure1,
       width=18, height=17, units="cm")


#posteriorTraces <- NULL
#for (oStr in outcome) {
#  posteriorTemp <- tidybayes::gather_draws(serologyModels$model[[oStr]],
#                                           ageSlope, ageSlopeSigma,
#                                           intercept, interceptSigma) %>%
#    dplyr::mutate(., type=oStr)
#  posteriorTraces <- rbind(posteriorTraces, posteriorTemp)
#}
#
#serologyTrace <- dplyr::mutate(posteriorTraces, .chain=factor(.chain))  %>%
#  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
#  geom_line() +
#  facet_grid(type~.variable, scales="free_y") +
#  theme_bw()
#
#pairs(serologyModels$model[["Deaths"]], pars=c("ageSlope", "ageSlopeSigma", "intercept",
#                                "interceptSigma"))
#pairs(serologyModels$model[["Deaths"]], pars=c("locationSlope"))

##############
# Extract serology data summary statistics
##############

# Get param summary statistics
mainParams <- c("ageSlope", "ageSlopeSigma", "intercept", "interceptSigma")
modelSummary <- NULL
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  modelParams <- extract_model_params_norm(serologyModels$model[[oStr]],
    xCenter=serologyModels$meanAge, xSd=serologyModels$sdAge)
  modelParams$outcome <- oStr
  modelSummary <- rbind(modelSummary, modelParams)
}

allPars <- summary(serologyModels$model[["Deaths"]], pars=mainParams)[["summary"]]
allPars <- signif(allPars[, c("mean", "2.5%", "97.5%")], 3)

# Compare ages
ratiosAges <- data.frame()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  agesComparison <- proportion_samples(model=serologyModels$model[[oStr]],
           ageVec=c(22.5, 72.5),
           meanAge=serologyModels$meanAge,
           sdAge=serologyModels$sdAge)
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
           ageVec=c(23),
           meanAge=serologyModels$meanAge,
           sdAge=serologyModels$sdAge)
  ratiosTemp <- agesComparison$samples %>%
    dplyr::select(., -age) %>%
    dplyr::filter(., ageInd==1) %>%
    dplyr::mutate(., propRatio=proportion/(1/80000))
  tempDf <- data.frame(meanRatio=mean(ratiosTemp$propRatio),
                       lowerRatio=quantile(ratiosTemp$propRatio, p=0.025),
                       upperRatio=quantile(ratiosTemp$propRatio, p=0.975))
  vaccineRatio <- rbind(vaccineRatio, tempDf)
}


exportFit <- dplyr::mutate(outcomeFitDf, Percentage=outcomeProp*100,
                           Percentage_L=outcome_L*100,
                           Percentage_H=outcome_H*100) %>%
  dplyr::select(., -outcomeProp, -outcome_L, -outcome_H)

write.csv(exportFit, file=serologyCsvName, row.names=FALSE)

estimateString <- as.character(signif(exportFit$Percentage, 3))
lowerString <- as.character(signif(exportFit$Percentage_L, 3))
upperString <- as.character(signif(exportFit$Percentage_H, 3))
fullString <- paste(estimateString, " (", lowerString, "-",
  upperString, ")", sep="")
tidyFitDf <- data.frame(meanAge=exportFit$meanAge,
                        outcome=exportFit$Outcome_type,
                        interval=fullString)
write.csv(tidyFitDf,
          file="../data/processed_data/3_serology_fits_corrected_tidy.csv",
          row.names=FALSE)

############################
############################
###
### Plot hospital lethality data fit
###
############################
############################

mortalityData <- read.csv("../data/collected_data/hospitalized_patient_studies.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.)
mortalityModels <- readRDS("../data/processed_data/3_hospital_mortality_fit.RDS")

# get fitted line
ageVec <- seq(2.5, 90, 5)
letType <- c("Hospitalized", "ICU")
mortalityPosterior <- list()
mortalityFitDf <- NULL
mortalitySamplesDf <- tibble()
for (no in c(1:length(letType))) {
  oStr <- letType[no]
  # extract model fit results
  stdAgeVec <- (ageVec-mortalityModels$meanAge[[oStr]])/mortalityModels$sdAge[[oStr]]
  mortalityPosterior[[oStr]] <- proportion_samples(model=mortalityModels$model[[oStr]],
                                                  ageVec=stdAgeVec)
  tempFitDf <- data.frame(meanAge=ageVec,
                          outcomeProp=mortalityPosterior[[oStr]]$prop_mean,
                          outcome_L=mortalityPosterior[[oStr]]$prop_L,
                          outcome_H=mortalityPosterior[[oStr]]$prop_H,
                          Type=oStr)
  mortalityFitDf <- rbind(mortalityFitDf, tempFitDf)
  tempSamplesDf <- mortalityPosterior[[oStr]]$samples
  tempSamplesDf$meanAge <- rep(ageVec, max(tempSamplesDf$sample))
  tempSamplesDf$Type <- oStr
  mortalitySamplesDf <- as_tibble(rbind(tempSamplesDf, mortalitySamplesDf))
}

letType2 <- c("Hospital", "ICU")
mortalityData$Type <- factor(mortalityData$Type, levels=letType, labels=letType2)
mortalityData$Location <- factor(mortalityData$Location)
mortalityFitDf$Type <- factor(mortalityFitDf$Type, levels=letType, labels=letType2)
mortalityFitDf <- dplyr::mutate(mortalityFitDf, outcomeProp=outcomeProp*100,
                                outcome_L=outcome_L*100, outcome_H=outcome_H*100)

mortalityPlot <- mortalityData %>%
  ggplot(., aes(x=meanAge, y=Deaths/Patients*100, color=Location, facet=Type)) +
  geom_point(size=locPointSize) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  facet_rep_wrap(Type~., nrow=2, strip.position="top") +
  geom_line(data=mortalityFitDf, aes(y=outcomeProp),
            color="black", linetype="solid", size=regLineSize) +
  geom_ribbon(data=mortalityFitDf,
              aes(x=meanAge,ymin=outcome_L, ymax=outcome_H),
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

ggsave("../data/plots/3_hospital_mortality_regression.png", mortalityPlot,
       width=10, height=14, units="cm")


### Extract and export in table format the estimates of mortality data
write.csv(mortalityFitDf, file="../data/processed_data/3_hospital_mortality_fit.csv",
          row.names=FALSE)

estimateString <- as.character(signif(mortalityFitDf$outcomeProp, 2))
lowerString <- as.character(signif(mortalityFitDf$outcome_L, 2))
upperString <- as.character(signif(mortalityFitDf$outcome_H, 2))
fullString <- paste(estimateString, " (", lowerString, "-",
  upperString, ")", sep="")
tidyFitDf <- data.frame(meanAge=mortalityFitDf$meanAge,
                        outcome=mortalityFitDf$Type,
                        interval=fullString)

write.csv(tidyFitDf,
          file="../data/processed_data/3_hospital_mortality_fit_tidy.csv",
          row.names=FALSE)



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

driscollInd <- outcomePropLit$Study == "Driscoll"
outcomePropLit$Study[driscollInd] <- "O'Driscoll"

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
#  geom_segment(aes(x=ageMeanMyo, xend=ageMeanMyo, y=myocarditisProp[1]*100,
#                   yend=myocarditisProp[2]*100), color="black", size=1) +
#  geom_segment(aes(x=ageMeanMyo+2, xend=myoTextX-20, y=100/4000, yend=myoTextY*1.5),
#               arrow=arrow(length=unit(0.2, "cm")), color="black",
#               inherit.aes=FALSE, size=0.3) +
#  geom_text(aes(x=myoTextX-5, y=myoTextY*0.5), label="Israel vaccine\nmyocarditis rate",
#            color="black", inherit.aes=FALSE, size=3.1) +
  scale_color_brewer(palette="Set1") +
  scale_y_continuous(trans='log10', labels=scaleFun,
                     breaks=10^c(-3, -2, -1, 0, 1, 2),
                     limits=10^c(-3.9, 2)) +
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


figure2 <- ggarrange(plotlist=list(mortalityPlot, literaturePropPlot),
          ncol=2, widths=c(0.25, 0.74), labels=c("(A)", "(B)"))

ggsave("../data/plots/figure2.png", figure2, width=20, height=10, units="cm")


# Get param summary statistics
outcome2 <- c("Hospitalized", "ICU")
mainParams <- c("ageSlope", "ageSlopeSigma", "intercept", "interceptSigma")
modelSummary <- NULL
for (no in c(1:length(outcome2))) {
  oStr <- outcome2[no]
  modelParams <- extract_model_params_norm(mortalityModels$model[[oStr]],
    xCenter=mortalityModels$meanAge[[oStr]], xSd=mortalityModels$sdAge[[oStr]])
  modelParams$outcome <- oStr
  modelSummary <- rbind(modelSummary, modelParams)
}

allPars <- summary(mortalityModels$model[["ICU"]], pars=mainParams)[["summary"]]
allPars <- signif(allPars[, c("mean", "2.5%", "97.5%")], 3)

