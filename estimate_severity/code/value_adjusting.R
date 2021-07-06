library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(wpp2019)
source("./functions_auxiliary.R")
data(popM)
data(popF)

regLineSize=1.2
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1

# function used for plotting log axis
scaleFun <- function(x) sprintf("%1g", x)


###############
# Load estimated lethality fit
###############
serologyModel <- readRDS("../data/processed_data/3_serology_fit_Deaths.RDS")

##################################
# Load in-hospital lethality data
##################################
lethalityData <- read.csv("../data/collected_data/hospitalized_patient_studies.csv") %>%
  as_tibble(.)
lethalityModels <- readRDS("../data/processed_data/4_hospital_lethality_fit.RDS")

ageVec <- seq(0, 90, 2.5)
outcome <- c("Hospitalized", "ICU")
lethalityPosterior <- data.frame()
for (no in c(1:length(outcome))) {
  oStr <- outcome[no]
  stdAgeVec <- (ageVec-lethalityModels$meanAge[[oStr]])/lethalityModels$sdAge[[oStr]]
  posteriorTemp <- proportion_samples(model=lethalityModels$model[[oStr]], ageVec=stdAgeVec)
  lethalityPosterior <- rbind(lethalityPosterior,
                              data.frame(meanAge=ageVec,
                                           prop_mean=posteriorTemp$prop_mean,
                                           prop_L=posteriorTemp$prop_L,
                                           prop_H=posteriorTemp$prop_H,
                                           Outcome_type=oStr))
}

######################################
# Load country data and separate into dataframes for hospitalization and ICU
######################################
countryData <- read.csv("../data/collected_data/locations_serology_data.csv")

countryDataMod <- countryData

#################################################
#################################################
# Tidy Data
#################################################
#################################################

######################
# Substitute Netherlands severe cases for hospitalizations
######################
hospDataNL <- read.csv("../downloaded_data/netherlands/COVID-19_casus_landelijk.csv",
                       stringsAsFactors=FALSE, sep=";") %>%
  as_tibble(.) %>%
  dplyr::filter(., (Hospital_admission=="Yes" | Deceased=="Yes") &
                #(lubridate::date(Date_statistics)<="2020-09-01")) %>%
                (lubridate::date(Date_statistics)<="2020-05-11")) %>%
  group_by(., Agegroup) %>%
  summarize(., nHosp=sum((Hospital_admission=="Yes")),
            nDead=sum(Deceased=="Yes"),
            nDeadHosp=sum((Hospital_admission=="Yes" & Deceased=="Yes")),
            nSevere=sum((Hospital_admission=="Yes" | Deceased=="Yes"))) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Agegroup)) %>%
  dplyr::filter(., Agegroup!="<50")

countryDataMod[countryDataMod$Location=="Netherlands",][["Hospitalized"]] <-
  hospDataNL$nHosp

######################
# Add estimated deaths for ICU death ranges in England
######################
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/latest#deaths-registered-by-age-group
icnarcAge <- c("16-39", "40-49", "50-59", "60-69", "70-79", "80+")
icuDeathsEngland <- c(119, 290, 918, 1372, 1069, 176)
totICUEng <- 10287
propICUEng <- c(2.2, 5.8, 13.7, 27.7, 29.7, 18.0, 2.9)
ICU_England <- round(propICUEng/100*totICUEng)

ageDeathsEngland <- c("0-14", "15-44", "45-64", "65-74", "75-84", "85+")
deathsTotEngland <- c(6, 536, 4800, 7392, 16226, 21179)

# redistribute deaths to match ICNARC bins
# https://www.populationpyramid.net/united-kingdom/2020/
splitVecs <- list("15-44"=list(seq(15, 19, 1), seq(20, 39, 1), seq(40, 44, 1)),
                  "46-64"=list(seq(45, 49, 1), seq(50, 59, 1), seq(60, 64, 1)),
                  "65-74"=list(seq(65, 69, 1), seq(70, 74, 1)),
                  "75-84"=list(seq(75, 79, 1), seq(80, 84, 1)))
reducedDeathVec <- deathsTotEngland[-c(1, 6)]
splitDeaths <- list()
for (s in c(1:length(splitVecs))) {
  sumMort <- NULL
  for (ss in c(1:length(splitVecs[[s]]))) {
    mort <- proportion_samples(model=serologyModel$model,
                                   ageVec=splitVecs[[s]][[ss]],
                                   meanAge=serologyModel$meanAge,
                                   sdAge=serologyModel$sdAge)
    sumMort[ss] <- sum(mort$prop_mean)
  }
  if (s==4) {sumMort <- sumMort*c(0.6, 0.4)}
  splitDeaths[[s]] <- round(reducedDeathVec[s]*sumMort/sum(sumMort))
}

# redistributed deaths
ageDeathsEngland_red <- c("15-39", "40-49", "50-59", "60-69", "70-79", "80+")
totalDeathsEngland_red <- c(splitDeaths[[1]][1]+splitDeaths[[1]][2],
                            splitDeaths[[1]][3]+splitDeaths[[2]][1],
                            splitDeaths[[2]][2],
                            splitDeaths[[2]][3]+splitDeaths[[3]][1],
                            splitDeaths[[3]][2]+round(splitDeaths[[4]][1]),
                            round(splitDeaths[[4]][2])+deathsTotEngland[6])

englandICUind <- with(countryDataMod, Location=="England" & !is.na(ICU) & (meanAge>40))
countryDataMod[englandICUind,][["Deaths"]] <- totalDeathsEngland_red[2:6]


#######################
# Substitute critical and severe for hospital and ICU in New Zealand
#######################
age_NZ <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")
# absolute numbers are reported
Hospitalized_NZ <- c(1, 3, 7, 18, 17, 29, 25, 18, 8, 4)
ICU_NZ <-    c(0, 1, 1, 2, 3, 4, 4, 4, 0, 0)
deaths_NZ <- c(0, 0, 0, 0, 0, 3, 3, 7, 8, 5)

countryDataMod[countryDataMod$Location=="New_Zealand",][["Hospitalized"]] <- Hospitalized_NZ
countryDataMod[countryDataMod$Location=="New_Zealand",][["ICU"]] <- ICU_NZ
countryDataMod[countryDataMod$Location=="New_Zealand",][["Deaths"]] <- deaths_NZ

#######################
# Substitute critical and severe for hospital and ICU in Iceland
#######################
oohDeaths_Iceland <- c(0, 0, 0, 0, 0, 0, 0, 0, 2)
ooiDeaths_Iceland <- c(0, 0, 0, 0, 0, 0, 0, 1, 4)

countryDataMod[countryDataMod$Location=="Iceland",][["Hospitalized"]] <- 
  countryDataMod[countryDataMod$Location=="Iceland",][["Hospitalized"]] - oohDeaths_Iceland
countryDataMod[countryDataMod$Location=="Iceland",][["ICU"]] <-
  countryDataMod[countryDataMod$Location=="Iceland",][["ICU"]] - ooiDeaths_Iceland


#################################################
#################################################
# Plot observed vs expected ratio
#################################################
#################################################

#######################################
#### Enlongate format for plotting
#######################################
countryLong <- tidyr::pivot_longer(countryDataMod, cols=c("Hospitalized", "ICU"),
                                   names_to="Outcome_type",
                                   values_to="Outcome_count") %>%
  dplyr::filter(., !is.na(Deaths) & !is.na(Outcome_count))

deathData <- dplyr::filter(countryDataMod, !is.na(Deaths))

public_vs_expected_ratio <- countryLong %>%
  ggplot(., aes(x=meanAge, y=Deaths/Outcome_count*100, color=Location)) +
  geom_point(size=locPointSize) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  geom_line(data=lethalityPosterior, aes(x=meanAge, y=prop_mean*100),
            color="black", linetype="solid",
            size=regLineSize, inherit.aes=FALSE) +
  geom_ribbon(data=lethalityPosterior,
              aes(x=meanAge, ymin=prop_L*100, ymax=prop_H*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  facet_grid(.~Outcome_type) +
  theme_bw() +
  xlab("Age") +
  ylab(expression(frac("Deaths", "Hospitalization")%*%100))
#  ylab("Deaths/Hospitalization")

ggsave("../data/plots/6_public_vs_expected_death_ratio.png",
       public_vs_expected_ratio, width=30, height=10, units="cm")


###############################
# Calculate ratio of ratios, to estimate % hospital deaths
###############################
ratioDf <- dplyr::mutate(countryDataMod, deathHospRatio=Deaths/Hospitalized,
  deathICURatio=Deaths/ICU)

# Get estimated hospital/ICU mortality for each age
mortHosp_samp <-
  proportion_samples(model=lethalityModels$model[["Hospitalized"]],
                     ageVec=ratioDf$meanAge,
                     meanAge=lethalityModels$meanAge[["Hospitalized"]],
                     sdAge=lethalityModels$sdAge[["Hospitalized"]])
mortICU_samp <-
  proportion_samples(model=lethalityModels$model[["ICU"]],
                     ageVec=ratioDf$meanAge,
                     meanAge=lethalityModels$meanAge[["ICU"]],
                     sdAge=lethalityModels$sdAge[["ICU"]])

ratioDf$hospMort <- mortHosp_samp$prop_mean
ratioDf$icuMort <- mortICU_samp$prop_mean

ratioDf <- ratioDf %>%
  dplyr::filter(., Location!="South_Korea") %>%
  dplyr::mutate(., hospRatioRatio=hospMort/deathHospRatio,
              icuRatioRatio=icuMort/deathICURatio) %>%
  tidyr::pivot_longer(., cols=c("hospRatioRatio", "icuRatioRatio"),
  names_to="Outcome_type", values_to="ROR") %>%
  dplyr::filter(., !is.na(ROR) & !is.infinite(ROR))
ratioDf$Outcome_type <- factor(ratioDf$Outcome_type)
levels(ratioDf$Outcome_type) <- c("Hospitalized", "ICU")

ratio_of_ratios <- ratioDf %>%
  ggplot(., aes(x=meanAge, y=ROR*100, color=Location)) +
  geom_point(size=locPointSize) +
  geom_line(alpha=locLineAlpha, size=locLineSize) +
  scale_y_continuous(trans='log10', labels=scaleFun) +
  facet_grid(.~Outcome_type) +
  theme_bw() +
  xlab("Age") +
  ylab("Ratio of ratios")

ggsave("../data/plots/6_ratio_of_ratios.png",
       ratio_of_ratios, width=30, height=10, units="cm")


#################################################
#################################################
# Compare to actual data
#################################################
#################################################

###############
#### ICU deaths
###############

###############
# New Zealand
###############
age_NZ <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")
totalDeaths_NZ <- c(0, 0, 0, 0, 0, 2, 3, 7, 8, 5)
icuDeaths_NZ <- c(0, 0, 0, 0, 0, 2, 1, 1, 0, 0)
nzICUDeathDf <- data.frame(Age=age_NZ, icuDeaths=icuDeaths_NZ,
                        totalDeaths=totalDeaths_NZ,
                        Location="New_Zealand",
                        meanAge=mid_bin_age(age_NZ))

##############
# England
##############
englandICUDeathDf <- data.frame(Age=icnarcAge,
                                icuDeaths=icuDeathsEngland,
                                totalDeaths=totalDeathsEngland_red,
                                Location="England",
                                meanAge=mid_bin_age(icnarcAge))

icuDeathDf <- rbind(nzICUDeathDf, englandICUDeathDf) %>%
  dplyr::mutate(., icuDeathRatio=icuDeaths/totalDeaths,
    Type="Observation") %>%
  dplyr::select(Age, Location, meanAge, icuDeathRatio, Type)

#### Append the estimate from the ratio of ratios to dataframe
estimatedICURatio <- ratioDf %>%
  filter(., Location %in% c("New_Zealand", "England") & Outcome_type=="ICU") %>%
  dplyr::select(., Age, Location, meanAge, ROR) %>%
  dplyr::rename(., icuDeathRatio=ROR) %>%
  dplyr::mutate(., Type="Estimate")
icuDeathDf <- rbind(icuDeathDf, estimatedICURatio)


###############
#### Hospital deaths
###############

##############
# England
##############
# England data in https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/latest#deaths-registered-by-place-of-occurrence
# says that 70% of deaths occured in the hospital

ageHospitalDeaths <- c("0-19", "20-39", "40-59", "60-79", "80+")
hospitalDeathsEngland <- c(20, 208, 2231, 10945, 15389)

# redistributed deaths
ageDeathsEngland_red2 <- c("0-19", "20-39", "40-59", "60-79", "80+")
totalDeathsEngland_red2 <- c(deathsTotEngland[1]+splitDeaths[[1]][1],
                             splitDeaths[[1]][2],
                             sum(totalDeathsEngland_red[c(2:3)]),
                             sum(totalDeathsEngland_red[c(4:5)]),
                             totalDeathsEngland_red[6])

englandHospDeathDf <- data.frame(Age=ageHospitalDeaths,
                                 hospDeaths=hospitalDeathsEngland,
                                 totalDeaths=totalDeathsEngland_red2,
                                 Location="England",
                                 meanAge=mid_bin_age(ageHospitalDeaths))

######################
# Netherlands
######################
hospDataNL <- read.csv("../downloaded_data/netherlands/COVID-19_casus_landelijk.csv",
                       stringsAsFactors=FALSE, sep=";") %>%
  as_tibble(.) %>%
  dplyr::filter(., (Hospital_admission=="Yes" | Deceased=="Yes") &
                #(lubridate::date(Date_statistics)<="2020-09-01")) %>%
                (lubridate::date(Date_statistics)<="2020-04-11")) %>%
  group_by(., Agegroup) %>%
  summarize(., nHosp=sum((Hospital_admission=="Yes")),
            nDead=sum(Deceased=="Yes"),
            nDeadHosp=sum((Hospital_admission=="Yes" & Deceased=="Yes")),
            nSevere=sum((Hospital_admission=="Yes" | Deceased=="Yes"))) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Agegroup)) %>%
  dplyr::filter(., Agegroup != "<50" & nDead > 0)


netherlandsHospDeathDf <- data.frame(Age=hospDataNL$Agegroup,
                                 hospDeaths=hospDataNL$nDeadHosp,
                                 totalDeaths=hospDataNL$nDead,
                                 Location="Netherlands",
                                 meanAge=mid_bin_age(hospDataNL$Agegroup))

hospDeathDf <- rbind(englandHospDeathDf, netherlandsHospDeathDf) %>%
  dplyr::mutate(., hospDeathRatio=hospDeaths/totalDeaths,
                Type="Observation") %>%
  dplyr::select(Age, Location, meanAge, hospDeathRatio, Type)

#### Add ROR estimation of % of deaths in hospital
estimatedHospRatio <- ratioDf %>%
  filter(., Location %in% c("England", "Netherlands") &
         Outcome_type=="Hospitalized") %>%
  dplyr::select(., Age, Location, meanAge, ROR) %>%
  dplyr::rename(., hospDeathRatio=ROR) %>%
  dplyr::mutate(., Type="Estimate")

hospDeathDf <- rbind(hospDeathDf, estimatedHospRatio)


#### Plot death ratios
icuDeathRatioPlot <- icuDeathDf %>%
  dplyr::filter(., !is.na(icuDeathRatio)) %>%
  ggplot(., aes(x=meanAge, y=icuDeathRatio, color=Location,
                linetype=Type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  #ylim(c(0,1.2)) +
  ylab(expression(frac("ICU deaths", "Total deaths"))) +
  xlab("Age")

hospDeathRatioPlot <- hospDeathDf %>%
  dplyr::filter(., !is.na(hospDeathRatio)) %>%
  ggplot(., aes(x=meanAge, y=hospDeathRatio, color=Location,
                linetype=Type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  #ylim(c(0,1.2)) +
  ylab(expression(frac("Hospital deaths", "Total deaths"))) +
  xlab("Age")

deathRatioPlot <- ggarrange(plotlist=list(hospDeathRatioPlot, icuDeathRatioPlot))



#####################
#####################
#
# Generate corrected data set
#
#####################
#####################

unique(countryData$Location)

#########
# Iceland no change needed due to good data
#########
iceland <- filter(countryData, Location=="Iceland")

#########
# New Zealand no change needed due to good data
#########
newZealand <- filter(countryData, Location=="New_Zealand")

#########
# South Korea are also OK, with nested data and
# everything computed over a well defined population
# with known outcome
#########
southKorea <- filter(countryData, Location=="South_Korea")

#########
# Spain
#########
spain <- filter(countryData, Location=="Spain") %>%
  dplyr::mutate(., hospDeathRatio=Deaths/Hospitalized,
                icuDeathRatio=Deaths/ICU)

#ooh = out of hospital
#ooi = out of icu
mortHospSpain <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=spain$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)

spainHospRatioSamp <- matrix(nrow=nrow(spain), ncol=0)
for (s in unique(mortHospSpain$samples$sample)) {
  sampleMort <- dplyr::filter(mortHospSpain$samples, sample==s)
  oohProp <- (1-pmin(1, sampleMort$proportion/spain$hospDeathRatio))
  spainHospRatioSamp <- cbind(spainHospRatioSamp, oohProp)
}

meanHospRatioSpain <- rowMeans(spainHospRatioSamp)
quantsHospRatioSpain <- matrixStats::rowQuantiles(spainHospRatioSamp, probs=c(0.025, 0.975))
newHospSpain <- round(spain$Hosp + spain$Deaths*meanHospRatioSpain)

#ooh = out of hospital
#ooi = out of icu
mortICUSpain <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=spain$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)

spainICURatioSamp <- matrix(nrow=nrow(spain), ncol=0)
for (s in unique(mortICUSpain$samples$sample)) {
  sampleMort <- dplyr::filter(mortICUSpain$samples, sample==s)
  ooiProp <- (1-pmin(1, sampleMort$proportion/spain$icuDeathRatio))
  spainICURatioSamp <- cbind(spainICURatioSamp, ooiProp)
}
meanICURatioSpain <- rowMeans(spainICURatioSamp)
quantsICURatioSpain <- matrixStats::rowQuantiles(spainICURatioSamp, probs=c(0.025, 0.975))
newICUSpain <- round(spain$ICU + spain$Deaths*meanICURatioSpain)

spain$Hospitalized <- newHospSpain
spain$ICU <- newICUSpain


############
# Ireland
############
ireland <- filter(countryData, Location=="Ireland") %>%
  dplyr::mutate(., hospDeathRatio=Deaths/Hospitalized,
                icuDeathRatio=Deaths/ICU)

mortHospIreland <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=ireland$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
irelandHospRatioSamp <- matrix(nrow=nrow(ireland), ncol=0)
for (s in unique(mortHospIreland$samples$sample)) {
  sampleMort <- dplyr::filter(mortHospIreland$samples, sample==s)
  oohProp <- (1-pmin(1, sampleMort$proportion/ireland$hospDeathRatio))
  irelandHospRatioSamp <- cbind(irelandHospRatioSamp, oohProp)
}
meanHospRatioIreland <- rowMeans(irelandHospRatioSamp)
quantsHospRatioIreland <- matrixStats::rowQuantiles(irelandHospRatioSamp, probs=c(0.025, 0.975))
newHospIreland <- round(ireland$Hosp + ireland$Deaths*meanHospRatioIreland)

mortICUIreland <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=ireland$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
irelandICURatioSamp <- matrix(nrow=nrow(ireland), ncol=0)
for (s in unique(mortICUIreland$samples$sample)) {
  sampleMort <- dplyr::filter(mortICUIreland$samples, sample==s)
  ooiProp <- (1-pmin(1, sampleMort$proportion/ireland$icuDeathRatio))
  irelandICURatioSamp <- cbind(irelandICURatioSamp, ooiProp)
}
meanICURatioIreland <- rowMeans(irelandICURatioSamp)
quantsICURatioIreland <- matrixStats::rowQuantiles(irelandICURatioSamp, probs=c(0.025, 0.975))
newICUIreland <- round(ireland$ICU + ireland$Deaths*meanICURatioIreland)

ireland$Hospitalized <- newHospIreland
ireland$ICU <- newICUIreland


############
# Sweden
############
sweden <- filter(countryData, Location=="Sweden") %>%
  dplyr::mutate(., hospDeathRatio=Deaths/Hospitalized,
                icuDeathRatio=Deaths/ICU)

mortICUSweden <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=sweden$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
swedenICURatioSamp <- matrix(nrow=nrow(sweden), ncol=0)
for (s in unique(mortICUSweden$samples$sample)) {
  sampleMort <- dplyr::filter(mortICUSweden$samples, sample==s)
  ooiProp <- (1-pmin(1, sampleMort$proportion/sweden$icuDeathRatio))
  swedenICURatioSamp <- cbind(swedenICURatioSamp, ooiProp)
}
meanICURatioSweden <- rowMeans(swedenICURatioSamp)
quantsICURatioSweden <- matrixStats::rowQuantiles(swedenICURatioSamp, probs=c(0.025, 0.975))
newICUSweden <- round(sweden$ICU + sweden$Deaths*meanICURatioSweden)

sweden$ICU <- newICUSweden


############
# France
############
# Epidemiologique report specifies hospital deaths and deaths
# in care homes. Just add those care Home deaths as new severe
# cases in the oldest age strata
careHomeDeaths <- 4368
franceDeathsHosp <- c(12, 49, 132, 503, 1039, 1665, 3684)
oohDeaths <- c(0, 0, 0, 0, 0, round(careHomeDeaths*(1/3)), round(careHomeDeaths*(2/3)))
franceDeathsTot <- franceDeathsHosp + oohDeaths
# No death data to use
france <- filter(countryData, Location=="Ile_de_France")
newHospFrance <- france$Hospitalized + oohDeaths

france$Hospitalized <- newHospFrance
france$Deaths <- franceDeathsTot


############
# England
############
england <- dplyr::filter(countryData, Location=="England")

englandICU <- dplyr::filter(england, !is.na(ICU))
ooiDeathsEngland <- totalDeathsEngland_red - icuDeathsEngland
newICUEngland <- englandICU$ICU[-1] + ooiDeathsEngland

hospAgesEngland <- c("0-17", "18-64", "65+")
ageDeathsEngland <- c("0-14", "15-44", "45-64", "65-74", "75-84", "85+")
deathsTotEngland <- c(6, 536, 4800, 7392, 16226, 21179)
deathsTotEngland_red <- c(deathsTotEngland[1],
                          sum(deathsTotEngland[2:3]),
                          sum(deathsTotEngland[4:6]))

ageHospitalDeaths <- c("0-19", "20-39", "40-59", "60-79", "80+")
hospitalDeathsEngland <- c(20, 208, 2231, 10945, 15389)
hospDeathsEngland_red <- round(c(hospitalDeathsEngland[1],
                           sum(hospitalDeathsEngland[2:3])+hospitalDeathsEngland[4]/4,
                           hospitalDeathsEngland[4]*3/4+hospitalDeathsEngland[5]))

oohDeathsEngland <- pmax(0, deathsTotEngland_red-hospDeathsEngland_red)
englandHosp <- dplyr::filter(england, !is.na(Hospitalized))
newHospEngland <- englandHosp$Hospitalized + oohDeathsEngland

englandICU$ICU[-1] <- newICUEngland
englandHosp$Hospitalized <- newHospEngland

england <- rbind(englandICU, englandHosp)

############
# Netherlands
############
# Netherlands data already includes deaths outside of hospital
netherlands <- filter(countryData, Location=="Netherlands")


############
# Atlanta
############
atlanta <- filter(countryData, Location=="Atlanta") %>%
  dplyr::mutate(., hospDeathRatio=Deaths/Hospitalized,
                icuDeathRatio=Deaths/ICU)

mortHospAtlanta <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=atlanta$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
atlantaHospRatioSamp <- matrix(nrow=nrow(atlanta), ncol=0)
for (s in unique(mortHospAtlanta$samples$sample)) {
  sampleMort <- dplyr::filter(mortHospAtlanta$samples, sample==s)
  oohProp <- (1-pmin(1, sampleMort$proportion/atlanta$hospDeathRatio))
  atlantaHospRatioSamp <- cbind(atlantaHospRatioSamp, oohProp)
}
meanHospRatioAtlanta <- rowMeans(atlantaHospRatioSamp)
quantsHospRatioAtlanta <- matrixStats::rowQuantiles(atlantaHospRatioSamp, probs=c(0.025, 0.975))
newHospAtlanta <- round(atlanta$Hosp + atlanta$Deaths*meanHospRatioAtlanta)

mortICUAtlanta <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=atlanta$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
atlantaICURatioSamp <- matrix(nrow=nrow(atlanta), ncol=0)
for (s in unique(mortICUAtlanta$samples$sample)) {
  sampleMort <- dplyr::filter(mortICUAtlanta$samples, sample==s)
  ooiProp <- (1-pmin(1, sampleMort$proportion/atlanta$icuDeathRatio))
  atlantaICURatioSamp <- cbind(atlantaICURatioSamp, ooiProp)
}
meanICURatioAtlanta <- rowMeans(atlantaICURatioSamp)
quantsICURatioAtlanta <- matrixStats::rowQuantiles(atlantaICURatioSamp, probs=c(0.025, 0.975))
newICUAtlanta <- round(atlanta$ICU + atlanta$Deaths*meanICURatioAtlanta)

atlanta$Hospitalized <- newHospAtlanta
atlanta$ICU <- newICUAtlanta


############
# NYC
############
nyc <- filter(countryData, Location=="NYC") %>%
  dplyr::mutate(., hospDeathRatio=Deaths/Hospitalized,
                icuDeathRatio=Deaths/ICU)

nursingDeaths_probableCOVID <- 1226
hospitalDeaths_probableCOVID <- 2719

mortHospNYC <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=nyc$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
nycHospRatioSamp <- matrix(nrow=nrow(nyc), ncol=0)
for (s in unique(mortHospNYC$samples$sample)) {
  sampleMort <- dplyr::filter(mortHospNYC$samples, sample==s)
  oohProp <- (1-pmin(1, sampleMort$proportion/nyc$hospDeathRatio))
  nycHospRatioSamp <- cbind(nycHospRatioSamp, oohProp)
}
meanHospRatioNYC <- rowMeans(nycHospRatioSamp)
quantsHospRatioNYC <- matrixStats::rowQuantiles(nycHospRatioSamp, probs=c(0.025, 0.975))
newHospNYC <- round(nyc$Hosp + nyc$Deaths*meanHospRatioNYC)

nyc$Hospitalized <- newHospNYC


############
# Ontario
############
ontario <- filter(countryData, Location=="Ontario") %>%
  dplyr::mutate(., hospDeathRatio=Deaths/Hospitalized,
                icuDeathRatio=Deaths/ICU)

mortHospOntario <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=ontario$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
ontarioHospRatioSamp <- matrix(nrow=nrow(ontario), ncol=0)
for (s in unique(mortHospOntario$samples$sample)) {
  sampleMort <- dplyr::filter(mortHospOntario$samples, sample==s)
  oohProp <- (1-pmin(1, sampleMort$proportion/ontario$hospDeathRatio))
  ontarioHospRatioSamp <- cbind(ontarioHospRatioSamp, oohProp)
}
meanHospRatioOntario <- rowMeans(ontarioHospRatioSamp)
quantsHospRatioOntario <- matrixStats::rowQuantiles(ontarioHospRatioSamp, probs=c(0.025, 0.975))
newHospOntario <- round(ontario$Hosp + ontario$Deaths*meanHospRatioOntario)

mortICUOntario <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=ontario$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
ontarioICURatioSamp <- matrix(nrow=nrow(ontario), ncol=0)
for (s in unique(mortICUOntario$samples$sample)) {
  sampleMort <- dplyr::filter(mortICUOntario$samples, sample==s)
  ooiProp <- (1-pmin(1, sampleMort$proportion/ontario$icuDeathRatio))
  ontarioICURatioSamp <- cbind(ontarioICURatioSamp, ooiProp)
}
meanICURatioOntario <- rowMeans(ontarioICURatioSamp)
quantsICURatioOntario <- matrixStats::rowQuantiles(ontarioICURatioSamp, probs=c(0.025, 0.975))
newICUOntario <- round(ontario$ICU + ontario$Deaths*meanICURatioOntario)

ontario$Hospitalized <- newHospOntario
ontario$ICU <- newICUOntario

##############
# Put corrected locations together
##############
locationsList <- list(iceland, newZealand, southKorea, spain, ireland,
                      ireland, sweden, france, england, netherlands,
                      atlanta, nyc, ontario)

correctedDf <- data.frame()
for (l in c(1:length(locationsList))) {
  if ("hospDeathRatio" %in% names(locationsList[[l]])) {
    locationsList[[l]] <- dplyr::select(locationsList[[l]], -hospDeathRatio,
                                        -icuDeathRatio)
  }
  correctedDf <- rbind(correctedDf, locationsList[[l]])
}

write.csv(correctedDf, "../data/collected_data/locations_serology_data_corrected.csv",
          row.names=FALSE)

