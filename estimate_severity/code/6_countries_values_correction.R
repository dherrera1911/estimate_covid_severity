library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(wpp2019)
source("./functions_auxiliary.R")
data(popM)
data(popF)

##################################
# Load in-hospital lethality data
##################################
lethalityData <- read.csv("../data/collected_data/hospitalized_patient_studies.csv") %>%
  as_tibble(.)
lethalityModels <- readRDS("../data/processed_data/4_hospital_lethality_fit.RDS")

##################################
# Load serology fit, which is used to redistirbute some death data
##################################
serologyModel <- readRDS("../data/processed_data/3_serology_fit_Deaths.RDS")

######################################
# Load country data and separate into dataframes for hospitalization and ICU
######################################
countryData <- read.csv("../data/collected_data/locations_serology_data.csv")

countryDataMod <- countryData

#####################
#####################
#
# Generate corrected data set
#
#####################
#####################
#ooh = out of hospital
#ooi = out of icu

unique(countryData$Location)

#########
# Iceland no change needed due to good data
#########
iceland <- filter(countryData, Location=="Iceland")
iceland$oohDeaths <- NA
iceland$ooiDeaths <- NA

#########
# New Zealand no change needed due to good data
#########
newZealand <- filter(countryData, Location=="New_Zealand")
newZealand$oohDeaths <- NA
newZealand$ooiDeaths <- NA

#########
# South Korea are also OK, with nested data and
# everything computed over a well defined population
# with known outcome
#########
southKorea <- filter(countryData, Location=="South_Korea")
southKorea$oohDeaths <- NA
southKorea$ooiDeaths <- NA

#########
# Spain
#########
spain <- filter(countryData, Location=="Spain")

# Hospital
mortHospSpain <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=spain$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohSpain <- ooh_deaths_estimation(mortalitySamples=mortHospSpain,
                                  hospitalized=spain$Hospitalized,
                                  deaths=spain$Deaths)
newHospSpain <- spain$Hosp + oohSpain$mean

# ICU
mortICUSpain <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=spain$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
ooiSpain <- ooh_deaths_estimation(mortalitySamples=mortICUSpain,
                                  hospitalized=spain$ICU,
                                  deaths=spain$Deaths)
newICUSpain <- spain$ICU + ooiSpain$mean

spain$Hospitalized <- newHospSpain
spain$ICU <- newICUSpain
spain$oohDeaths <- oohSpain$mean
spain$ooiDeaths <- ooiSpain$mean


############
# Ireland
############
ireland <- filter(countryData, Location=="Ireland")

mortHospIreland <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=ireland$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohIreland <- ooh_deaths_estimation(mortalitySamples=mortHospIreland,
                                  hospitalized=ireland$Hospitalized,
                                  deaths=ireland$Deaths)
newHospIreland <- ireland$Hosp + oohIreland$mean



mortICUIreland <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=ireland$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
ooiIreland <- ooh_deaths_estimation(mortalitySamples=mortICUIreland,
                                  hospitalized=ireland$ICU,
                                  deaths=ireland$Deaths)
newICUIreland <- ireland$ICU + ooiIreland$mean

ireland$Hospitalized <- newHospIreland
ireland$ICU <- newICUIreland
ireland$oohDeaths <- oohIreland$mean
ireland$ooiDeaths <- ooiIreland$mean


############
# Sweden
############
sweden <- filter(countryData, Location=="Sweden")

mortICUSweden <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=sweden$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)

ooiSweden <- ooh_deaths_estimation(mortalitySamples=mortICUSweden,
                                  hospitalized=sweden$ICU,
                                  deaths=sweden$Deaths)
newICUSweden <- sweden$ICU + ooiSweden$mean

sweden$Hospitalized <- NA
sweden$ICU <- newICUSweden
sweden$oohDeaths <- NA
sweden$ooiDeaths <- ooiSweden$mean


############
# Atlanta
############
atlanta <- filter(countryData, Location=="Atlanta")

# Hospital
mortHospAtlanta <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=atlanta$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohAtlanta <- ooh_deaths_estimation(mortalitySamples=mortHospAtlanta,
                                  hospitalized=atlanta$Hospitalized,
                                  deaths=atlanta$Deaths)
newHospAtlanta <- atlanta$Hosp + oohAtlanta$mean

#ICU
mortICUAtlanta <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=atlanta$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)

ooiAtlanta <- ooh_deaths_estimation(mortalitySamples=mortICUAtlanta,
                                  hospitalized=atlanta$ICU,
                                  deaths=atlanta$Deaths)
newICUAtlanta <- atlanta$ICU + ooiAtlanta$mean

atlanta$Hospitalized <- newHospAtlanta
atlanta$ICU <- newICUAtlanta
atlanta$oohDeaths <- oohAtlanta$mean
atlanta$ooiDeaths <- ooiAtlanta$mean


############
# NYC
############
nyc <- filter(countryData, Location=="NYC")

nursingDeaths_probableCOVID <- 1226
hospitalDeaths_probableCOVID <- 2719

mortHospNYC <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=nyc$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohNYC <- ooh_deaths_estimation(mortalitySamples=mortHospNYC,
                                  hospitalized=nyc$Hospitalized,
                                  deaths=nyc$Deaths)
newHospNYC <- nyc$Hosp + oohNYC$mean

nyc$Hospitalized <- newHospNYC
nyc$oohDeaths <- oohNYC$mean
nyc$ooiDeaths <- NA


############
# Ontario
############
ontario <- filter(countryData, Location=="Ontario")

mortHospOntario <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=ontario$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohOntario <- ooh_deaths_estimation(mortalitySamples=mortHospOntario,
                                  hospitalized=ontario$Hospitalized,
                                  deaths=ontario$Deaths)
newHospOntario <- ontario$Hosp + oohOntario$mean


mortICUOntario <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=ontario$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
ooiOntario <- ooh_deaths_estimation(mortalitySamples=mortICUOntario,
                                  hospitalized=ontario$ICU,
                                  deaths=ontario$Deaths)
newICUOntario <- ontario$ICU + ooiOntario$mean

ontario$Hospitalized <- newHospOntario
ontario$ICU <- newICUOntario
ontario$oohDeaths <- oohOntario$mean
ontario$ooiDeaths <- ooiOntario$mean


############
# Geneva
############
geneva <- filter(countryData, Location=="Geneva")

mortHospGeneva <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=geneva$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohGeneva <- ooh_deaths_estimation(mortalitySamples=mortHospGeneva,
                                  hospitalized=geneva$Hospitalized,
                                  deaths=geneva$Deaths)
newHospGeneva <- geneva$Hosp + oohGeneva$mean

geneva$Hospitalized <- newHospGeneva
geneva$oohDeaths <- oohGeneva$mean
geneva$ooiDeaths <- NA


############
# Belgium
############
belgium <- filter(countryData, Location=="Belgium")

mortHospBelgium <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=belgium$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohBelgium <- ooh_deaths_estimation(mortalitySamples=mortHospBelgium,
                                  hospitalized=belgium$Hospitalized,
                                  deaths=belgium$Deaths)
newHospBelgium <- belgium$Hosp + oohBelgium$mean


mortICUBelgium <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=belgium$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
ooiBelgium <- ooh_deaths_estimation(mortalitySamples=mortICUBelgium,
                                  hospitalized=belgium$ICU,
                                  deaths=belgium$Deaths)
newICUBelgium <- belgium$ICU + ooiBelgium$mean


belgium$Hospitalized <- newHospBelgium
belgium$ICU <- newICUBelgium
belgium$oohDeaths <- oohBelgium$mean
belgium$ooiDeaths <- ooiBelgium$mean


############
# France
############
france <- filter(countryData, Location=="Ile_de_France")
# Epidemiologique report specifies hospital deaths and deaths
# in care homes. Just add those care Home deaths as new severe
# cases in the oldest age strata
careHomeDeaths <- 4368
franceDeathsHosp <- c(12, 49, 132, 503, 1039, 1665, 3684)
deathsCareProp <- franceDeathsHosp[6:7]/sum(franceDeathsHosp[6:7])
oohDeaths <- c(0, 0, 0, 0, 0, round(careHomeDeaths*deathsCareProp[1]),
               round(careHomeDeaths*deathsCareProp[2]))

franceDeathsTot <- franceDeathsHosp + oohDeaths
newHospFrance <- france$Hospitalized + oohDeaths

france$Hospitalized <- newHospFrance
france$Deaths <- franceDeathsTot
france$oohDeaths <- oohDeaths
france$ooiDeaths <- NA


############
# Netherlands
############
# Netherlands data already includes deaths outside of hospital
netherlands <- filter(countryData, Location=="Netherlands")

# get netherlands ooh deaths
hospDataNL <- read.csv("../downloaded_data/netherlands/COVID-19_casus_landelijk.csv",
                       stringsAsFactors=FALSE, sep=";") %>%
  as_tibble(.) %>%
  dplyr::filter(., (Hospital_admission=="Yes" | Deceased=="Yes") &
                (lubridate::date(Date_statistics)<="2020-05-11")) %>%
  group_by(., Agegroup) %>%
  summarize(., oohDeaths=sum(Hospital_admission=="No" & Deceased=="Yes")) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Agegroup)) %>%
  dplyr::filter(., Agegroup!="<50")

netherlands$oohDeaths <- hospDataNL$oohDeaths
netherlands$ooiDeaths <- NA


############
# England
############
england <- dplyr::filter(countryData, Location=="England")

#ICU
# Note that ICU death bins are different from ICU bins
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/latest#deaths-registered-by-age-group
# ICU deaths bins

ageDeathsEngland <- c("0-14", "15-44", "45-64", "65-74", "75-84", "85+")
deathsTotEngland <- c(6, 536, 4800, 7392, 16226, 21179)
reducedDeathVec <- deathsTotEngland[-c(1, 6)]
# redistribute deaths to match ICNARC bins
# https://www.populationpyramid.net/united-kingdom/2020/
# Redistribute according to relative ICU risk within each
# bin fitted to the non-corrected data
splitVecs <- list("15-44"=list(seq(15, 30, 1), seq(30, 39, 1), seq(40, 44, 1)),
                  "46-64"=list(seq(45, 49, 1), seq(50, 59, 1), seq(60, 64, 1)),
                  "65-74"=list(seq(65, 69, 1), seq(70, 74, 1)),
                  "75-84"=list(seq(75, 79, 1), seq(80, 84, 1)))
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
  if (s==4) {
    sumMort <- sumMort*c(0.6, 0.4)
  }
  splitDeaths[[s]] <- round(reducedDeathVec[s]*sumMort/sum(sumMort))
}

ageDeathsEngland_red <- c("15-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
totalDeathsEngland_red <- c(splitDeaths[[1]][1],
                            splitDeaths[[1]][2],
                            splitDeaths[[1]][3]+splitDeaths[[2]][1],
                            splitDeaths[[2]][2],
                            splitDeaths[[2]][3]+splitDeaths[[3]][1],
                            splitDeaths[[3]][2]+round(splitDeaths[[4]][1]),
                            round(splitDeaths[[4]][2])+deathsTotEngland[6])

englandICU <- dplyr::filter(england, !is.na(ICU))
englandICU$Deaths <- totalDeathsEngland_red


mortICUEngland <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=englandICU$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
ooiEngland <- ooh_deaths_estimation(mortalitySamples=mortICUEngland,
                                  hospitalized=englandICU$ICU,
                                  deaths=englandICU$Deaths)
newICUEngland <- englandICU$ICU + ooiEngland$mean

englandICU$Hospitalized <- NA
englandICU$ooiDeaths <- ooiEngland$mean
englandICU$ICU <- newICUEngland
englandICU$oohDeaths <- NA



# redistribute England deaths (only one bin) for Hopsital age bins
ageDeathsEngland <- c("0-14", "15-44", "45-64", "65-74", "75-84", "85+")
deathsTotEngland <- c(6, 536, 4800, 7392, 16226, 21179)
# redistribute deaths to match Hospital bins
# bin fitted to the non-corrected data
ageVec <- seq(15, 44, 1)
mort <- proportion_samples(model=serologyModel$model,
                               ageVec=ageVec,
                               meanAge=serologyModel$meanAge,
                               sdAge=serologyModel$sdAge)

# redistribute the deaths of the second bin according to
# the within bin mortality. Will assign minimal deaths
# to younger bin, and the procedure won't change that much
sumMort <- c(sum(mort$prop_mean[1:3]), sum(mort$prop_mean[4:30]))
splitDeathsBin2 <- round(deathsTotEngland[2]*sumMort/sum(sumMort))

# hosp ages 0-17, 18-64, 65-84, 85+
englandHospDeaths <- c(deathsTotEngland[1]+splitDeathsBin2[1],
                       splitDeathsBin2[2]+deathsTotEngland[3],
                       sum(deathsTotEngland[4:5]),
                       deathsTotEngland[6])

englandHosp <- dplyr::filter(england, !is.na(Hospitalized))
englandHosp$Deaths <- englandHospDeaths


mortHospEngland <- proportion_samples(model=lethalityModels$model$Hospitalized,
                               ageVec=englandHosp$meanAge,
                               meanAge=lethalityModels$meanAge$Hospitalized,
                               sdAge=lethalityModels$sdAge$Hospitalized)
oohEngland <- ooh_deaths_estimation(mortalitySamples=mortHospEngland,
                                  hospitalized=englandHosp$Hospitalized,
                                  deaths=englandHosp$Deaths)
newHospEngland <- englandHosp$Hospitalized + oohEngland$mean

englandHosp$Hospitalized <- newHospEngland
englandHosp$oohDeaths <- oohEngland$mean
englandHosp$ICU <- NA
englandHosp$ooiDeaths <- NA


englandDeaths <- dplyr::filter(england, !is.na(Deaths))
englandDeaths$ooiDeaths <- NA
englandDeaths$oohDeaths <- NA


england <- rbind(englandICU, englandHosp, englandDeaths)


##############
# Put corrected locations together
##############
locationsList <- list(iceland, newZealand, southKorea, spain, ireland,
                      ireland, sweden, france, england, netherlands,
                      atlanta, nyc, ontario, geneva, belgium)

correctedDf <- data.frame()
for (l in c(1:length(locationsList))) {
  correctedDf <- rbind(correctedDf, locationsList[[l]])
#  print(names(locationsList[[l]]))
}

write.csv(correctedDf, "../data/collected_data/locations_serology_data_corrected.csv",
          row.names=FALSE)

