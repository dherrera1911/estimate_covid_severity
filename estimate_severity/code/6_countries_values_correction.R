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

######################################
# Load country data and separate into dataframes for hospitalization and ICU
######################################
countryData <- read.csv("../data/collected_data/locations_serology_data.csv")

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
# Iceland no change needed due to ooh and ooi data
#########
iceland <- filter(countryData, Location=="Iceland")
oohDeaths_Iceland <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)
ooiDeaths_Iceland <- c(0, 0, 0, 0, 0, 0, 0, 1, 3)

iceland$oohDeaths <- oohDeaths_Iceland
iceland$ooiDeaths <- ooiDeaths_Iceland
iceland$oohDeaths_source <- "Data"
iceland$ooiDeaths_source <- "Data"

#########
# New Zealand no change needed due to ooh and ooi data
#########
newZealand <- filter(countryData, Location=="New_Zealand")
oohDeaths_NZ <- c(0, 0, 0, 0, 0, 0, 2, 6, 8, 5)
ooiDeaths_NZ <- c(0, 0, 0, 0, 0, 0, 2, 6, 8, 5)
newZealand$oohDeaths <- oohDeaths_NZ
newZealand$ooiDeaths <- ooiDeaths_NZ
newZealand$oohDeaths_source <- "Data"
newZealand$ooiDeaths_source <- "Data"

#########
# South Korea are also OK. Although the disagregated data is not
# presented, the severe and critical cases count considers deaths
#########
southKorea <- filter(countryData, Location=="South_Korea")
southKorea$oohDeaths <- NA
southKorea$ooiDeaths <- NA
southKorea$oohDeaths_source <- "Data"
southKorea$ooiDeaths_source <- "Data"

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
spain$oohDeaths_source <- "Estimated"
spain$ooiDeaths_source <- "Estimated"


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
ireland$oohDeaths_source <- "Estimated"
ireland$ooiDeaths_source <- "Estimated"


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
sweden$oohDeaths_source <- NA
sweden$ooiDeaths_source <- "Estimated"


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
atlanta$oohDeaths_source <- "Estimated"
atlanta$ooiDeaths_source <- "Estimated"


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
nyc$oohDeaths_source <- "Estimated"
nyc$ooiDeaths_source <- NA


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
ontario$oohDeaths_source <- "Estimated"
ontario$ooiDeaths_source<- "Estimated"


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
geneva$oohDeaths_source <- "Estimated"
geneva$ooiDeaths_source <- NA


############
# Belgium 
############
# Severe cases information is already corrected with real data
belgium <- filter(countryData, Location=="Belgium")

mortICUBelgium <- proportion_samples(model=lethalityModels$model$ICU,
                               ageVec=belgium$meanAge,
                               meanAge=lethalityModels$meanAge$ICU,
                               sdAge=lethalityModels$sdAge$ICU)
ooiBelgium <- ooh_deaths_estimation(mortalitySamples=mortICUBelgium,
                                  hospitalized=belgium$ICU,
                                  deaths=belgium$Deaths)
newICUBelgium <- belgium$ICU + ooiBelgium$mean


# ooh deaths Belgium
nursingHomeDeaths_Belgium <- c(0, 0, 58, 508, 952, 4251)
nh_in_hospital <- 1276
nh_out_hospital <- 4494
oohDeaths_Belgium <- round(nursingHomeDeaths_Belgium*nh_out_hospital/
                                                      (nh_in_hospital+nh_out_hospital))
belgium$ICU <- newICUBelgium
belgium$oohDeaths <- oohDeaths_Belgium
belgium$ooiDeaths <- ooiBelgium$mean
belgium$oohDeaths_source <- oohDeaths_Belgium
belgium$ooiDeaths_source <- "Estimated"


############
# France
############
france <- filter(countryData, Location=="Ile_de_France")
# Epidemiologique report specifies age-stratified hospital deaths
# and total deaths in care homes. Since there is no
# age-stratified data on total deaths, the simple correction we are
# applying to other countries can't be used here. We need
# to allocate care-home deaths to get total deaths

# We will allocate the care home deaths in France with the same
# age distribution as in Belgium
careHomeDeathsTot_France <- 4405 # deaths in care home
oohDeaths_France <- round(careHomeDeathsTot_France*oohDeaths_Belgium/sum(oohDeaths_Belgium))
oohDeaths_France <- c(0, oohDeaths_France)

france$oohDeaths <- oohDeaths_France
france$ooiDeaths <- NA
france$oohDeaths_source <- "Data"
france$ooiDeaths_source <- NA


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
netherlands$oohDeaths_source <- "Data"
netherlands$ooiDeaths_source <- NA


############
# England
############
# England is already corrected with real data. Get oohDeaths and ooiDeaths
# to put in the dataframe, but leave severe and critical N alone

# ICU deaths
# ICU and in ICU deaths from:
# "ICNARC report on COVID-19 in critical care 03 July 2020"
# You need to look on google "ICNARC report on COVID-19 in critical
# care 03 July 2020"
age_England_ICU <- c("16-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
totICUEng <- 10287
propICUEng <- c(2.2, 5.8, 13.7, 27.7, 29.7, 18.0, 2.9)
ICU_England <- round(propICUEng/100*totICUEng)
# Note that ICU death bins are different from ICU bins
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/latest#deaths-registered-by-age-group
icnarcDeathAge <- c("16-39", "40-49", "50-59", "60-69", "70-79", "80+")
icuDeathsEngland <- c(119, 290, 918, 1372, 1069, 176)
# redistribute ICU deaths to match ICU patients bins
icuDeathsEngland_y <- round(icuDeathsEngland[1]*ICU_England[1:2]/sum(ICU_England[1:2]))
icuDeathsEngland <- c(icuDeathsEngland_y, icuDeathsEngland[2:6])

# Total deaths England
ageDeathsEngland <- c("0-14", "15-44", "45-64", "65-74", "75-84", "85+")
deathsTotEngland <- c(6, 536, 4800, 7392, 16226, 21179)
reducedDeathVec <- deathsTotEngland[-c(1, 6)]

# Redistribute deaths according to relative risk of death for the ages
# within each bin (taken from Levin). Using Levin within-bin
# balances many relevant factors, since these are expected
# to be relatively homogeneous within bins
levinIFR <- function(age){return(10^(-3.27+0.0524*age))}

# redistribute deaths to match ICNARC bins
# https://www.populationpyramid.net/united-kingdom/2020/
# Redistribute according to relative ICU risk within each
# bin fitted to the non-corrected data
splitVecs <- list("15-44"=list(seq(15, 29, 1), seq(30, 39, 1), seq(40, 44, 1)),
                  "46-64"=list(seq(45, 49, 1), seq(50, 59, 1), seq(60, 64, 1)),
                  "65-74"=list(seq(65, 69, 1), seq(70, 74, 1)),
                  "75-84"=list(seq(75, 79, 1), seq(80, 84, 1)))

# use Levin IFR to distribute deaths
splitDeaths <- list()
for (s in c(1:length(splitVecs))) {
  sumMort <- NULL
  for (ss in c(1:length(splitVecs[[s]]))) {
    mort <- levinIFR(splitVecs[[s]][[ss]])
    sumMort[ss] <- sum(mort)
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

ooiDeaths_England <- totalDeathsEngland_red - icuDeathsEngland

england1 <- dplyr::filter(countryData, Location=="England" & !is.na(ICU))
england1$Deaths <- totalDeathsEngland_red
england1$ooiDeaths <- ooiDeaths_England
england1$oohDeaths <- NA
england1$ooiDeaths_source <- "Data"
england1$oohDeaths_source <- NA 


####
# Out of hospital deaths England
####
# Redistribute England deaths for Hospital age bins (only requires redistirbuting
# one bin)
ageDeathsEngland <- c("0-14", "15-44", "45-64", "65-74", "75-84", "85+")
deathsTotEngland <- c(6, 536, 4800, 7392, 16226, 21179)
# redistribute deaths to match Hospital bins
# bin fitted to the non-corrected data
ageVec <- seq(15, 44, 1)
mort <- levinIFR(ageVec)

# redistribute the deaths of the second bin according to
# the within bin mortality. Will assign minimal deaths
# to younger bin, and the procedure won't change that much
sumMort <- c(sum(mort[1:3]), sum(mort[4:30]))
splitDeathsBin2 <- round(deathsTotEngland[2]*sumMort/sum(sumMort))

# hosp ages 0-17, 18-64, 65-84, 85+
totalDeathsEngland_red2 <- c(deathsTotEngland[1]+splitDeathsBin2[1],
                       splitDeathsBin2[2]+deathsTotEngland[3],
                       sum(deathsTotEngland[4:5]),
                       deathsTotEngland[6])


# In hospital death data from:
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/
ageHospitalDeaths <- c("0-19", "20-39", "40-59", "60-79", "80+")
hospitalDeathsEngland <- c(20, 208, 2231, 10945, 15389)

# redistribute hospital deaths to same age bins as hospitalizations
# redistribute just uniformly by the number of years in each new bin
hospitalDeathsEngland_red <- c(floor(hospitalDeathsEngland[1]/2),
  ceiling(hospitalDeathsEngland[1]/2)+sum(hospitalDeathsEngland[2:3])+
    round(hospitalDeathsEngland[4]/4),
  round(hospitalDeathsEngland[4]*3/4)+round(hospitalDeathsEngland[5]/2),
  round(hospitalDeathsEngland[5]/2))


# get ooh deaths and severe cases
oohDeaths_England <- totalDeathsEngland_red2 - hospitalDeathsEngland_red

england2 <- dplyr::filter(countryData, Location=="England" & !is.na(Hospitalized))
england2$Deaths <- totalDeathsEngland_red2
england2$ooiDeaths <- NA
england2$oohDeaths <- oohDeaths_England
england2$ooiDeaths_source <- NA
england2$oohDeaths_source <- "Data" 


# Extract death england data, to complete the england dataset
england3 <- dplyr::filter(countryData, Location=="England" & !is.na(Deaths))
england3$ooiDeaths <- NA
england3$oohDeaths <- NA
england3$ooiDeaths_source <- NA
england3$oohDeaths_source <- NA 

england <- rbind(england1, england2, england3)

##############
# Put corrected locations together
##############
locationsList <- list(iceland, newZealand, southKorea, spain,
                      ireland, sweden, france, england, netherlands,
                      atlanta, nyc, ontario, geneva, belgium)

correctedDf <- data.frame()
for (l in c(1:length(locationsList))) {
  correctedDf <- rbind(correctedDf, locationsList[[l]])
#  print(names(locationsList[[l]]))
}

write.csv(correctedDf, "../data/collected_data/locations_serology_data_corrected.csv",
          row.names=FALSE)

