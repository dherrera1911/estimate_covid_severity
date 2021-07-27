####################################
####################################
# 
# 
# This script gets the cumulative reported deaths
# in all locations used, around the date where the outcome
# data was collected.
# 
# 
####################################
####################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(COVID19)
source("./functions_auxiliary.R")

countryData <- read.csv("../data/collected_data/locations_serology_data.csv",
                                                stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::select(., Location, EndPointOutcome, EndPointCases) %>%
  unique(.) %>%
  dplyr::filter(., !(Location=="England" & EndPointOutcome!="2020-07-13"))
countryData$Location[countryData$Location=="South_Korea"] <- "Korea, South"
countryData$Location <- str_replace_all(countryData$Location, "_", " ")

countries <- c("Iceland", "New Zealand", "KOR", "Spain", "Ireland",
"Sweden", "Netherlands", "Belgium")

countriesDeaths <- covid19(country=countries, level=1) %>%
  dplyr::mutate(., Location=administrative_area_level_1) %>%
  dplyr::select(., date, deaths, Location)

englandDeaths <- covid19(country=c("United Kingdom"), level=2) %>%
  dplyr::filter(., !is.na(key) & !(key %in% c("N92000002", "S92000003", "W92000004"))) %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths, na.rm=TRUE)) %>%
  dplyr::mutate(., Location="England")

ilDeFranceDeaths <- covid19(country="France", level=2) %>%
  dplyr::filter(., (!is.na(key)) & key=="REG-11") %>%
  dplyr::mutate(., Location="Ile de France") %>%
  dplyr::select(., date, deaths, Location)

nyDeaths <- covid19(country="USA", level=3) %>%
  dplyr::filter(., administrative_area_level_3=="New York City") %>%
  dplyr::mutate(., Location="NYC") %>%
  dplyr::select(., date, deaths, Location)

atlantaDeaths <- covid19(country="USA", level=3) %>%
  dplyr::filter(., administrative_area_level_3=="Fulton" &
                administrative_area_level_2=="Georgia") %>%
  dplyr::mutate(., Location="Atlanta") %>%
  dplyr::select(., date, deaths, Location)

genevaDeaths <- covid19(country="CHE", level=2) %>%
  dplyr::filter(., key_alpha_2=="GE") %>%
  dplyr::mutate(., Location="Geneva") %>%
  dplyr::select(., date, deaths, Location)

ontarioDeaths <- covid19(country="Canada", level=2) %>%
  dplyr::filter(., administrative_area_level_2=="Ontario") %>%
  dplyr::mutate(., Location="Ontario") %>%
  dplyr::select(., date, deaths, Location)

# put all together
deathsDf <- rbind(countriesDeaths, englandDeaths, ilDeFranceDeaths,
                  nyDeaths, atlantaDeaths, genevaDeaths, ontarioDeaths)

locationList <- str_replace_all(unique(countryData$Location), "_", " ")
centeredDf <- NULL
for (loc in c(1:length(locationList))) {
  locationName <- locationList[loc]
  locationDataDate <- lubridate::date(with(countryData, EndPointCases[Location==locationName]))
  locDeathDf <- dplyr::filter(deathsDf, Location==locationName)
  daysSinceData <- lubridate::interval(locationDataDate, locDeathDf$date) %/%
    lubridate::days(1)
  locDeathDf$daysSinceData <- daysSinceData
  centeredDf <- rbind(centeredDf, locDeathDf)
}

centeredDf <- dplyr::filter(centeredDf, (daysSinceData>=-60) & (daysSinceData<=60))

icelandDf <- dplyr::filter(centeredDf, Location=="Iceland")
icelandDeaths <- c(rep(0, 38), 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                   2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 7, 8, 8,
                   8, 8, 8, 8, 9, 9, 9, rep(10, 46))

icelandDf$deaths <- icelandDeaths

write.csv(centeredDf, file="../data/collected_data/death_dynamics_countries.csv",
          row.names=FALSE)

