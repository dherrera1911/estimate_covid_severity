#################################
#################################
#  
#  This script puts together the information
#  reported in several studies, and by official
#  organisms, of estimates of SARS-CoV-2 infections
#  and age-stratified hospitalizations, ICU admissions
#  and deaths for different locations. The sources are
#  indicated in the script and in the manuscript
#  
#  Daniel Herrera-Esposito gathered the data from the
#  cited sources, and wrote the code.
#  03/07/2021
#  
#  Direct any questions to dherrera@fcien.edu.uy
#  
#  
#################################
#################################
#################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(wpp2019)
library(rriskDistributions)
source("./functions_auxiliary.R")
data(popM)
data(popF)


################
# Iceland
################

# Test and outcome data obtained by mail from the authors of
# "Clinical spectrum of coronavirus disease 2019 in Iceland: population based cohort study"
# BMJ. E Eythorsson et al.
icelandTestAges <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                         "60-69", "70-79", "80+")
Tested_Iceland <- c(38, 113, 251, 222, 300, 279, 195, 64, 18)
Hospitalized_Iceland <- c(1, 0, 4, 7, 8, 17, 25, 19, 13)
ICU_Iceland <- c(0, 0, 0, 1, 2, 6, 10, 5, 0)
deaths_Iceland <- c(0, 0, 0, 0, 0, 0, 2, 3, 3)

# data obtained: 
oohDeaths_Iceland_o65 <- 1 
ooiDeaths_Iceland_o65 <- 4
# all ooh and ooi deaths occured in +65
# evidently, 3 ooi deaths occurred in 80+. There's a remaining ooi death
# that we add to the 70-79 range
# There's one out of hospital death, which we add to the oldest age
oohDeaths_Iceland <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)
ooiDeaths_Iceland <- c(0, 0, 0, 0, 0, 0, 0, 1, 3)

severeCases_Iceland <- Hospitalized_Iceland + oohDeaths_Iceland
criticalCases_Iceland <- ICU_Iceland + ooiDeaths_Iceland

# Seroprevalence data was obtained from
# "Spread of SARS-CoV-2 in the Icelandic Population", NEJM, Gudbjartsson et al
icelandSeroprevAges <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                         "60-69", "70-79", "80+")
icelandSeroprev <- c(0, 0.32, 0.90, 1.01, 1.44, 0.82, 0.49, 0.28, 0)
icelandSeroprevL <- c(0, 0.1, 0.52, 0.63, 0.98, 0.47, 0.19, 0, 0)
icelandSeroprevH <- c(0.27, 0.74, 1.42, 1.53, 2.04, 1.34, 0.99, 1.24, 2.42)

population_Iceland <- extract_country_population(popM=popM, popF=popF,
                                                 countryName="Iceland",
                                                 ageBins=icelandSeroprevAges)

# Estimate the number of cases from seroprevalence, and
# from that, the ratio between testing and seroprevalence cases
seroprevCasesIceland <- round(population_Iceland*icelandSeroprev/100)
seroprevCasesIceland_L <- round(population_Iceland*icelandSeroprevL/100)
seroprevCasesIceland_H <- round(population_Iceland*icelandSeroprevH/100)

seroprevRatioIceland <- pmax(1, seroprevCasesIceland/Tested_Iceland)
seroprevRatioIceland_L <- pmax(1, seroprevCasesIceland_L/Tested_Iceland)
seroprevRatioIceland_H <- pmax(1, seroprevCasesIceland_H/Tested_Iceland)

caseRatioIcelandDf <- data.frame(age=icelandSeroprevAges,
   testedCases=Tested_Iceland,
   ratioSeroprev=pmax(1, seroprevRatioIceland),
   ratioSeroprevL=pmax(1, seroprevRatioIceland_L),
   ratioSeroprevH=pmax(1, seroprevRatioIceland_H))


testedPrevalenceIceland <- Tested_Iceland/population_Iceland*100
indSeroprevL <- which(icelandSeroprevL < testedPrevalenceIceland)
icelandSeroprevL_2 <- icelandSeroprevL
icelandSeroprevL_2[indSeroprevL] <- testedPrevalenceIceland[indSeroprevL]
indSeroprev <- which(icelandSeroprev < testedPrevalenceIceland)
icelandSeroprev_2 <- icelandSeroprev
icelandSeroprev_2[indSeroprev] <- testedPrevalenceIceland[indSeroprev]

# absolute numbers are reported
prevalenceL_Iceland2 <- Tested_Iceland/population_Iceland*100

Iceland <- data.frame(Age=icelandSeroprevAges,
                      Population=population_Iceland,
                      Prevalence=icelandSeroprev_2,
                      PrevalenceL=icelandSeroprevL_2,
                      PrevalenceH=icelandSeroprevH,
                      Severe=severeCases_Iceland,
                      Critical=criticalCases_Iceland,
                      Deaths=deaths_Iceland,
                      Type="Seroprevalence",
                      Location="Iceland",
                      EndPointOutcome="2020-04-04",
                      EndPointCases="2020-04-04")

################
# New Zealand
################
# Test and outcome demographics obtained through a combination of
# the official website, accessed January 3 2021
# https://www.health.govt.nz/our-work/diseases-and-conditions/covid-19-novel-coronavirus/covid-19-data-and-statistics/covid-19-case-demographics
# And by mail to the queries from the site 
# https://www.health.govt.nz/our-work/diseases-and-conditions/covid-19-novel-coronavirus/covid-19-data-and-statistics/covid-19-case-demographics#age-gender
age_NZ <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")
# absolute numbers are reported
Tested_NZ <- c(92, 196, 517, 395, 308, 319, 222, 91, 32, 9)
Severe_NZ <- c(1, 3, 7, 16, 16, 29, 27, 27, 16, 9)
Critical_NZ <-    c(0, 1, 1, 1, 3, 4, 6, 10, 8, 5)
deaths_NZ <- c(0, 0, 0, 0, 0, 2, 3, 7,  8, 5)

population_NZ <- extract_country_population(popM=popM, popF=popF,
                                            countryName="New Zealand",
                                            ageBins=age_NZ)

prevalenceRaw_NZ <- Tested_NZ/population_NZ*100
prevalence_NZ <- prevalenceRaw_NZ*c(seroprevRatioIceland, seroprevRatioIceland[9])
prevalenceL_NZ <- prevalenceRaw_NZ*c(seroprevRatioIceland_L, seroprevRatioIceland_L[9])
prevalenceH_NZ <- prevalenceRaw_NZ*c(seroprevRatioIceland_H, seroprevRatioIceland_H[9])

NewZealand <- data.frame(Age=age_NZ,
                         Population=population_NZ,
                         Prevalence=prevalence_NZ,
                         PrevalenceL=prevalenceL_NZ,
                         PrevalenceH=prevalenceH_NZ,
                         Severe=Severe_NZ,
                         Critical=Critical_NZ,
                         Deaths=deaths_NZ,
                         Type="Testing",
                         Location="New_Zealand",
                         EndPointOutcome="2021-01-13",
                         EndPointCases="2021-01-13")

######################
# South korea
######################
# analysis on severe and critical disease extracted from:
# https://www.thelancet.com/journals/lanwpc/article/PIIS2666-6065(20)30061-4/fulltext
# Has clearly nested statistics, can be useful
age_Korea <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
            "60-69", "70-79", "80+")
Tested_Korea <- c(89, 397, 2174, 780, 1037, 1490, 1007, 517, 312)
Critical_Korea <- c(0, 0, 2, 4, 2, 26, 50, 82, 117)
Severe_Korea <- c(0, 1, 19, 14, 32, 128, 174, 199, 189)
# Deaths from Severe COVID-19 Illness: Risk Factors and Its Burden on Critical Care Resources
deaths_Korea2 <- c(0, 0, 0, 2, 1, 14, 36, 67, 107)

population_Korea <- extract_country_population(popM=popM, popF=popF,
                                            countryName="Republic of Korea",
                                            ageBins=age_Korea)

prevalenceRaw_Korea <- Tested_Korea/population_Korea*100
prevalence_Korea <- prevalenceRaw_Korea*c(seroprevRatioIceland)
prevalenceL_Korea <- prevalenceRaw_Korea*c(seroprevRatioIceland_L)
prevalenceH_Korea <- prevalenceRaw_Korea*c(seroprevRatioIceland_H)
 
Korea <- data.frame(Age=age_Korea,
                    Population=population_Korea,
                    Prevalence=prevalence_Korea,
                    PrevalenceL=prevalenceL_Korea,
                    PrevalenceH=prevalenceH_Korea,
                    Severe=Severe_Korea,
                    Critical=Critical_Korea,
                    Deaths=deaths_Korea2,
                    Type="Testing",
                    Location="South_Korea",
                    EndPointOutcome="2020-04-30",
                    EndPointCases="2020-04-30")

################
# Spain
################
# Spain file with case data from
# https://cnecovid.isciii.es/covid19/#documentaci%C3%B3n-y-datos

age_Spain <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")
spainDailyOutcome <- read.csv("../data/downloaded_datasets/spain/casos_hosp_uci_def_sexo_edad_provres.csv",
                    stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., fecha=lubridate::date(fecha))

spainDeaths <- dplyr::filter(spainDailyOutcome, fecha <= "2020-05-11" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., Deaths=sum(num_def)) %>%
  dplyr::mutate(., proportion=Deaths, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., Deaths=proportion)

spainICU <- dplyr::filter(spainDailyOutcome, fecha <= "2020-05-11" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., ICU=sum(num_uci)) %>%
  dplyr::mutate(., proportion=ICU, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., ICU=proportion)

spainHosp <- dplyr::filter(spainDailyOutcome, fecha <= "2020-05-11" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., Hospitalized=sum(num_hosp)) %>%
  dplyr::mutate(., proportion=Hospitalized, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., Hospitalized=proportion)

Hospitalized_Spain <- spainHosp$Hospitalized
ICU_Spain <- spainICU$ICU
deaths_Spain <- spainDeaths$Deaths

# Spain seroprevalence Pollán et al
# Spain seroprevalence from:
# Prevalence of SARS-CoV-2 in Spain (ENE-COVID): a nationwide,
# population-based seroepidemiological study
# I take adjusted seroprevalence from Levin et al
seroprev_Spain <- c(2.4, 2.9, 3.4, 3.0, 3.6, 3.7, 3.4, 3.1)
seroprevL_Spain <- c(1.5, 2.3, 2.7, 2.4, 3.0, 3.1, 2.8, 2.5)
seroprevH_Spain <- c(3.4, 3.5, 4.1, 3.6, 4.5, 4.7, 4.5, 5.0)

population_Spain <- extract_country_population(popM=popM, popF=popF,
                                               countryName="Spain",
                                               ageBins=age_Spain)

Spain <- data.frame(Age=age_Spain,
                    Population=population_Spain,
                    Prevalence=seroprev_Spain,
                    PrevalenceL=seroprevL_Spain,
                    PrevalenceH=seroprevH_Spain,
                    Severe=Hospitalized_Spain,
                    Critical=ICU_Spain,
                    Deaths=deaths_Spain,
                    Type="Seroprevalence",
                    Location="Spain",
                    EndPointOutcome="2020-05-11",
                    EndPointCases="2020-05-11")


######################
# Ireland
######################
# Outcome data from: July 16th
# https://www.hpsc.ie/a-z/respiratory/coronavirus/novelcoronavirus/casesinireland/epidemiologyofcovid-19inireland/july2020/
age_Ireland_outcome <- c("15-24", "25-34", "35-44", "45-54", "55-64")
deaths_Ireland <- c(1, 5, 12, 23, 65)
ICU_Ireland <- c(5, 15, 36, 91, 127)
Hospitalized_Ireland <- c(75, 198, 274, 448, 497)

# Seroprevalence data from:
# https://www.hpsc.ie/a-z/respiratory/coronavirus/novelcoronavirus/scopi/SCOPI%20report%20preliminary%20results%20final%20version.pdf
# Note that seroprevalence values are shifted 5 years in the final
# dataframe, to match outcome data
age_Ireland_sero <- c("20-29", "30-39", "40-49", "50-59", "60-69")
seroprev_Ireland <- c(2.3, 1.4, 1.8, 1.5, 1.7)
seroprevL_Ireland <- c(0.8, 0.4, 0.7, 0.5, 0.6)
seroprevH_Ireland <- c(5.1, 3.5, 3.7, 3.6, 3.8)

# Apply test characteristics correction to seroprevalence data
sensitivity <- 0.93
specificity <- 1
factor <- 1/(sensitivity + specificity - 1)
seroprev_Ireland <- (seroprev_Ireland/100 + specificity - 1)*factor*100
seroprevL_Ireland <- seroprevL_Ireland*factor
seroprevH_Ireland <- seroprevH_Ireland*factor

population_Ireland <- extract_country_population(popM=popM, popF=popF,
                                                 countryName="Ireland",
                                                 ageBins=c("0-14", age_Ireland_outcome, "65+"))
population_Ireland <- population_Ireland[-c(1, length(population_Ireland))]


Ireland <- data.frame(Age=age_Ireland_outcome,
                      Population=population_Ireland,
                      Prevalence=seroprev_Ireland,
                      PrevalenceL=seroprevL_Ireland,
                      PrevalenceH=seroprevH_Ireland,
                      Severe=Hospitalized_Ireland,
                      Critical=ICU_Ireland,
                      Deaths=deaths_Ireland,
                      Type="Seroprevalence",
                      Location="Ireland",
                      EndPointOutcome="2020-07-16",
                      EndPointCases="2020-07-16")


######################
# Sweden
######################

# ICU data and death ICU data from:
# https://www.icuregswe.org/en/data--results/covid-19-in-swedish-intensive-care/
# On the information panel, look for the "Period" menu on the left
# and for Vardresultat to get deaths
# https://portal.icuregswe.org/siri/report/corona.alderkon
age_Sweden_ICU <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")
ICU_Sweden <- c(5, 7, 73, 91, 231, 527, 599, 386, 78, 1)
ICU_deaths_Sweden <- c(0, 0, 5, 6, 29, 90, 146, 146, 35, 1)

# Age stratified COVID-19 deaths from:
# https://www.socialstyrelsen.se/statistik-och-data/statistik/statistik-om-covid-19/statistik-over-antal-avlidna-i-covid-19/
# Extracted up to week 21
# Section Avlidna med covid-19 veckovis
age_Sweden_deaths <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")
deathsSweden <- c(0, 0, 12, 12, 40, 133, 314, 987, 1894, 1161)

ooiDeaths_Sweden <- deathsSweden - ICU_deaths_Sweden
criticalSweden <- ICU_Sweden + ooiDeaths_Sweden

# Seroprevalence from Levin et al, citation:
# 154. Sweden Public Health Authority. COVID-19 Report for Week 24 - COVID-19 veckorapport vecka 24. 2020.
# Dates of seroprev April 27-May24
age_Sweden_seroprev <- c("0-19", "20-49", "50-69", "70+")
seroprev_Sweden <- c(5.7, 6.5, 4.8, 3.1)
seroprevL_Sweden <- c(4.5, 5.2, 3.6, 2.1)
seroprevH_Sweden <- c(7.0, 7.8, 6.0, 4.1)

# adjust to ICU ages, using the large bin seroprevalence for sub-bins
# (e.g. use seroprevalence of 20-49 for 20-29, 30-39, 40-49
seroprev_Sweden_ICU <- c(rep(seroprev_Sweden[1], 2),
                     rep(seroprev_Sweden[2], 3),
                     rep(seroprev_Sweden[3], 2),
                     rep(seroprev_Sweden[4], 3))
seroprevL_Sweden_ICU <- c(rep(seroprevL_Sweden[1], 2),
                     rep(seroprevL_Sweden[2], 3),
                     rep(seroprevL_Sweden[3], 2),
                     rep(seroprevL_Sweden[4], 3))
seroprevH_Sweden_ICU <- c(rep(seroprevH_Sweden[1], 2),
                     rep(seroprevH_Sweden[2], 3),
                     rep(seroprevH_Sweden[3], 2),
                     rep(seroprevH_Sweden[4], 3))


population_Sweden <- extract_country_population(popM=popM, popF=popF,
                                                countryName="Sweden",
                                                ageBins=age_Sweden_ICU)

Sweden <- data.frame(Age=age_Sweden_ICU,
                     Population=population_Sweden,
                     Prevalence=seroprev_Sweden_ICU,
                     PrevalenceL=seroprevL_Sweden_ICU,
                     PrevalenceH=seroprevH_Sweden_ICU,
                     Severe=NA,
                     Critical=criticalSweden,
                     Deaths=deathsSweden,
                     Type="Seroprevalence_convenience",
                     Location="Sweden",
                     EndPointOutcome="2020-05-24",
                     EndPointCases="2020-05-24")

######################
# Belgium
######################
# Hospital data
# https://covid-19.sciensano.be/fr/covid-19-situation-epidemiologique
age_Belgium_outcome <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")
Hospitalized_proportion_Belgium <- c(1.2, 0.5, 1.8, 4.1, 8.3, 15.4,
                                     17.7, 21.0, 23.4, 6.7)/100
totalHospBelgium <- 16061
Hospitalized_Belgium <- round(Hospitalized_proportion_Belgium*totalHospBelgium)

# Death data:
# Belgian Covid-19 Mortality, Excess Deaths, Number of
# Deaths per Million, and Infection Fatality Rates 
totalDeaths <- 8521
hospitalDeaths <- 4041
age_Belgium_deaths <- c("0-24", "25-44", "45-64", "65-74", "75-84", "85+")
deaths_Belgium <- c(1, 30, 409, 1061, 2144, 5087)

# reshape hospitalizations to fit deaths. Just divide them equally
# within each bin, since bins are small and so should be ~ homogeneous
hB <- Hospitalized_Belgium
Hospitalized_Belgium2 <- round(c(sum(hB[1:2])+hB[3]/2,
                           hB[3]/2+hB[4]+hB[5]/2,
                           hB[5]/2+hB[6]+hB[7]/2,
                           hB[7]/2+hB[8]/2,
                           hB[8]/2+hB[9]/2,
                           hB[9]/2+hB[10]))

# ratio between ICU and hospitalizations for Belgium is
# availale up to 14 June, together with ICU age distribution
# Extrapolate to used date
propHospToICU <- 1696/17628
totICU <- round(propHospToICU*totalHospBelgium)
ICU_proportion_Belgium <- c(0.1, 0.5, 1.1, 3.3, 8.3, 20.9, 25.6, 26.9, 12.9, 0.4)/100
ICU_Belgium <- round(totICU*ICU_proportion_Belgium)
iB <- ICU_Belgium
ICU_Belgium2 <- round(c(sum(iB[1:2])+iB[3]/2,
                           iB[3]/2+iB[4]+iB[5]/2,
                           iB[5]/2+iB[6]+iB[7]/2,
                           iB[7]/2+iB[8]/2,
                           iB[8]/2+iB[9]/2,
                           iB[9]/2+iB[10]))

# data on nursing home deaths is also available, use those to
# correct for ooh deaths. From 
# Belgian Covid-19 Mortality, Excess Deaths, Number of
# Deaths per Million, and Infection Fatality Rates. Molenberghs et al
nursingHomeDeaths_Belgium <- c(0, 0, 58, 508, 952, 4251)
nh_in_hospital <- 1276
nh_out_hospital <- 4494
oohDeaths_Belgium <- round(nursingHomeDeaths_Belgium*nh_out_hospital/
                           (nh_in_hospital+nh_out_hospital))
severeCases_Belgium <- Hospitalized_Belgium2 + oohDeaths_Belgium

# Seroprevalence data:
# Seroprevalence of IgG antibodies against SARS coronavirus 2 in
# Belgium – a serial prospective cross-sectional nationwide study of residual samples
seroprev_Belgium <- c(6.7, 6.6, 6.9, 4.6, 7.8, 14.7)
seroprevL_Belgium <- c(4.7, 4.7, 5.2, 2.6, 4.7, 9.9)
seroprevH_Belgium <- c(9.6, 9.2, 9.2, 8.0, 13.0, 21.8)

population_Belgium <- c(3228894, 2956684, 3080528, 1147009, 690685, 326659)

Belgium <- data.frame(Age=age_Belgium_deaths,
                  Population=population_Belgium,
                  Prevalence=seroprev_Belgium,
                  PrevalenceL=seroprevL_Belgium,
                  PrevalenceH=seroprevH_Belgium,
                  Severe=severeCases_Belgium,
                  Critical=ICU_Belgium2,
                  Deaths=deaths_Belgium,
                  Type="Seroprevalence_convenience",
                  Location="Belgium",
                  EndPointOutcome="2020-05-08",
                  EndPointCases="2020-04-26")


######################
# France
######################
# Hospital data from https://www.santepubliquefrance.fr/regions/ile-de-france/documents/bulletin-regional/2020/covid-19-point-epidemiologique-en-ile-de-france-du-28-mai-2020
age_France_outcome <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")

ilDeFrance_inHosp <- c(24, 16, 82, 151, 313, 669, 1025, 1296, 1520, 944)
ilDeFrance_inICU <- c(5, 4, 16, 18, 45, 128, 161, 142, 23, 5)
ilDeFrance_hospRecov <- c(274, 166, 758, 1754, 2800, 4302, 4748, 4188, 3861, 1835)
ilDeFrance_hospDead <- c(2, 3, 12, 49, 132, 503, 1039, 1665, 2272, 1412)
ilDeFrance_hospTot <- ilDeFrance_inHosp + ilDeFrance_hospDead + ilDeFrance_hospRecov

Hospitalized_France <- c(ilDeFrance_hospTot[3:8], sum(ilDeFrance_hospTot[9:10]))
age_France_outcome <- c(age_France_outcome[3:8], "80+")

# Epidemiologique report specifies age-stratified hospital deaths
# and total deaths in care homes. Since there is no
# age-stratified data on total deaths, the simple correction we are
# applying to other countries (in another script) can't be used here. We need
# to allocate care-home deaths to get total deaths
# We will allocate the care home deaths in France with the same
# age distribution as in Belgium
careHomeDeathsTot_France <- 4405 # deaths in care home
oohDeaths_France <- round(careHomeDeathsTot_France*oohDeaths_Belgium/sum(oohDeaths_Belgium))
oohDeaths_France <- c(0, oohDeaths_France)

# get total deaths and total severe
ilDeFrance_hospDead_2 <- c(ilDeFrance_hospDead[3:8], sum(ilDeFrance_hospDead[9:10]))
ilDeFrance_totalDeath <- ilDeFrance_hospDead_2 + oohDeaths_France
severeCases_ilDeFrance <- Hospitalized_France + oohDeaths_France

# Seroprevalence from Carrat et al. Seroprevalence of SARS-CoV-2 among
# adults in three regions of France following the lockdown and
# associated risk factors: a multicohort study Figure 2
# Seroprevalence endpoint = June 23, beginPoint=May 4
age_France_Seroprev <- c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                          "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
                          "80-84", "85+")
seroprev_France <- c(8.0, 8.5, 16.0, 13.8, 17.1, 14.2, 8.4, 6.1, 8.0, 8.3, 
                        5.3, 5.7, 4.5, 5.7)
seroprevL_France <- c(0.9, 4.2, 12.1, 10.3, 13.5, 11.0, 5.7, 3.2, 4.9, 4.7,
                         3.6, 3.2, 0.9, 0.8)
seroprevH_France <- c(16.1, 12.2, 19.6, 16.8, 20.2, 17.3, 10.5, 8.3, 10.5, 11.0,
                         6.86, 8.3, 8.0, 12.5)

# Average adjacent seroprevalences, to fit demography data
seroprev_France <- colMeans(matrix(seroprev_France, nrow=2, ncol=7))
seroprevL_France <- colMeans(matrix(seroprevL_France, nrow=2, ncol=7))
seroprevH_France <- colMeans(matrix(seroprevH_France, nrow=2, ncol=7))

# Demography https://www.citypopulation.de/en/france/reg/admin/R11__%C3%AEle_de_france/
age_France_pop <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")
ilDeFrancePop <- c(1589856, 1568051, 1667319, 1762912, 1667546, 1537661,
                   1178725, 815274, 418062, 121023)
ilDeFrancePop <- c(ilDeFrancePop[3:8], sum(ilDeFrancePop[9:10]))

Ile_de_France <- data.frame(Age=age_France_outcome,
                         Population=ilDeFrancePop,
                         Prevalence=seroprev_France,
                         PrevalenceL=seroprevL_France,
                         PrevalenceH=seroprevH_France,
                         Severe=severeCases_ilDeFrance,
                         Critical=NA,
                         Deaths=ilDeFrance_totalDeath,
                         Type="Seroprevalence",
                         Location="Ile_de_France",
                         EndPointOutcome="2020-05-26",
                         EndPointCases="2020-06-23")


######################
# England
######################
# population data from Levin et al
age_England_Pop_Levin <- c("0-17", "18-24", "25-34", "35-44", "45-54", "55-64",
                       "65-74", "75+")
population_England_Levin <- c(12023568, 4746616, 7609363, 7147939, 7623273, 6782486,
                        5576066, 4777650)
age_England_Pop <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                     "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                     "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
population_England <- c(3496750, 3135711, 3258677, 2079229, 5267401,
                        3836609, 3683915, 3732161, 4099089, 4100526,
                        3601694, 3183915, 3377162, 2674161, 2178672,
                        1777547, 1338005, 1254688)

# Seroprevalence from:
# 1) SARS-CoV-2 antibody prevalence in England following the first peak of the pandemic
# (supplementary)
# 2) Seroprevalence of SARS-CoV-2 antibodies in children: a prospective multicentre cohort study
age_England_seroprev <- c("0-17", "18-24", "25-34", "35-44", "45-54", "55-64",
                       "65-74", "75+")
seroprev_England <- c(9.2, 7.9, 7.8, 6.1, 6.4, 5.9, 3.2, 3.3)
seroprevL_England <- c(6.2, 7.3, 7.4, 5.7, 6.0, 5.5, 2.8, 2.9)
seroprevH_England <- c(12.2, 8.5, 8.3, 6.6, 6.9, 6.4, 3.6, 3.8)

##########
# England ICU
##########
# ICU and in ICU deaths from:
# "ICNARC report on COVID-19 in critical care 03 July 2020"
# You need to look on google "ICNARC report on COVID-19 in critical
# care 03 July 2020"
age_England_ICU <- c("16-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
totICUEng <- 10287
propICUEng <- c(2.2, 5.8, 13.7, 27.7, 29.7, 18.0, 2.9)
ICU_England <- round(propICUEng/100*totICUEng)

# ICU deaths
icnarcDeathAge <- c("16-39", "40-49", "50-59", "60-69", "70-79", "80+")
icuDeathsEngland <- c(119, 290, 918, 1372, 1069, 176)
# redistribute ICU deaths to match ICU patients bins
icuDeathsEngland_y <- round(icuDeathsEngland[1]*ICU_England[1:2]/sum(ICU_England[1:2]))
icuDeathsEngland <- c(icuDeathsEngland_y, icuDeathsEngland[2:6])

# Total death data from:
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregisteredweeklyinenglandandwalesprovisional/weekending18june2021
# Figure 4
ageDeathsEngland <- c("0-14", "15-44", "45-64", "65-74", "75-84", "85+")
deathsTotEngland <- c(6, 536, 4800, 7392, 16226, 21179)
# redistribute total deaths to match ICU bins
reducedDeathVec <- deathsTotEngland[-c(1,6)]

# Redistribute deaths according to relative risk of death for the ages
# within each bin (taken from Levin). Using Levin within-bin
# balances many relevant factors, since these are expected
# to be relatively homogeneous within bins
levinIFR <- function(age){return(10^(-3.27+0.0524*age))}

splitVecs <- list("15-44"=list(seq(15, 29, 1), seq(30, 39, 1), seq(40, 44, 1)),
                  "46-64"=list(seq(45, 49, 1), seq(50, 59, 1), seq(60, 64, 1)),
                  "65-74"=list(seq(65, 69, 1), seq(70, 74, 1)),
                  "75-84"=list(seq(75, 79, 1), seq(80, 84, 1)))

splitDeaths <- list()
# use Levin IFR to distribute deaths
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
criticalEngland <- ICU_England + ooiDeaths_England

# redistribute population to match ICU bins
population_ICU_England <- c(sum(population_England[4:6]),
                            sum(population_England[7:8]),
                            sum(population_England[9:10]),
                            sum(population_England[11:12]),
                            sum(population_England[13:14]),
                            sum(population_England[15:16]),
                            sum(population_England[17:18]))

England1 <- data.frame(Age=age_England_ICU,
                       Population=population_ICU_England,
                       Prevalence=seroprev_England[-1],
                       PrevalenceL=seroprevL_England[-1],
                       PrevalenceH=seroprevH_England[-1],
                       Severe=NA,
                       Critical=criticalEngland,
                       Deaths=NA,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2020-07-03",
                       EndPointCases="2020-07-13")

##########
# England Hosp
##########
# Hospital data from:
# https://coronavirus.data.gov.uk/details/download
# Hospital 2 https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
# In hospital death data from:
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/

age_England_Hosp <- c("0-17", "18-64", "65-84", "85+")
Hospitalized_England <- c(1224, 35024, 46212, 24476)

# Match total deaths to hospitalization bins
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
severeCases_England <- oohDeaths_England + Hospitalized_England

# Estimate population for the hospitalization bins
population_Hosp_England <- c(population_England_Levin[1],
                             sum(population_England_Levin[2:6]),
                             sum(population_England[14:17]),
                             population_England[18])

# Estimate seroprevalence for the hospitalization bins,
# by estimating cases in sub-bins and rearranging them
age_England_seroprev2 <- c("0-17", "18-24", "25-34", "35-44", "45-54",
                           "55-64", "65-84", "85+")
population_seroprevHosp_England <- c(population_England_Levin[1:6],
                                     sum(population_England[14:17]),
                                     population_England[18])

cases_England <- round(population_seroprevHosp_England*seroprev_England)
casesL_England <- round(population_seroprevHosp_England*seroprevL_England)
casesH_England <- round(population_seroprevHosp_England*seroprevH_England)

seroprev_England2 <- c(cases_England[1], sum(cases_England[2:6]),
                    cases_England[7:8])/population_Hosp_England
seroprevL_England2 <- c(casesL_England[1], sum(casesL_England[2:6]),
                    casesL_England[7:8])/population_Hosp_England
seroprevH_England2 <- c(casesH_England[1], sum(casesH_England[2:6]),
                    casesH_England[7:8])/population_Hosp_England

England2 <- data.frame(Age=age_England_Hosp,
                       Population=population_Hosp_England,
                       Prevalence=seroprev_England2,
                       PrevalenceL=seroprevL_England2,
                       PrevalenceH=seroprevH_England2,
                       Severe=severeCases_England,
                       Critical=NA,
                       Deaths=NA,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2020-07-13",
                       EndPointCases="2020-06-23")

##########
# England Deaths
##########

age_England_death <- c("0-1", "1-14", "15-44", "45-64", "65-74", "75-84", "85+")
deaths_England_Tot_Jul17 <- c(2, 4, 543, 4857, 7492, 16436, 21466)


age_England_death_hosp <- c("0-19", "20-39", "40-59", "60-79", "80+")
deaths_England_Hosp_Jul08 <- c(20, 210, 2247, 11007, 15509)

# I could not find the online data cited in Levin, but the total
# deaths from Levin and from the ONS are off by only 2%, and
# Levin's data match the bins of the serology study
age_England_death_Levin <- c("0-17", "18-24", "25-34", "35-44", "45-54",
                             "55-64", "65-74", "75+")
deaths_England_Levin <- c(11, 30, 131, 394, 1348, 3605, 7631, 38629)

England3 <- data.frame(Age=age_England_death_Levin,
                       Population=population_England_Levin,
                       Prevalence=seroprev_England,
                       PrevalenceL=seroprevL_England,
                       PrevalenceH=seroprevH_England,
                       Severe=NA,
                       Critical=NA,
                       Deaths=deaths_England_Levin,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2020-07-29",
                       EndPointCases="2020-07-13")

England <- rbind(England1, England2, England3)


######################
# Netherlands
######################

# Hospital data
# https://data.rivm.nl/geonetwork/srv/dut/catalog.search#/metadata/2c4357c8-76e4-4662-9574-1deb8a73f724?tab=general

hospDataNL <- read.csv("../data/downloaded_datasets/netherlands/COVID-19_casus_landelijk.csv",
                       stringsAsFactors=FALSE, sep=";") %>%
  as_tibble(.) %>%
  dplyr::filter(., (Hospital_admission=="Yes" | Deceased=="Yes") &
                (lubridate::date(Date_statistics)<="2020-05-11")) %>%
  group_by(., Agegroup) %>%
  summarize(., nHosp=sum((Hospital_admission=="Yes")),
            nDead=sum(Deceased=="Yes"),
            nSevere=sum((Hospital_admission=="Yes" | Deceased=="Yes"))) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Agegroup)) %>%
  dplyr::filter(., Agegroup!="<50")

# Deaths for people under 50 are aggregated in the previous
# dataset, but are disagregated here:
# https://www.rivm.nl/documenten/epidemiologische-situatie-covid-19-in-nederland-11-mei-2020
deathsTot_U50 <- c(0, 1, 3, 10, 25)
hospDataNL$nDead[which(hospDataNL$meanAge<50)] <- deathsTot_U50

age_Netherlands_outcome <- hospDataNL$Agegroup
deaths_Netherlands <- hospDataNL$nDead
Hospitalized_Netherlands <- hospDataNL$nHosp
Severe_Netherlands <- hospDataNL$nSevere

# Seroprevalence estimate first week april:
# Nationwide seroprevalence of SARS-CoV-2 andidentification of risk factors in the general populationof the Netherlands during thefirst epidemic wave
# Newer seroprevalence https://www.rivm.nl/en/pienter-corona-study/results
# Seroprevalence 2 Associations between measures of social distancing and SARS-CoV-2 seropositivity: a nationwide population-based study in the Netherlands
# But data is not extractable

specificity <- 1
sensitivity <- 0.844
age_Netherlands_seroprev <- c("2-12", "13-17", "18-24", "25-39", "40-49",
  "50-59", "60-69", "70+")
nTests <- c(457, 129, 226, 696, 429, 485, 377, 301)
nPositive <- c(4, 1, 12, 24, 11, 8, 7, 7)
rawPrevalence <- nPositive/nTests
rawPrevalenceCI <- binomial_confint(nTests, nPositive)
adjustedPrevalence <- (rawPrevalence/100+specificity-1)/(sensitivity+specificity-1)*100
adjustedPrevalenceL <- rawPrevalenceCI$lower/(sensitivity+specificity-1)
adjustedPrevalenceH <- rawPrevalenceCI$upper/(sensitivity+specificity-1)

# repeat last age to match age vecs
adjustedPrevalence <- c(adjustedPrevalence, rep(adjustedPrevalence[8],2))
adjustedPrevalenceL <- c(adjustedPrevalenceL, rep(adjustedPrevalenceL[8],2))
adjustedPrevalenceH <- c(adjustedPrevalenceH, rep(adjustedPrevalenceH[8],2))

population_Netherlands <- extract_country_population(popM=popM, popF=popF,
                                                countryName="Netherlands",
                                                ageBins=age_Netherlands_outcome)

Netherlands <- data.frame(Age=age_Netherlands_outcome,
                          Population=population_Netherlands,
                          Prevalence=adjustedPrevalence*100,
                          PrevalenceL=adjustedPrevalenceL*100,
                          PrevalenceH=adjustedPrevalenceH*100,
                          Severe=Severe_Netherlands,
                          Critical=NA,
                          Deaths=deaths_Netherlands,
                          Type="Seroprevalence",
                          Location="Netherlands",
                          EndPointOutcome="2020-05-11",
                          EndPointCases="2020-05-11")

######################
# Georgia, USA
######################

# Outcome age distribution "Characteristics and Risk Factors for Hospitalization
# and Mortality among Persons with COVID-19 in Atlanta Metropolitan Area"
# Chishinga et al. Paper data are until 31 May. 
age_Atlanta_outcome <- c("0-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75+")
Hospitalized_Atlanta <- c(18, 62, 79, 115, 172, 180, 260)
ICU_Atlanta <- c(4, 13, 14, 33, 44, 57, 74)
deaths_Atlanta <- c(2, 3, 2, 16, 33, 69, 170)

# Hospital data https://www.fultoncountyga.gov/covid-19/epidemiology-reports
casesAtlantaFulton <- 1398 + round(0.425*399) # atlanta cases + portion of unknown
Hospitalized_AtlantaFulton <- round(0.179 * casesAtlantaFulton)
ICU_AtlantaFulton <- round(0.049 * casesAtlantaFulton)
deaths_AtlantaFulton <- round(0.040 * casesAtlantaFulton)

# use distribution of Chishinga with absolute numbers from governmnet report
Hospitalized_Atlanta <- round(Hospitalized_Atlanta*Hospitalized_AtlantaFulton/sum(Hospitalized_Atlanta))
ICU_Atlanta <- round(ICU_Atlanta*ICU_AtlantaFulton/sum(ICU_Atlanta))
deaths_Atlanta <- round(deaths_Atlanta*deaths_AtlantaFulton/sum(deaths_Atlanta))

# redistribute to match serology ages
Hospitalized_Atlanta <- c(sum(Hospitalized_Atlanta[2:3]),
                          sum(Hospitalized_Atlanta[4:5]),
                          sum(Hospitalized_Atlanta[6:7]))
ICU_Atlanta <- c(sum(ICU_Atlanta[2:3]),
                 sum(ICU_Atlanta[4:5]),
                 sum(ICU_Atlanta[6:7]))
deaths_Atlanta <- c(sum(deaths_Atlanta[2:3]),
                    sum(deaths_Atlanta[4:5]),
                    sum(deaths_Atlanta[6:7]))

# Seroprev: Estimated Community Seroprevalence of SARS-CoV-2 Antibodies —
# Two Georgia Counties, Biggs et al
# remove ranges 0-17 which had 0 positive serology
# samples following Levin et al
# Note: Seroprevalence estimate is from Atlanta at Fulton + Kalb counties,
# but since outcomes are only Fulton, we use the populaiton of Fulton alone to
# estimate number of infections
age_Atlanta_seroprev <- c("18-49", "50-64", "65+")
seroprev_Atlanta <- c(3.3, 4.9, 0.7)
seroprevL_Atlanta <- c(1.6, 1.8, 0.1)
seroprevH_Atlanta <- c(6.4, 12.9, 4.5)

### Correct for test characteristics
sensitivity <- 0.932
specificity <- 0.99
factor <- 1/(sensitivity + specificity - 1)
seroprev_Atlanta <- (seroprev_Atlanta/100 + specificity - 1)*factor*100
seroprevL_Atlanta <- seroprevL_Atlanta*factor
seroprevH_Atlanta <- seroprevH_Atlanta*factor

# Population of Atlanta
# https://data.census.gov/cedsci/table?tid=ACSST5Y2019.S0101&g=1600000US1304000
atlantaPopAge <- c("0-5", "5-9", "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                   "60-64", "65-69", "70-74", "75-79", "80-84", "85+")

atlantaPop <- c(26577, 27559, 22290, 32902, 48054, 55760, 46168,
                       36429, 31238, 29916, 27719, 26051, 21839, 18172,
                       14829, 9696, 6821, 6780)
# 10% of Atlanta is in DeKalb, which is not included in Fulton
# country numbers. adjust population
atlantaPop <- round(atlantaPop*0.9)

atlantaPop <- c(sum(atlantaPop[5:10])+round(atlantaPop[4]*2/5),
                sum(atlantaPop[11:13]),
                sum(atlantaPop[14:18]))

Atlanta <- data.frame(Age=age_Atlanta_seroprev,
                      Population=atlantaPop,
                      Prevalence=seroprev_Atlanta,
                      PrevalenceL=seroprevL_Atlanta,
                      PrevalenceH=seroprevH_Atlanta,
                      Severe=Hospitalized_Atlanta,
                      Critical=ICU_Atlanta,
                      Deaths=deaths_Atlanta,
                      Type="Seroprevalence",
                      Location="Atlanta",
                      EndPointOutcome="2020-05-04",
                      EndPointCases="2020-05-03")

# Remove oldest age, which has invalid prevalence after
# adjusting for test characteristics
Atlanta <- Atlanta[-3,]

######################
# New York City, USA
######################
# New York hospital and death data at 28 April
# https://www1.nyc.gov/site/doh/covid/covid-19-data-archive.page
age_NYC_outcome <- c("0-17", "18-44", "45-64", "65-74", "75+")
Hospitalized_NYC <- c(281, 5872, 14366, 9404, 11391)
deaths_NYC <- c(6, 497, 2742, 3027, 6013)

# Demography https://www.baruch.cuny.edu/nycdata/population-geography/pop-demography.htm
newYorkDemAgeVec <- c("0-5", "5-9", "10-14", "15-19",
                      "20-24", "25-29", "30-34", "35-39", "40-44",
                      "45-49", "50-54", "55-59", "60-64", "65-69",
                      "70-74", "75-79", "80-84", "85+")
population_NYC <- c(553277, 496622, 467016, 466963, 588268, 804436,
                728985, 625351, 550081, 553115, 553489, 530749,
                464246, 388657, 265894, 199912, 139369, 161243)

population_NYC <- c(sum(population_NYC[1:4]), sum(population_NYC[5:9]),
                sum(population_NYC[10:13]), sum(population_NYC[14:15]),
                sum(population_NYC[16:18]))

# Seroprev data from Rosenberg et al Cumulative incidence and diagnosis
# of SARS-CoV-2 infection in New York
age_NYC_seroprev <- c("18-34", "35-44", "45-54", "55+")
seroprev_NYC_source <- c(21.8, 23.4, 26.5, 21.5)
seroprevL_NYC_source <- c(19.2, 20.6, 23.8, 19.6)
seroprevH_NYC_source <- c(24.4, 26.2, 29.2, 23.5)

# interpolated and extrapolated seroprev to match outcome bins
seroprev_NYC <- c(21.8, 22.6, 26.5, 21.5, 21.5)
seroprevL_NYC <- c(19.2, 19.9, 23.8, 19.6, 19.6)
seroprevH_NYC <- c(24.4, 23.5, 29.2, 23.5, 23.5)

NYC <- data.frame(Age=age_NYC_outcome,
                  Population=population_NYC,
                  Prevalence=seroprev_NYC,
                  PrevalenceL=seroprevL_NYC,
                  PrevalenceH=seroprevH_NYC,
                  Severe=Hospitalized_NYC,
                  Critical=NA,
                  Deaths=deaths_NYC,
                  Type="Seroprevalence",
                  Location="NYC",
                  EndPointOutcome="2020-04-28",
                  EndPointCases="2020-04-28")

######################
# Canada
######################
# Toronto data with stratified outcomes:
# https://public.tableau.com/app/profile/tphseu/viz/EpidemiologicalSummaryofCOVID-19Cases/EpiSummary
age_Toronto_outcome <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80-89", "90+")

ICU_Toronto_sinceAug1 <- c(11, 13, 41, 88, 163, 308, 493, 393, 190, 35)
ICU_Toronto_allTime <- c(12, 15, 54, 100, 202, 405, 605, 483, 233, 50)
ICU_Toronto <- ICU_Toronto_allTime - ICU_Toronto_sinceAug1

Hospitalized_Toronto_sinceAug1 <- c(58, 62, 325, 476, 646, 1062, 1238, 1351, 1448, 695)
Hospitalized_Toronto_allTime <- c(65, 67, 362, 548, 758, 1269, 1482, 1627, 1824, 877)
Hospitalized_Toronto <- Hospitalized_Toronto_allTime - Hospitalized_Toronto_sinceAug1 + ICU_Toronto

deaths_Toronto_sinceAug1 <- c(0, 0, 9, 23, 42, 133, 291, 441, 738, 523)
deaths_Toronto_allTime <- c(1, 0, 10, 24, 51, 172, 410, 675, 1207, 950)
deaths_Toronto <- deaths_Toronto_allTime - deaths_Toronto_sinceAug1

# Ontario death data, and non-stratified hospital/ICU data:
# https://covid-19.ontario.ca/covid-19-epidemiologic-summaries-public-health-ontario
# https://covid-19.ontario.ca/data
# https://data.ontario.ca/en/dataset/covid-19-cases-in-hospital-and-icu-by-ontario-health-region
# https://health-infobase.canada.ca/covid-19/epidemiological-summary-covid-19-cases.html
totalHospOntario <- 4672
totalICUOntario <- 1000
totalDeathsOntario <- 2777
ICU_Ontario <- round(ICU_Toronto*totalICUOntario/sum(ICU_Toronto))
ICU_Ontario <- c(ICU_Ontario[-c(9,10)], sum(ICU_Ontario[c(9,10)]))
Hospitalized_Ontario <- round(Hospitalized_Toronto*totalHospOntario/sum(Hospitalized_Toronto))
Hospitalized_Ontario <- c(Hospitalized_Ontario[-c(9,10)], sum(Hospitalized_Ontario[c(9,10)]))

deaths_Ontario_coarse_ages <- c("0-19", "20-39", "40-59", "60-79", "80+")
deaths_Ontario_coarse <- c(1, 11, 117, 744, 1904)

# Redistribute Ontarios death with ontario proportions within bins
withinBinList <- list(c(1,2), c(3,4), c(5,6), c(7,8))
ratiosToronto <- list()
for (i in c(1:length(withinBinList))) {
  binDeaths <- deaths_Toronto[withinBinList[[i]]]
  ratiosToronto[[i]] <- binDeaths/sum(binDeaths)
}

deaths_Ontario <- round(c(deaths_Ontario_coarse[1]*ratiosToronto[[1]],
                    deaths_Ontario_coarse[2]*(ratiosToronto[[2]]+c(-0.1, 0.1)),
                    deaths_Ontario_coarse[3]*ratiosToronto[[3]],
                    deaths_Ontario_coarse[4]*ratiosToronto[[4]],
                    deaths_Ontario_coarse[5]))

# Seroprevalence from July seroprevalence report
age_Ontario_seroprev <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                "60-69", "70-79", "80+")
seroprev_Ontario_male_H <- c(7.1, 4.2, 3.0, 2.8, 2.1, 1.3, 1.7, 2.0, 2.7)
seroprev_Ontario_female_H <- c(4.2, 1.1, 2.8, 2.3, 2.1, 2.1, 0.9, 1.8, 1.9)
seroprev_Ontario_male_L <- c(0.20, 0.6, 0.5, 0.5, 0.2, 0, 0.1, 0.1, 0)
seroprev_Ontario_female_L <- c(0, 0, 0.5, 0.4, 0.1, 0, 0, 0, 0)
seroprev_Ontario_male <- c(2.1, 2.4, 1.8, 1.7, 1.1, 0.2, 0.9, 0.7, 0.8)
seroprev_Ontario_female <- c(0.0, 0.2, 1.7, 1.4, 0.9, 0.7, 0.2, 0.5, 0.4)

population_Ontario_male <- c(760313, 840140, 1099167, 1020536, 901607,
                             1011444, 853563, 528877, 263801)
population_Ontario_female <- c(725357, 805212, 1017927, 1014260, 951829,
                               1030621, 911642, 605684, 392034)

population_Ontario <- population_Ontario_male + population_Ontario_female
propMale_Ontario <- population_Ontario_male/(population_Ontario)
propFeale_Ontario <- population_Ontario_female/(population_Ontario)

seroprev_Ontario <- seroprev_Ontario_male*propMale_Ontario +
  seroprev_Ontario_female*propFeale_Ontario
seroprevL_Ontario <- seroprev_Ontario_male_L*propMale_Ontario +
  seroprev_Ontario_female_L*propFeale_Ontario
seroprevH_Ontario <- seroprev_Ontario_male_H*propMale_Ontario +
  seroprev_Ontario_female_H*propFeale_Ontario

Ontario <- data.frame(Age=age_Ontario_seroprev,
                  Population=population_Ontario,
                  Prevalence=seroprev_Ontario,
                  PrevalenceL=seroprevL_Ontario,
                  PrevalenceH=seroprevH_Ontario,
                  Severe=Hospitalized_Ontario,
                  Critical=ICU_Ontario,
                  Deaths=deaths_Ontario,
                  Type="Seroprevalence_convenience",
                  Location="Ontario",
                  EndPointOutcome="2020-07-31",
                  EndPointCases="2020-07-31")

######################
# Switzerland
######################
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataGE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2020-19" &
                altersklasse_covid19!="Unbekannt")

deathDataGE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2020-19" &
                altersklasse_covid19!="Unbekannt")

Hospitalized_Geneva <- hospDataGE$sumTotal
deaths_Geneva <- deathDataGE$sumTotal

# rearrange hospitalizations and deaths
Hospitalized_Geneva <- c(Hospitalized_Geneva[2],
                         sum(Hospitalized_Geneva[3:5]),
                         Hospitalized_Geneva[6]+Hospitalized_Geneva[7]/2,
                         Hospitalized_Geneva[7]/2+sum(Hospitalized_Geneva[8:9]))
deaths_Geneva <- c(deaths_Geneva[2],
                         sum(deaths_Geneva[3:5]),
                         deaths_Geneva[6]+deaths_Geneva[7]/2,
                         deaths_Geneva[7]/2+sum(deaths_Geneva[8:9]))

# Seroprevalence data from:
# Serology-informed estimates of SARS-CoV-2 infection
# fatality risk in Geneva, Switzerland
age_Geneva_seroprev <- c("5-9", "10-19", "20-49", "50-64", "65+")
infected_Geneva <- c(1200, 6100, 28800, 10300, 5700)
infectedL_Geneva <- c(400, 3900, 21400, 7200, 3200)
infectedH_Geneva <- c(2400, 8800, 37700, 13900, 8800)
population_Geneva <- c(26466, 53180, 219440, 98528, 83574)
seroprev_Geneva <- infected_Geneva/population_Geneva*100
seroprevL_Geneva <- infectedL_Geneva/population_Geneva*100
seroprevH_Geneva <- infectedH_Geneva/population_Geneva*100

Geneva <- data.frame(Age=age_Geneva_seroprev[-1],
                  Population=population_Geneva[-1],
                  Prevalence=seroprev_Geneva[-1],
                  PrevalenceL=seroprevL_Geneva[-1],
                  PrevalenceH=seroprevH_Geneva[-1],
                  Severe=Hospitalized_Geneva,
                  Critical=NA,
                  Deaths=deaths_Geneva,
                  Type="Seroprevalence",
                  Location="Geneva",
                  EndPointOutcome="2020-05-10",
                  EndPointCases="2020-05-06")

####################################
# Put countries together and export
####################################
countriesDf <- dplyr::bind_rows(Iceland, NewZealand, Korea, Spain, Ireland,
  Sweden, Ile_de_France, England, Netherlands, Atlanta, NYC, Ontario,
  Geneva, Belgium) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Age))

write.csv(countriesDf, "../data/collected_data/locations_serology_data.csv",
          row.names=FALSE)


