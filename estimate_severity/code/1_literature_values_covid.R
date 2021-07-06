############################################
############################################
# This script puts several literature values for
# different covid variables into a tidy format
# and exports them into files.
#
# Written by Daniel Herrera, November 2020
# Contact at dherrera@fcien.edu.uy
############################################
############################################

library(dplyr)
library(tidyr)
source("./functions_auxiliary.R")


##########################################
##########################################
#### Writing down of literature values
##########################################
##########################################

##############################
##############################
##############################
# 1 ) Age segregated IFR and IHR
##############################
##############################
##############################

##############
# IFR
##############

# O'Driscoll et. al. Meta-analysis
ageODriscoll <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
            "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
            "70-74", "75-79", "80+") 
ifrDriscollVec <- c(0.003, 0.001, 0.001, 0.003, 0.006, 0.013, 0.024, 0.040,
         0.075, 0.121, 0.207, 0.323, 0.456, 1.075, 1.674,
         3.203, 8.292)
ifrDriscollLow <- c(0.002, 0.000, 0.001, 0.002, 0.005, 0.011, 0.021,
                    0.034, 0.064, 0.104, 0.177, 0.277, 0.392, 0.921,
                    1.435, 2.744, 7.105)
ifrDriscollHigh <- c(0.004, 0.001, 0.001, 0.003, 0.008, 0.015, 0.028,
                    0.047, 0.087, 0.140, 0.239, 0.373, 0.527, 1.244,
                    1.937, 3.705, 9.593)
ifrDriscoll <- data.frame(Age=ageODriscoll, Proportion=ifrDriscollVec,
                          Proportion_L=ifrDriscollLow,
                          Proportion_H=ifrDriscollHigh,
                          Study="Driscoll", Type="IFR")

# Brazeau et al. Meta-analysis
ageBrazeau <- c("0-9", "10-14", "15-19", "20-24", "25-29", "30-34",
            "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
            "70-74", "75-79", "80-84", "85-89", "90+") 
ifrBrazeauVec <- c(0.01, 0.01, 0.02, 0.03, 0.04, 0.06, 0.10, 0.16, 0.24,
            0.38, 0.60, 0.94, 1.47, 2.31, 3.61, 5.66, 8.86, 17.37)
ifrBrazeauLow <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.01, 0.02,
            0.03, 0.05, 0.10, 0.18, 0.35, 0.65, 1.21, 2.23, 4.06, 9.7)
ifrBrazeauHigh <- c(0.06, 0.11, 0.18, 0.3 , 0.46, 0.71, 1.03, 1.47, 2.03,
             2.74, 3.64, 4.79, 6.27, 8.21, 10.8, 14.3, 19.3, 31.12)
ifrBrazeau <- data.frame(Age=ageBrazeau, Proportion=ifrBrazeauVec,
                         Proportion_L=ifrBrazeauLow,
                         Proportion_H=ifrBrazeauHigh,
                         Study="Brazeau", Type="IFR")

# Levin et al meta-analysis. These values are not the main
# presented, but the ones in the supplementary section with
# a more fine-grained age stratification
ageLevin <- c("0-9", "10-19", "20-29", "30-39", "40-49",
            "50-59", "60-69", "70-79", "80+")
ifrLevinVec <- c(0.001, 0.003, 0.011, 0.035, 0.116, 0.384, 1.27, 4.19, 15.61)
ifrLevinLow <- c(0.00007, 0.002, 0.009, 0.030, 0.101, 0.335, 1.09, 3.45, 12.2)
ifrLevinHigh <- c(0.0013, 0.004, 0.013, 0.042, 0.134, 0.441, 1.49, 5.10, 20.0)
ifrLevin <- data.frame(Age=ageLevin, Proportion=ifrLevinVec,
                       Proportion_L=ifrLevinLow,
                       Proportion_H=ifrLevinHigh,
                       Study="Levin", Type="IFR")

##############
# IHR
##############

# Age stratified Verity et. al.
ageLevin <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69",
            "70-79", "80+")
cfrVerityVec <- c(0.0026, 0.0148, 0.06, 0.146, 0.295, 1.25, 3.99,
                  8.61, 13.4)
ifrVerityVec <- c(0.00161, 0.00659, 0.0309, 0.0844, 0.161, 0.595,
                  1.93, 4.29, 7.80)
severeVerityVec <- c(0, 0.0408, 1.04, 3.43, 4.25, 8.16, 11.8, 16.6, 18.4)
severeVerityL <- c(0, 0.0243, 0.622, 2.04, 2.53, 4.86, 7.01, 9.87, 11.0)
severeVerityH <- c(0, 0.0832, 2.13, 7.0, 8.68, 16.7, 24.0, 33.8, 37.6)

ihrVerity <- data.frame(Age=ageLevin, Proportion=severeVerityVec,
                       Proportion_L=severeVerityL,
                       Proportion_H=severeVerityH,
                       Study="Verity", Type="IHR")

# Salje et al.
ageSalje <- c("0-20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
severeSalje <- c(0.1, 0.5, 1.1, 1.4, 2.9, 5.8, 9.3, 26.2)
severeSaljeL <- c(0.08, 0.3, 0.6, 0.8, 1.6, 3.3, 5.2, 14.8)
severeSaljeH <- c(0.2, 0.8, 1.7, 2.3, 4.7, 9.5, 15.1, 42.7)
severe2ICU <- c(22.2, 11.6, 15.9, 22.2, 27.6, 30.8, 24.9, 5.6)
severe2Death <- c(0.6, 1.1, 1.9, 3.3, 6.5, 12.6, 21.0, 31.6)

ihrSalje <- data.frame(Age=ageSalje, Proportion=severeSalje,
                      Proportion_L=severeSaljeL,
                      Proportion_H=severeSaljeH,
                      Study="Salje", Type="IHR")

# critical = severeSalje * severe2ICU/100
# criticalL = severeSaljeL * severe2ICU/100
# criticalH = severeSaljeH * severe2ICU/100
# 
# icrSalje <- data.frame(Age=agesSalje, Proportion=severeSaljeVec,
#                       Proportion_L=severeSaljeL,
#                       Proportion_H=severeSaljeH,
#                       Study="Salje", Type="IHR")


# For 80+ substitute for death prop, since they don't seem to go into ICU much
#dfSalje$critical[nrow(dfSalje)] <- tail(severeSaljeVec,1) * tail(severe2Death,1)/100
#dfSalje$criticalL[nrow(dfSalje)] <- tail(severeSaljeL,1) * tail(severe2Death,1)/100
#dfSalje$criticalH[nrow(dfSalje)] <- tail(severeSaljeH,1) * tail(severe2Death,1)/100
#
#veritySevereDf <- dplyr::mutate(dfVerity, study="Verity") %>%
#  dplyr::select(., age, severe, severeL, severeH, study)
#saljeSevereDf <- dplyr::mutate(dfSalje, study = "Salje") %>%
#  dplyr::select(., age, severe, severeL, severeH, study)
#
#literatureSevereDf <- rbind(veritySevereDf, saljeSevereDf)
#write.csv(literatureSevereDf, "../data/0_percentage_severe_literature.csv",
#          row.names=FALSE)
#


literatureDf <- rbind(ifrDriscoll, ifrBrazeau, ifrLevin, ihrVerity,
                      ihrSalje) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Age))

write.csv(literatureDf, "../data/collected_data/literature_rates_estimations.csv",
          row.names=FALSE)


###############################
###############################
###############################
# 2) Studies of patients in hospital or ICU
###############################
###############################
###############################

############
############
# ICU mortality
############
############

# overall letality of some studies
criticalFatalityDf <- data.frame(fatalityICU = c(62, 46.4, 44.3),
                                 study = c("Xu", "Veneces", "ICNARC"))

###############
# Letality by age in UK, ICNARC report for 9 September
###############

age_ICNARC <- c("16-39", "40-49", "50-59", "60-69", "70-79", "80+")
# absolute numbers are reported
discharged_ICNARC <- c(734, 1131, 1982, 1689, 786, 137)
deaths_ICNARC <- c(131, 311, 991, 1467, 1145, 191)
ICU_ICNARC <- discharged_ICNARC + deaths_ICNARC
criticalFatalityICNARC <- data.frame(Age=age_ICNARC,
                                     Patients=ICU_ICNARC,
                                     Deaths=deaths_ICNARC,
                                     Study="ICNARC",
                                     Type="ICU",
                                     Location="UK",
                                     EndPoint="2020-09-07")

###############
# Critical patients NY, Cummings et al
###############
age_NYC <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79",
             "80-89", "90+")
# absolute numbers are reported
ICU_NYC <- c(8, 19, 28, 52, 69, 52, 23, 6)
deaths_NYC <- c(0, 5, 6, 18, 22, 28, 17, 5)
criticalFatalityCummings <- data.frame(Age=age_NYC,
                                     Patients=ICU_NYC,
                                     Deaths=deaths_NYC,
                                     Study="Cummings",
                                     Type="ICU",
                                     Location="NYC, USA",
                                     EndPoint="2020-04-28")

##############
# Critical care, France, Belgium, Switzerland, REVA
##############
# Clinical characteristics and day-90 outcomes of 4244 critically ill
# adults with COVID-19: a prospective cohort study
age_REVA <- c("16-39", "40-59", "60-74", "75+")
# absolute numbers are reported
ICU_REVA <- c(220, 1480, 1972, 572)
deaths_REVA <- c(29, 302, 611, 306)
criticalFatalityREVA <- data.frame(Age=age_REVA,
                                     Patients=ICU_REVA,
                                     Deaths=deaths_REVA,
                                     Study="REVA Network",
                                     Type="ICU",
                                     Location="France, Belgium, Switzerland",
                                     EndPoint="2020-05-04")

##############
# Ranzani, Brazil
##############
# Characterisation of the first 250 000 hospital admissions for
# COVID-19 in Brazil: a retrospective analysis of nationwide data
age_Ranzani <- c("20-39", "40-49", "50-59", "60-69", "70-79", "80+")
# absolute numbers are reported
ICU_Ranzani <- c(7512, 9478, 14034, 18058, 16848, 13757)
deaths_ICU_Ranzani <- c(2225, 3503, 6732, 11372, 12257, 10913)
criticalFatalityRanzani <- data.frame(Age=age_Ranzani,
                                     Patients=ICU_Ranzani,
                                     Deaths=deaths_ICU_Ranzani,
                                     Study="Ranzani",
                                     Type="ICU",
                                     Location="Brazil",
                                     EndPoint="2020-08-15")

##############
# Oliveria, Florida
##############
# ICU outcomes and survival in patients with severe COVID-19 in the
# largest health care system in central Florida
age_Oliveira <- c("18-54", "55-64", "65-74", "75+")
# absolute numbers are reported
ICU_Oliveira <- c(41, 36, 34, 20)
deaths_ICU_Oliveira <- c(3, 5, 8, 10)
criticalFatalityOliveira <- data.frame(Age=age_Oliveira,
                                     Patients=ICU_Oliveira,
                                     Deaths=deaths_ICU_Oliveira,
                                     Study="Oliveira",
                                     Type="ICU",
                                     Location="Florida, USA",
                                     EndPoint="2020-05-18")



##############
# Pediatric Critical care, González-Dambrauskas et al
##############
# Since the age of individual patients is given, we disagregated
# them by age
age_Gonzalez <- c("0-10", "11-18")
# absolute numbers are reported
ICU_Gonzalez <- c(12, 5)
deaths_Gonzalez <- c(1, 0)
criticalFatalityGonzalez <- data.frame(Age=age_Gonzalez,
                                     Patients=ICU_Gonzalez,
                                     Deaths=deaths_Gonzalez,
                                     Study="Gonzalez-Dambrauskas",
                                     Type="ICU",
                                     Location="Intercontinental",
                                     EndPoint="2020-04-23")


##############
# Pediatric Critical care, Brazil, Prata-Barbosa et al
##############
age_Prata <- c("0-1", "1-3", "3-5", "5-12", "12-18")
# absolute numbers are reported
ICU_Prata <- c(19, 18, 7, 19, 14)
deaths_Prata <- c(0, 1, 0, 0, 1)
criticalFatalityPrata <- data.frame(Age=age_Prata,
                                     Patients=ICU_Prata,
                                     Deaths=deaths_Prata,
                                     Study="Prata-Barbosa",
                                     Type="ICU",
                                     Location="Brazil_2",
                                     EndPoint="2020-05-31")

############
############
# Hospital mortality
############
############

# Globlal letality reported in different studies
severeFatalityDf <- data.frame(letalityHosp = c(24.9, 21, 25.7, 22, 18.1, 4.4),
                                 study = c("Westblade", "Richardson", "RECOVERY",
                                           "Karagiannidis", "Salje", "Xia"))

###############
# Richardson et al. NY
###############

# Note, "hospitalized" patients here are those that were either
# recovered or dead by the endpoint. Other hospitalized patients were
# still in the hospital by the end of the study
# Length of stay were similar for those discharged alive and those
# who died, but there may be some bias in the data from the
# early cutting point

ages_Richardson <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69",
             "70-79", "80-89", "90+")
deathsMale <- c(0, 0, 3, 6, 19, 40, 56, 91, 94, 28)
deathsFemale <- c(0, 0, 1, 2, 3, 13, 28, 54, 76, 39)
hospitalizedMale <- c(13, 1, 42, 130, 233, 327, 300, 254, 155, 44)
hospitalizedFemale <- c(13, 7, 55, 81, 119, 188, 233, 197, 158, 84)
deaths_Richardson <- deathsFemale + deathsMale
hospitalized_Richardson <- hospitalizedFemale + hospitalizedMale
notDischarged <- c(7, 9, 52, 142, 319, 594, 771, 697, 369, 106)

severeFatalityRichardson <- data.frame(Age=ages_Richardson,
                                       Patients=hospitalized_Richardson,
                                       Deaths= deaths_Richardson,
                                       Study="Richardson",
                                       Type="Hospitalized",
                                       Location="NYC, USA",
                                       EndPoint="2020-04-04")

###############
# Karagiannidis
###############
age_Karagiannidis <- c("18-59", "60-69", "70-79", "80+")
patientsNoVentilator <- c(2474, 1239, 1623, 2958)
patientsVentilator <- c(422, 382, 535, 388)
mortalityNoVentilator <- c(0.7, 5.4, 14.6, 33.8)
mortalityVentilator <- c(27.7, 45.5, 62.6, 72.2)
hospitalized_Karagiannidis <- patientsVentilator + patientsNoVentilator
deaths_Karagiannidis <- round((patientsVentilator*mortalityVentilator +
                    patientsNoVentilator*mortalityNoVentilator)/100)

severeFatalityKaragiannidis <- data.frame(Age=age_Karagiannidis,
                                          Patients=hospitalized_Karagiannidis,
                                          Deaths= deaths_Karagiannidis,
                                          Study="Karagiannidis",
                                          Type="Hospitalized",
                                          Location="Germany",
                                          EndPoint="2020-04-19")

###############
# Salje
###############
age_Salje <- c("0-20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
letality_Salje <- c(0.6, 1.1, 1.9, 3.3, 6.5, 12.6, 21.0, 31.6)
letalityL_Salje <- c(0.2, 0.7, 1.5, 2.9, 6.0, 12.0, 20.3, 30.9)
letalityH_Salje <- c(1.3, 1.6, 2.3, 3.8, 7.0, 13.2, 21.7, 32.4)
# Estimate counts from the Salje intervals
saljeVar <- ((letalityL_Salje-letalityH_Salje)/100/4)^2
saljeNum <- (letality_Salje/100)*(1-letality_Salje/100)
hospitalized_Salje <- round(saljeNum/saljeVar)
deaths_Salje <- round(hospitalized_Salje*letality_Salje/100)
severeFatalitySalje <- data.frame(Age=age_Salje,
                                  Patients=hospitalized_Salje,
                                  Deaths=deaths_Salje,
                                  Study="Salje",
                                  Type="Hospitalized",
                                  Location="France",
                                  EndPoint="2020-05-13")


###############
# Docherty, UK
###############
age_Docherty <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                  "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                  "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+")
discharged_Docherty <- c(83, 21, 12, 24, 40, 53, 99, 159, 222, 309, 415, 459,
                477, 453, 496, 477, 466, 309, 137)
ongoingCare <- c(29, 5, 10, 5, 13, 27, 39, 71, 106, 195, 287, 380, 363,
                 395, 479, 506, 497, 414, 267)
deaths_Docherty <- c(1, 0, 0, 1, 1, 5, 5, 11, 19, 37, 63, 124, 184, 243, 451,
          563, 653, 552, 356)
resolved_Docherty <- deaths_Docherty + discharged_Docherty
severeFatalityDocherty <- data.frame(Age=age_Docherty,
                                  Patients=resolved_Docherty,
                                  Deaths=deaths_Docherty,
                                  Study="Docherty",
                                  Type="Hospitalized",
                                  Location="UK",
                                  EndPoint="2020-04-19")

###############
# Berenguer, Spain
###############
age_Berenguer <- c("0-10", "11-20", "21-30", "31-40", "41-50", "51-60",
                  "61-70", "71-80", "81-90", "91+")
alive_Berenguer <- c(13, 18, 89, 210, 373, 483, 624, 675, 355, 61)
deaths_Berenguer <- c(2, 0, 2, 5, 18, 68, 167, 357, 389, 122)
patients_Berenguer <- alive_Berenguer + deaths_Berenguer
severeFatalityBerenguer <- data.frame(Age=age_Berenguer,
                                  Patients=patients_Berenguer,
                                  Deaths=deaths_Berenguer,
                                  Study="Berenguer",
                                  Type="Hospitalized",
                                  Location="Spain",
                                  EndPoint="2020-04-17")

###############
# Maquilon, Chile
###############
age_Maquilon <- c("1-18", "19-39", "40-49", "50-59", "60-69", "70+")
alive_Maquilon <- c(14, 141, 97, 87, 81, 28)
deaths_Maquilon <- c(0, 3, 2, 4, 8, 37)
patients_Maquilon <- alive_Maquilon + deaths_Maquilon
severeFatalityMaquilon <- data.frame(Age=age_Maquilon,
                                  Patients=patients_Maquilon,
                                  Deaths=deaths_Maquilon,
                                  Study="Maquilon",
                                  Type="Hospitalized",
                                  Location="Chile",
                                  EndPoint="2020-06-04")


##############
# Ranzani, Brazil
##############
# Characterisation of the first 250 000 hospital admissions for
# COVID-19 in Brazil: a retrospective analysis of nationwide data
age_Ranzani <- c("20-39", "40-49", "50-59", "60-69", "70-79", "80+")
# absolute numbers are reported
hospitalized_Ranzani <- c(30603, 33968, 43376, 48270, 41434, 34385)
deaths_Ranzani <- c(3780, 6162, 11818, 20317, 22651, 22787)
severeFatalityRanzani <- data.frame(Age=age_Ranzani,
                                     Patients=hospitalized_Ranzani,
                                     Deaths=deaths_Ranzani,
                                     Study="Ranzani",
                                     Type="Hospitalized",
                                     Location="Brazil",
                                     EndPoint="2020-08-15")


################
# Netherlands Data
################
hospDataNL <- read.csv("../downloaded_data/netherlands/COVID-19_casus_landelijk.csv",
                       stringsAsFactors=FALSE, sep=";") %>%
  as_tibble(.) %>%
  dplyr::filter(., (Hospital_admission=="Yes") &
                (lubridate::date(Date_statistics)<="2020-05-11")) %>%
  group_by(., Agegroup) %>%
  summarize(., nHosp=n(), nDead=sum(Deceased=="Yes")) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Agegroup)) %>%
  dplyr::filter(., meanAge>50)


age_Netherlands <- hospDataNL$Agegroup
deaths_Netherlands <- hospDataNL$nDead
Hospitalized_Netherlands <- hospDataNL$nHosp

severeFatalityNetherlands <- data.frame(Age=age_Netherlands,
                                  Patients=Hospitalized_Netherlands,
                                  Deaths=deaths_Netherlands,
                                  Study="Public_data",
                                  Type="Hospitalized",
                                  Location="Netherlands",
                                  EndPoint="2020-05-11")

#####################
# Put together studies on hospitalized populations
#####################
controledStudies <- rbind(criticalFatalityICNARC, criticalFatalityCummings,
                          criticalFatalityREVA,
                          criticalFatalityRanzani,
                          criticalFatalityGonzalez, criticalFatalityPrata,
                          criticalFatalityOliveira,
                          severeFatalityRichardson, severeFatalityKaragiannidis,
                          severeFatalitySalje, severeFatalityNetherlands,
                          severeFatalityRanzani, severeFatalityDocherty,
                          severeFatalityBerenguer, severeFatalityMaquilon) %>%
  dplyr::mutate(., meanAge=mid_bin_age(Age))

write.csv(controledStudies, "../data/collected_data/hospitalized_patient_studies.csv",
          row.names=FALSE)


