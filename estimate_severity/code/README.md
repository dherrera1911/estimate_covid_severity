# Code for "Age-specific rate of severe and critical SARS-CoV-2 infections estimated with multi-country seroprevalence studies"

This is the code used to generate the analysis and figures for the paper
"Age-specific rate of severe and critical SARS-CoV-2 infections estimated with multi-country seroprevalence studies",
Daniel Herrera-Esposito, Gustavo de los Campos *BMC Infectious Diseases* vol 22.
https://link.springer.com/article/10.1186/s12879-022-07262-0

### 1_literature_values_covid.R

This script contains the values of different age-stratified IFR and IHR
estimates extracted from the literature, as well as the values of different
age-stratified hospital lethality reports. It puts the two kind of
data in table format and exports them to the csv files
*/data/collected_data/literature_rates_estimations.csv* and
*/data/collected_data/hospitalized_patient_studies.csv*
respectively.


### 2_countries_values.R

This script contains data gathered about age-stratified seroprevalence
estimates or testing reports for several countries, and the corresponding
age-stratified numbers of cumulative hospitalized, critical care and
dead patients. It also performs some pre-processing, mainly to
match the age bins of different data sources, as well as to
combine different data (e.g. total hospitalization numbers
on the one side with the age distribution in hospitalizations
on the other). It puts these data into the file 
*/data/collected_data/locations_serology_data.csv*

### 3_estimate_hospital_mortality.R

Fits a Bayesian regression model to the proportion of
deaths among hospitalized patients (in general hospital
or in ICU). Uses the model defined in
*3_hospital_mortality.stan*, and the data
from */data/collected_data/hospitalized_patient_studies.csv*.
The fitted models are saved in
*data/processed_data/3_hospital_mortality_fit.RDS*.


### 4_countries_values_correction.R

Apply the correction for out-of-hospital (ooh) and
out-of-ICU (ooi) deaths to the number of severe and
critical cases extracted from the original data
sources. Takes as input the raw data file,
*/data/collected_data/locations_serology_data.csv*
and applies the corresponding corrections
using the fitted in-hospital and in-ICU
mortalities found in 
*data/processed_data/3_hospital_mortality_fit.RDS*.


### 5_estimate_outcome_rate_literature.R

Estimates the rate of severe and critical outcomes
through the indirect method, by taking the ratio between
the IFR values reported in Levin et al,
O'Driscoll et al, and Brazeau et al and the
estimated in-hospital and in-ICU mortality 
fitted in *3_estimate_hospital_mortality.R*.


### 6_estimate_serology_outcome_rates.R

Fit a Bayesian regression model to the proportion of outcomes
(severe disease, critical disease, deaths), taken over
the number of estimated infections from seroprevalence data
(or from testing in some cases). The model is defined
in the file *6_estimate_serology_outcome_rates.stan*,
and uses the data from *locations_serology_data.csv*.
The fitted models are saved in
*data/processed_data/6_serology_fits.RDS*.


### 7_get_changes_in_deaths.R

Obtain the dynamics of death data for some controls.

### 8_children_point_estimate.R

Do the point estimate of children severity discussed
in the supplementary.

### plot_main_figures.R

Obtain the main figures shown in the manuscript

### plot_supplementary_figures.R

Obtain the figures shown in the supplementary materials

