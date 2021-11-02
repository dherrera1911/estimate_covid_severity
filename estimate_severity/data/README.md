# Data dictionary

## Collected_data

The *collected_data* directory contains the minimally
preprocessed data that we use as input.

### literature_rates_estimations.csv

Contains estimates of age-stratified outcome rates
among infected individuals from previous studies in the literature.
In particular, the three IFR estimates are meta-analyses based on
serology data.

| Variable | Description |
| -------- | ----------- |
| *Age* | Age stratum |
| *Proportion* | Estimated proportion of disease outcome among infected individuals for the age stratum |
| *Proportion_L* | Lower bound of reported 95 CI (if CI not reported, predictive interval is used instead) |
| *Proportion_H* | Same as above, but for the upper bound |
| *Study* | Variable identifying the study reporting the estimates |
| *Type* | Variable indicating type of outcome. Takes value of IFR for estimates of proportion of death, and IHR for estimate of proportion of severe disease
| *meanAge* | Mean age in the age stratum |


### hospitalized_patient_studies.csv

Contains data reported on the outcomes of patients hospitalized or
in critical care, in studies from different locations. 

| Variable | Description |
| -------- | ----------- |
| *Age* | Age range of the group of patients |
| *Patients* | Number of patients included in the group |
| *Deaths* | Number of deaths among the group patients |
| *Study* | Study from which data were extracted |
| *Type* | Type of medical atention received by patient group: *Hospital* or *ICU* |
| *Location* | Country or region where data were collected |
| *EndPoint* | Final date of the study (may vary between last endpoint date considered, or last date of hospital admission) |
| *meanAge* | Mean age in the age stratum |


### locations_serology_data.csv

Contains the age-stratified estimations of number of infections in
different locations, obtained from serology studies or from extensive
testing and tracing, together with reported data on age-stratified
counts of COVID disease outcomes.
Slight data processing was applied in some cases to match the age
stratifications of seroprevalence and disease outcome data sources.

| Variable | Description |
| -------- | ----------- |
| *Age* | Age range of population |
| *Population* | Number of people in this location and age range |
| *Prevalence* | Estimated seroprevalence |
| *PrevalenceL* | Lower 2.5% confidence bound on the estimated seroprevalence |
| *PrevalenceH* | Same as above, but for the upper 97.5% confidence bound |
| *Severe* | Number of individuals that developed severe disease in this location and age range|
| *Critical* | Number individuals that developed critical disease this location and age range |
| *Deaths* | Number deaths reported in this location and age range |
| *Type* | Method used to obtain the estimate of cases: *Testing*, *Seroprevalence* or *Seroprevalence_convenience* |
| *Location* | Location where the reported data corresponds to |
| *EndPointOutcome* | Last day of the period in which outcomes were counted |
| *EndPointCases* | Last day of the period in which seroprevalence data were gathered |
| *meanAge* | Mean age in the age stratum |



### locations_serology_data_corrected.csv

Contains the serology data, with correction for out-of-hospital and
out-of-ICU deaths in places where data on these quantities was not available.
It has all the same columns as locations_serology_data.csv, and 4 extra columns:

| Variable | Description |
| -------- | ----------- |
| *oohDeaths* | The number of out-of-hospital deaths for this location and age range |
| *ooiDeaths* | The number of out-of-ICU deaths for this location and age range |
| *oohDeaths_source* | Indicates if oohDeaths were obtained from data or estimated |
| *ooiDeaths_source* | Indicates if ooiDeaths were obtained from data or estimated |


### death_dynamics_countries.csv

Dynamics of deaths around data collection for used countries

| Variable | Description |
| -------- | ----------- |
| *id* | Location identifier |
| *date* | Date of the reported deaths | 
| *deaths* | Cumulative deaths at this date | 
| *Location* | Location name | 
| *daysSinceData* | Days since data collection date | 


