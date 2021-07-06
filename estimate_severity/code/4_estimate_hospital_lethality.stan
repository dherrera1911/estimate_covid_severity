data {
  int<lower=0> N;                       // number of observations
  int<lower=0> K;                       // number of locations
  int<lower=1, upper=K> location[N];    // location index vector
  vector[N] ageVec;                     // predictor
  int cases[N];                         // Case count
  int outcomes[N];                      // outcome count
}

parameters {
  // Countries outcome fit
  real<lower=0> ageSlope;             // mean slope
  real<lower=0> ageSlopeSigma;        // sd of slope across locations
  vector<lower=0>[K] locationSlope;   // vector with slope of each location
  real intercept;                     // mean intercept
  real<lower=0> interceptSigma;       // sd of the intercept
  vector[K] locationIntercept;        // intercept of each location 
}

transformed parameters{
  vector<lower=0, upper=1>[N] outcomeRate;  // variable to contain % outcome
  outcomeRate = inv_logit(locationIntercept[location] + locationSlope[location] .* ageVec);  // logistic regression model 
}

model {
  ageSlope ~ normal(0, 1); 
  ageSlopeSigma ~ exponential(0.5);
  locationSlope ~ normal(ageSlope, ageSlopeSigma);
  intercept ~ normal(-3, 1);
  interceptSigma ~ exponential(0.5);
  locationIntercept ~ normal(intercept, interceptSigma);
  outcomes ~ binomial(cases, outcomeRate);
}

