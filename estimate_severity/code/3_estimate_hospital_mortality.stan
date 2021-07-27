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

model {
  ageSlope ~ normal(0, 100); 
  //ageSlopeSigma ~ scaled_inv_chi_square(5, 1);
  ageSlopeSigma ~ gamma(4, 4);
  locationSlope ~ normal(ageSlope, ageSlopeSigma);
  intercept ~ normal(0, 100);
  //interceptSigma ~ scaled_inv_chi_square(5, 2);
  interceptSigma ~ gamma(4, 4);
  locationIntercept ~ normal(intercept, interceptSigma);
  outcomes ~ binomial_logit(cases, locationIntercept[location] + locationSlope[location] .* ageVec);
}

