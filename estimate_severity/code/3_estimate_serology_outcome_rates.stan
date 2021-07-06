functions{
  real bin_lpmf(int x, real n, real theta){
    real ans  = lchoose(n, x) + x*log(theta) + (n-x)*log1m(theta);
    return(ans);
  }
}

data {
  int<lower=0> N;                       // number of observations
  int<lower=0> K;                       // number of locations
  int<lower=1, upper=K> location[N];    // location index vector
  vector[N] ageVec;                     // predictor
  vector<lower=0>[N] population;                    // population count
  vector<lower=0>[N] seroprevShape;              // Prevalence gamma shape
  vector<lower=0>[N] seroprevRate;           // Prevalence gamma rate
  int<lower=0> outcomes[N];                      // outcome count
  //vector<lower=0, upper=1>[N] lowerBoundPrev;
}

parameters {
  // Countries outcome fit
  real ageSlope;                      // mean slope
  real<lower=0> ageSlopeSigma;        // sd of slope across locations
  vector[K] locationSlope;            // vector with slope of each location
  real intercept;                     // mean intercept
  real<lower=0> interceptSigma;       // sd of the intercept
  vector[K] locationIntercept;        // intercept of each location 
  vector<lower=0, upper=1>[N] prevalence_raw; 
}

transformed parameters{
  vector<lower=0, upper=1>[N] outcomeRate;  // variable to contain % outcome
  //vector<lower=0, upper=1>[N] prevalence;
  outcomeRate = inv_logit(locationIntercept[location] + locationSlope[location] .* ageVec);  // logistic regression model 
  //prevalence = prevalence_raw + lowerBoundPrev;
}

model {
  ageSlope ~ normal(2, 1); // normal(0, 0.1);
  print("Is this thing on?");
  ageSlopeSigma ~ exponential(0.5);
  locationSlope ~ normal(ageSlope, ageSlopeSigma);
  intercept ~ normal(-6, 2);
  interceptSigma ~ exponential(0.5);
  locationIntercept ~ normal(intercept, interceptSigma);
  prevalence_raw ~ gamma(seroprevShape, seroprevRate);
  for (n in 1:N) {
    outcomes[n] ~ bin(population[n] * prevalence_raw[n], outcomeRate[n]);
  }
}

