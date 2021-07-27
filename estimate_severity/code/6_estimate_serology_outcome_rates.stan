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
  vector<lower=0>[N] seroprevShape1;           // Prevalence beta shape1
  vector<lower=0>[N] seroprevShape2;           // Prevalence beta shape2
  int<lower=0> outcomes[N];                      // outcome count
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
  outcomeRate = inv_logit(locationIntercept[location] + locationSlope[location] .* ageVec);  // logistic regression model 
}

model {
  ageSlope ~ normal(0, 100);
  //ageSlopeSigma ~ scaled_inv_chi_square(1, 1);
  ageSlopeSigma ~ gamma(4, 4);
  locationSlope ~ normal(ageSlope, ageSlopeSigma);
  intercept ~ normal(0, 100);
  //interceptSigma ~ scaled_inv_chi_square(1, 1);
  interceptSigma ~ gamma(4, 4);
  locationIntercept ~ normal(intercept, interceptSigma);
  prevalence_raw ~ beta(seroprevShape1, seroprevShape2);
  for (n in 1:N) {
    outcomes[n] ~ bin(population[n] * prevalence_raw[n], outcomeRate[n]);
  }
}

