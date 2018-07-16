data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=1> Nclusters3; // number of levels for group 3 for random intercepts
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each observation to group 2
  int<lower=1, upper=Nclusters3> NgroupMap3[Ntotal]; // mapping variable to map each observation to group 3
  //int<lower=1> Ncol; // total number of columns in model matrix
  vector[Ntotal] X; // vector of covariate for slope
  int y[Ntotal]; // response variable
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
  real intercept_mean;
  real intercept_sd;
  real slope_mean;
  real slope_sd;
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  real intercept; // population intercept 
  real slope; // population slope
  real<lower=0> sigmaRan1; // random effect standard deviation for group 1
  real<lower=0> sigmaRan2; // random effect standard deviation for group 2
  real<lower=0> sigmaRanSlope1; // random effect slope standard deviation for group 3
  real<lower=1, upper=1000> iSize; // size parameter for the nb distribution
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
  vector[Nclusters2] rGroupsJitter2; // number of random jitters for each level of cluster/group 2
  vector[Nclusters3] rGroupsSlope3; // number of random slope jitters for each level of cluster/group 3
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ntotal] newSlope; // expanded slope vector after adding the random slope jitters
  newSlope = slope + rGroupsSlope3[NgroupMap3];
  // calculate the new fitted values after multiplication of predictor with slope coefficient
  for (i in 1:Ntotal){
    mu[i] = X[i] * newSlope[i];
  }
  // log link inverse i.e. exp
  mu = mu + intercept + rGroupsJitter1[NgroupMap1] + rGroupsJitter2[NgroupMap2];
  mu = exp(mu);
}
model {
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  sigmaRan2 ~ gamma(gammaShape, gammaRate);
  sigmaRanSlope1 ~ gamma(gammaShape, gammaRate);
  intercept ~ normal(intercept_mean, intercept_sd);
  slope ~ normal(slope_mean, slope_sd);
  // random effects sample
  rGroupsJitter1 ~ normal(0, sigmaRan1);
  rGroupsJitter2 ~ normal(0, sigmaRan2);
  rGroupsSlope3 ~ normal(0, sigmaRanSlope1);
  // likelihood function
  y ~ neg_binomial_2(mu, iSize);
}
