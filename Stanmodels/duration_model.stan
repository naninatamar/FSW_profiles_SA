
// 
data {
  int<lower=0> N; // Total number of observations (number of mean SW durations)
  int<lower=1> L; // number of different groups (study populations)
  array[N] real<lower=0> y; // outcome (mean SW durations)
  array[N] int<lower=1, upper=L> ll; // group indicator (study population)
  array[N] real<lower=0> s; // studysizes
  array[N] real x; // model matrix (study year)
  int<lower=0> N_new; // number of new data for predictions
  array[N_new] real x_new; // new data for predctions
  }

parameters{
  real theta_mu; // hyperparameter (pop-mean [random intercept])
  real<lower=0> theta_tau; // hyperparemeter (variance random effect [random intercept])
  real gamma; // parameter (overall slope)
  array[L] real theta; // parameter (group level intercept)
}

transformed parameters{
  array[N] real lambda; // group-level rate parameter of exponential distribution
  array[N] real lambda_overall; // population-level rate parameter of exp distr (excluding random effect variability)
  for (n in 1:N){
    lambda[n] = 1/(exp(theta[ll[n]] + x[n]*gamma)); // expected SW duration 
    lambda_overall[n] =  1/(exp(theta_mu + x[n]*gamma)); // expected SW duration excluding random effect variability
    }
}

model {
  array[N] real Alpha;
  array[N] real Beta;
  theta_mu ~ normal(0,10); // weakly informative hyper-prior distribution
  gamma ~ normal(0,10); // weakly informative prior distribution
  theta_tau ~ exponential(01); // weakly informative hyper-prior distribution
  for (l in 1:L){
      theta[l] ~ normal(theta_mu, theta_tau); 
    }
  for (n in 1:N){
    Alpha[n] = s[n]; 
    Beta[n] = 1/exp(theta[ll[n]] + x[n]*gamma)*s[n];
    target += gamma_lpdf(y[n] | Alpha[n], Beta[n]); 
  }
}

generated quantities{
  vector[N] log_lik;
  vector[N_new] expmeandur_pred;
  vector[N_new] meandur_pred;
  vector[N_new] dur_pred;
  for (n in 1:N){
    log_lik[n] = gamma_lpdf(y[n] | s[n], 1/exp(theta[ll[n]] + x[n]*gamma)*s[n]);
  }
  for (h in 1:N_new){
    expmeandur_pred[h] = exp(theta_mu + x_new[h]*gamma);  ## predictions of expected SW duration without random effect variability
    meandur_pred[h] = exp(normal_rng(theta_mu, theta_tau) + x_new[h]*gamma); ## predictions of expected SW duration including random effect variability
    dur_pred[h] = exponential_rng(1/meandur_pred[h]); ## predictions of individual SW durations
  }
}
