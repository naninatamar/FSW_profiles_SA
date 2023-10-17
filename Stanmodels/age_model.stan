
// 
data {
  int<lower=0> N; // Total number of observations (number of mean FSW ages)
  int<lower=1> L; // number of different groups (study populations)
  array[N] real<lower=0> y; // outcome (mean FSW ages)
  array[N] int<lower=1, upper=L> ll; // group indicator (study population)
  array[N] real<lower=0> s; // studysize
  array[N] real x; // model matrix (study year)
  int<lower=0> c; // cut-off for assumed gamma distribution of individual ages
  // real mean_log_coeff_var; // estimated mean of logartithm(coefficient of variation)
  // real<lower=0> sd_log_coeff_var; // estimated standard deviation of logarithm(coefficient of variation)
  int<lower=0> N_sd; // number of reported standard deviations
  array[N_sd] real<lower=0> sd_y; // standard deviation of age (if reported)
  array[N_sd] int<lower=0> ind_sd; // indicator of reported standard deviations of age
  int<lower=0> N_new; // number of new data for predictions
  array[N_new] real x_new; // new data for predctions
  }

transformed data{
  array[N] real<lower=c> y_trans;
  // array[N] real<lower=0> coeff_var;
  array[N_sd] real log_coeff_var;
  real mean_log_coeff_var; 
  real<lower=0> sd_log_coeff_var;
  // real<lower=0> sd_coeff_var;
  for (i in 1:N){
    y_trans[i] = y[i]-c;
  }
  for (j in 1:N_sd){
    log_coeff_var[j] = log(sd_y[j]/y_trans[ind_sd[j]]); ## logarithm of coeff of variation for those study populations that reported SD of age
  }
  mean_log_coeff_var = mean(log_coeff_var); ## mean among study populations (used for informative prior distribution)
  sd_log_coeff_var = sd(log_coeff_var); ## sd among study populations (used for informative prior distribution)
  }
  
parameters{
  array[L] real a; // parameter (group level intercept)
  real b; // paramter (overall slope)
  real log_c_var; // parameter (coefficient of variation)
  real a_mu; // hyperparameter (pop-mean [random intercept])
  real<lower=0> a_tau; //hyperparameter (variance random effect [random intercept])
}

transformed parameters{
  real alpha; // shape of assumed gamma distribution (for individual FSW ages)
  array[N] real beta; // rate of assumed gamma distribution (including random effect variability) 
  array[N] real beta_overall; // rate of assumed gamma distribtuion (excludion random effect variabilty)
  alpha = 1/exp(log_c_var)^2;
  for (n in 1:N){
    beta[n] = 1/(exp(log_c_var)^2 *exp(a[ll[n]] + x[n]*b));
    beta_overall[n] =  1/(exp(log_c_var)^2 *exp(a_mu + x[n]*b));
    }
}

model {
  array[N] real Alpha; // shape of gamma distribution (for mean FSW ages)
  array[N] real Beta; // rate of gamma distribution (for mean FSW ages)
    a_mu ~ normal(0,10); // weakly informative hyper-prior distribution
    b ~ normal(0,10); // weakly informative hyper-prior distribution
    a_tau ~ exponential(01);
    log_c_var ~ normal(mean_log_coeff_var, sd_log_coeff_var); ## informative prior distribution
    for (l in 1:L){
      a[l] ~ normal(a_mu, a_tau); 
    }
  for (n in 1:N){
    Alpha[n] = s[n]* alpha;
    Beta[n] = s[n]* beta[n];
    target += gamma_lpdf(y_trans[n] | Alpha[n], Beta[n]);
    }
  }
generated quantities{
  vector[N] log_lik;
  vector[N_new] expmeanage_pred;
  vector[N_new] meanage_pred;
  vector[N_new] age_pred;
  vector[N_new] beta_pred;
  vector[N_new] expbeta_pred;
  for (n in 1:N){
    log_lik[n] = gamma_lpdf(y_trans[n] | s[n]* alpha, s[n]* beta[n]);
  }
  for (h in 1:N_new){
    expmeanage_pred[h] = exp(a_mu + x_new[h]*b); //  predictions of expected FSW age excluding random effect variability
    expbeta_pred[h] = 1/(exp(log_c_var)^2 * exp(a_mu + x_new[h]*b));  // predictions of rate parameter excluding random effect variability
    beta_pred[h] = 1/(exp(log_c_var)^2 * exp(normal_rng(a_mu, a_tau) + x_new[h]*b)); // predictions of the rate parameter including random effect variability
    meanage_pred[h] = exp(normal_rng(a_mu, a_tau) + x_new[h]*b); // predicitons of expected FSW age including random effect variability
    age_pred[h] = gamma_rng(1/(exp(log_c_var)^2),1/(exp(log_c_var)^2)*1/meanage_pred[h]); // predictions of individual level FSW ages
  }
    for (k in 1:N_new){
    expmeanage_pred[k] = expmeanage_pred[k]+c;
    meanage_pred[k] = meanage_pred[k] + c;
    age_pred[k] = age_pred[k]  + c;
  }
}

