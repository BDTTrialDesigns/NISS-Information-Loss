data {
  
  int<lower=1> M; //number of trials
  
  real mean_stat[M]; //treatment effect 
  real<lower=0> sd_stat[M]; //sd of treatment effect
  real<lower=0> sigma_nu; //sd of sd of treatment effect across trials
  
}
parameters {
  
  real mu_b_treat; // mean change due to treatment at 24 (from 12) 
  real<lower=0> nu; // sd of beta12  and beta24
  vector[M] eta; //params errors
  
}
transformed parameters {

  vector[M] b_treat = mu_b_treat + nu*eta; //treatment effct borrows information across trials
  
}
model 
{
  target += normal_lpdf(mu_b_treat | 0, 10);
  //mu_b_treat ~ normal(0,10);
  // nu ~ normal(0,sigma_nu);
  target += normal_lpdf(nu | 0,sigma_nu) - 1 * normal_lccdf(0 | 0,sigma_nu);
  //eta ~normal(0,1);
  target += normal_lpdf(eta| 0,1);
  //mean_stat ~normal(b_treat,sd_stat);
  target += normal_lpdf(mean_stat| b_treat,sd_stat);
  
}

