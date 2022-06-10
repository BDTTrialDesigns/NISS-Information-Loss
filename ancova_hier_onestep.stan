data {
  
  int<lower=1> M; //number of trials

  real<lower=0> sigma_nu; //sd of sd of treatment effect across trials
  
  int<lower=1> N[M];  // number of observations
  vector[N[1]] Y;  // response variable main trial
  int<lower=1> K;  // number of population-level effects
  matrix[N[1], K] X;  // population-level design matrix
  
  vector[N[2]] Y0;  // response variable
  matrix[N[2], K] X0;  // population-level design matrix

  matrix[K, 2] priorMS;

}
transformed data {
  int Kc = K - 1;
  matrix[N[1], Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  matrix[N[2], Kc] Xc0;  // centered version of X0 without an intercept
  vector[Kc] means_X0;  // column means of X0 before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
    means_X0[i - 1] = mean(X0[, i]);
    Xc0[, i - 1] = X0[, i] - means_X0[i - 1];
  }
}
parameters {
  vector<lower=0>[M] sigmaj; //residual sd
  
  vector[Kc-1] b;  // population-level effects
  vector[Kc-1] b0;  // population-level effects previous trial
  // temporary intercept for centered predictors
  real Intercept;
  real Intercept0;
  
  vector[M] eta; //params errors
  real mu_b_treat; // mean effect of treatment at 24 
  real<lower=0> nu; // sd of effect of treatment at 24 
  
}
transformed parameters {
  
  vector[M] b_treat = mu_b_treat + nu*eta; //treatment effct borrows information across trials
}
model 
{
  vector[K-1] btot;
  vector[K-1] btot0;
  btot[1] = b_treat[1];
  btot0[1] = b_treat[2];
  btot[2] = b[1];
  btot[3] = b[2];
  btot0[2] = b0[1];
  btot0[3] = b0[2];

// priors including all constants
  target += normal_lpdf(b[1] | priorMS[2,1],priorMS[2,2]);
  target += normal_lpdf(b[2] | priorMS[3,1],priorMS[3,2]);
  target += normal_lpdf(Intercept | priorMS[4,1]+dot_product(means_X, sub_col(priorMS,1,1,3)),priorMS[4,2]);
  target += student_t_lpdf(sigmaj | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);

  target += normal_lpdf(b0[1] | priorMS[2,1],priorMS[2,2]);
  target += normal_lpdf(b0[2] | priorMS[3,1],priorMS[3,2]);
  target += normal_lpdf(Intercept0 | priorMS[4,1]+dot_product(means_X0, sub_col(priorMS,1,1,3)),priorMS[4,2]);

  target += normal_lpdf(mu_b_treat | priorMS[1,1],priorMS[1,2]);
  target += normal_lpdf(nu | 0,sigma_nu) 
  - 1 * normal_lccdf(0 | 0,sigma_nu);;
  target += normal_lpdf(eta| 0,1);
  
  

  // likelihood including all constants

  target += normal_id_glm_lpdf(Y | Xc, Intercept, btot, sigmaj[1]);
  target += normal_id_glm_lpdf(Y0 | Xc0, Intercept0, btot0, sigmaj[2]);
}
//generated quantities {
  // actual population-level intercept
//  real b_Intercept = Intercept - dot_product(means_X, btot);
//  real b_Intercept0 = Intercept0 - dot_product(means_X0, btot0);
//}

