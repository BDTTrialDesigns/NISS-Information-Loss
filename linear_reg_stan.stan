// script modified from the brsm package
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  real<lower=0, upper=1> alpha; // exponent of the likelihood
  matrix[K, 2] priorMS; // mean and sigma for population level parameter
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  // temporary intercept for centered predictors
  real Intercept;
  real<lower=0> sigma;  // residual SD
}
transformed parameters {
}
model {
  // priors including all constants
  target += normal_lpdf(b[1] | priorMS[1,1],priorMS[1,2]);
  target += normal_lpdf(b[2] | priorMS[2,1],priorMS[2,2]);
  target += normal_lpdf(b[3] | priorMS[3,1],priorMS[3,2]);
  target += normal_lpdf(Intercept | priorMS[4,1]+dot_product(means_X, sub_col(priorMS,1,1,3)),priorMS[4,2]);
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  // likelihood including all constants
  target += alpha*normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
