// generated with brms 2.17.0
functions {
  /* Efficient computation of the R2D2 prior
   * Args:
   *   z: standardized population-level coefficients
   *   phi: local weight parameters
   *   tau2: global scale parameter
   * Returns:
   *   population-level coefficients following the R2D2 prior
   */
  vector R2D2(vector z, vector phi, real tau2) {
    return z .* sqrt(phi * tau2);
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // concentration vector of the D2 prior
  vector<lower=0>[K-1] R2D2_cons_D2; //SR
  // data for the R2D2 prior
  real<lower=0> R2D2_mean_R2;  // mean of the R2 prior
  real<lower=0> R2D2_prec_R2;  // precision of the R2 prior
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kn = 1; //number of parameters with normal prior - SR
  int Kc = K - 1;
  vector<lower=0>[Kc-Kn] R2D2_cons_D2_updated;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  R2D2_cons_D2_updated = rep_vector(R2D2_cons_D2[1], Kc-Kn); //update dimension of concentraition parameter - SR
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  // local parameters for the R2D2 prior
  vector[Kc-Kn] zb; //SR
  simplex[Kc-Kn] R2D2_phi; //SR
  real Intercept;  // temporary intercept for centered predictors
  vector[Kn] b_n; //population effects with normal prior - SR
  // R2D2 shrinkage parameters
  real<lower=0,upper=1> R2D2_R2;  // R2 parameter
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  vector[Kc-Kn] b_r2d2;  // population-level effects with r2d2 prior SR
  vector[Kc] b; // all population-level effects
  real R2D2_tau2;  // global R2D2 scale parameter
  real lprior = 0;  // prior contributions to the log posterior
  R2D2_tau2 = sigma^2 * R2D2_R2 / (1 - R2D2_R2);
  // compute actual regression coefficients
  b_r2d2 = R2D2(zb, R2D2_phi, R2D2_tau2); //SR
  b = append_row(b_n,b_r2d2); //SR
  lprior += normal_lpdf(Intercept | 0, 1);
  lprior += normal_lpdf(b_n | 0, 1); //SR
  lprior += beta_lpdf(R2D2_R2 | R2D2_mean_R2 * R2D2_prec_R2, (1 - R2D2_mean_R2) * R2D2_prec_R2);
  lprior += exponential_lpdf(sigma | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zb);
  target += dirichlet_lpdf(R2D2_phi | R2D2_cons_D2_updated);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
