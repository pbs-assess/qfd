data {
  int<lower=1> N;
  vector[N] C;
  vector[N] log_Uobs;
  real<lower=0> init_depletion;
  real log_K_mean;
  real log_q_mean;
  real log_r_mean;
  real log_K_sd;
  real log_q_sd;
  real log_r_sd;
}
parameters {
  real log_K;
  real<lower=0> sigma_pro;
  real<lower=0> sigma_obs;
  real log_q;
  real log_r;
  vector[N] eps_raw; // N(0, 1) process deviations
}
transformed parameters {
  vector[N] B;
  vector[N] log_U;
  real K;
  real r;
  vector[N] eps; // N(0, sigma_pro) process deviations

  K = exp(log_K);
  r = exp(log_r);

  // 'non-centered' parameterization
  // reshapes posterior to help sampler
  eps = eps_raw * sigma_pro; // implies eps ~ N(0, sigma_pro);

  B[1] = init_depletion * K;
  for (i in 2:N) {
    B[i] = B[i - 1] + r * B[i - 1] * (1 - (B[i - 1] / K)) - C[i - 1];
    B[i] = B[i] * exp(eps[i]); // add process error
    if (B[i] < 0.001) B[i] = 0.001;
  }
  log_U = log_q + log(B); // create log index
}
model {
  eps_raw ~ std_normal();  // part of process errors 'non-centered' parameterization
  sigma_pro ~ normal(0, 0.5);
  sigma_obs ~ normal(0, 0.5);
  log_K ~ normal(log_K_mean, log_K_sd);
  log_q ~ normal(log_q_mean, log_q_sd);
  log_r ~ normal(log_r_mean, log_r_sd);
  log_Uobs ~ normal(log_U, sigma_obs); // observation error likelihood
}
