data {
  
  // Poll data
  int N;
  int id_cand;
  array[N] int tot_eff;
  array[N] real<lower=0,upper=1> vshare_raw;
  array[N] real sum_share;
  array[N] int rounding_ind;
  array[N] int r_0;
  int N_0;
  array[N] int r_1;
  int N_1;
  array[N] int r_2;
  int N_2;
  array[N] int r_3;
  int N_3;
  array[N] int r_4;
  int N_4;
  array[N] int r_5;
  int N_5;
  array[N] int id_poll;
  int P;
  array[N] int id_date;
  array[N] int id_house;
  int F;
  matrix[N,3] X;
  
  // Splines
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int num_basis;
  int D;
  matrix[num_basis,D] S;
  
}
parameters {
  
  // Splines
  real alpha0;
  array[num_basis] real alpha_raw;
  real<lower=0> tau_alpha;
  
  // Covariates
  array[P] real mu;
  array[F] real lambda;
  array[3] real beta;
  array[2] real nu;
  real<lower=0> tau_mu;
  real<lower=0> tau_lambda;
  
  // Rounding error
  array[N_1] real<lower=-0.005,upper=0.005> epsilon1;
  array[N_2] real<lower=-0.01,upper=0.01> epsilon2;
  array[N_3] real<lower=0.001,upper=0.015> epsilon3;
  array[N_4] real<lower=0.001,upper=0.01> epsilon4;
  array[N_5] real<lower=0.001,upper=0.005> epsilon5;
  
}
transformed parameters {
  
  array[num_basis] real alpha;
  array[N] real<lower=0,upper=1> theta;
  array[N] real<lower=-.015,upper=.015> epsilon;
  array[N] real<lower=0,upper=1> vshare;
  
  // Spline coefficients, specified as a random walk
  // to avoid overfit
  alpha[1] = alpha_raw[1];
  for (i in 2:num_basis)
    alpha[i] = alpha[i-1] + alpha_raw[i] * tau_alpha;
  
  
  for (i in 1:N) {
    
    // Theta
    theta[i] = inv_logit(alpha0 * id_date[i] + to_row_vector(alpha) * S[,id_date[i]] + // Spline
                                tau_mu * mu[id_poll[i]] + // Poll effect
                                tau_lambda * lambda[id_house[i]] + // House effect
                                X[i,1] * (beta[1] + nu[1] * (id_date[i] - 1)) + // Population definition
                                X[i,2] * (beta[2] + nu[2] * (id_date[i] - 1)) +
                                X[i,3] * beta[3]);
    
    
    // Rounding error
    if (rounding_ind[i] == 0) {
      epsilon[i] = epsilon1[r_0[i]];
    } else if (rounding_ind[i] == 1) {
      epsilon[i] = epsilon1[r_1[i]];
    } else if (rounding_ind[i] == 2) {
      epsilon[i] = epsilon2[r_2[i]];
    } else if (rounding_ind[i] == 3) {
      epsilon[i] = epsilon3[r_3[i]];
    } else if (rounding_ind[i] == 4) {
      epsilon[i] = epsilon4[r_4[i]];
    } else {
      epsilon[i] = epsilon5[r_5[i]];
    }
    vshare[i] = (vshare_raw[i] + epsilon[i]) / sum_share[i];
    
    
  }
  
}
model {
  
  
  // Priors
  
    // Splines
    alpha0 ~ normal(0, 1);
    alpha_raw ~ normal(0, 2);
    tau_alpha ~ student_t(3, 0, 1);
    
    // Poll effect
    mu ~ normal(0, 1);
    tau_mu ~ student_t(3, 0, 1);
    
    // House effect
    lambda ~ normal(0, 1);
    tau_lambda ~ student_t(3, 0, 1);
    
    // Population definition and poll type effects
    beta ~ normal(0, 1);
    nu ~ normal(0, 1);
  
  
  // Likelihood
  for (i in 1:N)
    target += beta_proportion_lpdf(vshare[i] | theta[i], tot_eff[i]);
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  array[D] real<lower=0,upper=1> prob;
  for (d in 1:D)
      prob[d] = inv_logit(alpha0 * d + to_row_vector(alpha) * S[,d]);
      
}
