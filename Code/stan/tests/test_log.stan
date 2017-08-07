// saved as test_log.stan
functions {
  real[] logistic(real t, real[] B, real[] theta, real[] x_r, int[] x_i){
    real dBdt[1];
    dBdt[1] = theta[1]*B[1]*(1 - (B[1]/x_r[1]));
    return dBdt;
  }
}

data {
  int<lower=1> N[1];
  int<lower=1> L; // Num steps
  real B0[1]; //init B state
  real<lower=0> z[L,N[1]]; //measures of B
  real t0; //init value of t
  real ts[L]; //values of t
  real<lower=0> K[1];
}
transformed data {
  int n = N[1];
  real x_r[1] = K;
  int x_i[0];
}
parameters {
  real<lower=0,upper=1> theta[1];
  // vector[1] sigma;
  real<lower=0> sigma;
}
model {
  real y_hat[L,1];
  real z_hat[L,n];
  // theta ~ normal(0,1);
  // sigma ~ cauchy(0,2.5);
  y_hat = integrate_ode_rk45(logistic, B0, t0, ts, theta, x_r, x_i);
  for (t in 1:L) {
    for (i in 1:n) {
      z_hat[t,i] = y_hat[t,1];
      z[t,i] ~ normal(z_hat[t,i], sigma);
    }
  }
}
