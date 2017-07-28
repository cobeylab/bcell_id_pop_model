// saved as 5log_test.stan
functions {
  real[] logistic(real t, real[] B, real[] theta, real[] x_r, int[] x_i){
    real dBdt[x_i[1]];
    for (i in 1:x_i[1]) {
      dBdt[i] = theta[i]*B[i]*(1 - (B[i]/x_r[i]));
    }
    return dBdt;
  }
}
data {
  int<lower=1> N[1];
  int<lower=1> L; // Num steps
  real B0[5]; //init B state
  real<lower=0> z[L,N[1]]; //measures of B
  real t0; //init value of t
  real ts[L]; //values of t
  real<lower=0> K[5];
}
transformed data {
  int n = N[1];
  real x_r[5] = K;
  int x_i[1] = N;
}
parameters {
  real<lower=0> theta[5];
  // vector[1] sigma;
  real<lower=0> sigma[5];
}
model {
  real y_hat[L,n];
  real z_hat[L,n];
  for (i in 1:n) {
    theta[i] ~ uniform(0,1);
    sigma[i] ~ cauchy(0,2.5);
  }
  y_hat = integrate_ode_rk45(logistic, B0, t0, ts, theta, x_r, x_i);
  for (t in 1:L) {
    for (i in 1:n) {
      z_hat[t,i] = y_hat[t,i];
      z[t,i] ~ normal(z_hat[t,i], sigma[i]);
    }
  }
}
