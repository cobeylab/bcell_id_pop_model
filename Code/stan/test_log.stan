// saved as test_log.stan
functions {
  real[] logistic(real t, real[] B, real[] theta, real[] x_r, int[] x_i){
    real dBdt[1];
    dBdt[1] = theta[1]*B[1]*(1 - (B[1]/x_r[1]));
    return dBdt;
  }
}
data {
  int<lower=1> Time; // Num steps
  real B0[1]; //init B state
  real<lower=0> z[Time]; //measures of B
  real t0; //init value of t
  real ts[Time]; //values of t
  real K[1];
}
transformed data {
  real x_r[1] = K;
  int x_i[0];
}
parameters {
  real theta[1];
  vector<lower=0>[1] sigma;
}
model {
  real y_hat[Time,1];
  real z_hat[Time];
  theta ~ normal(0,1);
  sigma ~ cauchy(0,2.5);
  y_hat = integrate_ode_rk45(logistic, B0, t0, ts, theta, x_r, x_i);
  for (t in 1:Time) {
    z_hat[t] = y_hat[t,1];
    z[t] ~ normal(z_hat[t], sigma);
  }
}
