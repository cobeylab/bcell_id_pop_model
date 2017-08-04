// saved as model1.stan

functions {
  real[] model1(real t, real[] B, real[] theta, real[] x_r, int[] x_i){
    real dXdt[2];
    if (t < theta[4])
      dXdt[1] = (((x_r[1]*theta[3])/x_r[3])*B[1])*(1 - (B[1]/(theta[1]*x_r[4])));
      dXdt[2] = theta[2]*B[1] - theta[5]*B[2];
    else
      dXdt[1] = (((x_r[2]*theta[3])/x_r[3])*B[1])*(1 - (B[1]/(theta[1]*x_r[4])));
      dXdt[2] = theta[2]*B[1] - theta[5]*B[2];
    return dXdt
  }
}

data {
  int<lower=1> n; //number of observations
  real<lower=0> B0[1]; //initial state of B
  real z[n,2]; //values of each observation
  real t0; //time at B0
  real ts[n]; //time at each observation
  real<lower=0> Ags; //conc. of antigen before t_peak
  real<lower=0> Agf; //conc. after t_peak
  real<lower=0> AC50; //AC50 score
  real<lower=0> bi; //bcells bound at Ags
}

transformed data {
  real x_r[4];
  int x_i[0];
  x_r[1] = Ags;
  x_r[2] = Agf;
  x_r[3] = AC50;
  x_r[4] = bi;
}

parameters {
  real theta[5]; // M, gamma, tau, t_peak, mu
  real<lower=0> sigma[2];
}

model {
  real z_hat[n,2];
  z_hat = integrate_ode_rk45(model1, B0, t0, ts, theta, x_r, x_i);
  for (i in 1:n){
    z[i,1] ~ normal(z_hat[i,1],sigma[1]);
    z[i,2] ~ normal(z_hat[i,2],sigma[2]);
  }

}
