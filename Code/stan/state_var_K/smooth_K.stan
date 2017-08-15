// saved as state_var_K.stan

functions {
  real[] model1(real t, real[] B, real[] theta, real[] x_r, int[] x_i){
    real dXdt[2];
    dXdt[1] = (((theta[2])/x_r[2])*B[1])*(1 - (B[1]/(B[2])));
    dXdt[2] = x_r[4]*t;
    return dXdt;
  }
}

data {
  int<lower=1> n; //number of observations
  real<lower=0> B0; //initial state of B
  real z[n,2]; //values of each observation
  real t0; //time at B0
  real ts[n]; //time at each observation
  real<lower=0> Ags; //conc. of antigen for max specific binding
  real<lower=0> AC50; //AC50 score
  real<lower=0> bi; //bcells bound at Ags
  real omega; //K waning parameter
}

transformed data {
  real x_r[4];
  int x_i[0];
  x_r[1] = Ags;
  x_r[2] = AC50;
  x_r[3] = bi;
  x_r[4] = omega;
}

parameters {
  real m;
  real<lower=0> tau;
  real<lower=0> sigma[2];
}

transformed parameters {
  real theta[2]; // M, tau
  theta[1] = 12500 + 1000*m;
  theta[2] = tau;
}

model {
  real z_hat[n,2];
  real B[2];
  B[1] = B0;
  B[2] = theta[1]*bi;
  sigma[1] ~ student_t(3,0,1);
  sigma[2] ~ student_t(3,0,1);
  m ~ normal(0,1);
  tau ~ normal(0.7,0.5);
  z_hat = integrate_ode_rk45(model1, B, t0, ts, theta, x_r, x_i);
  for (i in 1:n){
    z[i,1] ~ normal(z_hat[i,1],sigma[1]);
    z[i,2] ~ normal(z_hat[i,2],sigma[2]);
  }

}
