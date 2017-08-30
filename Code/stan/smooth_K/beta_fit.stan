functions{
  real[] single_beta(real t, real[] B, real[] theta, real[] x_r, int[] x_i){
  real dXdt[4];
  dXdt[1] = (x_r[4]/x_r[1])*B[1]*(1-((B[1] + theta[1]*B[2])/B[3]));
  dXdt[2] = (x_r[5]/x_r[2])*B[2]*(1-((B[2] + theta[1]*B[1])/B[4]));
  dXdt[3] = x_r[6]*(t);
  dXdt[4] = x_r[6]*(t);
  return dXdt;
  }
}

data {
  int<lower=1> n; //number of observations
  real<lower=0> B0[2]; //initial states of B
  real z[n,4]; //values of each observation
  real t0; //time at B0
  real ts[n]; //time at each observation
  real<lower=0> AC50[2]; //AC50 score
  real<lower=0> rho; //max bcells
  real<lower=0> tau[2]; //tau values for each B
  real<lower=0> t_peak; //time at which K begins decreasing
  real omega; //K waning parameter
}

transformed data {
  real x_r[7];
  int x_i[0];
  x_r[1] = AC50[1];
  x_r[2] = AC50[2];
  x_r[3] = rho;
  x_r[4] = tau[1];
  x_r[5] = tau[2];
  x_r[6] = omega;
  x_r[7] = t_peak;
}

parameters {
  real theta[1];
  real<lower=0> sigma[2];
}

model {
  real z_hat[n,4];
  real B[4];
  B[1] = B0[1];
  B[2] = B0[2];
  B[3] = rho/AC50[1];
  B[4] = rho/AC50[2];
  sigma[1] ~ student_t(3,0,1);
  sigma[2] ~ student_t(3,0,1);
  theta[1] ~ normal(0,1);
  z_hat = integrate_ode_rk45(single_beta, B, t0, ts, theta, x_r, x_i);
  for (i in 1:n){
    z[i,1] ~ normal(z_hat[i,1],sigma[1]*500);
    z[i,2] ~ normal(z_hat[i,2],sigma[2]*500);
  }

}
