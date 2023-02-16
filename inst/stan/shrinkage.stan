data {

  int N; // Num of genes
  int M; // Num of samples
  int K; // Num of cell types

  int<lower = 1> cell_type[M]; // Cell type label

  int<lower = 0> counts[N, M];

  int<lower = 1> OC[M];

  int<lower = 0, upper = 1> run_estimation; // estimate model vs. generate prior predictive

  simplex[N] p;

  real delta_sd;
  real background_sd;

}


parameters {

  simplex[N] r[K]; // Native count rate

  vector<lower = 0, upper = 0.5>[M] background[N];

  vector<lower = 0, upper = 1>[N] delta[M];

  vector<lower = 0, upper = 1>[M] delta_mean;
  vector<lower = 0, upper = 0.5>[N] background_mean;

  real<lower = 0> tau_a;
  real<lower = 0> tau_b;

}


model {

  tau_a ~ normal(0, delta_sd);
  tau_b ~ normal(0, background_sd);



  for(m in 1:M){
    delta[m] ~ normal(delta_mean[m], tau_a);
  }


  // background_mean ~ normal(0.05, 1e-1);

  for(n in 1:N){
    background[n] ~ normal(background_mean[n], tau_b);
  }


  // Likelihood
  if (run_estimation == 1){

    for(n in 1:N){
      for(m in 1:M){

        counts[n, m] ~ poisson(OC[m]*(1-background[n, m])*delta[m, n]*p[n] +
                               OC[m]*(1-background[n, m])*(1-delta[m, n])*r[cell_type[m],n] +
                               OC[m]*background[n, m]);
      }
    }

  }

}
