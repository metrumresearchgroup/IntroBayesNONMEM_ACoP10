data{
  int<lower = 1> nObs;
  real<lower = 0> c24[nObs];
  //  vector[nObs] fxa24;
  real fxa24[nObs];
}

parameters{
  real<lower = 0, upper = 100> emax;
  real<lower = 0> ec50;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
 }

transformed parameters{
  //  vector[nObs] fxa24Hat;
  real fxa24Hat[nObs];

  for(i in 1:nObs){
    fxa24Hat[i] = emax * c24[i]^gamma / (ec50^gamma + c24[i]^gamma);
  }
}

model{
  emax ~ uniform(0, 100);
  ec50 ~ normal(0, 250);
  gamma ~ normal(0, 5);
  sigma ~ cauchy(0, 10);

  fxa24 ~ normal(fxa24Hat, sigma);
}

generated quantities{
  real fxa24Pred[nObs];
  
  for(i in 1:nObs){
    fxa24Pred[i] = normal_rng(fxa24Hat[i], sigma);
  }
}

