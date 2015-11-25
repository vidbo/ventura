# resource selection function with asymmetric effects

data {
   int<lower=1> M;
   int<lower=1> D;
   int<lower=1> Dd;
   int<lower=1> U; 
   int<lower=1> S;
   int<lower=1> Y;
   int<lower=1> m_yr[M];
   int<lower=1> m_unit[M];
   int<lower=1> m_site[M];
   real<lower=0> m_len[M];
   int<lower=1> d_yr[D];
   int<lower=1> d_pass[D];
   int<lower=0> d_fry[D];
   int<lower=0> d_juv[D];
   int<lower=1> d_unit[D];
   int<lower=1> d_site[D];
   real<lower=0> d_len[D];
   int<lower=1> d_event[D];
   int<lower=1> Nmax;
}

transformed data {
  real log_unif;
  log_unif <- -log(D);  # uniform prior for D dive counts
}

parameters {
  vector<lower=0.01, upper=0.99>[Dd] d_p;
  vector<lower=0, upper=1000>[Dd] rho[2]; # density of fry (1) and juv (2), using non-marginalized method
  # second set for comparison of two estimation methods
  simplex[Dd] d_p2[2];
  vector<lower=0, upper=1000>[Dd] rho2[2]; # density of fry (1) and juv (2), using non-marginalized method
}
transformed parameters {
//   simplex[I] thetu[J];
//   vector[I] ratio[J];
  vector[D] lp_fry;   # log-probability of fry counts
  lp_fry <- rep_vector(log_unif, D);
  for(d in 1:D)
    for(n in d_fry[d]:Nmax)  # marginalizing over n should involve sums of non-logs, since n's are mutually exclusive, not independent
      lp_fry[d] <-  lp_fry[d] 
                + poisson_log(n, rho[1, d_event[d]] * d_len[d]) 
                + binomial_log(d_fry[d], n, d_p[d_event[d]])  ;

}

model {


  for(d in 1:D) {
    d_fry[d] ~ poisson(rho2[1, d_event[d]] * d_len[d] * d_p2[1, d_event[d]]);
    d_juv[d] ~ poisson(rho2[2, d_event[d]] * d_len[d] * d_p2[1, d_event[d]]);
  }

  increment_log_prob(log_sum_exp(lp_fry)); # these are independent and should be a sum of logs (product of non-logs)
  
}


generated quantities {
//   real lratio[J,I];
//   for(j in 1:J) {
//     for(i in 1:I) {
//       lratio[j,i] <- log10(ratio[j,i]);
//     }
//   }
}

  
  
  
  
  
  
  
  
  
  
  