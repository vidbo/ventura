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
   int<lower=1> d_yr[Dd];
   int<lower=1> d_pass[D];
   int<lower=0> d_fry[D];
   int<lower=0> d_juv[D];
   int<lower=1> d_unit[Dd];
   int<lower=1> d_site[Dd];
   real<lower=0> d_len[Dd];
   int<lower=1> d_event[D];
   int<lower=1> d_pos[Dd];
   int<lower=1> d_reps[Dd];
}

transformed data {
  # set the maximum possible count for each snorkel survey
  int<lower=0> d_fry_max[Dd];
  int<lower=0> d_juv_max[Dd];
  int<lower=0> d_fry_min[Dd];
  int<lower=0> d_juv_min[Dd];

  for(d in 1:Dd) {
    d_fry_min[d] <- min(segment(d_fry, d_pos[Dd], d_reps[Dd]));
    d_juv_min[d] <- min(segment(d_juv, d_pos[Dd], d_reps[Dd]));
    d_fry_max[d] <- max(20+d_fry_min[d], d_fry_min[d] * 2);
    d_juv_max[d] <- max(20+d_juv_min[d], d_juv_min[d] * 2);
  }

}

parameters {
  real<lower=-3, upper=3> mu;
  real<lower=0.001, upper=5> sigma;
  vector[Dd] theta;
  #vector[Dd] d_p;
  vector<lower=-7, upper=7>[Dd] lrho[2]; # log-density of fry (1) and juv (2), using marginalized method
}
transformed parameters {
  vector[Dd] rho[2]; # density of fry (1) and juv (2), using marginalized method
  for(i in 1:2) rho[i] <- exp(lrho[i]);
  #for(d in 1:Dd) d_p[d] <- inv_logit(theta[d]);

}

model {
  # hierarchical dive observation prob, per dive event

  # dive counts follow N-mixture model of Royle 2004
  for(d in 1:Dd) {
    theta[d] ~ normal(mu, sigma);
    for(n in d_fry_min[d]:d_fry_max[d]) {  # marginalizing over n should involve sums of anti-logs, since n's are mutually exclusive, not independent
       segment(d_fry, d_pos[Dd], d_reps[Dd]) ~ binomial(n, inv_logit(theta[d]));  # likelihood of replicated counts for each n
       n ~ poisson(rho[1, d] * d_len[d]) ;                        # likelihood of n
    }
    for(n in d_juv_min[d]:d_juv_max[d]) {  # marginalizing over n should involve sums of anti-logs, since n's are mutually exclusive, not independent
       segment(d_juv, d_pos[Dd], d_reps[Dd]) ~ binomial(n, inv_logit(theta[d]));  # likelihood of replicated counts for each n
       n ~ poisson(rho[2, d] * d_len[d]) ;                        # likelihood of n
    }
  }

}


generated quantities {
//   real lratio[J,I];
//   for(j in 1:J) {
//     for(i in 1:I) {
//       lratio[j,i] <- log10(ratio[j,i]);
//     }
//   }
}

  
  
  
  
  
  
  
  
  
  
  