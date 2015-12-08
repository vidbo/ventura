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
    d_fry_min[d] <- max(segment(d_fry, d_pos[d], d_reps[d]));
    d_juv_min[d] <- max(segment(d_juv, d_pos[d], d_reps[d]));
    d_fry_max[d] <- max(20+d_fry_min[d], d_fry_min[d] * 2);
    d_juv_max[d] <- max(20+d_juv_min[d], d_juv_min[d] * 2);
  }

}

parameters {
  real<lower=-3, upper=3> mu;
  real<lower=0.001, upper=2> sigma;
  vector[Dd] theta;
  #real<lower=0.01, upper=0.99> d_p;
  vector<lower=-7, upper=7>[Dd] lrho[2]; # log-density of fry (1) and juv (2), using marginalized method
}
transformed parameters {
  vector[Dd] rho[2]; # density of fry (1) and juv (2), using marginalized method
  real LLd;
  vector[Dd] d_p;

  for(i in 1:2) rho[i] <- exp(lrho[i]);
  #for(d in 1:Dd) d_p[d] <- inv_logit(theta[d]);

  # hierarchical dive observation prob, per dive event
  # dive counts follow N-mixture model of Royle 2004
  {
  real LL[Dd];
  for(d in 1:Dd) {
    int j;
    int  R[d_reps[d]];                         # replicate dive counts at this dive event
    real LLn[d_fry_max[d]-d_fry_min[d]+1];     # log-likelihood for each n for fry
    real LLn2[d_juv_max[d]-d_juv_min[d]+1];    # log-likelihood for each n for juv
    d_p[d] <- inv_logit(theta[d]);
    # fry counts
    R <- segment(d_fry, d_pos[d], d_reps[d]);  # extract fry dive counts for this event
    for(n in d_fry_min[d]:d_fry_max[d]) {      # marginalizing over mutually-exclusive n's
       j <- n - d_fry_min[d] + 1;
       LLn[j] <- 0;
       for(r in 1:d_reps[d]) LLn[j] <- LLn[j] + binomial_log(R[r], n, d_p[d]); # likelihood of replicated counts for each n
       LLn[j] <- LLn[j] + poisson_log(n, rho[1, d] * d_len[d]) ;      # likelihood of n
    }
    LL[d] <- log_sum_exp(LLn);                 # sum of (non-log) likelihoods for each discrete n
    
    # juvenile counts
    R <- segment(d_juv, d_pos[d], d_reps[d]);  # extract juv dive counts for this event
    for(n in d_juv_min[d]:d_juv_max[d]) {      # marginalizing over mutually-exclusive n's
       j <- n - d_juv_min[d] + 1;
       LLn2[j] <- 0;
       for(r in 1:d_reps[d]) LLn2[j] <- LLn2[j] + binomial_log(R[r], n, d_p[d]); # likelihood of replicated counts for each n
       LLn2[j] <- LLn2[j] + poisson_log(n, rho[2, d] * d_len[d]) ;      # likelihood of n
    }
    LL[d] <- LL[d] + log_sum_exp(LLn2);         # sum of (non-log) likelihoods for each discrete n
  }
  LLd <- sum(LL);
  }


}

model {
  for(d in 1:Dd) {
    theta[d] ~ normal(mu, sigma);
  }
  increment_log_prob(LLd);   # log-likelihood for dive counts conditional on rho
}


generated quantities {
//   real lratio[J,I];
//   for(j in 1:J) {
//     for(i in 1:I) {
//       lratio[j,i] <- log10(ratio[j,i]);
//     }
//   }
}

  
  
  
  
  
  
  
  
  
  
  