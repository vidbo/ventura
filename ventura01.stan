# resource selection function with asymmetric effects

data {
   int<lower=1> M;
   int<lower=1> D;
   int<lower=1> Dd;
   int<lower=1> E;
   int<lower=1> Ee;
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
   int<lower=1> e_yr[Ee];
   int<lower=1> e_pass[E];
   int<lower=0> e_fry[E];
   int<lower=0> e_juv[E];
   int<lower=1> e_unit[Ee];
   int<lower=1> e_site[Ee];
   real<lower=0> e_len[Ee];
   int<lower=1> e_event[E];
   int<lower=1> e_pos[Ee];
   int<lower=1> e_reps[Ee];
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
    d_fry_max[d] <- max(20+d_fry_min[d], d_fry_min[d] * 10);
    d_juv_max[d] <- max(20+d_juv_min[d], d_juv_min[d] * 10);
  }

}

parameters {
  # parameters for dive counts
  real<lower=-4, upper=4> d_mu;         # linear parameters for obs prob.
  real<lower=-4, upper=4> e_mu;         # linear parameters for electrofishing prob.
  real<lower=0, upper=4> d_sigma_s;    # linear parameters for obs prob.
  real<lower=0, upper=4> e_sigma_s;    # linear parameters for electrofishing prob.
  real<lower=0, upper=4> d_sigma_y;    # linear parameters for obs prob.
  real<lower=0, upper=4> e_sigma_y;    # linear parameters for electrofishing prob.
  real d_epsilon_s[S];                  # linear parameters for obs prob.
  real e_epsilon_s[S];                  # linear parameters for electrofishing prob.
  real d_epsilon_y[Y];                  # linear parameters for obs prob.
  real e_epsilon_y[Y];                  # linear parameters for electrofishing prob.
  vector<lower=-7, upper=7>[Dd] d_lrho[2]; # log-density of fry (1) and juv (2), using marginalized method
  vector<lower=-7, upper=7>[Ee] e_lrho[2]; # log-density of fry (1) and juv (2), using marginalized method
}
transformed parameters {
  vector[Dd] d_rho[2]; # density of fry (1) and juv (2), fro dive-count units
  vector[Ee] e_rho[2]; # density of fry (1) and juv (2), from electrofishing units
  real LLd;
  real LLe;
  real d_p[S,Y];    # dive-count observation probabilities
  real e_p[S,Y];     # electrofishing capture probabilities

  ###### Estimation of density from Dive Counts #######
  for(i in 1:2) d_rho[i] <- exp(d_lrho[i]);
  # obs. prob estimated separately for each site-year
  for(s in 1:S) {for(y in 1:Y){ d_p[s,y] <-inv_logit(d_mu + d_epsilon_s[s] + d_epsilon_y[y]); }}
  # dive counts follow N-mixture model of Royle 2004
  {
  real LL[Dd];
  for(d in 1:Dd) {
    int j;
    int  R[d_reps[d]];                         # replicate dive counts at this dive event
    real LLn[d_fry_max[d]-d_fry_min[d]+1];     # log-likelihood for each n for fry
    real LLn2[d_juv_max[d]-d_juv_min[d]+1];    # log-likelihood for each n for juv

    # fry counts
    R <- segment(d_fry, d_pos[d], d_reps[d]);  # extract fry dive counts for this event
    for(n in d_fry_min[d]:d_fry_max[d]) {      # marginalizing over mutually-exclusive n's
       j <- n - d_fry_min[d] + 1;
       LLn[j] <- 0;
       for(r in 1:d_reps[d]) LLn[j] <- LLn[j] + binomial_log(R[r], n, d_p[d_site[d], d_yr[d]]); # likelihood of replicated counts for each n
       LLn[j] <- LLn[j] + poisson_log(n, d_rho[1, d] * d_len[d]) ;      # likelihood of n
    }
    LL[d] <- log_sum_exp(LLn);                 # sum of (non-log) likelihoods for each discrete n
    
    # juvenile counts
    R <- segment(d_juv, d_pos[d], d_reps[d]);  # extract juv dive counts for this event
    for(n in d_juv_min[d]:d_juv_max[d]) {      # marginalizing over mutually-exclusive n's
       j <- n - d_juv_min[d] + 1;
       LLn2[j] <- 0;
       for(r in 1:d_reps[d]) LLn2[j] <- LLn2[j] + binomial_log(R[r], n, d_p[d_site[d], d_yr[d]]); # likelihood of replicated counts for each n
       LLn2[j] <- LLn2[j] + poisson_log(n, d_rho[2, d] * d_len[d]) ;      # likelihood of n
    }
    LL[d] <- LL[d] + log_sum_exp(LLn2);         # sum of (non-log) likelihoods for each discrete n
  }
  LLd <- sum(LL);
  }

  ###### Estimation of density from Depletion Electrofishing #######
  for(i in 1:2) e_rho[i] <- exp(e_lrho[i]);
  # obs. prob estimated separately for each site-year
  for(s in 1:S) {for(y in 1:Y){ e_p[s,y] <-inv_logit(e_mu + e_epsilon_s[s] + e_epsilon_y[y]); }}
  # electrofishing counts follow standard depletion estimator

  {
  real LL[Ee];
  for(e in 1:Ee) {
    int j;
    int  R[e_reps[e]];                         # replicate depletion counts at this EF event

    # fry counts
    R <- segment(e_fry, e_pos[e], e_reps[e]);  # extract fry EF counts for this event
    LL[e] <- 0;
    for(r in 1:e_reps[e]) LL[e] <- LL[e] + poisson_log(R[r], e_rho[1, e] * e_len[e] * pow((1-e_p[e_site[e], e_yr[e]]), (r-1)) * e_p[e_site[e], e_yr[e]]) ;

    # juvenile counts
    R <- segment(e_juv, e_pos[e], e_reps[e]);  # extract juv EF counts for this event
    for(r in 1:e_reps[e]) LL[e] <- LL[e] + poisson_log(R[r], e_rho[2, e] * e_len[e] * pow((1-e_p[e_site[e], e_yr[e]]), (r-1)) * e_p[e_site[e], e_yr[e]]) ;
  }
  LLe <- sum(LL);
  }






}

model {
  for(s in 1:S) {
    d_epsilon_s[s] ~ normal(d_mu, d_sigma_s);
    e_epsilon_s[s] ~ normal(e_mu, e_sigma_s);
  }
  for(y in 1:Y) {
    d_epsilon_y[y] ~ normal(d_mu, d_sigma_y);
    e_epsilon_y[y] ~ normal(e_mu, e_sigma_y);
  }
  increment_log_prob(LLd);   # log-likelihood for dive counts conditional on rho
  increment_log_prob(LLe);   # log-likelihood for electrofishing counts conditional on rho
}


generated quantities {
//   real lratio[J,I];
//   for(j in 1:J) {
//     for(i in 1:I) {
//       lratio[j,i] <- log10(ratio[j,i]);
//     }
//   }
}

  
  
  
  
  
  
  
  
  
  
  