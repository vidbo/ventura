# resource selection function with asymmetric effects

data {
   int<lower=1> M;
   int<lower=1> D;
   int<lower=1> U; 
   int<lower=1> S;
   int<lower=1> Y;
   int<lower=1> m_yr;
   int<lower=1> m_unit;
   int<lower=1> m_site;
   real<lower=0> m_len;
   int<lower=1> d_yr;
   int<lower=1> d_pass;
   int<lower=0> d_fry;
   int<lower=0> d_juv;
   int<lower=1> d_unit;,
   int<lower=1> d_site;

}

#transformed data {
#}

parameters {
  simplex[I] theta;
  real<lower=0, upper=100> mu[J];
  real<lower=0, upper=100> sigma[J];
  real<lower=0, upper=100> L[J];  // maximum value of logistic curve
}
transformed parameters {
  simplex[I] thetu[J];
  vector[I] ratio[J];
  vector[I] ustar[J];
  for(i in 1:I) {
    ratio[3,i] <-  L[3]/(1+exp(-sigma[3]*(cut1[i]-mu[3] )));    // * (logistic_cdf(cut1[i+1], mu[3], sigma[3]) - logistic_cdf(cut1[i], mu[3], sigma[3]))  ;
    ratio[2,i] <-  L[2]/(1+exp(-sigma[2]*(cut1[i]-mu[2] )));    // * (logistic_cdf(cut1[i+1], mu[2], sigma[2]) - logistic_cdf(cut1[i], mu[2], sigma[2]))  ;
    ratio[1,i] <-  L[1]/(1+exp(-sigma[1]*(cut1[i]-mu[1] )));    // * (logistic_cdf(cut1[i+1], mu[1], sigma[1]) - logistic_cdf(cut1[i], mu[1], sigma[1]))  ;
  }
  for(j in 1:J) {
    for(i in 1:I) {
      ustar[j,i] <- ratio[j,i] * theta[i];
  }  }
  for(j in 1:J) {
    for(i in 1:I) {
      thetu[j,i] <- ustar[j,i] / sum(ustar[j]);
    }
  }
}

model {

  avail ~ multinomial(theta);
  for(j in 1:J) {
    used[j] ~ multinomial(thetu[j]);
  }
  
}


generated quantities {
  real lratio[J,I];
  for(j in 1:J) {
    for(i in 1:I) {
      lratio[j,i] <- log10(ratio[j,i]);
    }
  }
}

  
  
  
  
  
  
  
  
  
  
  