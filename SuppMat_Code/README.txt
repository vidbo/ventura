This zip file contains three files:

1. README.txt - descriptive of files and details of how to run the R code

2. mv_model_sim.r - contains R code to simulate a dataset, and fit the multivariate poisson/negative-binomial model, as well at the Royle (2004) N-mixture model

3. mv_model.r - contains the R functions required to optimise the N-mixture model, using either the multivariate Poisson/negative-binomial formulation or Royle (2004) formulation


The R code in mv_model_sim.r contains multiple stages:

1. Specify the model (poisson/negative binomial, covariates)

2. Set up parameters for simulation

3. Create a simulation dataset

4. Generate starting values

5. Fit the N-mixture model using the multivariate formulation with optim() 

