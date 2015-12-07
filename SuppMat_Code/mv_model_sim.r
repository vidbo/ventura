# This file contains code to simulate an example dataset, and fit the N-mixture model using either formulation
 
# Set working directory, contatining mv_model.r 
 setwd(...)
# Load file containing functions to optimise the model
source(file="mv_model.r")

#=====================================================
# Specify model 
#=====================================================
# Specify mixing distribution to be Poisson "pois", or negative-binomial "nb"
	model <- "pois"
# Specify whether lambda is constant spatially "const" or varies with respect to a covariate "covar"
	lambda.cov <- "const" 
# Specify whether p is constant "const" or varying with respect to a covariate "covar" (which can vary with respect to space and/or time)
	p.cov <- "const"
  
#=====================================================
# Set up parameters values for simulation
#=====================================================
# Specify lambda
	lambda.int <- log(5)
	#lambda.slope <- 1 # If lambda.cov = "covar"
# Specify p
	p.int <- logit(0.25)
	#p.slope <- 1 # If p.cov = "covar"
# Number of sites
	nsites <- 20
# Number of sampling occasions
    nv <- 3
# Dispersion parameter
	#alpha <- 1 # If model = "nb"
# Simulate covariates if required
	#lambda.var <- runif(nsites,0,1)
	#p.var <- matrix(runif(nsites*nv,0,1),nrow=nsites,ncol=nv)
# Define lambda
	switch(lambda.cov,
		covar = {lambda <- exp(lambda.int +lambda.slope*lambda.var)},
		const = {lambda <- exp(lambda.int)})
# Define p
	switch(p.cov,
		covar = {p <- expit(p.int + p.slope*p.var)},
		const = {p <- matrix(expit(p.int),nrow = nsites, ncol = nv)})
	
#=====================================================
# Create a simulation dataset
#=====================================================		
# Simulate true abundance (N)	
	switch(model,
		pois = {true_ab <- rpois(nsites,lambda)},
		nb = {true_ab <- rnbinom(nsites, mu = lambda, size = alpha)})
# Simulate observed abundance (y matrix)
	yMat <- matrix(rbinom(n = nsites*nv, size = true_ab , prob = p), nrow = nsites, ncol = nv)
	
#=====================================================
# Generate starting values
#=====================================================		
	# Evaluate the diagnostic(s)
	m1 <- mean(yMat)
	m2 <- mean(yMat^2)
	m12 <- covstarm(yMat)+m1^2
	diag1 <- covstarm(yMat)
	diag2 <- m1 - m2 + m12
	# Evaluate MOM estimates 
	ptilde <- (m1-m2+m12)/m1
	ltilde <- m1/ptilde	
	stilde <- (m12-m1^2)*ptilde^2
	switch(model,
		pois = {if(diag1 <= 0){
				cat("Lambda is infinite, simulate an alternative dataset","\n")
			} else {
				# Evaluate MOM estimates and use as starting values	
				p.int.st <- logit(ptilde)
				lambda.int.st <- log(ltilde)
			}
			alpha.st <- NULL
			},
		nb = {if(ptilde >= 0 & ltilde >= 0){
					p.int.st <- logit(ptilde)
					lambda.int.st <- log(ltilde)
					alpha.st <- log(stilde)}
			  else { # If any negative-binomial covariance diagnostics are negative, use random starting values
					p.int.st <- logit(runif(1,.2,.8))
					lambda.int.st <- log(mean(apply(yMat,1,max))/expit(p.int.st))
					alpha.st <- log(sample(1:5,1))
			  }
			})
	 	
  switch(lambda.cov,
		covar={lambda.slope.st <- 0},
		const={lambda.slope.st <- NULL})
  switch(p.cov,
		covar={p.slope.st <- 0},
		const={p.slope.st <- NULL})

#=====================================================
# Model fitting
#=====================================================  
# Fit multivariate Poisson/negative-binomial model		
	fit_mv <- optim(c(lambda.int.st,lambda.slope.st,p.int.st,p.slope.st,alpha.st), logl.mvp,, yMat, hessian = TRUE, control = list(reltol = 1e-12))
# Output from the model
	switch(lambda.cov,
		const = {lambda.out <- exp(fit_mv$par[1])
				 outnum <- 1},
		covar = {lambda.out <- exp(fit_mv$par[1] + fit_mv$par[2]*lambda.var)
				 outnum <- 2})
	switch(p.cov,
		const = {p.out <- expit(fit_mv$par[outnum + 1])
				 outnum <- outnum + 1},
		covar = {p.out <- expit(fit_mv$par[outnum + 1] + fit_mv$par[outnum + 2]*p.var)
				 outnum <- outnum + 2})
	if(model=="nb"){alpha.out <- exp(fit_mv$par[outnum + 1])}
	
#=====================================================
# Fit Royle (2004) N-mixture model for comparison
# Note we have not included code for covariates here but this is accessible elsewhere, e.g. using unmarked						
	fit_nm <- try(optim(par=c(p.int.st, lambda.int.st,alpha.st), fn = ll_func2, hessian = TRUE,control = list(reltol = 1e-12)),silent = TRUE)
# Output from the model		
	expit(fit_nm$par[1])
	exp(fit_nm$par[2])
	if(model == "nb") exp(fit_nm$par[3])
    
 
		
		