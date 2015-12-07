# This file contains the functions required to optimise the N-mixture model, using the multivariate Poisson/negative-binomial 
# formulation, as well as the standard Royle (2004) for any number of sampling occasions
 

# Logit link functions for constraining p
  logit <- function(x){ log(x/(1-x)) }
  expit <- function(x){ exp(x)/(1+exp(x)) }
 
# Function to evaluate the first covariance diagnostic
  covstarm <- function(yMat){
				p1 <- ct <- 0
				for(i in 1:(nv-1)){
					for(j in (i+1):nv){
						p1 <- p1 + yMat[,i]*yMat[,j]
						ct <- ct+1}}
			
				mean(p1)/ct-mean(yMat)^2
				}
  
# Log-likelihood for multivariate Poisson or negative-binomial for nv visits
#====================================
 logl.mvp <- function(param, x)
#====================================
  {
  # Read in parameters   
  # Lambda can be constant or varying with a spatial covariate
  switch(lambda.cov,
	covar = {alpha0 <- param[1]; alpha1 <- param[2]
	parno <- 2;  lambda <- exp(alpha0 + alpha1*lambda.var)},
	const = {alpha0 <- param[1]; alpha1 <- 0
	parno <- 1; lambda <- rep(exp(alpha0),nrow(x))})
  # p can be constant or varying with a covariate
  switch(p.cov,
	covar = {beta0 <- param[parno+1]; beta1 <- param[parno+2];parno <- parno+2; p <-expit(beta0 + beta1*p.var)},
	const = {beta0 <- param[parno+1]; beta1 <- 0;parno <- parno+1;p <- matrix(expit(beta0), nrow = nsites, ncol = nv)})
  # Dispersion parameter for negative binomial case
  switch(model,
	pois = {alpha <- NULL},
	nb = {alpha <- exp(param[parno+1])})
  # Evaluate the likelihood at each site	
  nobs <- dim(x)[1]
  prb <- numeric(nobs)
  for (i in 1: nobs) {
		prb[i] <- mv_func(x[i,], lambda[i], p[i,], alpha)
		}
  -sum(log(prb))    
  }

# Evaluate the multivariate or negative binomial Poisson for nv visits
#==============================
 mv_func <- function(x, lambda, p, alpha)
#==============================
  {
  # This function evaluations the multivariate distribution by generating
  # all possible values of the latent variables
  # (and many other values that are not possible) and eliminating
  # the unwanted ones in a step-by-step manner. 
  if(model == "nb") beta <- alpha/lambda	
  # Function determines which rows of gg are possible 	
  keep_func <- function(gg){
  	# Each column of w represents one of the visits n_it
  	w <- matrix(NA, nrow = nrow(gg), ncol = nv)
  	for(i in 1:nv){
  		tgg <- gg[, which(apply(poss, 1, function(x){sum(x==i, na.rm=TRUE)})==0)+1]
  		if(is.data.frame(tgg)) tgg <- apply(tgg, 1, sum)
  		w[,i] <- x[i] - gg[,1] - tgg
		}
  	# Determine which rows of gg retain positive n_it
  	cnd <- apply(w, 1, function(x){sum(x < 0) == 0})
  	return(cnd)
  }	
  
  # Matrix of possibilities corresponding to occasions not included for each X_i,s
  # nv = number of visits
  poss <- matrix(c(1:nv, rep(NA,nv*max(1,nv-2))), ncol = max(2,nv-1), nrow = nv)
  if(nv > 3){
  	for(i in 2:max((nv-2),2)){
  		poss <- rbind(poss, rbind(cbind(t(combn(1:nv,i)),matrix(NA, nrow = nrow(t(combn(1:nv,i))), ncol = nv-i-1))))
		}
	}	
  # Evaluate all possible X_i,1:nv assuming all other X_i,s are zero
  gg <- merge(0:min(x), data.frame(setNames(replicate(nrow(poss), 0, simplify = F), paste("x", 2:(nrow(poss)+1), sep = ""))))
  colnames(gg) <- paste("x", 1:(nrow(poss)+1), sep = "")
  gg <- gg[keep_func(gg),]
  
  if(nv > 2){
	for(d in 1:nrow(poss)){
		name <- paste("x", d+1, sep = "")
		temp <- data.frame(0:min(x[-poss[d,][!is.na(poss[d,])]]))
		colnames(temp) <- name
		gg <- merge(unique(gg[,!(names(gg) %in% name)]), temp)
		gg <- gg[,paste("x", 1:(nrow(poss)+1), sep = "")]
		gg <- gg[keep_func(gg),]
		}
	}
  # Evalute w, where each row represents n_it consisting of the sum of the latent X_i,s
  w <- matrix(NA, nrow = nrow(gg), ncol = nv)
  for(i in 1:nv){
	tgg <- gg[, which(apply(poss,1,function(x){sum(x == i, na.rm = TRUE)}) == 0)+1]
	if(is.data.frame(tgg)) tgg <- apply(tgg, 1, sum)
	w[,i] <- x[i] - gg[,1] - tgg
	}
  # Evaluate the joint probability function, where for each possible value of x_i,s we find the product of the probabilities for the latent variables X_i,s where E(X_i,s) = lambda*p^(|s|)*q^(T-|s|)
  q <- 1-p
  switch(model,
	pois = {
	    tmp <- dpois(gg[,1], lambda * prod(p))
		if(nv > 2){
			for(i in 2:(nrow(poss)+1)){
				qvals <- unique(poss[i-1,][!is.na(poss[i-1,])])
				pvals <- (1:nv)[!(1:nv %in% qvals)]
				tmp <- tmp * dpois(gg[,i],lambda*prod(p[pvals]) * prod(q[qvals]))
				}
			}
		for(i in 1:nv){
			tmp <- tmp * dpois(w[,i], lambda * p[i] * prod(q[-i]))
			}},
	nb = {
		sumx <- rowSums(gg) + apply(w, 1, sum)
		tmp <- prod(p)^gg[,1] * beta^alpha / gamma(alpha) * gamma(sumx + alpha)  / (1 - prod(q) + beta)^(sumx + alpha) /apply(factorial(gg), 1, prod) / apply(factorial(w), 1, prod)
		if(nv > 2){
			for(i in 2:(nrow(poss)+1)){
				qvals <- unique(poss[i-1,][!is.na(poss[i-1,])])
				pvals <- (1:nv)[!(1:nv %in% qvals)]
				tmp <- tmp*(prod(p[pvals])*prod(q[qvals]))^gg[,i]
				}
			}
		for(i in 1:nv){
			tmp <- tmp * (p[i] * prod(q[-i]))^w[,i]
			}
		})	
  sum(tmp)
  }

# N-mixture model (for comparison)
#==============================
ll_func2 <- function(param) 
#==============================
  {
  # Evaluate the likelihood for the Royle (2004) N-mixture model
  p <- expit(param[1])
  lambda <- exp(param[2])
  if(model=="nb"){alpha <- exp(param[3])}
  nsites <- length(nrepls)
  logLike <- rep(NA, nsites)
  for(i in 1:nsites){
	y <- yMat[i, 1:nrepls[i]]
	Nmin <- max(y)	
	switch(model,
		pois = {Ninfty <- qpois(1-1e-8, lambda)},
		nb = {Ninfty <- qnbinom(1-1e-8, mu = lambda, size = alpha)})	
	siteSum <- 0
	N <- Nmin:Ninfty
	switch(model,
		pois = {logSum <- log(dpois(N, lambda = lambda))},
		nb = {logSum <- log(dnbinom(N, mu = lambda, size = alpha))})
	for(j in 1:nrepls[i]){
		logSum <- logSum + dbinom(y[j], size = N, prob = p, log = TRUE)
		}
	siteSum <- siteSum + sum(exp(logSum))	
	logLike[i] <- log(siteSum)
	} 
  -1*sum(logLike) 
  }
  
  
 
  