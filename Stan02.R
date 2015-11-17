

library("xlsx")

load(file="VenturaData.RData")








datdens <- read.xlsx(file="Seco_OmyDensAdj.xlsx", sheetName="OmyDensAdj")
datdens2 <- datdens
datdens2[,3:20] <- 10^datdens2[,3:20]  # convert from base-10 logarithms to densities

# Analysis of depth associations from 2 Spina papers

# relate the two tables
dvect <- paste0(dat$quantity, "_", dat$study, "_", dat$age)
svect <- paste0(scal$quantity, "_", scal$study, "_", scal$age)
idx <- match(dvect, svect)
scal2 <- scal[idx,] # rows of scal2 match rows of dat

# compute the actual counts and probs from the raw numbers (measurements taken from a graph)
dat$N <- scal2$n

ref <- aggregate(dat$rawFreq, by=list(dvect), FUN=sum)
idx <- match(dvect, ref[,1])
ref2 <- ref[idx,]

dat$p <- dat$rawFreq / ref2[,2]  # probabilities computed from measurements taken off graph
dat$freq <- round(dat$p * dat$N) # estimated counts will have rounding error

dat$freq

# check the magnitude of rounding error
tst <- aggregate(dat$freq, by=list(dvect), FUN=sum)
idx <- match(dvect, tst[,1])
cbind(dat$N, tst[idx,2])
# a bit of error, but not much


# compute the conditional probabilities
datfish <- dat[dat$age!="available",]
dataval <- dat[dat$age=="available",]

dfishvect <- paste0(datfish$quantity, "_", datfish$study, "_", datfish$x.axis)
davalvect <- paste0(dataval$quantity, "_", dataval$study, "_", dataval$x.axis)
idx <- match(dfishvect, davalvect)
dataval <- dataval[idx,]

datfish$p_aval <- dataval$p
datfish$freq_aval <- dataval$freq
datfish$relfreq <- datfish$p / datfish$p_aval


jdx <- datfish$study==1 & datfish$quantity=="depth" & datfish$age==0
kdx <- datfish$study==2 & datfish$quantity=="depth" & datfish$age==0

plot(datfish$x.axis[jdx], datfish$relfreq[jdx], type="b")
points(datfish$x.axis[kdx], datfish$relfreq[kdx], type="b", col="red")


jdx <- datfish$study==1 & datfish$quantity=="depth" & datfish$age==1
kdx <- datfish$study==2 & datfish$quantity=="depth" & datfish$age==1

plot(datfish$x.axis[jdx], datfish$relfreq[jdx], type="b", ylim=c(0,15))
points(datfish$x.axis[kdx], datfish$relfreq[kdx], type="b", col="red")

jdx <- datfish$study==1 & datfish$quantity=="depth" & datfish$age==2
kdx <- datfish$study==2 & datfish$quantity=="depth" & datfish$age==2

plot(datfish$x.axis[jdx], datfish$relfreq[jdx], type="b", ylim=c(0,15))
points(datfish$x.axis[kdx], datfish$relfreq[kdx], type="b", col="red")



jdx <- datfish$study==1 & datfish$quantity=="depth" & datfish$age==0
kdx <- datfish$study==1 & datfish$quantity=="depth" & datfish$age==1
ldx <- datfish$study==1 & datfish$quantity=="depth" & datfish$age==2

plot(  datfish$x.axis[jdx], datfish$relfreq[jdx]/sum(datfish$relfreq[jdx], na.rm=T), type="b", ylim=c(0,1))
points(datfish$x.axis[kdx], datfish$relfreq[kdx]/sum(datfish$relfreq[kdx], na.rm=T), type="b", col="red")
points(datfish$x.axis[ldx], datfish$relfreq[ldx]/sum(datfish$relfreq[ldx], na.rm=T), type="b", col="blue")

tot <- datfish$relfreq[jdx]/sum(datfish$relfreq[jdx], na.rm=T) + datfish$relfreq[kdx]/sum(datfish$relfreq[kdx], na.rm=T) + datfish$relfreq[ldx]/sum(datfish$relfreq[ldx], na.rm=T)
points(datfish$x.axis[ldx],tot, type="b",lwd=2)

plot(  datfish$relfreq[jdx]/sum(datfish$relfreq[kdx], na.rm=T), datfish$relfreq[kdx]/sum(datfish$relfreq[jdx], na.rm=T), type="b", col="red")
points(datfish$relfreq[jdx]/sum(datfish$relfreq[ldx], na.rm=T), datfish$relfreq[ldx]/sum(datfish$relfreq[jdx], na.rm=T), type="b", col="blue")



#####



jdx <- datfish$study==2 & datfish$quantity=="depth" & datfish$age==0
kdx <- datfish$study==2 & datfish$quantity=="depth" & datfish$age==1
ldx <- datfish$study==2 & datfish$quantity=="depth" & datfish$age==2

plot(  datfish$x.axis[jdx], datfish$relfreq[jdx]/sum(datfish$relfreq[jdx], na.rm=T), type="b", ylim=c(0,1))
points(datfish$x.axis[kdx], datfish$relfreq[kdx]/sum(datfish$relfreq[kdx], na.rm=T), type="b", col="red")
points(datfish$x.axis[ldx], datfish$relfreq[ldx]/sum(datfish$relfreq[ldx], na.rm=T), type="b", col="blue")

tot <- datfish$relfreq[jdx]/sum(datfish$relfreq[jdx], na.rm=T) + datfish$relfreq[kdx]/sum(datfish$relfreq[kdx], na.rm=T) + datfish$relfreq[ldx]/sum(datfish$relfreq[ldx], na.rm=T)
points(datfish$x.axis[ldx],tot, type="b",lwd=2)



# library(mgcv)
# 
# # recode for presence versus available data
# dat$y <- as.integer(! dat$age=="available")
# 
# jdx <- dat$study==1 & dat$quantity=="depth" & (dat$age=="0" | dat$age=="available")
# dat[jdx,]
# 
# 
# #dat$freq <- dat$freq*100
# 
# b <- gam(y ~ s(x.axis), family=poisson(), data=dat[jdx,], weights=freq)
# 
# summary(b)
# 
# plot(b, trans=exp)
# 
# plot(b,pages=1,residuals=TRUE, trans=function(x) {1/(1+exp(-x))})
# 
# kdx <- datfish$study==1 & datfish$quantity=="depth" & datfish$age==1
# points(datfish$x.axis[kdx], datfish$relfreq[kdx], type="b", col="red")
# 
# 
# 
# gam.check(b)


library(rstan)
library(reshape2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


mdat <- melt(dat, id.vars=c("study", "quantity", "age", "x.axis"), measure.vars="freq")

datm <- acast(mdat, study + quantity + age ~ x.axis, value.var="value")
datm

# study 1
dat1 <- list(
   avail = as.vector(datm[4,3:17]),
   used  = datm[1:3,3:17],
   I = 15, J = 3,
   cut1 = c(7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5))
dat1

fit1 <- stan(file="ratiofunc.stan", data=dat1,iter=1000, chains=4)

# study 2
dat2 <- list(
  avail = as.vector(datm[12,2:18]),
  used  = datm[9:11,2:18],
  I = 17, J = 3 )
dat2

fit2 <- stan(file="ratio.stan", data=dat2,iter=1000, chains=4)



library(ggplot2)

# study 1
rat <- extract(fit1, pars="ratio")$ratio   # MCMC samples with warmup removed, chains merged and permuted
dim(rat)
rat <- apply(rat, 2:3, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
dim(rat)
age.x <- 0:2       # 3 age classes
depth.x  <- 5*(2:16)  # depth classes

# study 2
rat2 <- extract(fit2, pars="ratio")$ratio   # MCMC samples with warmup removed, chains merged and permuted
dim(rat2)
rat2 <- apply(rat2, 2:3, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
dim(rat2)
depth2.x  <- 5*(1:17)  # depth classes


# dataframes for 3 age classes
rdat <- data.frame(study=1, age=0, Depth=depth.x, t(rat[,1,]))
rdat <- rbind(rdat, data.frame(study=1, age=1, Depth=depth.x, t(rat[,2,])))
rdat <- rbind(rdat, data.frame(study=1, age=2, Depth=depth.x, t(rat[,3,])))

# same, but for study 2
rdat <- rbind(rdat, data.frame(study=2, age=0, Depth=depth2.x, t(rat2[,1,])))
rdat <- rbind(rdat, data.frame(study=2, age=1, Depth=depth2.x, t(rat2[,2,])))
rdat <- rbind(rdat, data.frame(study=2, age=2, Depth=depth2.x, t(rat2[,3,])))



h <- ggplot(rdat, aes(x=Depth)) + theme_bw()  + ylab("Ratio of Used-to-Available")
h <- h + geom_ribbon(aes(ymin=X25., ymax=X75.), fill="light grey") + geom_line(aes(y=X50.))
h + facet_grid(age ~ study, scales = "free_y", labeller = label_both)





# analysis of densities from Seco 2006 survey

plot(datdens2$OmyS2mean, datdens2$OmyL2mean)
plot(sort(datdens2$OmyS2mean))

plot(datdens2$OmyS1mean, datdens2$OmyM1mean)



