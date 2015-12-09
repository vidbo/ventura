

library("xlsx")

load(file="VenturaData.RData")

library(rstan)
library(reshape2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# temporarily reduce sample size
# diveDat <- diveDat[diveDat$yr < 2008,]
# EFDat <- EFDat[EFDat$yr < 2008,]
# mapDat <- mapDat[mapDat$yr < 2008,]


# study 1
dat1 <- list(
   M = length(mapDat$siteN),
   D = length(diveDat$d_event),
   Dd= max(diveDat$d_event),
   E = length(EFDat$e_event),
   Ee= max(EFDat$e_event),
   #U = as.vector(tapply(mapDat$unith, INDEX=mapDat$siteN, FUN=max, na.rm=TRUE)), ## max unit id in each site
   U = max(mapDat$unit),
   S = max(mapDat$siteN),
   Y = max(mapDat$yr)-2005,
   m_yr =   mapDat$yr-2005,
   m_unit = mapDat$unit,
   m_site = mapDat$siteN,
   m_len  = mapDat$Length,
   d_yr   = tapply(diveDat$yr-2005, diveDat$d_event, FUN=max),
   d_pass = diveDat$Dive_Pass,
   d_fry  = diveDat$Fry,
   d_juv  = diveDat$Juv,
   d_unit = tapply(diveDat$unit, diveDat$d_event, FUN=max),
   d_site = tapply(diveDat$siteN, diveDat$d_event, FUN=max),
   d_len  = tapply(diveDat$Length, diveDat$d_event, FUN=max),
   d_event= diveDat$d_event,
   d_pos  = match(unique(diveDat$d_event), diveDat$d_event),  # starting position in d_fry/d_juv of each d_event
   d_reps = tapply(diveDat$Dive_Pass, diveDat$d_event, FUN=max), # total number of dive passes for each dive event
   e_yr   = tapply(EFDat$yr-2005, EFDat$e_event, FUN=max),
   e_pass = EFDat$EF_Pass,
   e_fry  = EFDat$Fry,
   e_juv  = EFDat$Juv,
   e_unit = tapply(EFDat$unit, EFDat$e_event, FUN=max),
   e_site = tapply(EFDat$siteN, EFDat$e_event, FUN=max),
   e_len  = tapply(EFDat$Length, EFDat$e_event, FUN=max),
   e_event= EFDat$e_event,
   e_pos  = match(unique(EFDat$e_event), EFDat$e_event),  # starting position in e_fry/e_juv of each e_event
   e_reps = tapply(EFDat$EF_Pass, EFDat$e_event, FUN=max) # total number of EF passes for each e_event
)
dat1

fit1 <- stan(file="ventura01.stan", data=dat1,iter=2000, chains=2)

library(rv)
var1 <- as.rv(fit1)










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



