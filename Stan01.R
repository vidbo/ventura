

library("xlsx")

maps <- paste0("map.", sprintf("%02i", 6:12))
dives <- paste0("dive.", sprintf("%02i", 6:12))
EFs <- paste0("EF.", sprintf("%02i", c(6:7,10:12)))
LWts <- paste0("L-Wt.", sprintf("%02i", c(6:7,10:12)))

getventura <- function(sheetNames, colIndex=NULL, file="Matilija raw data 06-12.xlsx") {
  out <- as.data.frame(NULL)
  for(sheetName in sheetNames) {
    this <-read.xlsx(file=file, sheetName=sheetName, colIndex=colIndex, stringsAsFactors=FALSE)
    yr <- rep(2000+as.integer(substr(sheetName, nchar(sheetName)-1,15)), dim(this)[1])
    this <- data.frame(this, yr=yr)
    out <- rbind(out, this)
  }
  out[! is.na(out$Unit_ID),]
}

mapDat <- getventura(maps, colIndex=1:3)
diveDat<- getventura(dives, colIndex=1:9)
EFDat  <- getventura(EFs,   colIndex=1:9)
LWtDat <- getventura(LWts, colIndex=1:5)

parseUnit_ID <- function(dat) {
  obj <- as.character(dat$Unit_ID)
  a <- regexpr("_", obj)
  site <- (substr(obj, start=1, stop=a-1))
  rem  <- substr(obj, start=a+1, stop=15)
  a <- regexpr("_", rem)
  unit <- as.integer(substr(rem, start=a+1, stop=15))
  out <- data.frame(dat, site=as.character(site), unith=unit, stringsAsFactors=FALSE)
}

mapDat <- parseUnit_ID(mapDat)
diveDat <- parseUnit_ID(diveDat)
EFDat <- parseUnit_ID(EFDat)
LWtDat <- parseUnit_ID(LWtDat)

# some fields that are characters should be integers
diveDat$Dive_Pass <- as.integer(as.character(diveDat$Dive_Pass))
diveDat$Fry <- as.integer(as.character(diveDat$Fry))
EFDat$EF_Pass <- as.integer(as.character(EFDat$EF_Pass))
EFDat$Fry <- as.integer(as.character(EFDat$Fry))

# remove rows that do not have useful data
idx <- ! (is.na(diveDat$Dive_Pass) | is.na(diveDat$Fry) | is.na(diveDat$Juv))
diveDat <- diveDat[idx,]

idx <- ! (is.na(EFDat$EF_Pass) | is.na(EFDat$Fry) | is.na(EFDat$Juv))
EFDat <- EFDat[idx,]

# make a factor with common coding for all datasets
u <- mapDat$site
w <- mapDat$unith
v <- rep(1, length(mapDat$site))
u <- c(u, diveDat$site)
w <- c(w, diveDat$unith)
v <- c(v, rep(2, length(diveDat$site)))
u <- c(u, EFDat$site)
w <- c(w, EFDat$unith)
v <- c(v, rep(3, length(EFDat$site)))
u <- c(u, LWtDat$site)
w <- c(w, LWtDat$unith)
v <- c(v, rep(4, length(LWtDat$site)))
u <- factor(u)

# now cut the factor back up to individual datasets
mapDat$sitef <- u[v==1]
diveDat$sitef <- u[v==2]
EFDat$sitef <- u[v==3]
LWtDat$sitef <- u[v==4]

mapDat$siteN <- unclass(mapDat$sitef)
diveDat$siteN <- unclass(diveDat$sitef)
EFDat$siteN <- unclass(EFDat$sitef)
LWtDat$siteN <- unclass(LWtDat$sitef)


# create unique id numbers for each unit (may need to include year as well if units differ year to year)
units <- unique(data.frame(siteN=unclass(u), unith=w))
ord <- order(units$siteN, units$unith)
units <- units[ord,]
units$unit <- 1:dim(units)[1]  # new unique id number for units

# assign unique unit ids to datasets
idx <- match(paste0(as.character(mapDat$siteN), "_", as.character(mapDat$unith)) , 
             paste0(as.character(units$siteN),  "_", as.character(units$unith)))
mapDat$unit <- units$unit[idx]

idx <- match(paste0(as.character(diveDat$siteN), "_", as.character(diveDat$unith)) , 
             paste0(as.character(units$siteN),   "_", as.character(units$unith)))
diveDat$unit <- units$unit[idx]

idx <- match(paste0(as.character(EFDat$siteN), "_", as.character(EFDat$unith)) , 
             paste0(as.character(units$siteN), "_", as.character(units$unith)))
EFDat$unit <- units$unit[idx]

idx <- match(paste0(as.character(LWtDat$siteN), "_", as.character(LWtDat$unith)) , 
             paste0(as.character(units$siteN), "_", as.character(units$unith)))
LWtDat$unit <- units$unit[idx]

# save the final datasets
save(mapDat, diveDat, EFDat, LWtDat, file="VenturaData.RData")




