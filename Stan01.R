

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
  out <- data.frame(dat, site=as.character(site), unit=unit, stringsAsFactors=FALSE)
}

mapDat <- parseUnit_ID(mapDat)
diveDat <- parseUnit_ID(diveDat)
EFDat <- parseUnit_ID(EFDat)
LWtDat <- parseUnit_ID(LWtDat)

u <- mapDat$site
v <- rep(1, length(mapDat$site))
u <- c(u, diveDat$site)
v <- c(v, rep(2, length(diveDat$site)))
u <- c(u, EFDat$site)
v <- c(v, rep(3, length(EFDat$site)))
u <- c(u, LWtDat$site)
v <- c(v, rep(4, length(LWtDat$site)))
u <- factor(u)

mapDat$sitef <- u[v==1]
diveDat$sitef <- u[v==2]
EFDat$sitef <- u[v==3]
LWtDat$sitef <- u[v==4]

mapDat$siteN <- unclass(mapDat$sitef)
diveDat$siteN <- unclass(diveDat$sitef)
EFDat$siteN <- unclass(EFDat$sitef)
LWtDat$siteN <- unclass(LWtDat$sitef)

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


save(mapDat, diveDat, EFDat, LWtDat, file="VenturaData.RData")




