
library(ncdf4)
library(raster)
library(bigleaf)

setwd("C:/Users/kna016/Documents/Projects/WUE_treerings/sites_met")
sitelist <- read.csv("../isonet_tr_metadata.csv")

vars_out  <- c("wind","Tair","VPD","pressure","precip","PPFD")
years     <- c(1901:2003)
light_cor <- 24/16 # assuming average daylength of 16 hours, should be site-specific


### 1) get met files (as processed in extract_climate_vars_daily.ksh)
vars <- c("vgrd","ugrd","tmax","spfh","pres","pre","dswrf")
for (var in vars){
  assign(paste0("cru_",var),raster::brick(paste0("CRUJRA_daily/",var,"_1901_2003.nc")))
}




### 2) loop over sites and extract met data of matching pixel
for (i in 1:nrow(sitelist)){
  
  site <- as.character(sitelist[i,"Site.Abbreviation"])
  lat  <- sitelist[i,"Latitude"]
  lon  <- sitelist[i,"Longitude"]
  coords <- SpatialPoints(matrix(c(lon,lat),nrow=1),proj4string=CRS(projection(cru_pre)))
  
  res <- matrix(NA,ncol=length(vars_out),nrow=length(years),
                dimnames=list(years,vars_out))
  
  # extract all vars
  for (var in vars){
    
    r <- get(paste0('cru_',var))
    assign(paste0(var,"_vals"),as.vector(extract(r,coords)))
    
  }
  
  
  # some more processing
  res[,"wind"]     = sqrt(vgrd_vals^2 + ugrd_vals^2)  # ms-1
  res[,"Tair"]     = tmax_vals - 273.15  # degC
  res[,"pressure"] = pres_vals / 1000    # kPa 
  res[,"VPD"]      = bigleaf::q.to.VPD(spfh_vals,res[,"Tair"],res[,"pressure"])
  res[,"precip"]   = pre_vals   # mm precip in summer months
  res[,"PPFD"]     = bigleaf::Rg.to.PPFD(dswrf_vals) * light_cor   # umol m-2 s-1
  
  write.csv(res,file=paste0(site,"_yearly_met.csv"))
  
}