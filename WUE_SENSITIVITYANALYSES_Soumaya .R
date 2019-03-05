

# NOTATION: for linear Farquhar model <- _l
#           for Full Farquhar model <- _f

# define 2 functions for WUE: WUE_l $ WUE_c

setwd('C:/Users/kna016/Documents/Projects/WUE_treerings/')

# INPUT DATA --------------------------------------------------------------

## read input file for all data
#data_TR<- read.csv("isonet_frank.csv", header=T) #time series of tree-rings 13C
data_TR<- read.csv("isonet_tr.csv", header=T) #time series of tree-rings 13C


startyear <- 1900
endyear   <- 2003

## data from CMIP6/ Graven et al 2017
# CMIP6 <- read.csv("~/Documents/GCB_PAPER_NE13C/ATM_GHG/TableS1_170710_submitted.csv", header=T)
CMIP6 <- read.csv("TableS1_170710_submitted.csv", header=T)
# atm_13C <- ts(CMIP6$Global.delta13co2, start=1850, frequency = 1) # needs subsetting the data for the time window of obs
atm_13C <- CMIP6[floor(CMIP6[,"Date"]) >= startyear & floor(CMIP6[,"Date"]) <= endyear,"Global.delta13co2"]

## CO2 concentration data from SCRIPPS, spline ice core and merged obs downlanded June 5 2017
CO2_atm <- read.csv("spline_merged_ice_core_yearly_SB.csv", header=T)
# atm_C <- ts(data= CO2_atm$CO2[CO2_atm$Year>=1901], start=1901, frequency = 1) # needs subsetting the data for the time window of obs
atm_C <- CO2_atm[CO2_atm[,"Year"] >= startyear & CO2_atm[,"Year"] <= endyear,"CO2"]
  
  

leaf <- 2.1 # offset between leaf and stem
btv <- 0.8 # between tree variability based on idividual tree measurements
uncert <- 1.5+btv #uncertainty of terms [FRANK ET AL NCC-2016]
d13_plant <- ts(data=(data_TR-leaf), start=1901, end=2003) # upscale from tree ring to leaf



# WUE & CARBOXYLATION [Linear Model] -----------------------------------------------------

# The most commonly invoked values for b in higher plants range between 26 and 30 permil
# (Christeller et al., 1976; Wong et al., 1979; Farquhar et al., 1982; Roeske and O’Leary, 1984; 
# Guy et al., 1993; Lloyd and Farquhar, 1994; Suits et al., 2005)

# because the diffusion of carbon dioxide into plant leaves is a passive process and pCO2 cannot
# be concentrated within tissues relative to the atmosphere: b>= Δ as in EQ of Farquhar et al., 1989
# i.e.,ca–ci may not be less than zero and therefore ci/ca must not exceed 1.

# NOTE: IN THE LINEAR FARQUHAR MODEL  * b IS 27 PERMIL (FARQUHAR ET AL 1982-OECOLOGIA & FARQUHAR AND RICHARDS 1984)
#       IN THE CLASSIC FRQUHAR MODEL  * b IS SET AT 29 PERMIL IN SEIBT ET AL Oecologia-2008
#                                       [Roeske & O'Leary 1984, Guy et al. 1993]
#                                     * b IS SET AT 29 IN CERNUSAK ET AL NEW PHYTO-2013
#                                     * b IS SET AT 30 PERMIL IN KEELING ET AL PNAS-2017
#                                     * b is described as 30 PERMIL in FARQUHAR ET AL 1982 (Aust. J. Plant Physiol)

#  Here we use the linear model of net discrimination 
#  Eq 10, Farquhar et al. 1982
a=4.4
b-27
# wue_l


## main function: calculates ci/ca and WUE from isotope data with different models
WUE_isotopes <- function(ca=atm_C, dplant=d13_plant, datm=atm_13C, 
                         a     = 4.4,   # fractionation during CO2 diffusion through stomata
                         b_lin = 27,    # fractionation during carboxylation in the simple model 
                         b     = 30,    # fractionation during carboxylation
                         am    = 1.8,   # fractionation during internal CO2 transfer
                         f     = 12,    # fractionation during photorespiration
                         gam   = 43,    # gamma_star
                         An    = 9,     # net assimilation
                         gm    = 0.2,   # mesophyll conductance
                         model = c('linear','classical')){

  model <- match.arg(model)
  
  big_delta = (datm-dplant)/(1+dplant/1000) # Discrimination Farquhar 1982
  
  if (model == 'linear'){
  
    ciOca = (big_delta - a) / (b_lin - a)
  
  } else if (model == 'classical'){
    
    ciOca = ((b - am) * An/(gm * ca) + f*gam/ca - a + big_delta) / (b - a)  # Eq. 6 in Seibt et al. 2008 solved for Ci/Ca
    
  }

  ci    <- ca * ciOca
  caMci <- ca-ci
  iWUE  <- (ca-ci)/1.6
  
  
  return(data.frame(ci,ciOca,caMci,iWUE,big_delta))
  
}


## some tests with the function:
WUE_isotopes(ca=380,dplant=-29,datm=-8,model='linear')
WUE_isotopes(ca=380,dplant=-29,datm=-8,model='classical')              # with default arguments (see function definition)
WUE_isotopes(ca=380,dplant=-29,datm=-8,model='classical',gm=0.1,An=4)  # changing some arguments (see function definition)


# first test with real data
COL_CA_simple    = WUE_isotopes(ca=atm_C,dplant=data_TR[,"COL_CA"],datm=atm_13C,model='linear')
COL_CA_classical = WUE_isotopes(ca=atm_C,dplant=data_TR[,"COL_CA"],datm=atm_13C,model='classical')

# plot
par(mfrow=c(1,1),mar=c(4,5,1,1))
plot(COL_CA_simple[,'ciOca'] ~ c(startyear:endyear),type='l',ylim=c(0.4,0.7),xlab='Year',ylab='Ci/Ca',main='COL_CA')
points(COL_CA_classical[,'ciOca'] ~ c(startyear:endyear),type='l',col='green3')
legend('topleft',c('simple','classical'),col=c('black','green3'),lty=1,bty='n',y.intersp=1.5)

dev.copy2pdf(file="C:/Users/kna016/Documents/Projects/WUE_treerings/plots/COL_CA.pdf",
             width=6,height=4)


# next steps:
# - loop over all sites
# - determine values for An and gm








# FULL FARQUHAR MODEL -----------------------------------------------------

# Here we use the full Farquhar model to derive ci 
# We use Eq1 in Keeling ET AL PNAS-2017 & Eq 6 in SEIBT ET AL Oecologia-2008
# which is Eq B24 in Farquhar et al. 1982


#LIST OF PARAMETERS
#DEFAULTS VALUES ARE FROM SEIBT ET AL Oecologia-2008 & KEELING ET AL PNAS-2017

A=9 # leaf-level gross photosynthesis (μmol·m−2·s−1); SEIBT ET AL Oecologia-2008
gi=0.2 # mesophyll conductance (mol·m−2·s−1); SEIBT ET AL Oecologia-2008
Γ=43   # CO2 compensation point in the absence of day respiration; Keeling ET AL PNAS-2017.
f= 12  # the discrimination due to photorespiration (permil); Keeling ET AL PNAS-2017
       # FROM CERNUSAK EL AL NEW PHYTOLOGIST-2013 (AS CITED BY KEELING IN PNAS2017)
       # - Recent estimates for f range from c. 8 to 16 (Gillon & Griffiths, 1997;Lanigan et al., 2008; 
       #   Evans&von Caemmerer, 2013) with a value of 11 suggested on theoretical grounds (Tcherkez, 2006).
       # - SEIBT ET AL Oecologia-2008 uses 8 permil


am=1.8 # fractionation during the internal (mesophyll) CO2 transfer
#        SUM OF: 1-Fractionation during liquid diffusion 0.7 permil  O’Leary (1984)
#                2-Fractionation of CO2 entering solution 1.1 PERMIL Mook et al. (1974)

b=30 # fractionation during carboxylation as set in the Original Farquhar model
#Δ= this will be big delta as computed from the WUE function 

# UNCERTAINTIES OF 
# 1-uncertainties for classic model parameters
# VALUES ARE : 4 PEMIL FOR F; 3 PERMIL FOR A; 0.05 PERMIL FOR GI

uncert_f <- 4+3+0.05 


ci equation will be in WUE_c function
ci_c= ((Δ+(f*Γ/ca)+(((b1-am)*A)/(ca*gi))-a)*ca)/ (b1-a)  # JK: where does this one come from?


  
  
  




# ADDITIONAL DISC FROM KEELING --------------------------------------------

# Total change of CO2 concetration from 1901 to 2012 is 96.87 ppm
# sensitivity of mesophyl and photorespiration is 0.010 0.004 per 1 ppm of CO2
MP <- 0.010*96.87
UMP <- 0.004*96.87
DK <- 0.0144*96.87
UDK <- 0.007*96.87
# FROM 1901 TO 2012
# TOTAL CONTRIBUTION OF MESOPHYL AND PHOTORESPIRATION IS 0.9687 ± 0.38748 ‰
# TOTAL DISCRIMINATION 1.394928 ± 0.67809 ‰

