
#### Prevalence of top-down control in coastal planktonic food webs

# This script contains the complete statistical analysis, as well as the code 
# for all figures in the results section

####

library(plyr)      # To work with dataframes  
library(lubridate) # To manage the dates
library(stats)
library(gridExtra)
#library(lme4)     # Package for GLMMs
source('functions_main.R')

#### Phytoplankton data
stationf = read.csv("Phytoplankton_NLWKN.csv")

#### Mixotrophs data
dino.ess = read.csv("Mixotrophs_NLWKN.csv")
dino.ess = dino.ess[-which(dino.ess$bioC==0),] # Remove the taxa absent from the samples

#### Fish landings data
datah = read.csv("fishlandings.csv", skip=4, header=T) # Extracting the dataset
planktinames = c( "Hering" ) # Focus on Herring

dataf = fish.approx(datah, planktinames) # Herring

# There are only 2 months per year (non NaN) before 2013 on average
ind = which(dataf$year >= 2013) 
dataf = dataf[ind, ] # There are not much data before 2013

# Creating dates variable
dataf$date = ymd( paste( dataf$year, "-", dataf$month, "-15", sep = "" ) )

#### Zoo data
data = read.csv("Zooplankton_NLWKN.csv")
data$Date = ymd( data$Date )

# Removing lines where no OPS is available
data.ops = data[ -which( is.na(data$ops) ), ]
data$ops.alt = NA

## Environmental forcing data (Nutrients, Temperature, Salinity)
waddenEnvtot = read.table("Environment.csv", header = TRUE, sep=",")
waddenEnvtot = waddenEnvtot[ which( waddenEnvtot$StationID %in% unique(stationf$stationID) | waddenEnvtot$StationID %in% unique(data$Station) ), ] # Remove stations where no plankton information
statin = waddenEnvtot$StationID

#### Computing mean daily series of Phyto and Zoo ----

## Phytoplankton
# Computing the daily means of phytoplankton
phy.day.comp = function(stationf){
  un.dates = unique(stationf$USI) # Sampling each station (and day) once (to extract the total biomass on that day/station)
  
  statbio = ddply(stationf, .(USI), summarize, biovol.tot = sum(biovolume, na.rm = T))
  stationf = merge(statbio, stationf, by="USI")
  stationf$relbiov = stationf$biovolume / stationf$biovol.tot # Weight for mean esd
  
  statsum = ddply( stationf, .(USI), summarize, 
                   usisize = sum( esd * relbiov, na.rm = T), 
                   usibioc = sum( bioC, na.rm = T) ) # Stats per USI
  stationf = merge( stationf, statsum, by = "USI" )
  ind.un = c()
  for(iun in ( 1 : length(un.dates) )){
    iun = un.dates[iun]
    ind.un = c( ind.un, which(stationf$USI == iun)[1] )
  } # Sampling every USI once (to get the stats only once in the mean)
  
  phy.day = ddply( stationf[ind.un,], .(date), summarize, 
                   bioc = mean( log(usibioc +1), na.rm=T), # Mean Log BioC per day
                   size = mean( log(usisize), na.rm=T),    # Mean Log Size per day
                   sdc = sd( log(usibioc +1), na.rm=T),
                   sds = sd( log(usisize +1), na.rm=T) )  
  phy.day$yr = year(phy.day$date) ; phy.day$mth = month(phy.day$date)
  return(phy.day) } 
phy.day = phy.day.comp(stationf)
dfl.day = phy.day.comp(dino.ess)

## Zooplankton
data$USI = paste( data$Station, "_", data$Date, sep="" ) # Station codes for zooplankton

## Computing the daily averaged series of all zooplankton 

# Removing the NAs in OPS, as they will not be grouped between small and large feeders
# ops.alt is for the switching by adult copepods (between herbivory and carnivory)
zoo.day = ddply( data[ -which(is.na(data$ops)), ], .(USI, Date), summarize, 
                 dw = sum( DryW, na.rm=T),
                 esd = sum( DryW * mean.esd, na.rm = T ),
                 ops = sum( DryW * ops, na.rm = T ),
                 ops.alt = sum( DryW * ops.alt, na.rm = T ) ) # Sum Dry Weight per USI
zoo.day$esd = zoo.day$esd / zoo.day$dw # DW-weighted ESD
zoo.day$ops = zoo.day$ops / zoo.day$dw 
zoo.day$ops.alt = zoo.day$ops.alt / zoo.day$dw

zoo.day = ddply( zoo.day, .(Date), summarize, 
                 dw = mean( log(dw +1), na.rm=T),
                 esd = mean( log(esd), na.rm = T ),
                 ops = mean( log(ops), na.rm = T ),
                 ops.alt = mean( log(ops.alt), na.rm = T ),
                 sddw = sd( log(dw +1), na.rm=T),
                 sdesd = sd( log(esd), na.rm = T ) ) # Mean Dry Weight per day
zoo.day$yr = year(zoo.day$Date) ; zoo.day$mth = month(zoo.day$Date)
## ----

## Checking the variability between stations in the phytoplankton series
un.dates = unique(stationf$USI) # Sampling each station (and day) once (to extract the total biomass on that day/station)

statbio = ddply(stationf, .(USI), 
                summarize, biovol.tot = sum(biovolume, na.rm = T))
station = merge(statbio, stationf, by="USI")
station$relbiov = station$biovolume / station$biovol.tot # Weight for mean esd

statday = ddply( station, .(USI, date, stationID), summarize,
                 usisize = sum( esd * relbiov, na.rm = T),
                 usibioc = sum( bioC, na.rm = T) ) #  C/ESD per USI

# Checking a moving average on one week
statday$date = ymd( statday$date )
statday$mth = month( statday$date )

# Taking out the data in winter (not used in the analysis)
statday = statday[ which(statday$mth >= 3 & statday$mth <= 10), ]

spatvar = function( ind = unique(statday$stationID) ){
  statday = statday[ which(statday$stationID %in% ind), ]
  i = min(statday$date)
  rm.stat = data.frame()
  while(i < max(statday$date) ){
    ind = which( statday$date >= i & statday$date < i + days(7) ) # one week moving average
    
    rm.ind = data.frame( date = i + days(3),
                         Cl =  mean( log(statday$usibioc[ind] +1), na.rm = T),
                         C =  mean( statday$usibioc[ind], na.rm = T),
                         sdCl =  sd( log(statday$usibioc[ind] +1), na.rm = T),
                         sdC =  sd( statday$usibioc[ind], na.rm = T),
                         n = length( ind ) )
    rm.stat = rbind( rm.stat, rm.ind )
    i = i + days(7)
  }
  
  rm.stat = rm.stat[ -which( is.na( rm.stat$sdCl ) ), ]
  dspat = density( rm.stat$sdCl / rm.stat$Cl, bw="nrd0", kernel = "gaussian" )
  
  x11()
  plot( dspat$x, dspat$y, type = "l", lwd = 2, xaxs = "i", yaxs = "i",
        las = 1, cex.axis = 1.4,
        ylab = "Density", xlab = "CV of log phytoplankton C (1 week average)", cex.lab = 1.4)
}
spatvar()
## 

#### Checking the Herring Larvae feeding range using the literature ----

# Max copepod size, from Hufnagl and Peck (2011)
max.lherr = function( L ){
  pmax = 2200 / ( 1 + ( L / 14 )**(-2) ) 
  return( pmax / 1000 )
} # L = Larvae length (mm) --> Copepod length (mm)  by comparing to Fig.2

# Max Herring in the WS, from Maathuis et al, 2024
pmax = max.lherr( 27.7 * 10 ) 

# Majority of Herring below 12cm
pmean = max.lherr( 12 * 10 )

lherr = seq(0, 12*10, 0.5)
plot(lherr, max.lherr(lherr), type = "l", xlim = c(0, 40), las=1,
     xlab = "Herring larvae length, mm",
     ylab = "Copepod prey length, mm")

# Converting the Length to ESD, using width from Cohen and Lough, 1983
esd.estimate = function( L ){
  w1 = 0.31 * L + 0.05
  w2 = 0.23 * L + 0.09
  w3 = 0.40 * L + 0.04 # These are the three equations given
  
  w = mean( c(w1, w2, w3 ) )
  
  esd = ( ( L * ( w**2 ) ) ** (1/3) ) #/ 2
return( esd )} # L is copepod Length in mm, gives the ESD in mm

esd.herr = esd.estimate( pmean * 1000 ) # Microm
# Max ESD eaten by Herring

#### ----

#### Classes Using the Herring larvae OPS criteria ----
ops.lim.phy = 20 # Size classes of phytoplankton

# OPS limit for feeding on small and large phyto 
ops.lim.zoo = exp(log(ops.lim.phy))
esd.lim.zoo = 250

# Both zoo small and large eaten by Herring
#zoo.herr = data[ which(data$mean.esd <= esd.herr/2), ]

# Split small and large zooplankton while making sure they don't overlap
ind = which( data$ops < ops.lim.phy )
zoos = data[ ind, ] 
zoo = data[-ind,]

# Large zooplankton
ind = which( (zoo$ops < esd.lim.zoo & zoo$ops >= ops.lim.phy) & zoo$Taxa != "Chaetognath" & zoo$Taxa != "Gelatinous" )
zool = zoo[ ind, ] 
zoo = zoo[-ind,]

# Carnivorous zooplankton
ind = which( zoo$ops >= esd.lim.zoo | zoo$Taxa == "Chaetognath" | zoo$Taxa == "Gelatinous" )
zooc = zoo[ ind, ]  # Not eaten

zoohl = zooc[ which( log(zooc$ops) > 7 ), ] # They do not feed on large zoo, but bigger plankton
# Mostly Beroe sp.

zooc = zooc[-which( log(zooc$ops) > 7 ),] # Remove Beroe sp., as it feeds on the carnivorous zooplankton.

# Small phytoplankton
phys = stationf[which( stationf$esd <= ops.lim.phy ), ] 

# Large phytoplankton
phyl = stationf[which( stationf$esd > ops.lim.phy ), ] 

#### Splitting the dinoflagellates into two OPS groups as well ----
dfls = dino.ess[ which( dino.ess$ops <= ops.lim.phy ), ]
dfll = dino.ess[ which( dino.ess$ops > ops.lim.phy ), ]

# Heterotrophic or Mixotrophic dinoflagellates
dfls.day = phy.day.comp( dfls )[ c("date", "bioc", "size", "yr", "mth" ) ]
names(dfls.day) = c("date", "dflsc", "dflsesd", "yr", "mth") # Rename for a function later
# The large dinoflagellates group is almost empty
dfll.day = phy.day.comp( dfll )[ c("date", "bioc", "size", "yr", "mth" ) ]
names(dfll.day) = c("date", "dfllc", "dfllesd", "yr", "mth") # Rename for a function later
## ----

# Checking the dominance of large and small phytoplankton
phyl.day = ddply( phyl, .(USI), summarize, bioCtot = sum( bioC, na.rm = T ) )
phy.day = ddply( stationf, .(USI), summarize, bioCtot = sum( bioC, na.rm = T ) )
ind = which( phyl.day$USI %in% phy.day$USI ) ; ind2 = which( phy.day$USI %in% phyl.day$USI )
boxplot( phyl.day$bioCtot[ind] / phy.day$bioCtot[ind2], ylim = c(0, 1),
         ylab = "% of large phytoplankton biomass")
mean(phyl.day$bioCtot[ind] / phy.day$bioCtot[ind2]) # 73 % large phytoplankton over the small species

#### Computing the daily series --------------------------------------------------------------------
phys.day = phy.day.comp( phys )[ c("date", "bioc", "size", "yr", "mth" ) ]
names(phys.day) = c("date", "biocs", "sizes", "yr", "mth") # Rename for a function later
phyl.day = phy.day.comp( phyl )[ c("date", "bioc", "size", "yr", "mth" ) ]
names(phyl.day) = c("date", "biocl", "sizel", "yr", "mth") # Rename for a function later

phy.day = phy.day.comp( stationf ) # For all phytoplankton, used for bloom detection

# Phytoplankton growth (change rate on log bioC) # See Mieruch et al (2010)
phy.rate.day = phy.day$bioc[-1] - phy.day$bioc[-length(phy.day$bioc)]
date.seq = ymd( phy.day$date ) # Dates, to scale the change rate
day.scale = day( days(date.seq[-1]) ) - day( days(date.seq[-length(date.seq)]) )
phy.rate.day = phy.rate.day / day.scale
phy.rate.day = data.frame( rate = phy.rate.day, date = date.seq[-length(date.seq)] )
phy.rate.day$mth = month(phy.rate.day$date) # Adding months and year as variables
phy.rate.day$yr = year(phy.rate.day$date)

day.apx = day( days( max( phy.rate.day$date ) ) ) - day( days( min( phy.rate.day$date ) ) )
day.apx = min( phy.rate.day$date ) + days( 0:day.apx )
phy.rate.apx = approx( x = phy.rate.day$date, y = phy.rate.day$rate,
                       xout = day.apx, method = "linear" )
phy.rate.apx = data.frame( date = phy.rate.apx$x,
                           mth = month( phy.rate.apx$x ),
                           yr = year( phy.rate.apx$x ),
                           rate = phy.rate.apx$y ) # Linear approximation of the change rate

apx.phy.day = function(phy.day.date, phy.day.data){
  date.seq = ymd( phy.day.date )
  day.apx = day( days( max( date.seq ) ) ) - day( days( min( date.seq ) ) )
  day.apx = min( date.seq ) + days( 0:day.apx )
  phy.apx = approx( x = date.seq, y = phy.day.data,
                    xout = day.apx, method = "linear" )
  
  phy.apx = data.frame( date = phy.apx$x,
                        mth = month( phy.apx$x ),
                        yr = year( phy.apx$x ),
                        bioc = phy.apx$y )
return(phy.apx) }

phy.apx = apx.phy.day( phy.day$date, phy.day$bioc )
phys.apx = apx.phy.day( phys.day$date, phys.day$biocs )
phyl.apx = apx.phy.day( phyl.day$date, phyl.day$biocl ) # Linear approximations for missing dates

# Zooplankton
zoo.day.comp = function( zoo ){
  zoo$USI = paste( zoo$Station, "_", zoo$Date, sep="" )
  
  zoo.day = ddply( zoo, .(USI, Date), summarize, 
                   dw = sum( DryW, na.rm=T), 
                   esd = sum( DryW * mean.esd, na.rm = T ) ) # Sum Dry Weight per USI
  
  zoo.day$esd = zoo.day$esd / zoo.day$dw # DW-weighted ESD
  
  zoo.day = ddply( zoo.day, .(Date), summarize, 
                   dw = mean( log(dw +1), na.rm=T ),
                   esd = mean( log(esd), na.rm = T ) ) # Mean Dry Weight per day
  
  zoo.day$yr = year(zoo.day$Date) ; zoo.day$mth = month(zoo.day$Date)
return(zoo.day) }

zoos.day = zoo.day.comp( zoos )
zool.day = zoo.day.comp( zool )
zooc.day = zoo.day.comp( zooc )
#zoohl.day = zoo.day.comp( zoohl )

# Nutrients
nut.day = ddply( waddenEnvtot, .(date), summarize,
                 n = mean( as.numeric( nitrate ), na.rm = T ), 
                 p = mean( as.numeric( phos ), na.rm = T ),
                 spm = mean( as.numeric( suspended.particulates ), na.rm = T ),
                 sal = mean( as.numeric( salinity ), na.rm = T ),
                 si = mean( as.numeric( silicon ), na.rm = T ),
                 temperature = mean( as.numeric( temperature ), na.rm = T ),
                 I0 = mean( as.numeric( I0 ), na.rm = T ) )
nut.day$yr = year(nut.day$date) ; nut.day$mth = month(nut.day$date)

## ----

#### Bloom dates using the change rate of phytoplankton (approximated for daily data) ----
# We use the definition of Trombetta et al (2019), so 5 consecutive days of > 0 growth

mean( day(days( ymd(phy.day$date[-1]) )) - day(days( ymd(phy.day$date[-dim(phy.day)[1]]) )) )
# Mean of 5 days between samples

bloom.detec = function( rate.apx, bloom.count = 5 ){ # How many days of positive growth until a bloom begins
  rate.pos = rate.apx[ which( rate.apx$rate > 0 ), ] # Positive change rates
  
  # Checking the time difference
  day.diff = day( days(rate.pos$date[-1]) ) - day( days(rate.pos$date[-length(rate.pos$date)]) )
  
  rate.pos$diff = c(0, day.diff) # Days difference between positive rates
  
  # Checking for X consecutive >0 for blooms. Ends if <0 for more than X days.
  bl.ind = c()  # Indices
  bl.seas = c() # Seasons
  pos.count = 0
  for( ii in 1:length(rate.pos$date) ){
    idiff = rate.pos$diff[ii]
    
    if( idiff == 1 ){ pos.count = pos.count +1                     # Consecutive positive days
    }else if( idiff > 1 & pos.count < bloom.count ){ pos.count = 0 # Chain break and restart
    }else if( idiff > bloom.count ){ pos.count = 0                  # Bloom collapse
    }
    
    if( pos.count == bloom.count ){ 
      bl.ind = c( bl.ind, c( (ii - (bloom.count -1) ):ii ) ) # Adding bloom start dates
      
      # Checking spring or summer bloom, using June (6) as the limit
      if( month( rate.pos$date[ii- (bloom.count -1) ] ) < 6 ){ seas = "spring" }else{ seas = "autumn" } 
      # Start of the bloom date defines the season
      bl.seas = c( bl.seas, rep( seas, bloom.count ) )
      
    }else if( pos.count >= bloom.count ){ 
      bl.ind = c( bl.ind, ii )  # Adding bloom date
      bl.seas = c( bl.seas, seas )
    }
  }
  
  season.f = data.frame( date = rate.pos$date[bl.ind], season = bl.seas )
  season.f$mth = month( season.f$date) ; season.f$yr = year( season.f$date )
  
return( season.f ) }

season.f = bloom.detec( phy.rate.apx, bloom.count = 7 )
season.f = season.f

# Removing winter maxima
season.f = season.f[ -which(season.f$mth < 3 | season.f$mth > 10), ]

## ----

#### Phytoplankton Monthly series ------------------------------------------------------------------

# Small phy
physm = seriesphy(phys)
phys.clim = physm[[2]] ; physm = physm[[1]]  # Monthly data
physm$bioctrend = filter( log(physm$bioc), filter=filter.wind ) # Smoothed 
physm$szetrend = filter( log(physm$vol), filter=filter.wind )

# Large phy
phylm = seriesphy(phyl)
phyl.clim = phylm[[2]] ; phylm = phylm[[1]]  # Monthly data
phylm$bioctrend = filter( log(phylm$bioc), filter=filter.wind )
phylm$szetrend = filter( log(phylm$vol), filter=filter.wind )

# All phy
phyallm = seriesphy(stationf)
phy.clim = phyallm[[2]] ; phyallm = phyallm[[1]]  # Monthly data

#### Other plankton and environment ----
# Dinoflagellates heterotophs
dflsm = seriesphy(dfls) # Small dinoflagellates
dflsm.clim = dflsm[[2]] ; dflsm = dflsm[[1]]  # Monthly data
dflsm$bioctrend = filter( log(dflsm$bioc), filter=filter.wind )
dflsm$szetrend = filter( log(dflsm$vol), filter=filter.wind )

dfllm = seriesphy(dfll) # Large dinoflagellates
dfllm.clim = dfllm[[2]] ; dfllm = dfllm[[1]]  # Monthly data
dfllm$bioctrend = filter( log(dfllm$bioc), filter=filter.wind )
dfllm$szetrend = filter( log(dfllm$vol), filter=filter.wind )

# Environment
listenvw = envapprox(waddenEnvtot) 
listenvw = trendenv(listenvw)
t.series = listenvw[[1]] # Temperature
n.series = listenvw[[2]] # Nitrogen
p.series = listenvw[[3]] # Phosphorus

# Small zooplankton
zoosm = zooseries(zoos) 
zoos.clim = zoosm[[2]] # Climatology
zoosm = zoosm[[1]]
zoosm$trend = log( filter(zoosm$dataprx, filter=filter.wind) ) # Smoothed biomass

# Large zooplankton
zoolm = zooseries(zool) 
zool.clim = zoolm[[2]] # Climatology
zoolm = zoolm[[1]]
zoolm$trend = log( filter(zoolm$dataprx, filter=filter.wind) )

# Carn zooplankton
zoocm = zooseries(zooc) 
zooc.clim = zoocm[[2]] # Climatology
zoocm = zoocm[[1]]
zoocm$trend = log( filter(zoocm$dataprx, filter=filter.wind) )

# Highest TL zooplankton
# zoohlm = zooseries(zoohl) 
# zoohl.clim = zoohlm[[2]] # Climatology
# zoohlm = zoohlm[[1]]
# zoohlm$trend = log( filter(zoohlm$dataprx, filter=filter.wind) )

## Adding the dates for plotting later
physm$date = ymd(paste(physm$yr, "-", physm$mth, "-15", sep=""))
phylm$date = ymd(paste(phylm$yr, "-", phylm$mth, "-15", sep=""))

dflsm$date = ymd(paste(dflsm$yr, "-", dflsm$mth, "-15", sep=""))
dfllm$date = ymd(paste(dfllm$yr, "-", dfllm$mth, "-15", sep=""))

n.series$date = ymd(paste(n.series$yr, "-", n.series$mth, "-15", sep=""))
p.series$date = ymd(paste(p.series$yr, "-", p.series$mth, "-15", sep=""))

zoosm$date = ymd(paste(zoosm$yr, "-", zoosm$mth, "-15", sep="")) 
zoolm$date = ymd(paste(zoolm$yr, "-", zoolm$mth, "-15", sep=""))
zoocm$date = ymd(paste(zoocm$yr, "-", zoocm$mth, "-15", sep=""))
#zoohlm$date = ymd(paste(zoohlm$yr, "-", zoohlm$mth, "-15", sep=""))
## ----

#### Defining the blooming seasons according to the month with max temperature ----
springVaut.temp = function( season.f, t.series ){
  for( iy in unique( season.f$yr) ){
    
    seasi = season.f[ which( season.f$yr == iy), ]
    ti = t.series[ which( t.series$yr == iy ), ]
    
    max.t = ti$mth[ which.max(ti$dmn) ] # <= Spring ; > Autumn
    # The definition is when the bloom starts
    
    if( length(max.t) == 0 ){
      max.t = 7 # t.clim
    } # If no temperature data available, use the mean max month as a proxy
    max.t = 7 # t.clim
    
    day.diff = day( days(seasi$date[ -1 ]) ) - day( days(seasi$date[ -length(seasi$date) ]) ) 
    
    bl.winds = which(day.diff > 1) +1 # Indices of the different bloom windows
    iin = 1 # Start of the first bloom
    for( iind in bl.winds ){
      if( seasi$mth[iin] < max.t ){ seasi$season[ iin:(iind -1) ] = "spring"
      }else{ seasi$season[ iin:(iind -1) ] = "autumn" }
      
      if( iind == bl.winds[ length(bl.winds) ] ){ 
        seasi$season[ iind:length(seasi$season) ] = "autumn" } # End of year is autumn
      
      iin = iind # Next bloom window
    } 
    
    season.f$season[ which( season.f$yr == iy) ] = seasi$season
  }
return( season.f ) }

season.f = springVaut.temp( season.f, t.series )
#### ----

#### Light data, from Copernicus (CAMS) ----

# Convert SPM (pt) to mg.m-3 (for *1000)
nut.day$lm = nut.day$I0 * ( 1 - exp( -attenuation(nut.day$sal, nut.day$spm*1000) * 10 ) ) /
  attenuation(nut.day$sal, nut.day$spm*1000) / 10

rad = nut.day # Daily resolved data

# Computing climatology
rad$clim.use = 1
rad$clim.use[ which( is.na(rad$lm) ) ] = 0 # Not using months where NA

rad = ddply( rad, .(yr, mth), summarize, 
             l = mean(I0, na.rm = T), 
             lm = mean(lm, na.rm = T),
             nb = sum( mth * clim.use ) ) # Monthly
weight.mth = ddply( rad, .(mth), summarize, totnb = sum( nb ) )
rad = merge( rad, weight.mth, by = "mth" )
rad = rad[ order(rad$yr), ]
rad$nb = rad$nb / rad$totnb

l.clim = ddply( rad, .(mth), summarize, lm = sum( lm * nb, na.rm = T ) ) # W.m-2.h-1

#### Climatologies of the environment ----

t.clim = listenvw[[8]]
n.clim = listenvw[[9]]
p.clim = listenvw[[10]]

spm.series = listenvw[[6]] # Check of SPM climatology
spm.series$date = ymd(paste(spm.series$yr, "-", spm.series$mth, "-15", sep=""))
spm.clim = listenvw[[13]]

sal.series = listenvw[[5]] # Check of Salinity climatology
sal.series$date = ymd(paste(sal.series$yr, "-", sal.series$mth, "-15", sep=""))
sal.clim = listenvw[[12]]

si.series = listenvw[[7]] # Check of Silicates climatology
si.series$date = ymd(paste(si.series$yr, "-", si.series$mth, "-15", sep=""))
si.clim = listenvw[[14]]

##  ----

#### Figure to explain the correlation and TD / BU analysis ----
plot.ideal.rel = function(rel, xlab, ylab, colo, line.xsub=1, line.ysub=1){
  r.pt = seq(-1, 1, length.out = 20)      # Sampling random points to plot a relation
  r.y = rel * (r.pt + rnorm(length(r.pt), 0, 0.2))  # Adding a bit of variation to the points
  
  matplot( x = r.pt, y = r.y, pch = 19, cex = 2.2, col = colo,
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(-1, 1), xlim = c(-1, 1) )
  mod = lm( r.y ~ r.pt )
  abline( mod, lwd = 3, lty=3, col=rgb(0.4, 0.4, 0.4, 0.7) )

  mtext( side = 1, xlab, cex = 1.4, line = line.xsub )
  mtext( side = 2, ylab, cex = 1.4, line = line.ysub )
}

x11(width = 10, height = 8)
par( mfrow = c(2, 3), mar = c(4, 4, 3, 1) )

col.bu = rgb( 0.2, 0.2, 0.9, alpha = 0.8)
col.td = rgb( 0.9, 0.2, 0.2, alpha = 0.8)

matplot(0, 0, axes = F, xlab="", ylab="", type="n")

# Bottom-up, Annual
plot.ideal.rel(1, "Resource", "Consummer", col.bu)

# Top-down, Annual
plot.ideal.rel(1, expression( frac(1,"Predator") ), "Prey", col.td, line.xsub=3.5)

matplot(0, 0, axes = F, xlab="", ylab="", type="n")

# Bottom-up, Season
plot.ideal.rel(1, "Resource", 
               expression(Delta*"Consummer / " *Delta*"t"), col.bu)

# Top-down, Season
plot.ideal.rel(-1, "Predator",
               expression(Delta*"Prey / " *Delta*"t"), col.td)

# Text legends
par( fig = c(0, 1, 0, 1), new = T )
matplot( x = 0, y = 0, type = "n", axes = F, xlab = "", ylab = "", ylim = c(0, 1), xlim = c(0, 1) )

# Subtitles
text( x = 0.5, y = 1.07, "Bottom-up", cex = 2.5, font = 1, col = adjustcolor(col.bu, alpha.f=2), xpd = T )
text( x = 0.9, y = 1.07, "Top-down", cex = 2.5, font = 1, col = adjustcolor(col.td, alpha.f=2), xpd = T )

# Legend
text( x = 0.1, y = 1, "Annual", cex = 2.5, font = 2, xpd = T )
text(x=-0.07, y=0.92, "At steady state:", cex = 2, xpd = T, adj=0)
text(x=-0.07, y=0.85, "dX/dt = r.(1-X/K) = 0 if X = K", cex = 2, xpd = T, adj=0)
text(x=-0.07, y=0.78, expression("K" ~ "\u2191" *" if Resource" ~ "\u2191"), cex = 2, xpd = T, adj=0)
text(x=-0.07, y=0.71, expression("K" ~ "\u2191" *" if Predator" ~ "\u2193"), cex = 2, xpd = T, adj=0)

text( x = 0.1, y = 0.42, "Seasonal", cex = 2.5, font = 2, xpd = T )
text(x=-0.07, y=0.35, "dX/dt = Resource - Predation", cex = 2, xpd = T, adj=0)
text(x=-0.07, y=0.28, expression("dX/dt" ~ "\u2191" *" if Resource" ~ "\u2191"), cex = 2, xpd = T, adj=0)
text(x=-0.07, y=0.21, expression("dX/dt" ~ "\u2193" *" if Predator" ~ "\u2191"), cex = 2, xpd = T, adj=0)
#### ---- 

# Checking the mean windows for the seasons (blooms and low)----
date.sp = season.f$date[ which( season.f$season == "spring" ) ]
date.at = season.f$date[ which( season.f$season == "autumn" ) ]

mth.frame = data.frame()
for( yri in unique(season.f$yr) ){
  date.sp.i = date.sp[ which( year(date.sp) == yri ) ]
  date.at.i = date.at[ which( year(date.at) == yri ) ]
  
  min.sp = min( month(date.sp.i) ) ; max.sp = max( month(date.sp.i) )
  min.at = min( month(date.at.i) ) ; max.at = max( month(date.at.i) )
  
  mth.frame = rbind( mth.frame, data.frame( min.sp = min.sp, max.sp = max.sp,
                                            min.at = min.at, max.at = max.at ) )
}

min.sp = mean( mth.frame$min.sp, na.rm = T )
max.sp = mean( mth.frame$max.sp, na.rm = T )

min.at = mean( mth.frame$min.at, na.rm = T )
max.at = mean( mth.frame$max.at, na.rm = T )
# ----

# Plot smooth lines
smooth.lines = function(points, data, colo, span=0.05, points.p = T ){
  curve = c()
  xaxis = c()
  
  predni = seq( min(points), max(points), length.out = 1000) # Prediction points
  
  model = loess( data ~ points, span = span )
  smooth = predict( model, newdata = predni )
  
  
  matplot( predni, smooth, type = "l", col = colo, lwd = 3.5, add = T)

  if(points.p){
    matplot( points, data, pch=19, cex = 2, col = "white", add = T ) 
    matplot( points, data, pch=19, cex = 1.5, col = colo, add = T ) 
  }
} 

#### Plot of the climatologies of the plankton system ----
mth.lab = c( "Jan", "Feb", "Mar", "Apr", "May", "June", "Jul",
             "Aug", "Sep", "Oct", "Nov", "Dec" )

axis.scale = function( data, side, off = 0, ... ){
  scale.data = scale( data )
  axis.data = ( pretty( range(data, na.rm = T) ) - attr(scale.data, "scaled:center") ) /
    attr( scale.data, "scaled:scale" )
  axis( side, at = axis.data +off, labels = pretty( range(data, na.rm = T) ), ... )
}

x11(height = 12, width = 6)
par( mar = c( 2.5, 6, 2, 6 ) )

matplot( x= c(0, 1), y = c(0, 1), xlim = c(1, 12), ylim = c(-2, 18.5), 
         type = "n", xlab = "", ylab = "", axes = F )

# Polygons of blooms
alp = 0.2
col.sp = rgb(0.2, 0.8, 0.1, alpha = alp)
col.sw = rgb(0.8, 0.6, 0.4, alpha = alp)

x.pol.sp = c( min.sp, max.sp, max.sp, min.sp )
x.pol.sw = c( max.sp, min.at, min.at, max.sp )
y.pol = c( -10, -10, 30, 30 ) 

polygon( x = x.pol.sp, y = y.pol, col = col.sp, border = NA )
polygon( x = x.pol.sw, y = y.pol, col = col.sw, border = NA )

axis( 1, at = 0.5:12.5, labels = F, cex.axis = 1.5, lwd=2 )
axis( 1, at = 1:12, labels = mth.lab, cex.axis = 1.5, tck=F )

smooth.lines( n.clim$mth, scale( n.clim$dmn ), colo = "darkblue", span = 0.2 )

# Normalized SiO2 to plot on the same plot as NO3-
n.scale = scale( n.clim$dmn )
si.clim.plot = (si.clim$dmn - attr(n.scale, 'scaled:center'))/attr(n.scale, 'scaled:scale')
smooth.lines( si.clim$mth, si.clim.plot, colo = "gray60", span = 0.2 )
axis.scale( n.clim$dmn, side = 2, cex.axis = 1.5, las = 1 )
mtext( side = 2, expression( "NO"[3]^-"" ),
       line = 4.2, cex = 1.5, adj = 0.07, col = "darkblue" )
mtext( side = 2, 'and', line = 4.5, cex = 1.5, adj = 0.13, col = "black" )
mtext( side = 2, expression( "SiO"[2] ), 
       line = 4.2, cex = 1.5, adj = 0.19, col = "gray40" )
mtext( side = 2, expression( mu *"mol.l"^-1 ), 
       line = 2.5, cex = 1.5, adj = 0.12, col = "black" )
# mtext( side = 2, expression( "NO"[3]^-"" *", "*mu *"mol.l"^-1 ), 
#        line = 3.2, cex = 1.5, adj = 0.05, col = "darkblue" )

smooth.lines( p.clim$mth, scale( p.clim$dmn ), colo = "dodgerblue", span = 0.2 )
axis.scale( p.clim$dmn, side = 4, cex.axis = 1.5, las = 1 )
mtext( side = 4, expression( "PO"[4]^3^-"" *", "*mu *"mol.l"^-1 ), 
       line = 4.5, cex = 1.5, adj = 0.05, col = "dodgerblue" )

off = 3.8 # Offset for plotting the second row of scaled climatologies
smooth.lines( phys.clim$mth, scale( phys.clim$bmean ) +off, colo = "chartreuse2", span = 0.2 )
axis.scale( phys.clim$bmean, side = 2, off = off, cex.axis = 1.5, las = 1 )
mtext( side = 2, expression( "Small, "*mu *"gC.l"^-1 ), 
       line = 3.2, cex = 1.5, adj = 0.32, col = "chartreuse3" )

smooth.lines( phyl.clim$mth, scale( phyl.clim$bmean ) +off, colo = "darkolivegreen", span = 0.2 )
axis.scale( phyl.clim$bmean, side = 4, off = off, cex.axis = 1.5, las = 1 )
mtext( side = 4, expression( "Large, "*mu *"gC.l"^-1 ), 
       line = 4.5, cex = 1.5, adj = 0.32, col = "darkolivegreen" )

off = 8.8 # Offset for plotting the fourth
smooth.lines( dflsm.clim$mth, scale( dflsm.clim$bmean ) +off, colo = "black", span = 0.2 )
axis.scale( dflsm.clim$bmean, side = 2, off = off, cex.axis = 1.5, las = 1 )
mtext( side = 2, expression( mu *"gC.l"^-1 ), 
       line = 3.2, cex = 1.5, adj = 0.56, col = "black" )

off = 13.3 # Offset for plotting the fourth
smooth.lines( zoos.clim$mth, scale( zoos.clim$dw ) +off, colo = "orange", span = 0.3 )
axis.scale( zoos.clim$dw, side = 2, off = off, cex.axis = 1.5, las = 1 )
mtext( side = 2, expression( "Small, mg.m"^-3 ), 
       line = 3.2, cex = 1.5, adj = 0.77, col = "orange" )

smooth.lines( zool.clim$mth, scale( zool.clim$dw ) +off, colo = "darkorange4", span = 0.3 )
axis.scale( zool.clim$dw, side = 4, off = off, cex.axis = 1.5, las = 1 )
mtext( side = 4, expression( "Large, mg.m"^-3 ), 
       line = 4.5, cex = 1.5, adj = 0.77, col = "darkorange4" )

off = 17 # Offset for plotting the fifth
smooth.lines( zooc.clim$mth, scale( zooc.clim$dw ) +off, colo = "darkred", span = 0.3 )
axis.scale( zooc.clim$dw, side = 2, off = off, cex.axis = 1.5, las = 1 )
mtext( side = 2, expression( "Carn., mg.m"^-3 ), 
       line = 3.2, cex = 1.5, adj = 1, col = "darkred" )

# Text legends
text( x = 6, y = 1.8, "Nutrient", cex = 1.5 )
text( x = 7, y = 6.5, "Phytoplankton", cex = 1.5 )
text( x = 6.5, y = 11, "Dinoflagellates", cex = 1.5 )
text( x = 6.5, y = 15.4, "Zooplankton", cex = 1.5 )

col.sp = rgb(0.2, 0.8, 0.1)
col.sw = rgb(0.8, 0.6, 0.4)
text( x = 4.5, y = 19, "Spring", col = adjustcolor(col.sp, r=0.8, g=0.8, b=0.8, alpha.f=5), cex = 1.5, xpd = T )
text( x = 7.6, y = 19, "Summer", col = adjustcolor(col.sw, r=0.8, g=0.8, b=0.8, alpha.f=5), cex = 1.5, xpd = T )
#### ----

#######################################################
### Mean taxonomic groups in each zooplankton level ### -------------------------------------
#######################################################

dom.taxa = function(usi, taxa, biomass){
  datf = data.frame(USI = usi, Taxa = taxa, bio = biomass)
  df = ddply( datf, .(USI, Taxa), summarize, biotax = sum( bio, na.rm = T ) )
  
  df.day = ddply( df, .(USI), summarize, biotot = sum( biotax, na.rm = T ) )
  df.day = merge(df, df.day, by='USI')
  
  df.day$w = df.day$biotax / df.day$biotot # Weight per sample
  taxa.day = df.day
  
  tot.obs = length(unique(df$USI))         # Number of sample
  df.day$w = df.day$w / tot.obs            # Biomass fraction weighted by the frequency of observation
  
  taxa.fraction = ddply( df.day, .(Taxa), summarize, 
                         w = sum(w, na.rm = T) )
  taxa.fraction = taxa.fraction[ order(taxa.fraction$w), ]
  
  return( list(taxa.fraction, taxa.day) )
}

dom.phys = dom.taxa(phys$USI, phys$phylum, phys$bioC)[[1]]
dom.phyl = dom.taxa(phyl$USI, phyl$phylum, phyl$bioC)[[1]]
dom.zoos = dom.taxa(zoos$USI, zoos$Taxa, zoos$DryW)[[1]]
dom.zool = dom.taxa(zool$USI, zool$Taxa, zool$DryW)[[1]]
dom.zooc = dom.taxa(zooc$USI, zooc$Taxa, zooc$DryW)[[1]]

dom.phys$group = "Small phytoplankton"
dom.phyl$group = "Large phytoplankton"
dom.zoos$group = "Small zooplankton"
dom.zool$group = "Large zooplankton"
dom.zooc$group = "Carnivorous zooplankton"

dom.zoo.all = rbind(dom.phys, dom.phyl, dom.zoos, dom.zool, dom.zooc)
dom.zoo.all = dom.zoo.all[ c("group", "Taxa", "w") ] # Mean taxonomic composition of each trophic group.

#####################################
### Annual means and correlations ### ----------------------------------------------------
#####################################

## Annual means using monthly series ----
n.anu = ddply(n.series, .(yr), summarize, nm = mean( dmn, na.rm = T ) )
p.anu = ddply(p.series, .(yr), summarize, pm = mean( dmn, na.rm = T ) )
si.anu = ddply(si.series, .(yr), summarize, sim = mean( dmn, na.rm = T ) )
spm.anu = ddply(spm.series, .(yr), summarize, spm = mean( dmn, na.rm = T ) )
sal.anu = ddply(sal.series, .(yr), summarize, sal = mean( dmn, na.rm = T ) )
l.anu = ddply(rad, .(yr), summarize, l = mean( lm, na.rm = T ) )
t.anu = ddply(t.series, .(yr), summarize, tm = mean( dmn, na.rm = T ) )

phytos.anu = ddply(physm, .(yr), summarize, 
                  biocs = mean( log( bioc ), na.rm = T ),
                  #biocs = mean( bioc, na.rm = T ),
                  esd.phys = mean( log( (vol * 3/4/pi)**(1/3)*2 ) ) )
phytol.anu = ddply(phylm, .(yr), summarize, 
                   biocl = mean( log( bioc ), na.rm = T ),
                   #biocl = mean(bioc, na.rm = T ),
                   esd.phyl = mean( log( (vol * 3/4/pi)**(1/3)*2 ) ) )
phytall.anu = ddply(phyallm, .(yr), summarize, 
                    biocl = mean( log( bioc ), na.rm = T ),
                    esd.phyl = mean( log( (vol * 3/4/pi)**(1/3)*2 ) ) )

dfls.anu = ddply(dflsm, .(yr), summarize, 
                 dflcs = mean( log( bioc +1), na.rm = T ),
                 #dflcs = mean( bioc, na.rm = T ),
                 esd.dfls = mean( log( (vol * 3/4/pi)**(1/3)*2 ) ) )
dfll.anu = ddply(dfllm, .(yr), summarize, 
                 dflcl = mean( log( bioc +1), na.rm = T ),
                 #dflcl = mean( bioc, na.rm = T ),
                 esd.dfll = mean( log( (vol * 3/4/pi)**(1/3)*2 ) ) )

zoos.anu = ddply(zoosm, .(yr), summarize, 
                 dws = mean( log( dw +1), na.rm=T),
                 #dws = mean( dw, na.rm=T),
                 esds = mean( log( (sze * 3/4/pi)**(1/3)*2 ), na.rm = T) )
zool.anu = ddply(zoolm, .(yr), summarize, 
                 dwl = mean( log( dw +1), na.rm=T),
                 #dwl = mean( dw, na.rm=T),
                 esdl = mean( log( (sze * 3/4/pi)**(1/3)*2 ), na.rm = T) )
ind.remove = which(zoocm$yr > 2015) # Extreme point first year 
zooc.anu = ddply(zoocm[ind.remove,], .(yr), summarize, 
                 dwc = mean( log( dw +1), na.rm=T),
                 #dwc = mean( dw, na.rm=T),
                 esdc = mean( log( (sze * 3/4/pi)**(1/3)*2 ), na.rm = T) )

# ind.remove = which(zoohlm$yr > 2015) # Extreme point first year 
# zoohl.anu = ddply(zoohlm[ind.remove,], .(yr), summarize, 
#                   dwc = mean( log( dw +1), na.rm=T),
#                   esdc = mean( log( (sze * 3/4/pi)**(1/3)*2 ), na.rm = T) )

fish.anu = ddply( dataf, .(year), summarize, 
                  ft = mean( log(t +1), na.rm = T ) )
## ----

## Checking correlations on the annual scale ----
anu.data = phytos.anu # New data frame containing everything
anu.data$biocl = NA ; anu.data$biocl[which(anu.data$yr %in% phytol.anu$yr)] = phytol.anu$biocl
anu.data$esd.phyl = NA ; anu.data$esd.phyl[which(anu.data$yr %in% phytol.anu$yr)] = phytol.anu$esd.phyl

anu.data$bioca = NA ; anu.data$bioca[which(anu.data$yr %in% phytall.anu$yr)] = phytall.anu$biocl
anu.data$esd.phya = NA ; anu.data$esd.phya[which(anu.data$yr %in% phytall.anu$yr)] = phytall.anu$esd.phyl

anu.data$dfllc = NA ; anu.data$dfllc[which(anu.data$yr %in% dfll.anu$yr)] = dfll.anu$dflcl
anu.data$esd.dfll = NA ; anu.data$esd.dfll[which(anu.data$yr %in% dfll.anu$yr)] = dfll.anu$esd.dfll

anu.data$dflsc = NA ; anu.data$dflsc[which(anu.data$yr %in% dfls.anu$yr)] = dfls.anu$dflcs
anu.data$esd.dfls = NA ; anu.data$esd.dfls[which(anu.data$yr %in% dfls.anu$yr)] = dfls.anu$esd.dfls

anu.data$temp = NA ; anu.data$temp[which(anu.data$yr %in% t.anu$yr)] = t.anu$tm
anu.data$n = NA ; anu.data$n[which(anu.data$yr %in% n.anu$yr)] = n.anu$nm
anu.data$p = NA ; anu.data$p[which(anu.data$yr %in% p.anu$yr)] = p.anu$pm
anu.data$si = NA ; anu.data$si[which(anu.data$yr %in% si.anu$yr)] = si.anu$sim
anu.data$spm = NA ; anu.data$spm[which(anu.data$yr %in% spm.anu$yr)] = spm.anu$spm
anu.data$sal = NA ; anu.data$sal[which(anu.data$yr %in% sal.anu$yr)] = sal.anu$sal
anu.data$lm = NA ; anu.data$lm[which(anu.data$yr %in% l.anu$yr)] = l.anu$l # Last NaN disarded

anu.data$s.dw = NA ; anu.data$s.dw[which(anu.data$yr %in% zoos.anu$yr)] = zoos.anu$dws
anu.data$esds = NA ; anu.data$esds[which(anu.data$yr %in% zoos.anu$yr)] = zoos.anu$esds

anu.data$l.dw = NA ; anu.data$l.dw[which(anu.data$yr %in% zool.anu$yr)] = zool.anu$dwl
anu.data$esdl = NA ; anu.data$esdl[which(anu.data$yr %in% zool.anu$yr)] = zool.anu$esdl

anu.data$c.dw = NA ; anu.data$c.dw[which(anu.data$yr %in% zooc.anu$yr)] = zooc.anu$dwc
anu.data$esdc = NA ; anu.data$esdc[which(anu.data$yr %in% zooc.anu$yr)] = zooc.anu$esdc

#anu.data$hl.dw = NA ; anu.data$hl.dw[which(anu.data$yr %in% zoohl.anu$yr)] = zoohl.anu$dwc
#anu.data$esdhl = NA ; anu.data$esdhl[which(anu.data$yr %in% zoohl.anu$yr)] = zoohl.anu$esdc

anu.data$fish = NA ; anu.data$fish[which(anu.data$yr %in% fish.anu$year)] = fish.anu$ft[which(fish.anu$year %in% anu.data$yr)]

# Adding Fish planktivores with one year lag 
anu.data$flag = NA
anu.data$flag[which(anu.data$yr %in% ( fish.anu$year +1) )] = fish.anu$ft[which( (fish.anu$year +1) %in% anu.data$yr)]

nvar = ncol(anu.data) -1 
cor.frame.anu = matrix(nrow = nvar, ncol = nvar) ; names.test = names(anu.data)[-1]
rownames(cor.frame.anu) = names.test ; colnames(cor.frame.anu) = names.test

p.frame.anu = matrix(nrow = nvar, ncol = nvar) 
rownames(p.frame.anu) = names.test ; colnames(p.frame.anu) = names.test

for( ii in ( 1:length(names.test) ) ){
  for( ij in ( 1:length(names.test) ) ){
    v1 = anu.data[, ii+1] ; v2 = anu.data[, ij+1]
    ctc = cor.test( v1, v2, method = "spearman", use = "complete.obs" )
    cor.frame.anu[ii, ij] = round(ctc$estimate, 2)
    p.frame.anu[ii, ij] = ctc$p.value
  } } # Computing spearman correlations and p.values

# Correlations using the inverse
nvar = ncol(anu.data) -1 
corinv.frame.anu = matrix(nrow = nvar, ncol = nvar) ; names.test = names(anu.data)[-1]
rownames(corinv.frame.anu) = names.test ; colnames(corinv.frame.anu) = names.test

pinv.frame.anu = matrix(nrow = nvar, ncol = nvar) 
rownames(pinv.frame.anu) = names.test ; colnames(pinv.frame.anu) = names.test

for( ii in ( 1:length(names.test) ) ){
  for( ij in ( 1:length(names.test) ) ){
    v1 = anu.data[, ii+1] ; v2 = anu.data[, ij+1]
    ctc = cor.test( 1/v1, v2, method = "spearman", use = "complete.obs" )
    corinv.frame.anu[ii, ij] = round(ctc$estimate, 2)
    pinv.frame.anu[ii, ij] = ctc$p.value
  } } # Computing spearman correlations and p.values
## ----

# Correlations on annual means (spearman) are available in cor.frame.anu
# The p-values are in p.frame.anu

plot.table.anu = function( corrdt, corrinvdt, symetric = T ){
  corr = round( as.data.frame( corrdt[[1]] ), 2) # Extracting the R2 and p-values
  p = as.data.frame( corrdt[[2]] )
  
  corrinv = round( as.data.frame( corrinvdt[[1]] ), 2) # Extracting the R2 and p-values
  pinv = as.data.frame( corrinvdt[[2]] )
  
  ## Selecting the variables to extract (season dataframes)
  rown = c( "lm", "n", "p", "si", "biocs", "biocl",
            "dflsc", "dfllc", "s.dw", "l.dw", "c.dw", 'flag' ) # Rows to look at
  cn = rown
  
  if(!symetric){cn = rown}
  
  corr = corr[ rown, cn ]
  p = p[ rown, cn ]
  
  corrinv = corrinv[ rown, cn ]
  pinv = pinv[ rown, cn ]
  
  # Rename the variables
  rown = c( "Light", "NO3", "PO4", "SiO2", "S.phy", "L.phy",
            "S.dfl", "L.dfl", "S.zoo", "L.zoo", "C.zoo", "Herring larvae" )
  rownames(corr) = rown
  names(corr) = rown
  
  if(!symetric){names(corr) = rown}
  
  ## Prepare the table, and replace with 1/X for 
  corrf = corr
  pf = p
  
  # Prints only one side of the correlation table
  if( symetric ){
    ind = which( row( corrf ) <= col(corrf), arr.ind = T )
    corrf[ind] = "-"
    pf[ind] = 1}
  
  # Add predator 1/X
  corrf[2:4, 5:6] = corrinv[2:4, 5:6] # Nutrients and phytoplankton
  corrf[5, 7] = corrinv[5, 7]   # Dinoflagellates and Phytoplankton
  corrf[5:8, 9] = corrinv[5:8, 9]   # S zoo
  corrf[5:9, 10] = corrinv[5:9, 10] # L zoo
  corrf[7:10, 11] = corrinv[7:10, 11] # C zoo
  corrf[7:10, 12] = corrinv[7:10, 12] # Herr
  
  pf[2:4, 5:6] = pinv[2:4, 5:6] # Nutrients and phytoplankton
  pf[5, 7] = pinv[5, 7]   # Dinoflagellates and Phytoplankton
  pf[5:8, 9] = pinv[5:8, 9]   # S zoo
  pf[5:9, 10] = pinv[5:9, 10] # L zoo
  pf[7:10, 11] = pinv[7:10, 11] # C zoo
  pf[7:10, 12] = pinv[7:10, 12] # Herr
  
  ## Print the significance level on the table
  # p < 0.1
  ind = which(pf < 0.1 & pf > 0.05, arr.ind = T)
  corrf[ind] = paste( corrf[ind], "*", sep = "" )
  
  # p < 0.05
  ind = which(pf < 0.05 & pf > 0.01, arr.ind = T)
  corrf[ind] = paste( corrf[ind], "**", sep = "" )
  
  # p < 0.01
  ind = which(pf < 0.01, arr.ind = T)
  corrf[ind] = paste( corrf[ind], "***", sep = "" )
  
  corrf=corrf[-1,]
  
  ## Print the table
  x11(width = 12, height = 8)
  grid.table(corrf)
  
  return(corrf)
}

corr.anu.f = plot.table.anu( list(cor.frame.anu, p.frame.anu), list(corinv.frame.anu, pinv.frame.anu) )

### Annual time series ----
smooth.lines = function(points, data, colo, span=0.05, points.p = T ){
  curve = c()
  xaxis = c()
  
  predni = seq( min(points), max(points), length.out = 1000) # Prediction points
  
  model = loess( data ~ points, span = span )
  smooth = predict( model, newdata = predni )
  
  matplot( points, data, type = "l", lty = 1, col = colo, lwd = 2, add = T)
  
  if(points.p){
    matplot( points, data, pch=19, cex = 2, col = "white", add = T ) 
    matplot( points, data, pch=19, cex = 1.5, col = colo, add = T ) 
  }
  #return( list( points, smooth ) )
} # Plot smooth climatologies (plot only)

x11(width = 6, height = 10) # Plotting annual data ----

## Time series
par( mar = c(2, 3, 2, 2), mgp = c(3, 0.3, 0) ) 

matplot( x = c(0, 1), y = c(0, 1), ylim = c(-1, 18), xlim = c(2013, 2024), 
         xlab = "", ylab = "", axes = F )
y.axs = (2013:2024) ; y.labs = y.axs
y.labs[ seq(2, length(y.axs), 2) ] = NA
x.axis = abline( v = y.labs, col = "gray80", lty = 2, lwd = 0.5 ) # Time labels
axis( 1, at = y.axs, labels = y.labs, cex.axis = 1.5, tck = 0.02, hadj=0.1 )
axis( 3, at = y.axs, labels = y.labs, cex.axis = 1.5, tck = 0.02, hadj=0.1 )

stackl = c(0, 5, 8.4, 12, 16) # Y axis tick labels
for(ssi in stackl){
  if( ssi == 0 | ssi == max(stackl) ){
    labs = c(-1, NA, 1)
  }else{
    labs = rep(NA, 3)
  }
  y.axis = axis(2, at = c(-1, 0, 1) + ssi, labels = labs, cex.axis = 1.7, las = 1, font=1, tck = -0.01)
  abline( h = ssi, lty = 2, lwd = 0.5, col = "gray80" )
}

# Phytoplankton and nutrients
smooth.lines(anu.data$yr, scale(anu.data$n), colo = "darkblue", span = 0.2 )
smooth.lines(anu.data$yr, scale(anu.data$p), colo = "dodgerblue", span = 0.2 )
smooth.lines(anu.data$yr, scale(anu.data$si), colo = "gray60", span = 0.2 )

smooth.lines(anu.data$yr, scale(anu.data$biocs) +stackl[2], colo = "chartreuse1", span = 0.2 )
smooth.lines(anu.data$yr, scale(anu.data$biocl) +stackl[2], colo = "darkolivegreen", span = 0.2 )

# Dinoflagellates
smooth.lines(anu.data$yr, scale(anu.data$dflsc) +stackl[3], colo = "black", span = 0.4 )

# Zooplankton 
smooth.lines(anu.data$yr, scale(anu.data$s.dw) +stackl[4], colo = "orange", span = 0.4 )
smooth.lines(anu.data$yr, scale(anu.data$l.dw) +stackl[4], colo = "darkorange4", span = 0.4 )

# High trophic levels
smooth.lines(anu.data$yr, scale(anu.data$c.dw) +stackl[5], colo = "darkred", span = 0.4 )
smooth.lines(anu.data$yr, scale(anu.data$flag) +stackl[5], span = 0.3, colo = "gray50" ) # Larvae

# Text labels
mtext( side = 2, "Normalized Log Annual Concentration", cex = 1.5, line = 1.6 )

# Trophic labels
# text( x = 2016.5, y = 1.2, 
#       expression( "NO"[3]^-"" ), col = "darkblue", cex = 1.5, xpd = T )
# text( x = 2018.5, y = -0.2, 
#       expression( "PO"[4]^3-"" ), col = "dodgerblue", cex = 1.5, xpd = T )
text( x = 2014.2, y = 0.5, 
      expression( "NO"[3]^"-" ), col = "darkblue", cex = 1.5, xpd = T )
text( x = 2016.5, y = 1.6, 
      expression( "SiO"[2] ), col = "gray40", cex = 1.5, xpd = T )
text( x = 2015.8, y = -0.4, 
      expression( "PO"[4]^"3-" ), col = "dodgerblue", cex = 1.5, xpd = T )

text( x = 2021, y = 4, "Small phyto", col = "chartreuse3", cex = 1.5, xpd = T )
text( x = 2014, y = 6.2, "Large\nphytoplankton", col = "darkolivegreen",
      cex = 1.5, adj = 0, xpd = T )

text( x = 2015, y = 8.5, "Dinoflagellates", col = "black",
      cex = 1.5, adj = 0, xpd = T )

text( x = 2013, y = 11.5, "Large\nzooplankton", col = "darkorange4", cex = 1.5, 
      adj = 0, xpd = T )
text( x = 2018.5, y = 13.7, "Small zoo", col = "orange", cex = 1.5, xpd = T )

text( x = 2018.5, y = 18, "Carnivorous zoo", col = "darkred", cex = 1.5, xpd = T )
text( x = 2014.4, y = 17.8, "Herring\nlarvae", col = "gray30", cex = 1.5, xpd = T )
### ----

################################################################
## Plankton and their bottom / top regulations during seasons ## ---------------------------------------------
################################################################

#### New way to sample and compute the rates ####

# Data extraction function
impute.nan = function(date.seas, v.seas, clim.mth, clim.v, date.series, v.series){
  ind = which( is.na(v.seas) )
  
  if( length(ind) > 0 ){
    v.approx = approx(x = ymd(date.series), y = v.series, xout = ymd(date.seas))$y
    
    for(ii in ind){
      ind.clim = which( clim.mth == month(date.seas[ii]) )
      v.seas[ii] = ( clim.v[ind.clim] + v.approx[ii] ) / 2
    }
  }
  return(v.seas)
}

fill.nan.days = function(frame){
  datei = frame$date
  for( vi in names(frame) ){
    if(vi == "biocs"){
      vclim = phys.clim
      vday = phys.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$bmean), vday$date, vday$biocs)
    }
    if(vi == "biocl"){
      vclim = phyl.clim
      vday = phyl.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$bmean), vday$date, vday$biocl)
    }
    if(vi == "dflsc"){
      vclim = dflsm.clim
      vday = dfls.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$bmean), vday$date, vday$dflsc)
    }
    if(vi == "dfllc"){
      vclim = dfllm.clim
      vday = dfll.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$bmean), vday$date, vday$dfllc)
    }
    if(vi == "s.dw"){
      vclim = zoos.clim
      vday = zoos.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$dw+1), vday$Date, vday$dw)
    }
    if(vi == "l.dw"){
      vclim = zool.clim
      vday = zool.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$dw+1), vday$Date, vday$dw)
    }
    if(vi == "c.dw"){
      vclim = zooc.clim
      vday = zooc.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$dw+1), vday$Date, vday$dw)
    }
    if(vi == "hl.dw"){
      vclim = zoohl.clim
      vday = zoohl.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth,  log(vclim$dw+1), vday$Date, vday$dw)
    }
    if(vi == "n"){
      vclim = n.clim
      vday = nut.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth, vclim$dmn, vday$date, vday$n)
    }
    if(vi == "p"){
      vclim = p.clim
      vday = nut.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth, vclim$dmn, vday$date, vday$p)
    }
    if(vi == "lm"){
      vclim = l.clim
      vday = nut.day
      frame[vi] = impute.nan(datei, frame[[vi]], vclim$mth, vclim$lm, vday$date, vday$lm)
    }
  }
  return(frame)
}

season.extract.days = function(season.f, date.extend=0){
  
  season.f$yr = year( season.f$date ) #; low.dates = ret$date[s.dates] # Classic frame
  sp.sum = data.frame() ; at.sum = data.frame() # Placeholders
  sw.sum = data.frame()
  
  for( yri in unique(season.f$yr) ){
    t.i = season.f[ which( season.f$yr == yri ), ]
    
    # Using months as the sampling frequency is not the same for all series
    date.sp = t.i$date[ which( t.i$season == "spring" ) ] # Spring
    date.at = t.i$date[ which( t.i$season == "autumn" ) ] # Autumn

    # Summer period 
    if( length(date.sp)>0 & length(date.at)>0 ){
      dur.sw = day( days( min( date.at ) ) ) - day( days( max( date.sp ) ) ) -1
      date.sw = max(date.sp) + 1:dur.sw
    }else{
      date.sw=c()
    }
    
    ## Extend a bit the seasons, due to biweekly sampling
    if( length(date.sp)>0 ){ 
      date.sp = as_date( ( min(date.sp)-days(date.extend) ):( max(date.sp)+days(date.extend) ) )
    }
    if( length(date.at)>0 ){
      date.at = as_date( ( min(date.at)-days(date.extend) ):( max(date.at)+days(date.extend) ) )
    }
    if( length(date.sw)>0 ){
      date.sw = as_date( ( min(date.sw)-days(date.extend) ):( max(date.sw)+days(date.extend) ) )
    }
    
    ## Extracting variables
    # Dates for phytoplankton small
    phys.sp = phys.day[ which( phys.day$date %in% date.sp ), ]
    
    phys.at = phys.day[ which( phys.day$date %in% date.at ), ]
    
    phys.sw = phys.day[ which( phys.day$date %in% date.sw ), ]
    
    # Dates for phytoplankton large
    phyl.sp = phyl.day[ which( phyl.day$date %in% date.sp ), ]
    
    phyl.at = phyl.day[ which( phyl.day$date %in% date.at ), ]
    
    phyl.sw = phyl.day[ which( phyl.day$date %in% date.sw ), ]
    
    # Dates for dinoflagellates small feeders
    dfls.sp = dfls.day[ which( dfls.day$date %in% date.sp ), ]
    
    dfls.at = dfls.day[ which( dfls.day$date %in% date.at ), ]
    
    dfls.sw = dfls.day[ which( dfls.day$date %in% date.sw ), ]
    
    # Dates for dinoflagellates large feeders
    dfll.sp = dfll.day[ which( dfll.day$date %in% date.sp ), ]
    
    dfll.at = dfll.day[ which( dfll.day$date %in% date.at ), ]
    
    dfll.sw = dfll.day[ which( dfll.day$date %in% date.sw ), ]
    
    # Dates for nutrients
    nut.sp = nut.day[ which( nut.day$date %in% date.sp ), ]
    names(nut.sp) = c( "date", "n", "p", "spm", "sal", "si", "temp", "l", "yr", "mth", "lm" )
    
    nut.at = nut.day[ which( nut.day$date %in% date.at ), ]
    names(nut.at) = c( "date", "n", "p", "spm", "sal", "si", "temp", "l", "yr", "mth", "lm" )
    
    nut.sw = nut.day[ which( nut.day$date %in% date.sw ), ]
    names(nut.sw) = c( "date", "n", "p", "spm", "sal", "si", "temp", "l", "yr", "mth", "lm" )
    
    # Dates for zooplankton small
    zoos.sp = zoos.day[ which( zoos.day$Date %in% date.sp ), ]
    names(zoos.sp) = c( "date", "s.dw", "s.esd", "yr", "mth" )
    
    zoos.at = zoos.day[ which( zoos.day$Date %in% date.at ), ]
    names(zoos.at) = c( "date", "s.dw", "s.esd", "yr", "mth" )

    zoos.sw = zoos.day[ which( zoos.day$Date %in% date.sw ), ]
    names(zoos.sw) = c( "date", "s.dw", "s.esd", "yr", "mth" )
    
    # Large zooplankton
    zool.sp = zool.day[ which( zool.day$Date %in% date.sp ), ]
    names(zool.sp) = c( "date", "l.dw", "l.esd", "yr", "mth" )
    
    zool.at = zool.day[ which( zool.day$Date %in% date.at ), ]
    names(zool.at) = c( "date", "l.dw", "l.esd", "yr", "mth" )
    
    zool.sw = zool.day[ which( zool.day$Date %in% date.sw ), ]
    names(zool.sw) = c( "date", "l.dw", "l.esd", "yr", "mth" )
    
    # Carn zooplankton
    zooc.sp = zooc.day[ which( zooc.day$Date %in% date.sp ), ]
    names(zooc.sp) = c( "date", "c.dw", "c.esd", "yr", "mth" )
    
    zooc.at = zooc.day[ which( zooc.day$Date %in% date.at ), ]
    names(zooc.at) = c( "date", "c.dw", "c.esd", "yr", "mth" )
    
    zooc.sw = zooc.day[ which( zooc.day$Date %in% date.sw ), ]
    names(zooc.sw) = c( "date", "c.dw", "c.esd", "yr", "mth" )
    
    # HTL zooplankton
    # zoohl.sp = zoohl.day[ which( zoohl.day$Date %in% date.sp ), ]
    # names(zoohl.sp) = c( "date", "hl.dw", "hl.esd", "yr", "mth" )
    # 
    # zoohl.at = zoohl.day[ which( zoohl.day$Date %in% date.at ), ]
    # names(zoohl.at) = c( "date", "hl.dw", "hl.esd", "yr", "mth" )
    # 
    # zoohl.sw = zoohl.day[ which( zoohl.day$Date %in% date.sw ), ]
    # names(zoohl.sw) = c( "date", "hl.dw", "hl.esd", "yr", "mth" )
    
    # Creating temporal variables for storing in a dataframe
    mthsp = month(date.sp) ; mthat = month(date.at) ; mthsw = month(date.sw)
    
    # Aggregating everything in a common table
    # Spring ----
    if( length(date.sp) > 0 & !is.infinite(date.sp)[1] ){
      all.sp = phys.sp ; all.sp$yr = yri
      all.sp = merge( all.sp, phyl.sp[ c("biocl", "sizel", "date") ], by = "date", all.x = T )
      all.sp = merge( all.sp, dfls.sp[ c("dflsc", "dflsesd", "date") ], by = "date", all.x = T )
      all.sp = merge( all.sp, dfll.sp[ c("dfllc", "dfllesd", "date") ], by = "date", all.x = T )
      all.sp = merge( all.sp, nut.sp[ c("n", "p", "l", "lm", "temp", "si", "sal", "spm", "date") ], by = "date", all.x = T )
      all.sp = merge( all.sp, zoos.sp[ c("date", "s.dw", "s.esd") ], by = "date", all.x = T )
      all.sp = merge( all.sp, zooc.sp[ c("date", "c.dw", "c.esd") ], by = "date", all.x = T )
      all.sp = merge( all.sp, zool.sp[ c("date", "l.dw", "l.esd") ], by = "date", all.x = T )
      #all.sp = merge( all.sp, zoohl.sp[ c("date", "hl.dw", "hl.esd") ], by = "date", all.x = T )
      sp.sum = rbind( sp.sum, all.sp )
    }
    # Spring ----
    
    # Autumn ----
    if( length(date.at) > 0 & !is.infinite(date.at)[1] ){
      all.at = phys.at ; all.at$yr = yri
      all.at = merge( all.at, phyl.at[ c("biocl", "sizel", "date") ], by = "date", all.x = T )
      all.at = merge( all.at, dfls.at[ c("dflsc", "dflsesd", "date") ], by = "date", all.x = T )
      all.at = merge( all.at, dfll.at[ c("dfllc", "dfllesd", "date") ], by = "date", all.x = T )
      all.at = merge( all.at, nut.at[ c("n", "p", "l", "lm", "temp", "si", "sal", "spm", "date") ], by = "date", all.x = T )
      all.at = merge( all.at, zoos.at[ c("date", "s.dw", "s.esd") ], by = "date", all.x = T )
      all.at = merge( all.at, zooc.at[ c("date", "c.dw", "c.esd") ], by = "date", all.x = T )
      all.at = merge( all.at, zool.at[ c("date", "l.dw", "l.esd") ], by = "date", all.x = T )
      #all.at = merge( all.at, zoohl.at[ c("date", "hl.dw", "hl.esd") ], by = "date", all.x = T )
      at.sum = rbind( at.sum, all.at )
    }
    # Autumn ----
    
    # After spring ----
    if( length(date.sw) > 0 & all(!is.infinite(date.sw)) ){
      all.sw = phys.sw ; all.sw$yr = yri
      all.sw = merge( all.sw, phyl.sw[ c("biocl", "sizel", "date") ], by = "date", all.x = T )
      all.sw = merge( all.sw, dfls.sw[ c("dflsc", "dflsesd", "date") ], by = "date", all.x = T )
      all.sw = merge( all.sw, dfll.sw[ c("dfllc", "dfllesd", "date") ], by = "date", all.x = T )
      all.sw = merge( all.sw, nut.sw[ c("n", "p", "l", "lm", "temp", "si", "sal", "spm", "date") ], by = "date", all.x = T )
      all.sw = merge( all.sw, zoos.sw[ c("date", "s.dw", "s.esd") ], by = "date", all.x = T )
      all.sw = merge( all.sw, zooc.sw[ c("date", "c.dw", "c.esd") ], by = "date", all.x = T )
      all.sw = merge( all.sw, zool.sw[ c("date", "l.dw", "l.esd") ], by = "date", all.x = T )
      #all.sw = merge( all.sw, zoohl.sw[ c("date", "hl.dw", "hl.esd") ], by = "date", all.x = T )
      #all.sw$dur = NA ; all.sw$dur = dur.sw
      sw.sum = rbind( sw.sum, all.sw )
    }
    # After spring ----
  }
  
  # Removing winter months
  ind = which(sp.sum$mth<3 | sp.sum$mth>10)
  if(length(ind)>0){ 
  sp.sum = sp.sum[-ind,]}
  
  ind = which(at.sum$mth<3 | at.sum$mth>10)
  if(length(ind)>0){ 
    at.sum = at.sum[-ind,]}
  
  ind = which(sw.sum$mth<3 | sw.sum$mth>10)
  if(length(ind)>0){ 
    sw.sum = sw.sum[-ind,]}
  
  # Data imputation using 50% climatology and 50% regression
  sp.sum = fill.nan.days(sp.sum)
  at.sum = fill.nan.days(at.sum)
  sw.sum = fill.nan.days(sw.sum)
  
  return( list( sp.sum, at.sum, sw.sum ) ) }

sb = season.extract.days(season.f, date.extend=0) # Sampling the data for spring, after spring, and autumn
sp.sum = sb[[1]] ; at.sum = sb[[2]] ; sw.sum = sb[[3]]
# Spring           Autumn             After spring

# Add a nutrient saturation metric, as nitrogen and phosphorus are highly correlated
saturation.nutrients = function(esd, n, p){
  vol = (exp(esd)/2)**3 * 4/3*pi
  
  alpha.n = -0.84
  beta.n = 0.33

  alpha.p = -1.4
  beta.p = 0.41
  
  Kn = 10**alpha.n * vol**beta.n # Half saturation, in mu mol.L-1
  Kp = 10**alpha.p * vol**beta.p 
  
  saturation.n = n / (Kn + n)
  saturation.p = p / (Kp + p)
  
  saturation = data.frame( n=saturation.n, p=saturation.p )

  saturation = 1 / (1/saturation$n + 1/saturation$p)
  
  ind = is.infinite(saturation)
  if(length(ind)>0){
    saturation[ind] = NA
  }
  
  return(saturation)
}

sp.sum$sat.phys = saturation.nutrients(sp.sum$sizes, sp.sum$n, sp.sum$p)
sp.sum$sat.phyl = saturation.nutrients(sp.sum$sizel, sp.sum$n, sp.sum$p)
sp.sum$sat.dfls = saturation.nutrients(sp.sum$dflsesd, sp.sum$n, sp.sum$p)

at.sum$sat.phys = saturation.nutrients(at.sum$sizes, at.sum$n, at.sum$p)
at.sum$sat.phyl = saturation.nutrients(at.sum$sizel, at.sum$n, at.sum$p)
at.sum$sat.dfls = saturation.nutrients(at.sum$dflsesd, at.sum$n, at.sum$p)

sw.sum$sat.phys = saturation.nutrients(sw.sum$sizes, sw.sum$n, sw.sum$p)
sw.sum$sat.phyl = saturation.nutrients(sw.sum$sizel, sw.sum$n, sw.sum$p)
sw.sum$sat.dfls = saturation.nutrients(sw.sum$dflsesd, sw.sum$n, sw.sum$p)

# Testing the light saturation, following the model in Edwards et al (2015)
saturation.light = function(esd, l){
  vol = (exp(esd)/2)**3 * 4/3*pi
  
  # SMA model
  alpha = -0.62
  beta = -0.37
  
  # LM model
  #alpha = -1.36
  #beta = -0.13
  
  al = 10**alpha * vol**beta # Saturation, in mu quanta-1.m2.d-1
  lq = l * 24 / (6.626e-34 * 3e8 / 600e-9 * 6.022e23) * 1e6 # Conversion to mu mol quanta.d
  # Planck constant * light speed / wavelength (~ 600 mu m for visible ) * Avogadro's number
  
  max_growth = 10**0.7 * vol**(-0.24) # Max growth rate of phytoplankton
  
  saturation = 1 - exp(- al * lq / max_growth) # Poisson law of photon capture

  ind = is.infinite(saturation)
  if(length(ind)>0){
    saturation[ind] = NA
  }
  
  return(saturation)
}

sp.sum$satl.phys = saturation.light(sp.sum$sizes, sp.sum$lm)
sp.sum$satl.phyl = saturation.light(sp.sum$sizel, sp.sum$lm)
sp.sum$satl.dfls = saturation.light(sp.sum$dflsesd, sp.sum$lm)

at.sum$satl.phys = saturation.light(at.sum$sizes, at.sum$lm)
at.sum$satl.phyl = saturation.light(at.sum$sizel, at.sum$lm)
at.sum$satl.dfls = saturation.light(at.sum$dflsesd, at.sum$lm)

sw.sum$satl.phys = saturation.light(sw.sum$sizes, sw.sum$lm)
sw.sum$satl.phyl = saturation.light(sw.sum$sizel, sw.sum$lm)
sw.sum$satl.dfls = saturation.light(sw.sum$dflsesd, sw.sum$lm)

## 
summary(sp.sum[c('satl.phys', 'satl.phyl', 'satl.dfls')]) # Mean of light saturation
summary(at.sum[c('satl.phys', 'satl.phyl', 'satl.dfls')]) # Mean of light saturation
summary(sw.sum[c('satl.phys', 'satl.phyl', 'satl.dfls')]) # Mean of light saturation

mov.avg = function(datel, datal, wind, forecast){
  datel = ymd(datel)
  date0 = min(datel) + days(wind)
  datef = max(datel)
  
  date.avg = c()
  data.avg = c()
  
  move = wind
  if(forecast){
    move = 1
  }
  
  # Test
  date.to = datel[ which(datel < datef-days(wind) & datel > date0 ) ]
  
  #for(i in 1:length(date.to)){
  #  date0 = date.to[i]
  while( date0 < datef ){
    up.wind = date0+days(wind)
    low.wind = date0-days(wind)
    ind = which(datel < up.wind & datel >= low.wind)
    
    if( length(ind)>1 ){
      # Normal MV
      meani = mean(datal[ind], na.rm=T)
      datei = mean(datel[ind], na.rm=T)
      
    }else if( length(ind)==1 ){
      meani = datal[ind]
      datei = datel[ind]
    }
    
    if( length(ind)>0 ){
      data.avg = c(data.avg, meani)
      date.avg = c(date.avg, datei)
    }
    
    date0 = date0 + 1*days(move)
  }
  
  df = data.frame(date = date.avg, dat = data.avg)

  return(df)
}

seas.avg = function(frame, wind = 7, forecast=F){
  datef = frame$date
  
  ind.no = which( names(frame) == 'date' )
  frame = frame[,-ind.no]
  
  varl = names(frame)
  nvar = length(varl)
  
  avg.f = mov.avg(datef, frame[[varl[1]]], wind, forecast) 
  for(vi in 2:nvar){
    avg.f = cbind(avg.f, mov.avg(datef, frame[[varl[vi]]], wind, forecast)$dat )
  }
  names(avg.f)[-1] = varl
  
  avg.f$date = as_date(avg.f$date)
  
  # Remove duplicated dates, if any
  ind = which( duplicated( ymd(avg.f$date) ) )
  if( length(ind)>0 ){
    avg.f = avg.f[-ind,]
  }
  
  return(avg.f)
}

# Compute the rates and the means, with a moving average (or none)
sp.avg = seas.avg(sp.sum, wind=30, forecast=F)
at.avg = seas.avg(at.sum, wind=30, forecast=F)
sw.avg = seas.avg(sw.sum, wind=30, forecast=F)

sp.avg[ which(sp.avg == 0, arr.ind = T) ] = NA
at.avg[ which(at.avg == 0, arr.ind = T) ] = NA
sw.avg[ which(sw.avg == 0, arr.ind = T) ] = NA

## Test the importance of light
# Correlations with light
cor.test(sp.avg$lm, sp.avg$n, use="complete.obs", method="spearman")
cor.test(sp.avg$lm, sp.avg$p, use="complete.obs", method="spearman")

rate.frames = function(frame){ # Function to compute the cross-correlation matrix
  nvar = length( names(frame) )
  nyear = length( unique(frame$yr) )
  
  rate.all = data.frame()
  mean.all = data.frame()
  
  for( yri in unique(frame$yr)){
    fi = frame[ which(frame$yr == yri), ]
    
    date.f = ymd(fi$date)
    diff.days = as.vector( diff(date.f) ) #day(days(date.f))[-1] - day(days(date.f))[-length(date.f)]
    
    ratei = matrix( data = NA, nrow = nrow(fi)-1, ncol = nvar-1 )
    meani = matrix( data = NA, nrow = nrow(fi)-1, ncol = nvar-1 )
    datei = date.f[-nrow(fi)] + days(round(diff.days/2)) # Median date
    
    for( vi in 1:(nvar-1) ){
      ratei[,vi] = diff( fi[,vi+1] ) / diff.days
      
      meani[,vi] = (fi[-nrow(fi),vi+1] + fi[-1,vi+1])/2
      
      #ratei[,vi] = diff( exp(fi[,vi+1]) ) / exp(meani[,vi]) / diff.days
    }
    
    # Storing the data
    ratei = as.data.frame(ratei) ; names(ratei) = names(fi)[-1]
    meani = as.data.frame(meani) ; names(meani) = names(fi)[-1]
    
    ratei$date = datei ; meani$date = datei
    rate.all = rbind(rate.all, ratei)
    mean.all = rbind(mean.all, meani)
  }
  
  return( list(rate.all, mean.all) )}

sp.rate = rate.frames(sp.avg)
sw.rate = rate.frames(sw.avg)
at.rate = rate.frames(at.avg)

get.cor = function(v1, v2){
  # Get correlations and p.values
  p.lab = c( "p>0.1", "p<0.1", "p<0.05", "p<0.01" )
  p.test = c( 1, 0.1, 0.05, 0.01 )
  
  cor.st = cor.test( v1, v2, method = "spearman", use = "complete.obs" )
  p.st = p.lab[ max( which(cor.st$p.value <= p.test ) ) ]
  #cor.txt = paste( "r2:", round(cor.st$estimate**2, 2) , ", ", p.st, sep = "" )
  cor.txt = paste( "r2:", round(cor.st$estimate**2, 2) , "\n", p.st, sep = "" )
  return( list(cor.txt, cor.st$p.value, cor.st$estimate) )
}

subplot.cor = function(v2, v1, xlab, ylab, line.cor, control){
  cex.txt = 1.3
  
  # Colors
  col.bu = rgb( 0.3, 0.3, 0.8, alpha = 0.5)
  col.td = rgb( 0.9, 0.3, 0.3, alpha = 0.5)
  col.unc = rgb(0.3, 0.3, 0.3, 0.5)
  col.not = rgb(0.6, 0.6, 0.6, 0.4)
  
  cor.txt = get.cor(v1, v2)
  if( cor.txt[[2]] < 0.1 & control == "BU" & cor.txt[[3]] > 0 ){
    colo = col.bu
  }else if( cor.txt[[2]] < 0.1 & control == "TD" & cor.txt[[3]] < 0 ){
    colo = col.td
  }else if( cor.txt[[2]] < 0.1 ){
    colo = col.unc
  }else{
    colo = col.not
  }
  
  # Remove extremes (limit at 3.SD)
  v1.rge = c( mean(v1, na.rm=T) - 3*sd(v1, na.rm=T), mean(v1, na.rm=T) + 3*sd(v1, na.rm=T) )
  ind = which( (v1 >= v1.rge[1] & v1 <= v1.rge[2]) & v1 != 0 )
  v1=v1[ind]
  v2=v2[ind]
  
  x.lim=range(v2, na.rm=T)
  # y.lim=range(v1, na.rm=T)
  #x.lim=quantile(v2, probs=c(0.01, 0.99), na.rm=T)
  y.lim=quantile(v1, probs=c(0.01, 0.99), na.rm=T)
  
  matplot(c(0,1), c(0,1), type = "n", xlab="", ylab="", xaxt="n", yaxt="n",
          xlim=x.lim, ylim=y.lim)
  matplot(v2, v1, pch = 19, col = colo, cex = 2.5, add = T)
  abline( lm(v1~v2), col=adjustcolor(colo, alpha.f=2), lwd=3 )
  
  # Show the points outside the plotting range
  y.rge = diff(y.lim)*0.1 # Test if points are outside the 10% plot range
  y.abs.plot = diff(y.lim)*0.03
  
  par(xpd=NA)
  points(v2[v1 > y.lim[2]+y.rge], rep(y.lim[2]+y.abs.plot, sum(v1 > y.lim[2]+y.rge, na.rm=T)),
         pch = 24, bg="black", cex=1.5)
  points(v2[v1 < y.lim[1]-y.rge], rep(y.lim[1]-y.abs.plot, sum(v1 < y.lim[1]-y.rge, na.rm=T)),
         pch = 25, bg="black", cex=1.5)
  par(xpd=F)
  
  col.axis = rgb(0.3, 0.3, 0.3, 0.5)
  axis(1, at=pretty(x.lim), cex.axis=cex.txt*1.3, tck=0.02, col.axis=col.axis)
  axis(2, at=pretty(y.lim), labels=pretty(y.lim)*100, cex.axis=cex.txt*1.3, las = 1, tck=0.02, col.axis=col.axis)
  
  mtext(side = 1, xlab, line = 2.7, cex=cex.txt*0.92 )
  mtext(side = 2, ylab, line = 2.6, cex=cex.txt*0.95, adj=0.8 )
  
  col.cor = col.axis
  if( cor.txt[[2]]<0.1 ){
    col.cor = adjustcolor(col.cor, alpha.f = 2)
  } 
  
  if( cor.txt[[3]] > 0 ){
    mtext(side = 3, cor.txt[[1]], line = line.cor, cex=cex.txt-0.2, adj=0.05, font=2, col=col.cor )
  }else{
    mtext(side = 3, cor.txt[[1]], line = line.cor, cex=cex.txt-0.2, adj=0.95, font=2, col=col.cor )
  }
}

plot.rates = function(frame, line.cor=-5, season=NA){
  meanf = frame[[2]]
  ratef = frame[[1]]
  
  x11(height = 20, width = 15)
  par(mfrow=c(4,4), mar=c(4, 5, 1, 0.1), mgp=c(3,0.6,0))
  
  ## Small phy
  subplot.cor(meanf$sat.phys, ratef$biocs, expression("Nutrient co-limitation"),
              expression("Small phyto RCR, 10"^-2 *".d"^-1), line.cor, "BU")
  subplot.cor(meanf$s.dw, ratef$biocs, expression("Small zoo., mg.m"^-3),
              '', line.cor, "TD")
  subplot.cor(meanf$l.dw, ratef$biocs, expression("Large zoo., mg.m"^-3),
              '', line.cor, "TD")
  matplot(1,1, type="n", axes=F, xlab='', ylab='')
  
  ## Large phy
  subplot.cor(meanf$sat.phyl, ratef$biocl, expression("Nutrient co-limitation"),
              expression("Large phyto RCR, 10"^-2 *".d"^-1), line.cor, "BU")
  subplot.cor(meanf$l.dw, ratef$biocl, expression("Large zoo., mg.m"^-3), 
              '', line.cor, "TD")
  matplot(1,1, type="n", axes=F, xlab='', ylab='')
  matplot(1,1, type="n", axes=F, xlab='', ylab='')

  ## Small zoo
  subplot.cor(meanf$biocs, ratef$s.dw, expression("Small phyto, " *mu *"molC.L"^-1), 
              expression("Small zoo RCR, 10"^-2 *".d"^-1), line.cor, "BU")
  subplot.cor(meanf$l.dw, ratef$s.dw, expression("Large zoo., mg.m"^-3), 
              '', line.cor, "TD")
  subplot.cor(meanf$c.dw, ratef$s.dw, expression("Carn zoo., mg.m"^-3), 
              '', line.cor, "TD")
  matplot(1,1, type="n", axes=F, xlab='', ylab='')

  ## Large zoo
  subplot.cor(meanf$biocs, ratef$l.dw, expression("Small phyto, " *mu *"molC.L"^-1), 
              expression("Large zoo RCR, 10"^-2 *".d"^-1), line.cor, "BU")
  subplot.cor(meanf$biocl, ratef$l.dw, expression("Large phyto, " *mu *"molC.L"^-1), 
              '', line.cor, "BU")
  subplot.cor(meanf$s.dw, ratef$l.dw, expression("Small zoo., mg.m"^-3), 
              '', line.cor, "BU")
  subplot.cor(meanf$c.dw, ratef$l.dw, expression("Carn zoo., mg.m"^-3), 
              '', line.cor, "TD")

  ## Season label
  par(fig=c(0,1,0,1), new=T)
  matplot(1,1, type="n", axes=F, xlab='', ylab='', xlim=c(0,1), ylim=c(0,1))
  mtext(side = 3, season, line = -33, cex=1.4, adj=0.83, font=2 )
  
  col.bu = rgb( 0.3, 0.3, 0.8, alpha = 0.5)
  col.td = rgb( 0.9, 0.3, 0.3, alpha = 0.5)
  col.unc = rgb(0.3, 0.3, 0.3, 0.5)
  col.not = rgb(0.6, 0.6, 0.6, 0.3)
  
  legend(x=0.8, y=0.55, horiz=F, bty="n", 
         legend=c('Top-down regulation', 'Bottom-up regulation', 'Unclear regulation', 'Not significant'),
         pch=19, col=c(col.td, col.bu, col.unc, col.not), cex=2., pt.cex=4, 
         x.intersp=0.7, y.intersp=1.5)
}

plot.rates(sp.rate, line.cor=-3.5, season="Spring")
plot.rates(at.rate, line.cor=-3.5, season="Autumn")
plot.rates(sw.rate, line.cor=-3.5, season="Summer")

## Show the complete correlation table (in SI)
corr.rate = function(frame){
  meanf = frame[[2]]
  ratef = frame[[1]]
  
  var.tocheck = c("biocs", "biocl", "dflsc", "dfllc", 'n', 'p', 'lm', 'temp', 'si',
                  's.dw', 'l.dw', 'c.dw')
  meanf = meanf[var.tocheck]
  ratef = ratef[var.tocheck]
  
  nvar = ncol(ratef) 
  
  p.frame = matrix( data = NA, nrow = nvar, ncol = nvar )
  colnames(p.frame) = names(meanf) ; rownames(p.frame) = names(meanf)

  c.frame = matrix( data = NA, nrow = nvar, ncol = nvar )
  colnames(c.frame) = names(meanf) ; rownames(c.frame) = names(meanf)

  for( i in ( 1:nvar ) ){
    yi = ratef[,i]
    for( j in ( 1:nvar ) ){  # Computing Spearman's correlation by pairs
      yj = meanf[,j]

      ind = which( (!is.na(yi) & !is.na(yj)) )
      if( length(ind)>4 ){ # At least four points
        cr = cor.test( yi, yj, method="spearman", use = "complete.obs" )

        p.frame[i, j] = cr$p.value
        c.frame[i, j] = round( cr$estimate, 2 ) }
    }
  }

  return( list(c.frame, p.frame) )
}

plot.table = function( corrdt ){
  corr = round( as.data.frame( corrdt[[1]] ), 2) # Extracting the R2 and p-values
  p = as.data.frame( corrdt[[2]] )
  
  ## Selecting the variables to extract (season dataframes)
  cn = c( "lm", "n", "p", "si", "biocs", "biocl",
          "dflsc", "s.dw", "l.dw", "c.dw" )
  rown = c( "biocs", "biocl", "dflsc", "s.dw", "l.dw", "c.dw" ) # Rows to look at

  corr = corr[ rown, cn ]
  p = p[ rown, cn ]
  
  # Rename the variables
  names(corr) = c( "Light", "NO3", "PO4", "SiO2", "S.phy", "L.phy",
                   "S.dfl", "S.zoo", "L.zoo", "C.zoo" )
  rown = c( "S.phy", "L.phy", "S.dfl", "S.zoo", "L.zoo", 'C.zoo' )
  rownames(corr) = rown
  
  ## Print the significance level on the table
  
  # p < 0.1
  ind = which(p < 0.1 & p > 0.05, arr.ind = T)
  corr[ind] = paste( corr[ind], "*", sep = "" )
  
  # p < 0.05
  ind = which(p < 0.05 & p > 0.01, arr.ind = T)
  corr[ind] = paste( corr[ind], "**", sep = "" )
  
  # p < 0.01
  ind = which(p < 0.01, arr.ind = T)
  corr[ind] = paste( corr[ind], "***", sep = "" )
  
  ## Remove correlations that do not make sense
  corr[5:6, 1:4] = '-' # Nutrients and zooplankton
  corr[rown=="S.phy", names(corr) %in% c("S.phy", "L.phy", "C.zoo")] = '-'
  corr[rown=="L.phy", names(corr) %in% c("L.phy", "S.phy", "S.dfl", "S.zoo", "C.zoo")] = '-'
  corr[rown=="S.dfl", names(corr) %in% c("S.dfl", "L.phy", "C.zoo")] = '-'
  corr[rown=="S.zoo", names(corr) %in% c("S.zoo", "L.phy")] = '-'
  corr[rown=="L.zoo", names(corr) %in% c("L.zoo")] = '-'
  corr[rown=="C.zoo", names(corr) %in% c("L.phy", "S.phy", "S.dfl", "C.zoo")] = '-'
  
  ## Print the table
  x11(width = 8, height = 8)
  grid.table(corr)
  
  return(corr)
}

sp.rate.corr = corr.rate(sp.rate)
corr.sp = plot.table(sp.rate.corr) # Spring

at.rate.corr = corr.rate(at.rate)
corr.at = plot.table(at.rate.corr) # Autumn

sw.rate.corr = corr.rate(sw.rate)
corr.sw = plot.table(sw.rate.corr) # Summer

# Correlations of concentrations
corr.sp.bio = corr.rate( list(sp.rate[[2]], sp.rate[[2]]) )
corr.sw.bio = corr.rate( list(sw.rate[[2]], sw.rate[[2]]) )

## Print the dominant species in zooplankton
dom.spec = function(date.seas, usi.dt, date.dt, bio.dt, spec.dt, taxa.dt){
  ind = which( date.dt %in% date.seas )
  
  df = data.frame( usi = usi.dt[ind], date = date.dt[ind], spec = spec.dt[ind], taxa=taxa.dt[ind], bio = bio.dt[ind] )
  
  df.day = ddply( df, .(usi), summarize, biotot = sum( bio, na.rm = T ) )
  df.day = merge(df, df.day, by='usi')
  
  df.day$w = df.day$bio / df.day$biotot # Weight per sample
  tot.obs = length(unique(df$usi))      # Number of sample
  df.day$w = df.day$w / tot.obs         # Biomass fraction weighted by the frequency of observation
  
  spec.fraction = ddply( df.day, .(spec, taxa), summarize, 
                         w = sum(w, na.rm = T) )
  spec.fraction = spec.fraction[ order(spec.fraction$w), ]
  
  return(spec.fraction)
}

dom.seas = function(usi.dt, date.dt, bio.dt, spec.dt, taxa.dt){
  sp.dom = dom.spec(sp.sum$date, usi.dt, date.dt, bio.dt, spec.dt, taxa.dt)
  print('Spring')
  print( sp.dom[ (dim(sp.dom)[1]-3):dim(sp.dom)[1], ] )
  
  sw.dom = dom.spec(sw.sum$date, usi.dt, date.dt, bio.dt, spec.dt, taxa.dt)
  print('Summer')
  print( sw.dom[ (dim(sw.dom)[1]-3):dim(sw.dom)[1], ] )
  
  at.dom = dom.spec(at.sum$date, usi.dt, date.dt, bio.dt, spec.dt, taxa.dt)
  print('Autumn')
  print( at.dom[ (dim(at.dom)[1]-3):dim(at.dom)[1], ] )
}

dom.seas(zooc$USI, zooc$Date, zooc$DryW, zooc$Phylum, zooc$Taxa)
dom.seas(zoos$USI, zoos$Date, zoos$DryW, zoos$Phylum, zoos$Taxa)
dom.seas(zool$USI, zool$Date, zool$DryW, zool$Phylum, zool$Taxa)

ind = which( zooc$Date %in% sw.sum$date )
hist( log(zooc$mean.esd[ind]) )

# Small and large phytoplankton
dom.seas(dfls$USI, dfls$date, dfls$bioC, dfls$species, dfls$phylum)

phys.spec = dom.seas(phys$USI, phys$date, phys$bioC, phys$species, phys$phylum)
phyl.spec = dom.seas(phyl$USI, phyl$date, phyl$bioC, phyl$species, phyl$phylum)

################################################
## Fourier section ( interannual variations ) ## ---------------------------------------------
################################################

# Function to plot smooth spectra
specapprox = function(points, data, span=0.05, cut=T, max.x = max(points)){ # Data X and Y
  curve = c()
  xaxis = c()
  
  data = data[ which( points <= max.x ) ] # If it is to cut the series
  points = points[ which( points <= max.x ) ]
  
  if(cut){
    # Suppressing redundant points for plotting 
    # - Knots are the pits
    # - Only peaks and knots are kept
    ndetec = data[-1] - data[-length(data)] # Variation
    nvar = ndetec[-1] * ndetec[-length(ndetec)] # Variation
    ind.peak = which(nvar < 0) +1               # Variation change, peak or pit
    datai = c(data[1],  data[ind.peak] )
    pointsi = c( points[1], points[ind.peak]) # Get first point as well
  }else{
    pointsi = points
    datai = data
  }
  predni = seq(min(pointsi), max(pointsi), length.out=max(pointsi) * 200)
  model = loess(datai~pointsi, span=span)
  spec.smooth = predict(model, newdata = predni)

  curve = spec.smooth
  xaxis = predni
  return(list(xaxis, curve))
}

## Looking for a shift in the data
# Classic method (Student) : computing mean and var on two windows, are they different ?
# Student stat is 1.98
test.student = function(month.dat, year.dat, d.dat){
  
  ini = 12 #(6 months min)
  imax = length(month.dat)-ini
  varcx = c()
  mncx = c()
  
  varcy = c()
  mncy = c()
  
  for(icount in ini:imax){
    mncx = c(mncx, mean(d.dat[1:icount], na.rm=T))
    varcx = c(varcx, var(d.dat[1:icount], na.rm=T))
    
    mncy = c(mncy, mean(d.dat[(icount+1):(imax+ini)], na.rm=T))
    varcy = c(varcy, var(d.dat[(icount+1):(imax+ini)], na.rm=T))
  }
  
  D = mncx - mncy
  
  ncx = ini:imax
  ncy = rev(ncx)
  
  sd = sqrt(varcx / ncx + varcy / ncy) 
  
  St = c(rep(NaN, ini), D / sd, rep(NaN, ini))
  
  inds = which.max( abs( round(St, 2) ) )
  
  x11()
  par(mfrow=c(2,1))
  xticks=which(month.dat == 1)
  
  plot(d.dat, type="n", lwd=2, col="black", 
       xlab="Month", ylab="Variable", xaxt="n", cex.axis=1.5, cex.lab=1.6,
       main = paste("Student shift test"), cex.main=2)
  lines(d.dat, type="l", lwd=2, col="black")
  axis(1, at=xticks, labels = year.dat[xticks], cex.axis=1.5)
  par(new=T)
  plot(d.dat, type="l", lwd=3, xlab="", ylab="", axes=F)
  
  abline(v=inds, col="black", lty=2, lwd=1)
  
  plot(St, type="n", lwd=2, col="black", xlab="Months", ylab="Test statistic", xaxt="n", cex.axis=1.5, cex.lab=1.6)
  lines(St, type="l", lwd=2, col="black")
  axis(1, at=xticks, labels = year.dat[xticks], cex.axis=1.5)
  abline(1.98, 0, col="black", lty=2, lwd=2)
  abline(-1.98, 0, col="black", lty=2, lwd=2)
  
  abline(v=inds, col="black", lty=2, lwd=1)
  return(inds)}

# Function to normalize according to the shift
norm.stud = function( inds, d.dat ){
  # Spectrum after removing shift
  indlast = length(d.dat)
  
  mnx = mean(d.dat[1:inds], na.rm=T)
  varx = var(d.dat[1:inds], na.rm=T)
  mny = mean(d.dat[inds:indlast], na.rm=T)
  vary = var(d.dat[inds:indlast], na.rm=T)
  
  # Normalization according shift
  d.dat[1:inds] = (d.dat[1:inds] - mnx) / varx**0.5
  d.dat[(inds+1):indlast] = (d.dat[(inds+1):indlast] - mny) / vary**0.5
  return(d.dat)}

## List for Western stations only
Sizelist = list() # For storing the spectra and the data of each phylum
Bioclist = list()
Biovlist = list()

groupl = c("Small phyto", "Large phyto", "Dinophyta")

# Starting computation
for(s in 1:3){
  datam = list(physm, phylm, dflsm)[[s]]
  
  # Filter Moving average
  bioroll = exp(filter(log(datam$bioc), filter=filter.wind, sides = 2, method="convolution"))
  roller = exp(filter(log(datam$biovol), filter=filter.wind, sides = 2, method="convolution"))
  
  # Size as ESD (like in the rest of the analysis)
  datam$esd = (datam$vol *3/4/pi)**(1/3) *2
  sizeroll = exp(filter(log(datam$esd), filter=filter.wind, sides = 2, method="convolution"))
  
  plot( 1:length(sizeroll), sizeroll, type="l", log = "y" )
  
  # Storing data
  Sizelist = c( Sizelist, list(sizeroll) )
  Bioclist = c( Bioclist, list(bioroll) )
  Biovlist = c( Biovlist, list(roller) )
  
  # Size
  sizefft = spec.pgram(fast=T, log(sizeroll[- which( is.na(sizeroll) ) ]), demean=T, plot=F)
  Sizelist = c( Sizelist, list(sizefft) )
  
  #BioC
  biocfft = spec.pgram(fast=T, log(bioroll[- which( is.na(bioroll) ) ]), demean=T, plot=F)
  Bioclist = c( Bioclist, list(biocfft) )
  
  # Total Biovolume
  biovfft = spec.pgram(fast=T, log(roller[- which( is.na(roller) ) ]), demean=T, plot=F)
  Biovlist = c( Bioclist, list(biovfft) )
  
  x11()
  par(fig = c(0, 0.6, 0.45, 1), mar = c(4, 4, 4, 4) )
  xyr = which(datam$mth == 1) # Labels on years
  qt = quantile(datam$vol, probs=c(0.05, 0.95))
  plot(datam$esd, col="darkgray", xaxt="n", xlab="", ylab="", type="l", lwd=0.1,
       cex.axis=1.6, log="y")#, ylim=c(qt[1], qt[2]))
  lines(sizeroll, col="darkgreen", lwd=5)
  x_axis = axis(1, at=xyr, labels=F, tick = T)
  abline(v = x_axis, col="lightgray", lwd=2, lty=3)
  legend("bottomright", c("Monthly Data", "Year Mov. Av. filter"),
         col=c("darkgray", "lightgreen"), lty=c(1,1), inset=c(0,0.9), xpd=TRUE, horiz=TRUE, 
         bty="n", cex=1.8, seg.len = 0.8, lwd=c(2,8), text.width=c(20,30), x.intersp = 0.5)
  title(ylab=paste("ESD [µm]"), line=2.5, cex.lab=1.8)
  
  par(fig = c(0, 0.6, 0, 0.55), new=T )
  qt = quantile(datam$bioc, probs=c(0.05, 0.95))
  plot(datam$bioc, col="darkgray", xaxt="n", xlab="", ylab="", type="l", lwd=0.1,
       cex.axis=1.5, log="y")#, ylim=c(qt[1], qt[2]))
  axis(1, at=xyr, labels=datam$yr[xyr], cex.axis=1.5)
  lines(bioroll, col="green", lwd=5)
  title(xlab="Time", line=2.5, cex.lab=1.8)
  title(ylab=paste("Biomass [pgC]"), line=2.5, cex.lab=1.8)
  abline(v = x_axis, col="lightgray", lwd=2, lty=3)
  
  par(fig = c(0.6, 1, 0.45, 1), new=T )
  szmod = Sizelist[[s*2]]
  axises = specapprox(1/szmod$freq, szmod$spec)
  plot(axises[[1]], axises[[2]], ylim=c(0,max(axises[[2]])), xlim=c(8, floor(length(datam$mth)/3)), log="x",
       type="l", xlab="", ylab="", 
       xaxt='n', cex.axis=1.6, col="darkgreen", lwd=4)
  x_axis = axis(1, at=seq(0, floor(1/min(szmod$freq) ), 6), labels=F, tick=T)
  abline(v = x_axis, col="lightgray", lwd=2, lty=3)
  title(ylab="Spectrum intensity", line=2.5, cex.lab=1.8)
  
  par(fig = c(0.6, 1, 0, 0.55), new=T )
  szmod = Bioclist[[s*2]]
  axises = specapprox(1/szmod$freq, szmod$spec)
  plot(axises[[1]], axises[[2]], ylim=c(0,max(axises[[2]])), xlim=c(8, floor(length(datam$mth)/3)), log="x",
       type="l", xlab="", ylab="", 
       xaxt='n', cex.axis=1.6, col="green", lwd=4)
  axis(1, at=seq(0, floor(1/min(szmod$freq) ), 12), cex.axis=1.6)
  axis(1, at=seq(0, floor(1/min(szmod$freq) ), 6), label=F)
  abline(v = x_axis, col="lightgray", lwd=2, lty=3)
  title(ylab="Spectrum intensity", line=2.5, cex.lab=1.8)
  title(xlab="Frequency (Months)", line=2.5, cex.lab=1.8)
  
  ## Looking for a shift in the data (vol and bioC)
  # Classic method (Student) : computing mean and var on two windows, are they different ?
  # Student stat is 1.98
  datam$vol = log(datam$vol)
  datam$bioc = log(datam$bioc)
  datam$biovol = log(datam$biovol)
  
  indb = test.student(datam$mth, datam$yr, datam$bioc)
  inds = test.student(datam$mth, datam$yr, datam$vol)
  
  dev.off()
  
  if(groupl[s] == "Dinophyta"){
    # Spectrum after removing shift
    indmean = inds #mean(inds, indb)
    indlast = length(datam$mth)
    
    # Normalization according to regions (pre-post shift)
    datam.norm = datam
    datam.norm$vol = norm.stud(indmean, datam$vol)
    datam.norm$bioc = norm.stud(indmean, datam$bioc)
    datam.norm$biovol = norm.stud(indmean, datam$biovol)
    
    # Filter Moving average
    bioroll = filter( datam.norm$bioc, filter=filter.wind, sides = 2, method="convolution")
    sizeroll = filter( datam.norm$vol, filter=filter.wind, sides = 2, method="convolution")
    roller = filter( datam.norm$biovol, filter=filter.wind, sides = 2, method="convolution")  
    
    # Storing data
    Sizelist = c( Sizelist, list(sizeroll) )
    Bioclist = c( Bioclist, list(bioroll) )
    Biovlist = c( Biovlist, list(roller) )
    
    # Size
    sizefft = spec.pgram(fast=T, sizeroll[- which( is.na(sizeroll)  )], demean=T, plot=F)
    Sizelist = c( Sizelist, list(sizefft) )
    
    #BioC
    biocfft = spec.pgram(fast=T, bioroll[- which( is.na(bioroll) ) ], demean=T, plot=F)
    Bioclist = c( Bioclist, list(biocfft) )
    
    # Total Biovolume
    biovfft = spec.pgram(fast=T, roller[- which( is.na(roller) ) ], demean=T, plot=F)
    Biovlist = c( Bioclist, list(biovfft) )
    
    x11()
    par(fig = c(0, 0.6, 0.45, 1) )
    xyr = which(datam$mth == 1) # Labels on years
    qt = quantile(datam.norm$vol, probs=c(0.05, 0.95))
    plot(datam.norm$vol, col="darkgray", xaxt="n", xlab="", ylab="", type="l", lwd=0.1,
         cex.axis=1.6, ylim=c(qt[1], qt[2]))
    lines(sizeroll, col="darkgreen", lwd=5)
    x_axis = axis(1, at=xyr, labels=F, tick = T)
    abline(v = x_axis, col="lightgray", lwd=2, lty=3)
    legend("bottomright", c("Monthly Data", "Year Mov. Av. filter"),
           col=c("darkgray", "lightgreen"), lty=c(1,1), inset=c(0,0.9), xpd=TRUE, horiz=TRUE, 
           bty="n", cex=1.8, seg.len = 0.8, lwd=c(2,8), text.width=c(20,30), x.intersp = 0.5)
    title(ylab=paste(s, "ESD (normalized)"), line=2.5, cex.lab=1.8)
    
    par(fig = c(0, 0.6, 0, 0.55), new=T )
    qt = quantile(datam.norm$bioc, probs=c(0.05, 0.95))
    plot(datam.norm$bioc, col="darkgray", xaxt="n", xlab="", ylab="", type="l", lwd=0.1,
         cex.axis=1.6, ylim=c(qt[1], qt[2]))
    axis(1, at=xyr, labels=datam$yr[xyr], cex.axis=1.5)
    lines(bioroll, col="green", lwd=5)
    title(xlab="Time", line=2.5, cex.lab=1.8)
    title(ylab=paste(s, "Biomass (normalized)"), line=2.5, cex.lab=1.8)
    abline(v = x_axis, col="lightgray", lwd=2, lty=3)
    
    par(fig = c(0.6, 1, 0.45, 1), new=T )
    ifft = length(Sizelist) #which(specl == s)+1
    szmod = Sizelist[[ifft]]
    axises = specapprox(1/szmod$freq, szmod$spec)
    plot(axises[[1]], axises[[2]], ylim=c(0,max(axises[[2]])), xlim=c(8, floor(length(datam.norm$mth)/3)), log="x",
         type="l", xlab="", ylab="", 
         xaxt='n', cex.axis=1.6, col="darkgreen", lwd=4)
    x_axis = axis(1, at=seq(0, floor(1/min(szmod$freq) ), 6), labels=F, tick=T)
    abline(v = x_axis, col="lightgray", lwd=2, lty=3)
    title(ylab="Spectrum intensity", line=2.5, cex.lab=1.8)
    
    par(fig = c(0.6, 1, 0, 0.55), new=T )
    szmod = Bioclist[[ifft]]
    axises = specapprox(1/szmod$freq, szmod$spec)
    plot(axises[[1]], axises[[2]], ylim=c(0,max(axises[[2]])), xlim=c(8, floor(length(datam.norm$mth)/3)), log="x",
         type="l", xlab="", ylab="", 
         xaxt='n', cex.axis=1.6, col="green", lwd=4)
    axis(1, at=seq(0, floor(1/min(szmod$freq) ), 12), cex.axis=1.6)
    axis(1, at=seq(0, floor(1/min(szmod$freq) ), 6), labels=F)
    abline(v = x_axis, col="lightgray", lwd=2, lty=3)
    title(ylab="Spectrum intensity", line=2.5, cex.lab=1.8)
    title(xlab="Frequency (Months)", line=2.5, cex.lab=1.8)
    
  }
  graphics.off()
}

scaleplot = function(ax, xlimit, colo, lwd=5, stacker=0){ # Plots and scales automatically according to limit
  x = ax[[1]]
  y = ax[[2]]
  
  ind.scale = which(x <= xlimit+2)
  x = x[ind.scale] ; y = y[ind.scale]
  
  # Checking if negative values due to LOESS
  y[ which(y < 0) ] = 0
  
  matplot(x, y/max(y) + stacker, col = colo, lwd=lwd, type="l", add=T)
}
# Function for plotting the spectra
xmax = 36 ; linew = 4 # Period limit

#### Checking spectrum of Phytoplankton OK ----
i = 2
axz = Sizelist[[i]] ; axbc = Bioclist[[i]]

plot(1/axz$freq, axz$spec/max(axz$spec), type="l", xlim = c(0, 48), ylim = c(0,1))
axz = specapprox(1/axz$freq, axz$spec, cut=T, span = 0.035) # OK
scaleplot(axz , xmax, "sienna4", stacker = 0, lwd = linew) 

plot(1/axbc$freq, axbc$spec/max(axbc$spec), type="l", xlim = c(0, 48), ylim = c(0,1))
axbc = specapprox(1/axbc$freq, axbc$spec, cut=T, span = 0.04) # OK
scaleplot(axbc, xmax, "chartreuse3", stacker = 0, lwd = linew)

#### Checking NAO (log trend) ----
nao = read.table("nao.monthly.txt", sep="", skip=1)
naof = matrix(data.matrix(nao[,2:13]), ncol=1)
timeseq = rep(nao$V1, each = 12)
mthseq = rep(1:12, length(nao$V1))

naof = data.frame(yr = timeseq, mth = mthseq, nao = naof)
naof$trend = filter( naof$nao, filter=filter.wind, sides=2, method="convolution")
naof$date = ymd( paste(naof$yr, "-", naof$mth, "-01", sep="") )

#test.student(naof$mth, naof$yr, naof$nao)
naof$nao = norm.stud(293, naof$nao) # Shift at indice 293

ind.na = which(is.na(naof$trend))
specnao = spec.pgram(fast=T, naof$trend[-ind.na], demean=T, plot=F)

# Smoothing using a filter moving average
f3 = c(1, 2, 1) / 4
spa = filter( specnao$spec, filter = f3 )
plot(1/specnao$freq, specnao$spec, type="l", xlim = c(0, 48), ylim = c(0,2))
lines( 1/specnao$freq, spa, col = "red" )

# Make the filter visually smooth, but values identical to the filter
axnao = graphapprox(spa)

axnaox = 1/seq( min(specnao$freq), max(specnao$freq), length.out = length(axnao) )
plot(1/specnao$freq, specnao$spec, type="l", xlim = c(0, 48), ylim = c(0,2))
lines( 1/specnao$freq, spa, col = "blue" )
lines( axnaox, axnao, col = "red" )

axnao = list( axnaox, axnao ) # Stacking the data

#### Checking wind 2 (CDC) ----
wind22 = read.table("Wind_DWD.txt", sep=";", header=T)

# Removing non valid lines
ind.not.to = which(wind22$FK_TER == -999 | wind22$DK_TER == -999)
wind22 = wind22[-ind.not.to,] # This one is near Helgoland

# Daily, then monthly averages
wind22day = ddply(wind22, .(date), summarize, dir = mean(DK_TER, na.rm=T),
                  force = mean(FK_TER, na.rm=T))
wind22day$yr = year(wind22day$date) ; wind22day$mth = month(wind22day$date)

wind22 = ddply(wind22day, .(yr, mth), summarize, dirm = mean(dir, na.rm=T),
               fm = mean(force, na.rm=T))

wind22$date = ymd( paste(wind22$yr, "-", wind22$mth, "-01", sep="") )

# Checking all dates are in
wind22.t = data.approx(wind22$dirm, wind22$date)
wind22 = data.approx(wind22$fm, wind22$date)
wind22$dirm = wind22.t$dmn

# Trends
wind22$dirtrend = filter(wind22$dirm, filter=filter.wind)
wind22$ftrend = filter(wind22$dmn, filter=filter.wind)

# Spectral analysis
spec.verif = function(datatrend, span = 0.05, cut=T){
  ind.na = which( is.na(datatrend) )
  spec = spec.pgram(fast=F,  datatrend[-ind.na], demean=T, plot=F)
  
  axdata = specapprox(1/spec$freq, spec$spec, span=span, cut=cut)
  
  ind.plot = which( 1/spec$freq <= 50 )
  max = which.max( spec$spec[ind.plot]  )
  plot(1/spec$freq[ind.plot], spec$spec[ind.plot] / spec$spec[ind.plot][max], type="l", xlim=c(0, 50))
  
  ind.plot = which( axdata[[1]] <= 50 )
  max = which.max( axdata[[2]][ind.plot]  )
  lines(axdata[[1]], axdata[[2]] / axdata[[2]][max], col="red")
  return(axdata)
}

#test.student(wind22$mth, wind22$yr, log(wind22$dirtrend) )
# No real shift detected

axdir22 = spec.verif(wind22$dirtrend, span = 0.03, cut = T) # Both are good (or cut and 0.03)
axf22 = spec.verif(wind22$ftrend, span = 0.035, cut = F)

#### Spectra of environment ------------------------------------
temp = listenvw[[1]] ; N = listenvw[[2]] ; si = listenvw[[7]] ; sal = listenvw[[5]] 

temp$date = ymd( paste( temp$yr, "-", temp$mth, "-01", sep="" ) )
N$date = ymd( paste( N$yr, "-", N$mth, "-01", sep="" ) )
si$date = ymd( paste( si$yr, "-", si$mth, "-01", sep="" ) )
sal$date = ymd( paste( sal$yr, "-", sal$mth, "-01", sep="" ) )

ind.na = which(is.na(temp$trend))
spectemp = spec.pgram(fast=T, temp$trend[-ind.na], demean=T, plot=F) 
axtemp = specapprox(1/spectemp$freq, spectemp$spec, span=0.04, cut=F)
plot(1/spectemp$freq, spectemp$spec, type="l", xlim = c(0, 48), ylim = c(0,200))
lines( axtemp[[1]], axtemp[[2]], col = "red" ) 
plot(temp$date, temp$trend, type = "l")

ind.na = which(is.na(N$trend))
specn = spec.pgram(fast=T, N$trend[-ind.na], demean=T, plot=F) 
axn = specapprox(1/specn$freq, specn$spec, span=0.03, cut=F)
plot(1/specn$freq, specn$spec, type="l", xlim = c(0, 48), ylim = c(0,2000))
lines( axn[[1]], axn[[2]], col = "red" ) 
plot(N$date, N$trend, type = "l")

ind.na = which(is.na(si$trend))
specsi = spec.pgram(fast=T, si$trend[-ind.na], demean=T, plot=F) 
axsi = specapprox(1/specsi$freq, specsi$spec, span=0.03, cut=F)
plot(1/specsi$freq, specsi$spec, type="l", xlim = c(0, 48), ylim = c(0,600))
lines( axsi[[1]], axsi[[2]], col = "red" ) 
plot(si$date, si$trend, type = "l")

#inds = test.student(sal$mth, sal$yr, sal$trend)
sal$trend = norm.stud(71, sal$trend) # Shift at indice 73
ind.na = which(is.na(sal$trend))
specsal = spec.pgram(fast=T, sal$trend[-ind.na], demean=T, plot=F) 
axsal = specapprox(1/specsal$freq, specsal$spec, span=0.03, cut=F)
plot(1/specsal$freq, specsal$spec, type="l", xlim = c(0, 48), ylim = c(0,15))
lines( axsal[[1]], axsal[[2]], col = "red" ) 
plot(sal$date, sal$trend, type = "l")

# Plot (all combined) ----------------------------------------------------------

x11(width=14, height=20)

par( fig=c(0.5, 1, 0,1), mar=c(4, 1, 3, 5) ) # Spectra ---------------------
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim=c(10, xmax), ylim=c(0.15, 4.5), log="x" )
axis(1, at=seq(6, xmax, 6), labels=c(6, 12, 18, 24, NA, 36), cex.axis=1.8, font=1)
abline( v = c(12,24), lty=2, lwd=3, col=rgb(0.5, 0.5, 0.5, alpha = 0.4) )

title(xlab = "Period, months", cex.lab=1.8, font.lab=1, line=2.5)
mtext(side = 4, "Spectral intensity, scaled to maximum", cex=1.8, line=3)

y.ticks = c(0, NA, 1)
# Starting with NAO and Wind ---------------------------------------------
scaleplot(axnao, xmax, "dodgerblue3", lwd = linew) 
scaleplot(axdir22 , xmax, "gray20", lwd = linew) 
axis(4, y.ticks, cex.axis=1.7, font=1, las=1)

# Adding Salinity -------------------------------------------------
stacker = 1.2

scaleplot(axsal, xmax, "black", stacker = stacker, lwd = linew) 
axis(4, y.ticks+stacker, labels = F)
abline(h = stacker, col='black', lwd=1.5, lty=1)

# Adding S phy ------------------------------------------------------------------
stacker = 2.4 # 1.2
i = 4
axz = Sizelist[[i]] ; axbc = Bioclist[[i]]

axz = specapprox(1/axz$freq, axz$spec, cut=T, span = 0.05) # OK
axbc = specapprox(1/axbc$freq, axbc$spec, cut=T, span = 0.04) # OK

scaleplot(axz , xmax, "sienna4", stacker = stacker, lwd = linew) 
scaleplot(axbc, xmax, "chartreuse3", stacker = stacker, lwd = linew)
axis(4, y.ticks +stacker, labels=F)
abline(h = stacker, col='black', lwd=1.5, lty=1)

# Adding L phy ------------------------------------------------------------------
stacker = 3.6 # 2.4
i = 2
axz = Sizelist[[i]] ; axbc = Bioclist[[i]]

axz = specapprox(1/axz$freq, axz$spec, cut=T, span = 0.035) # OK
axbc = specapprox(1/axbc$freq, axbc$spec, cut=T, span = 0.04) # OK

scaleplot(axz , xmax, "sienna4", stacker = stacker, lwd = linew) 
scaleplot(axbc, xmax, "chartreuse3", stacker = stacker, lwd = linew)
axis(4, y.ticks +stacker, labels=y.ticks, cex.axis=1.7, font=1, las=1)
abline(h = stacker, col='black', lwd=1.5, lty=1)

par(fig=c(0, 0.5, 0, 0.32), mar=c(4, 5, 2, 3), new=T) # NAO and wind series --------------------
dateslimit = ymd( c( "2008-01-01", "2021-12-01" ) )

# Taking a sample of the TS, otherwise not readable
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim= dateslimit , ylim=c(-2.5, 3) )
axis(4, at=-1:1, labels=F) ; axis(2, at=-1:1, line=0, labels = c(-1, NA, 1), cex.axis=1.7, las=1)

x.dates = naof$date[ which( month(naof$date) == 1) ]
axis(1, at=x.dates, labels=F )
x.labs = ymd("1980-01-01")+years( seq(0, 40, 4) )
axis(1, at=x.labs, labels=year(x.labs), font=1, cex.axis=1.8, tck=-0.08, lwd.ticks = 2, 
     padj=0.6, hadj = 0.25 )
axis(1, at=x.labs+ years(1), labels=F, tck=-0.08, lwd.ticks = 2)

dir.scale = scale(wind22$dirtrend)
y = graphapprox(dir.scale[-which( is.na(dir.scale) )])
x = seq( min(wind22$date[-which( is.na(dir.scale) )]), max(wind22$date[-which( is.na(dir.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="gray30", add=T)

#y = graphapprox(naof$trend[-which( is.na(naof$trend) )])
#x = seq( min(naof$date[-which( is.na(naof$trend) )]), max(naof$date[-which( is.na(naof$trend) )]), length.out = length(y) )

# Taking only the NAO above 2006 (where the phyto starts) to plot, for better visualization
naofplot = naof[ which( naof$yr >= 2006  ), ]
y = graphapprox( scale( naofplot$trend[-which( is.na(naofplot$trend) )] ) )
x = seq( min(naofplot$date[-which( is.na(naofplot$trend) )]), max(naofplot$date[-which( is.na(naofplot$trend) )]), length.out = length(y) )

matplot(x, y, type="l", lwd = linew, col="dodgerblue3", add=T)

par(fig=c(0, 0.5, 0.23, 0.63), new=T) # Salinity series --------------------

sal.scale = scale(sal$trend)
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim= dateslimit, ylim=c(-2.5, 3.8) )
axis(4, at=-1:1, labels=F, line=0) ; axis(2, at=-1:1, labels=F, line=0)

axis(1, at=x.dates, labels=F )
axis(1, at=x.labs, labels=F, font=2, cex.axis=1.5, tck=-0.06, lwd.ticks = 2)
axis(1, at=x.labs+ years(1), labels=F, tck=-0.06, lwd.ticks = 2)

y = graphapprox(sal.scale[-which( is.na(sal.scale) )])
x = seq( min(sal$date[-which( is.na(sal.scale) )]), max(sal$date[-which( is.na(sal.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="black", add=T)

par(fig=c(0, 0.5, 0.45, 0.77), new=T) # S phys series --------------------

physm$esd = (physm$vol *3/4/pi)**(1/3) *2
physm$esdtrend = filter(physm$esd, filter=filter.wind)
bioc.scale = scale(physm$bioctrend)
esd.scale = scale(physm$esdtrend)

matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim= dateslimit, ylim=c(-2, 2.5) )
axis(4, at=-1:1, labels=F, line=0) ; axis(2, at=-1:1, labels=F, line=0)

axis(1, at=x.dates, labels=F )
axis(1, at=x.labs, labels=F, font=2, cex.axis=1.5, tck=-0.06, lwd.ticks = 2)
axis(1, at=x.labs+ years(1), labels=F, tck=-0.06, lwd.ticks = 2)

y = graphapprox(bioc.scale[-which( is.na(bioc.scale) )])
x = seq( min(physm$date[-which( is.na(esd.scale) )]), max(physm$date[-which( is.na(esd.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="chartreuse3", add=T)

y = graphapprox(esd.scale[-which( is.na(esd.scale) )])
matplot(x, y, type="l", lwd = linew, col="sienna4", add=T)

par(fig=c(0, 0.5, 0.67, 1), new=T) # L phy series --------------------

phylm$esd = (phylm$vol *3/4/pi)**(1/3) *2
phylm$esdtrend = filter(phylm$esd, filter=filter.wind)
bioc.scale = scale(phylm$bioctrend)
esd.scale = scale(phylm$esdtrend)

matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim= dateslimit, ylim=c(-2, 4) )
axis(4, at=-1:1, labels=F, line=0) ; axis(2, at=-1:1, labels=c(-1, NA, 1), line=0, cex.axis=1.7, las=1)

axis(1, at=x.dates, labels=F )
axis(1, at=x.labs, labels=F, font=2, cex.axis=1.5, tck=-0.06, lwd.ticks = 2)
axis(1, at=x.labs+ years(1), labels=F, tck=-0.06, lwd.ticks = 2)

y = graphapprox(bioc.scale[-which( is.na(bioc.scale) )])
x = seq( min(phylm$date[-which( is.na(bioc.scale) )]), max(phylm$date[-which( is.na(bioc.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="chartreuse3", add=T)

y = graphapprox(esd.scale[-which( is.na(esd.scale) )])
x = seq( min(phylm$date[-which( is.na(esd.scale) )]), max(phylm$date[-which( is.na(esd.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="sienna4", add=T)

# Adding Text ------------------------------------------------------------------
par(fig=c(0, 0.5, 0, 1), new=T)
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim=c(0, 1), ylim=c(0, 1) )

text(x = 0.7, y = 0.16, "NAO", cex=1.8, col="dodgerblue4", font=1)
text(x = 0.3, y = 0.16, "Wind direction", cex=1.8, col="black", font=1)

text(x = 0.6, y = 0.45, "Salinity", cex=1.8, col="black", font=1, xpd = T)

text(x = 0, y = 0.72, "Small phytoplankton", cex=1.8, col="black", font=1, adj=0)
text(x = 0.65, y = 0.7, "Size", cex=1.8, col="sienna4", font=1)
text(x = 0.57, y = 0.53, "Carbon", cex=1.8, col="chartreuse4", font=1)

text(x = 0.4, y = 0.95, "Large phytoplankton", cex=1.8, col="black", font=1, adj=0)
#text(x = 0.41, y = 1, "Size", cex=1.8, col="sienna4", font=1)
#text(x = 0.7, y = 1, "Carbon", cex=1.8, col="chartreuse4", font=1)

#text(x = 10, y = 4.2, "Diatoms \n / Dinoflagellates C", cex=1.8, col="black", font=1, xpd = T)
#text(x = 22, y = 4.05, "Grazing", cex=1.8, col="darkorange", font=1)

par(fig=c(0, 0.5, 0, 1), new=T) # Legend --------------------
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n")
mtext(side = 2, "Moving average, scaled", font = 1, cex = 1.8, line=2.5)

par( fig = c( 0, 1, 0, 1 ), mar = c( 0, 0, 0, 0 ), new = T )
matplot( -0.9, 1, type = "n", axes = F, xlab = "", ylab = "", 
         xlim = c( -1, 1 ), ylim = c( -1, 1 ), add = F, xpd = T)
text( x = -1, y = 1.03, "(a)", font = 2, cex = 1.8, xpd = T )
text( x = -0.63, y = 1.03, "Smoothed monthly time series", cex = 1.8, xpd = T )

text( x = 0.05, y = 1.03, "(b)", font = 2, cex = 1.8 )
text( x = 0.28, y = 1.03, "Fourier spectra", cex = 1.8 )

mtext(side = 1, adj=0.72, line=-38.5, "Small phytoplankton", font = 1, cex=1.8)
mtext(side = 1, adj=0.78, line=-49, "Large phytoplankton", font = 1, cex=1.8)
## ----

## Fourier spectrum of Dinoflagellates, Nitrate, Phosphate, and Temperature (Supplementary information) ----
i = 8
axz = Sizelist[[i]] ; axbc = Bioclist[[i]]

plot(1/axz$freq, axz$spec/max(axz$spec[1/axz$freq<36]), type="l", xlim = c(0, 48), ylim = c(0,1))
axz = specapprox(1/axz$freq, axz$spec, cut=T, span = 0.035) # OK
scaleplot(axz, 36, "sienna4", stacker = 0, lwd = 2) 
plot(dflsm$date, Sizelist[[i-1]], type = "l")

plot(1/axbc$freq, axbc$spec/max(axbc$spec[1/axbc$freq<36]), type="l", xlim = c(0, 48), ylim = c(0,1))
axbc = specapprox(1/axbc$freq, axbc$spec, cut=T, span = 0.04) # OK
scaleplot(axbc, 36, "chartreuse3", stacker = 0, lwd = 2)
plot(dflsm$date, Bioclist[[i-1]], type = "l")

# Extracting the series and the spectra
temp = listenvw[[1]] ; N = listenvw[[2]] ; si = listenvw[[7]] ; P = listenvw[[3]] 

temp$date = ymd( paste( temp$yr, "-", temp$mth, "-01", sep="" ) )
N$date = ymd( paste( N$yr, "-", N$mth, "-01", sep="" ) )
si$date = ymd( paste( si$yr, "-", si$mth, "-01", sep="" ) )
P$date = ymd( paste( P$yr, "-", P$mth, "-01", sep="" ) )

ind.na = which(is.na(temp$trend))
spectemp = spec.pgram(fast=T, temp$trend[-ind.na], demean=T, plot=F) 
axtemp = specapprox(1/spectemp$freq, spectemp$spec, span=0.03, cut=F)
plot(1/spectemp$freq, spectemp$spec, type="l", xlim = c(0, 48), ylim = c(0,200))
lines( axtemp[[1]], axtemp[[2]], col = "red" ) 
plot(temp$date, temp$trend, type = "l")

N$trend = norm.stud(114, N$trend) # Shift detected
ind.na = which(is.na(N$trend))
specn = spec.pgram(fast=T, N$trend[-ind.na], demean=T, plot=F) 
axn = specapprox(1/specn$freq, specn$spec, span=0.03, cut=F)
plot(1/specn$freq, specn$spec, type="l", xlim = c(0, 48), ylim = c(0,50))
lines( axn[[1]], axn[[2]], col = "red" ) 
plot(N$date, N$trend, type = "l")

ind.na = which(is.na(P$trend))
specp = spec.pgram(fast=T, P$trend[-ind.na], demean=T, plot=F) 
axp = specapprox(1/specp$freq, specp$spec, span=0.03, cut=F)
plot(1/specp$freq, specp$spec, type="l", xlim = c(0, 48), ylim = c(0,1))
lines( axp[[1]], axp[[2]], col = "red" ) 
plot(P$date, P$trend, type = "l")

# Plot 
x11(width=14, height=18)
xmax = 48 ; linew = 4

par( fig=c(0.5, 1, 0,1), mar=c(4, 1, 3, 5) ) # Spectra 
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim=c(10, xmax), ylim=c(0.05, 3.4), log="x" )
axis(1, at=seq(6, xmax, 6), labels=c(6, 12, 18, 24, NA, 36, NA, 48), cex.axis=1.8, font=1)

title(xlab = "Period (months)", cex.lab=1.8, font.lab=1, line=2.5)
mtext(side = 4, "Spectral intensity scaled to maximum", cex=1.8, line=3)

y.ticks = c(0, NA, 1)

# Spectra plots
scaleplot(axtemp , xmax, "firebrick3", lwd = linew) 
axis(4, y.ticks, cex.axis=1.7, font=1, las=1)

scaleplot(axn, xmax, "darkblue", lwd = linew, stacker = 1.2) 
scaleplot(axp , xmax, "dodgerblue", lwd = linew, stacker = 1.2)
abline(h = 1.2, col='black', lwd=1.5, lty=1)
axis(4, y.ticks +1.2, labels = F)

scaleplot(axbc , xmax, "chartreuse3", lwd = linew, stacker = 2.4)
scaleplot(axz , xmax, "sienna4", lwd = linew, stacker = 2.4)
abline(h = 2.4, col='black', lwd=1.5, lty=1)
axis(4, y.ticks + 2.4, labels = y.ticks, cex.axis=1.7, font=1, las=1)

par(fig=c(0, 0.5, 0, 0.4), mar=c(4, 5, 1, 3), new=T) # Plot series
dateslimit = ymd( c( "2008-01-01", "2021-12-01" ) )

# Taking a sample of the TS, otherwise not readable
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim= dateslimit , ylim=c(-2.5, 3) )

x.dates = N$date[ which( N$mth == 1) ]
axis(1, at=x.dates, labels=F )
x.labs = ymd("1980-01-01")+years( seq(0, 40, 4) )
axis(1, at=x.labs, labels=year(x.labs), font=1, cex.axis=1.8, tck=-0.08, lwd.ticks = 2, 
     padj=0.6, hadj = 0.25 )
axis(1, at=x.labs+ years(1), labels=F, tck=-0.08, lwd.ticks = 2)

t.scale = scale(temp$trend)
y = graphapprox(t.scale[-which( is.na(t.scale) )])
x = seq( min(temp$date[-which( is.na(t.scale) )]), max(temp$date[-which( is.na(t.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="firebrick3", add=T)
axis(4, at=-1:1, labels=F) ; axis(2, at=-1:1, line=0, labels = c(-1, NA, 1), cex.axis=1.7, las=1)

par(fig=c(0, 0.5, 0.3, 0.7), new=T)
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim= dateslimit , ylim=c(-2.5, 3.5) )

n.scale = scale(N$trend)
y = graphapprox(n.scale[-which( is.na(n.scale) )])
x = seq( min(N$date[-which( is.na(n.scale) )]), max(N$date[-which( is.na(n.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="darkblue", add=T)

p.scale = scale(P$trend)
y = graphapprox(p.scale[-which( is.na(p.scale) )])
x = seq( min(P$date[-which( is.na(p.scale) )]), max(P$date[-which( is.na(p.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="dodgerblue", add=T)
axis(4, at=-1:1, labels=F) ; axis(2, at=-1:1, line=0, labels = c(-1, NA, 1), cex.axis=1.7, las=1)

axis(1, at=x.dates, labels=F )
axis(1, at=x.labs, labels=F, font=2, cex.axis=1.5, tck=-0.06, lwd.ticks = 2)
axis(1, at=x.labs+ years(1), labels=F, tck=-0.06, lwd.ticks = 2)

par(fig=c(0, 0.5, 0.6, 1), new=T)
matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n", 
         xlim= dateslimit , ylim=c(-2, 2) )

sz.scale = Sizelist[[7]]
y = graphapprox(sz.scale[-which( is.na(sz.scale) )])
x = seq( min(dflsm$date[-which( is.na(sz.scale) )]), max(dflsm$date[-which( is.na(sz.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="sienna4", add=T)
bc.scale = Bioclist[[7]]
y = graphapprox(bc.scale[-which( is.na(bc.scale) )])
x = seq( min(dflsm$date[-which( is.na(bc.scale) )]), max(dflsm$date[-which( is.na(bc.scale) )]), length.out = length(y) )
matplot(x, y, type="l", lwd = linew, col="chartreuse3", add=T)
axis(4, at=-1:1, labels=F) ; axis(2, at=-1:1, line=0, labels = c(-1, NA, 1), cex.axis=1.7, las=1)

axis(1, at=x.dates, labels=F )
axis(1, at=x.labs, labels=F, font=2, cex.axis=1.5, tck=-0.06, lwd.ticks = 2)
axis(1, at=x.labs+ years(1), labels=F, tck=-0.06, lwd.ticks = 2)

par(fig=c(0, 0.5, 0, 1), new=T) # Legend 

matplot( c(1,1), c(0,1), axes=F, xlab="", ylab="", type="n")
mtext(side = 2, "Moving average (year)", font = 1, cex = 1.8, line=2.6)

par( fig = c( 0, 1, 0, 1 ), mar = c( 0, 0, 0, 0 ), new = T )
matplot( -0.9, 1, type = "n", axes = F, xlab = "", ylab = "", 
         xlim = c( -1, 1 ), ylim = c( -1, 1 ), add = F, xpd = T)
text( x = -1, y = 1.03, "(a)", font = 2, cex = 1.8, xpd = T )
text( x = -0.63, y = 1.03, "Smoothed monthly time series", cex = 1.8, xpd = T )

text( x = 0.05, y = 1.03, "(b)", font = 2, cex = 1.8 )
text( x = 0.28, y = 1.03, "Fourier spectra", cex = 1.8 )

mtext(side = 1, adj=0.09, line=-48, "Dinoflagellates", 
      font = 1, cex=1.8)
mtext(side = 1, adj=0.34, line=-38, "Size", col="sienna4",
      font = 1, cex=1.8)
mtext(side = 1, adj=0.3, line=-46.5, "Biomass", col="chartreuse3",
      font = 1, cex=1.8)
mtext(side = 1, adj=0.2, line=-33, expression("PO"[4] *""^3 *""^"-"), 
      col = "dodgerblue", font = 1, cex=1.8)
mtext(side = 1, adj=0.08, line=-33, expression("NO" *""^3 *""^"-"), 
      col = "darkblue", font = 1, cex=1.8)
mtext(side = 1, adj=0.08, line=-17, "Temperature", col="firebrick3", 
      font = 1, cex=1.8)
## ----
