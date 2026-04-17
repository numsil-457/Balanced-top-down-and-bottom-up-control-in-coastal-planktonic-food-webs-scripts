## This script contains functions of the main script for statistical analysis

#### General Functions ####
savefig = function( name, figpath = "" ){
  fname = paste(figpath, name, sep = "")
  dev.copy2pdf( file = fname )
  dev.off()
}

# Filters to smooth the monthly time series 
filter.wind = c( seq(1,7), rev(seq(1,6)) )
filter.wind = filter.wind / sum(filter.wind) # One year filter (on months)

## Approximation for missing dates in phytoplankton
# Replaced by 50% climatology and 50% linear approx
approxer = function(datam){
  
  # Climatology as in Barton et al, (2003)
  mth.weight = ddply(datam, .(mth), summarize, ntot = sum(nb, na.rm = T) ) # NaNs are not accounted
  
  datam = merge( datam, mth.weight, by = "mth" )
  datam = datam[order(datam$yr),] # Retidying after the merge
  
  datam$wi = datam$nb / datam$ntot
  climat = ddply(datam, .(mth), summarize, 
                 smean = sum( vol * wi, na.rm=T ), 
                 bmean = sum( bioc * wi, na.rm=T ), 
                 bvmean = sum( biovol * wi, na.rm=T ),
                 #size.sd = sqrt( sum( nb*(vol - sum( vol * wi, na.rm=T ))**2 ) / (sum(nb)-1) ),
                 bioc.sd = sqrt( sum( nb*( bioc - sum( bioc * wi, na.rm=T ))**2 ) / (sum(nb)-1) ),
                 bioc.qt90 = quantile( rep(bioc, times=nb), probs=0.9 ),
                 bioc.qt10 = quantile( rep(bioc, times=nb), probs=0.1 ) # Weighted quantiles
  )
                 #biovol.sd = sqrt( sum( nb*(biovol - sum( biovol * wi, na.rm=T ))**2 ) / (sum(nb)-1) ) )

  indlast = length(datam$yr)
  cmth = c() # Data vectors
  cyr=c()
  csze=c()
  cbio=c()
  cbiov=c()
  cnb = c()
  
  mini = datam$mth[1] ; yini = datam$yr[1]
  icheck=1
  while( yini <= datam$yr[indlast] & icheck <= length(datam$yr) ){
    if(mini != datam$mth[icheck]){
      csze=c(csze, NaN)
      cbio=c(cbio, NaN)
      cbiov=c(cbiov, NaN)
    }else{
      csze=c(csze, datam$vol[icheck])
      cbio=c(cbio, datam$bioc[icheck])
      cbiov=c(cbiov, datam$biovol[icheck])
      icheck=icheck+1
    }
    cmth=c(cmth, mini)
    cyr=c(cyr, yini)
    
    mini=mini+1
    if(mini>12){mini=1 ; yini=yini+1}
  }
  
  datam = data.frame(yr=cyr, mth=cmth, bioc=cbio, vol=csze, biovol=cbiov)
  
  # Approximation : Combining linear approx and climatology (50/50)
  
  indnan = which(is.nan(datam$bioc)) # Replace by approx
  
  final.approx=approx(1:length(datam$bioc), datam$bioc, xout=1:length(datam$bioc))
  datam$bioc = final.approx$y
  
  final.approx=approx(1:length(datam$vol), datam$vol, xout=1:length(datam$vol)) # If climatology not available
  datam$vol = final.approx$y
  
  final.approx=approx(1:length(datam$biovol), datam$biovol, xout=1:length(datam$biovol)) # If climatology not available
  datam$biovol = final.approx$y
  
  for(nan in indnan){
    datam$bioc[nan] = mean( c(climat$bmean[ which(climat$mth==datam$mth[nan]) ], datam$bioc[nan]), na.rm=T)
    datam$vol[nan] = mean( c(climat$smean[ which(climat$mth==datam$mth[nan]) ], datam$vol[nan]), na.rm=T)
    datam$biovol[nan] = mean( c(climat$bvmean[ which(climat$mth==datam$mth[nan]) ], datam$biovol[nan]), na.rm=T)
  }
  return( list(datam, climat) ) }

# Approximation of monthly means for a general variable
data.approx = function(datal, data.dates){ # Input monthly data series
  mindate = min(data.dates) ; maxdate = max(data.dates)
  
  cmth =c() # Data vectors
  cyr=c()
  cdata=c()
  
  mini = month(mindate) ; yini = year(mindate)
  icheck=1
  while(mini != month(maxdate) | yini < year(maxdate) ){ # TS complete ?
    if( mini == month(data.dates[icheck]) & yini == year(data.dates[icheck]) ){
      cdata=c(cdata, datal[icheck])
      icheck=icheck+1
    }else{
      cdata=c(cdata, NaN)
    }
    cmth=c(cmth, mini)
    cyr=c(cyr, yini)
    
    mini=mini+1
    if(mini>12){mini=1 ; yini=yini+1}
  }
  
  cmth = c(cmth, mini)
  cyr = c(cyr, yini)
  cdata = c(cdata, datal[length(datal)] ) # Last one was not included
  
  datam = data.frame(yr=cyr, mth=cmth, dmn = cdata)
  climat=ddply(datam, .(mth), summarize, cmn=mean(dmn, na.rm=T) )
  
  # Approximation : Combining linear approx and climatology (50/50)
  indnan = which(is.nan(datam$dmn)) # Replace by approx
  
  final.approx=approx(1:length(datam$dmn), datam$dmn, xout=1:length(datam$dmn))
  datam$dmn = final.approx$y
  
  for(nan in indnan){
    datam$dmn[nan] = mean( c(climat$cmn[ which(climat$mth==datam$mth[nan]) ], datam$dmn[nan]), na.rm=T)
  }
  
  datam$date = ymd( paste(datam$yr, "-", datam$mth, "-01", sep="") )
  return(datam)
}

# Makes plots smoother
graphapprox = function(data){         # Data series
  points = seq(1, 3, length.out = 30) # Poly regression on 3 points each, 10 points per month for smooth curve 
  
  curve = c()
  indice = 2
  predni = 1:3
  
  while(indice < length(data)){
    model = lm( data[ (indice-1):(indice+1) ] ~ poly( predni, 2 ))
    predi = predict(model, newdata = data.frame(predni = points))
    curve  = c(curve, predi[-30])
    indice = indice+2
  }
  curve = c(curve, predi[30]) # excluding ind 30 in loop to avoid repeating points
  
  # Adding last points
  if(indice == length(data)){
    indice.last = length(data) -1
    model = lm( data[ (indice.last-1):(indice.last+1) ] ~ poly( predni, 2 ))
    predi = predict(model, newdata = data.frame(predni = points))
    curve  = c(curve, predi[16:30]) # Adding data for the last point 
  }
  return( curve ) 
}
#### General Functions ####

## Compute phytoplankton monthly time series 
seriesphy = function(datas){ 
  dataday = ddply(datas, .(USI, date, stationID), summarize, 
                  biocday = sum(bioC, na.rm=T), # Sum biomass per day
                  biovday = sum(biovolume, na.rm=T) )
  
  datas$relbios = NA # Placeholder
  
  # Calculating the weights for the mean cell size
  for(d in unique(datas$USI)){
    ind1 = which(datas$USI == d)
    ind2 = which(dataday$USI == d)
    datas$relbios[ind1] = datas$bioC[ind1]/dataday$biocday[ind2] # Biomass weight for size
  }
  dataphyl = ddply(datas, .(USI, date), summarize, meanvol = sum(relbios*mean.size, na.rm=T))
  dataday$meanvol = dataphyl$meanvol
  
  ## Adding temporal variables to the dataset
  dataday$mth = month(dataday$date) ; dataday$yr = year(dataday$date)
  datasi = split(dataday, dataday$stationID)
  statin = names(datasi)
  
  datastat = data.frame()
  datastat = dataday # Storing the data
  
  datam = ddply(datastat, .(yr, mth), summarize,
                bioc = exp( mean(log(biocday), na.rm=T) ), 
                vol = exp( mean(log(meanvol), na.rm=T) ), 
                biovol = exp( mean(log(biovday), na.rm=T) ),
                nb = sum(mth) )
  
  datam$nb = datam$nb / datam$mth # Number of observations, to compute climatology as in Barton et al, (2003)
  datam.final = approxer(datam) # Filling missing months
  
  return( datam.final )
}

## Compute environmental data monthly time series
# Function to compute the climatology for every environmental variable
climat.env = function(data.env){
  data.env$clim.use = 1
  data.env$clim.use[ which( is.na(data.env$dmn) ) ] = 0 # If NA, ignored to compute the climatology
  mth.env = ddply( data.env, .(mth), summarize, tot.nb = sum( nb * clim.use ) )
  
  data.env = merge( data.env, mth.env, by = "mth" )
  data.env$wi = data.env$nb / data.env$tot.nb
  
  clim.comp = ddply( data.env, .(mth), summarize,
                     dmn.qt90 = quantile( rep(dmn, times=nb), probs=0.9 ),  # Weighted quantiles
                     dmn.qt10 = quantile( rep(dmn, times=nb), probs=0.1 ),
                     dmn = sum( dmn * wi, na.rm = T ) )
  return(clim.comp)} 

# Computes approximated series 
envapprox = function(waddenEnv){ 
  
  ## Creating the working dataframe of environmental data
  tempWadden = ddply(waddenEnv, .(date), summarize, mean.temp = mean( as.numeric(temperature), na.rm = T))
  NWadden = ddply(waddenEnv, .(date), summarize, mean.N = mean( as.numeric(nitrate), na.rm = T))
  PWadden = ddply(waddenEnv, .(date), summarize, mean.P = mean( as.numeric(phos), na.rm = T))
  phWadden = ddply(waddenEnv, .(date), summarize, mean.pH = mean( as.numeric(pH), na.rm = T))
  salWadden = ddply(waddenEnv, .(date), summarize, mean.sal = mean( as.numeric(salinity), na.rm = T))
  partWadden = ddply(waddenEnv, .(date), summarize, mean.part = mean( as.numeric(suspended.particulates), na.rm = T))
  siWadden = ddply(waddenEnv, .(date), summarize, mean.si = mean( as.numeric(silicon), na.rm = T))
  ##
  
  data=tempWadden # Monthly series of each
  data$mth = month(data$date) ; data$yr = year(data$date)
  # temp.m = ddply(data, .(yr, mth), summarize, 
  #                dmn = quantile(mean.temp, 0.5, na.rm=T), nb = sum(mth) )
  temp.m = ddply(data, .(yr, mth), summarize, 
                 dmn = mean(mean.temp, na.rm=T), nb = sum(mth) )
  
  data=NWadden
  data$mth = month(data$date) ; data$yr = year(data$date)
  # N.m = ddply(data, .(yr, mth), summarize, 
  #             dmn = quantile(mean.N, 0.5, na.rm=T), nb = sum(mth) )
  N.m = ddply(data, .(yr, mth), summarize, 
              dmn = mean(mean.N, na.rm=T), nb = sum(mth) )
  
  data=PWadden
  data$mth = month(data$date) ; data$yr = year(data$date)
  P.m = ddply(data, .(yr, mth), summarize, 
              dmn = mean(mean.P, na.rm=T), nb = sum(mth) )
  
  data=phWadden
  data$mth = month(data$date) ; data$yr = year(data$date)
  pH.m = ddply(data, .(yr, mth), summarize, 
               dmn = mean(mean.pH, na.rm=T), nb = sum(mth) )
  
  data=salWadden
  data$mth = month(data$date) ; data$yr = year(data$date)
  sal.m = ddply(data, .(yr, mth), 
                summarize, dmn = mean(mean.sal, na.rm=T), nb = sum(mth) )
  
  data=partWadden
  data$mth = month(data$date) ; data$yr = year(data$date)
  part.m = ddply(data, .(yr, mth), 
                 summarize, dmn = mean(mean.part, na.rm=T), nb = sum(mth) )
  
  data=siWadden
  data$mth = month(data$date) ; data$yr = year(data$date)
  si.m = ddply(data, .(yr, mth), 
               summarize, dmn = mean(mean.si, na.rm=T), nb = sum(mth) )
  
  ## Aproximation for missing dates
  
  # Replaced by the climatology of the month and linear approximation (50/50)
  for(ilist in 1:length(list(temp.m, N.m, P.m, pH.m, sal.m, part.m, si.m)) ){
    dataim = list(temp.m, N.m, P.m, pH.m, sal.m, part.m, si.m)[[ilist]]
    climat = climat.env(dataim)
    mthr = range(dataim$mth)
    yrr = range(dataim$yr)
    
    # Checking for missing months
    im = mthr[1] ; ilines=1
    dwseries=c()
    while(ilines <= length(dataim$yr)){
      if(im != dataim$mth[ilines]){
        dwseries = c(dwseries, NaN)
      }
      else{
        dwseries = c(dwseries, dataim$dmn[ilines])
        ilines = ilines+1
      }
      im=im+1
      if(im>mthr[2]){im=1}
    }
    dataim = data.frame(yr=rep(yrr[1]:yrr[2], each=mthr[2]-mthr[1]+1),
                        mth=rep(mthr[1]:mthr[2], yrr[2]-yrr[1]+1),
                        dmn=dwseries )
    dataim.approx = approx(1:length(dataim$dmn), dataim$dmn, xout = 1:length(dataim$dmn))
    
    indnan = which(is.nan(dataim$dmn)) # Replace by approx
    dataim$dmn = dataim.approx$y
    
    for(nan in indnan){
      dataim$dmn[nan] = mean( c(climat$dmn[ which(climat$mth==dataim$mth[nan]) ], dataim$dmn[nan]), na.rm=T)
    }
    if(ilist==1){temp.m=dataim ; t.clim = climat}
    if(ilist==2){N.m=dataim ; n.clim = climat}
    if(ilist==3){P.m=dataim ; p.clim = climat}
    if(ilist==4){pH.m=dataim ; pH.clim = climat}
    if(ilist==5){sal.m=dataim ; sal.clim = climat}
    if(ilist==6){part.m=dataim ; part.clim = climat}
    if(ilist==7){si.m=dataim ; si.clim = climat}
  }
  return( list(temp.m, N.m, P.m, pH.m, sal.m, part.m, si.m,
               t.clim, n.clim, p.clim, pH.clim, sal.clim, part.clim, si.clim) )}

# Computing moving average for all variables
trendenv = function(listenv){
  for(i in ( 1:( length(listenv) / 2 ) ) ){ 
    data = listenv[[i]]
    
    listenv[[i]]$trend = filter(data$dmn, filter=filter.wind, sides=2, method="convolution")
  }
  return(listenv)
}

## Compute fish data monthly time series

fish.approx = function(datah, names.fish){
  indh = which(datah$Species %in% names.fish) # Extracting the species
  datahi = datah[indh,]
  datahi = datahi[order(datahi$month),] # Ordering the dataset according to year and month
  datahi = datahi[order(datahi$year),]
  
  datahi = ddply( datahi, .(year, month), summarize, t = sum(t, na.rm = T),
                  k. = sum(k., na.rm = T), X..kg.1 = mean(X..kg.1, na.rm = T) )
  
  clim = ddply(datahi, .(month), summarize, 
               t = mean(t, na.rm=T), 
               k = mean(k., na.rm=T),
               kg = mean(X..kg.1, na.rm=T)) # Climatology
  
  ## Approx algorithm
  mthc = c()
  yrc = c()
  tc = c()
  kc = c()
  kgc = c()
  
  mc = c()
  yc = c()
  yini = datahi$year[1] ; mini = datahi$month[1]
  
  # Checking for missing months in the series
  iline = 1
  while(iline <= length(datahi$year)){                           # Preparing table
    if(mini == datahi$month[iline] & yini == datahi$year[iline]){ # Date known
      tc = c(tc, datahi$t[iline])
      kc = c(kc, datahi$k.[iline])
      kgc = c(kgc, datahi$X..kg.1[iline])
      iline = iline +1
    }else{
      tc = c(tc, NA)
      kc = c(kc, NA)
      kgc = c(kgc, NA)
    }
    mc=c(mc, mini)
    yc = c(yc, yini)
    mini = mini + 1
    if(mini > 12){mini = 1 ; yini = yini+1}
  }
  dataf = data.frame(month = mc, year = yc, t = tc, k = kc, kg = kgc) 
  ind.na = which(is.infinite(dataf$kg))
  dataf$kg[ind.na] = NA
  
  # Approximation using a 50% linear approx and a 50% climatology
  approx.frame = data.frame(year = dataf$year)
  apxd = approx(1:length(dataf$year), dataf$t, xout = 1:length(dataf$year))
  approx.frame = cbind(approx.frame, apxd$y)
  
  apxd = approx(1:length(dataf$year), dataf$k, xout = 1:length(dataf$year))
  approx.frame = cbind(approx.frame, apxd$y)
  
  apxd = approx(1:length(dataf$year), dataf$kg, xout = 1:length(dataf$year))
  approx.frame = cbind(approx.frame, apxd$y)
  
  names(approx.frame) = c("year", "t", "k", "kg") 
  approx.frame$month = dataf$month
  
  # Filling the NAs
  for(iline in 1:length(dataf$year)){ 
    dataline = dataf[iline,]
    
    ind.na = which(is.na(dataline))
    ind.mth = which(clim$month == approx.frame$month[iline])
    
    dataf[iline, ind.na] = apply( rbind( approx.frame[iline, ind.na-1], clim[ind.mth, ind.na-1]), 2, mean, na.rm=T )
  }
  
  return(dataf)
}

## Compute zooplankton data monthly time series

zooseries = function(data){ # Computing the monthly series and climatology
  stationffs = ddply(data, .(Station, Date), summarize, dw.tot = sum(DryW, na.rm=T))
  data = merge(data, stationffs, by=c("Station", "Date"))
  # Weights as biovolume for computing the mean cell size
  data$relbiov = data$DryW / data$dw.tot
  
  data$size.weight = data$mean.vol * data$relbiov # Computing the monthly means
  data = ddply(data, .(Date, Station), summarize, 
               DryW = mean(dw.tot, na.rm=T), 
               mean.size = sum(size.weight, na.rm=T))
  
  # Approximation for missing dates
  data$mth = month(data$Date) # Adding temporal variables
  data$yr = year(data$Date)
  dataim = ddply(data, .(yr, mth), summarize,
                 dw = exp( mean( log(DryW+1), na.rm=T) ),
                 size = exp( mean( log(mean.size), na.rm=T) ),
                 nb = sum( mth ) )
  # dataim = ddply(data, .(yr, mth), summarize,
  #                dw = mean( DryW, na.rm=T),
  #                size = exp( mean( log(mean.size), na.rm=T) ),
  #                nb = sum( mth ) )
  
  # Replaced by the climatology (as in Barton et al, 2003) of the month
  weight.mth = ddply( dataim, .(mth), summarize, tot.nb = sum( nb ) )
  dataim = merge( dataim, weight.mth,  by = "mth")
  dataim = dataim[ order(dataim$yr), ]
  dataim$wi = dataim$nb / dataim$tot.nb
  
  climat = ddply(dataim, .(mth), summarize, # Climatology
                 dw.qt90 = quantile( rep(dw, times=nb), probs=0.90 ),  # Weighted quantiles
                 dw.qt10 = quantile( rep(dw, times=nb), probs=0.10 ),
                 dw = sum( dw * wi, na.rm = T ),
                 sze = sum( size * wi, na.rm = T )
                  )
                 # dw.sd = sqrt( sum( nb*(dw - sum( dw * wi, na.rm=T ))**2 ) / (sum(nb)-1) ),
                 # sze.sd = sqrt( sum( nb*(sze - sum( sze * wi, na.rm=T ))**2 ) / (sum(nb)-1) ) )
  
  ## There is no SD, as there is only one point per month usually
  # climat = ddply(dataim, .(mth), summarize, # Climatology
  #                dw = exp( sum( log(dw) * wi, na.rm = T ) ), 
  #                sze = exp( sum( log(size) * wi, na.rm = T ) ),
  #                dw.sd = sqrt( sum( nb*( log(dw) - sum( log(dw) * wi, na.rm=T ))**2 ) / (sum(nb)-1) ) )
  #                #sze.sd = sqrt( sum( nb*(sze - sum( sze * wi, na.rm=T ))**2 ) / (sum(nb)-1) ) ) 
  # 
  # climat$dw.sd.min = exp( log(climat$dw) - climat$dw.sd )
  # climat$dw.sd.max = exp( log(climat$dw) + climat$dw.sd )
  
  # Checking for missing months
  mthr = range(climat$mth)
  yrr = range(dataim$yr)
  
  im = mthr[1]
  iyr = yrr[1]
  ilines=1 # Used for counting months
  
  dwseries=c()
  szeseries = c()
  while(iyr <= max(dataim$yr) & ilines <=length(dataim$yr)){
    if(im != dataim$mth[ilines] | iyr != dataim$yr[ilines]){
      dwseries = c(dwseries, NaN)
      szeseries = c(szeseries, NaN)
    }
    else{
      dwseries = c(dwseries, dataim$dw[ilines])
      szeseries = c(szeseries, dataim$size[ilines])
      ilines = ilines+1
    }
    im=im+1
    if(im>mthr[2]){im=mthr[1] ; iyr=iyr+1}
  } # The output is a time series with NaN where monthly data are missing
  
  # Filling the missing months with linear approximation first
  indnan = which(is.nan(dwseries))
  final.approx.dw = approx(1:length(dwseries), dwseries, xout = 1:length(dwseries))
  dwseries = final.approx.dw$y
  final.approx.sze = approx(1:length(szeseries), szeseries, xout = 1:length(szeseries))
  szeseries = final.approx.sze$y
  
  # Then computing a mean of the climatology and linear approximation
  mthl = rep(mthr[1]:mthr[2], yrr[2]-yrr[1]+1)[1:length(dwseries)]
  for(nan in indnan){
    dwseries[nan] = mean( c(climat$dw[which(climat$mth==mthl[nan])], dwseries[nan]), na.rm=T)  
    szeseries[nan] = mean( c(climat$sze[which(climat$mth==mthl[nan])], szeseries[nan]), na.rm=T)  
  }
  # Storing the monthly series in a dataframe
  dataim = data.frame(yr=rep(yrr[1]:yrr[2], each=mthr[2]-mthr[1]+1)[1:length(dwseries)],
                      mth=rep(mthr[1]:mthr[2], yrr[2]-yrr[1]+1)[1:length(dwseries)],
                      dw=dwseries, sze = szeseries)
  
  # Approximation (linear) for full year, including winter
  indlast = length(dataim$yr)
  cmth =c() # Data vectors
  cyr=c()
  cdw=c() ; csze = c()
  
  mini = dataim$mth[1] ; yini = dataim$yr[1]
  icheck=1
  while( mini != dataim$mth[indlast] | yini < dataim$yr[indlast]){
    if(mini != dataim$mth[icheck]){
      cdw=c(cdw, NaN)
      csze=c(csze, NaN)
    }else{
      cdw=c(cdw, dataim$dw[icheck])
      csze=c(csze, dataim$sze[icheck])
      icheck=icheck+1
    }
    cmth=c(cmth, mini)
    cyr=c(cyr, yini)
    
    mini=mini+1
    if(mini>12){mini=1 ; yini=yini+1}
  }
  if(icheck == indlast){
    cyr=c(cyr, dataim$yr[icheck])
    cmth=c(cmth, dataim$mth[icheck])
    cdw=c(cdw, dataim$dw[icheck]) ; csze=c(csze, dataim$sze[icheck])
  }
  
  # Storing the approximated dataset in a dataframe
  datam = data.frame(yr=cyr, mth=cmth, dw=cdw, sze = csze)
  
  # Marking the points using only linear approximation, if needed
  datam$aprx = datam$yr*0 +1
  indnan = which(is.nan(datam$dw))
  datam$aprx[indnan] = 0
  indlast = length(datam$aprx)
  datam.approx.dw = approx(1:indlast, datam$dw, xout=1:indlast) # Approximated series for missing months (Jan-Fev and Nov-Dec)
  datam$dataprx = datam.approx.dw$y
  
  datam.approx.sze = approx(1:indlast, datam$sze, xout=1:indlast)
  datam$szaprx = datam.approx.sze$y
  
  # datam[ c("dw", "sze", "dataprx", "szaprx") ] =
  #   exp( datam[ c("dw", "sze", "dataprx", "szaprx") ] )
  return( list( datam, climat ) ) }

#### Compute the light penetration time series

attenuation = function(sal, spm){
  # Parameters for attenuation model (Tian et al, 2009)
  Kw = 2.06      # m-1
  esal = 0.05714 # m-1.psu-1
  espm = 0.2e-4  # m2.mg-1
  
  Kb = Kw - esal * sal
  Kspm = espm * spm
  Kd = Kb + Kspm # Attenuation coeff.
  return(Kd)
}  

#### Read the NAO data from NOAA
read.nao = function(){
  # Read the NAO information from NOAA
  df = read_html('https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table')
  x = html_nodes(df, 'p')
  x = x %>% html_text()
  
  # Convert the html text to a dataframe
  x1 = strsplit(x, ' ')[[1]]
  x1 = strsplit(x, '\\s|\\n', fixed=F)[[1]]
  x1 = x1[-which(x1=='')]
  
  # Add NAs for the missing months of the current year, if any
  x1 = c(NA, x1) # The NA at the upper left of the table, empty box
  x1 = c(x1, rep(NA, (13-length(x1)%%13) ) ) # All the empty boxes left by the incomplete measurements for the current year
  
  xmat = matrix(x1, ncol = 13, byrow=T)
  nao = matrix(data.matrix( t(xmat[-1,2:13]) ), ncol=1)
  timeseq = rep(xmat[-1,1], each = 12)
  mthseq = rep(1:12, length( unique(timeseq) ))
  naof = data.frame(yr = as.numeric(timeseq),
                    mth = as.numeric(mthseq),
                    nao = as.numeric(nao))
  
  return(naof)
}
