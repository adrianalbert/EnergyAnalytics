require('solaR')
require('ggplot2')
require(gridExtra)

# zone elevation and azimuth breaks.
# elevation greater than the max here are assugned to the overheat zone z0
# elevation < 0 are assigned to the night zone
eBr = c(60,30,0)
aBr = c(seq(0,360,45))

nameZone = function(azim,elev,azimBreaks,elevBreaks) {
  if(elev <= 0) return('night')
  if(elev >= max(elevBreaks)) return('z0')
  for (i in 1:(length(azimBreaks)-1)) {
    for (j in 1:(length(elevBreaks)-1)) {
      aBounds = sort(c(azimBreaks[i],azimBreaks[i+1])) # ensure we know which is smallest
      eBounds = sort(c(elevBreaks[j],elevBreaks[j+1]))
      if( (azim >= aBounds[1] & aBounds[2] > azim) & 
            (elev >= eBounds[1] & eBounds[2] > elev) ) return(paste('z',i,'.',j,sep=''))
    }
  }
  return('?')
}

solarGeom = function(resData,lat=37.87,azimBreaks=seq(0,360,45),elevBreaks=c(45,0)) {
  # Berkeley 37.8717 N, 122.2728 W
  # calcSol computes the angles which describe the intradaily apparent movement of the Sun from the Earth
  # solObj is an S4 class with lots of solar info
  # recall that 'slots' in S4 objects are accessed via the @ operator
  solObj = calcSol(lat,sample="hour",BTi=resData$dates,EoT=T,keep.night=T)
  solarGeom = as.data.frameI(solObj)
  # the solI slot is a zoo object and we want the time stamps
  solarGeom$dates     = time(solObj@solI)
  solarGeom$elevation = solarGeom$AlS * 180/pi   # convery from radians
  solarGeom$azimuth   = solarGeom$AzS * 180/pi   # convery from radians
  solarGeom$azimuth   = solarGeom$azimuth %% 360 # no negative or > 360 degrees
  
  # get named zones based on breaks into a factor
  solarGeom$zone      = factor(apply(as.matrix(solarGeom[,c('azimuth','elevation')]),1,
                                   function(X) nameZone(X[1],X[2],azimBreaks,elevBreaks)))
  #solarGeom$AlS is the solar elevation
  #solarGeom$AzS is the solar asimuth
  class(solarGeom) <- c('solarGeom',class(solarGeom))
  return(solarGeom)
}

plot.solarGeom = function(solarGeom,azimBreaks=seq(0,360,45),elevBreaks=c(45,0)) {
  g = ggplot(data.frame(azimuth=c(0,360),elevation=rep(max(elevBreaks),2)),aes(y=elevation,x=azimuth)) + 
              geom_polygon(color="grey",size=1,fill=NA) + 
              coord_polar(start=pi) + 
              scale_x_continuous(breaks=c(seq(0,315,by=45)),limits=c(0,360), # Note that 0 is S in this case.
                  labels=paste(c('S','SW','W','NW','N','NE','E','SE'),seq(0,315,by=45)) ) + 
              scale_y_reverse(breaks=c(seq(0,90,by=15)),limits=c(90,0)) + 
              labs(title="Annual solar path for each of 24 hours of the day") + 
              geom_text(x=0,y=-90,label='z0',color='grey') # no reverse axis for text make negative instead
  
  for (i in 1:(length(azimBreaks)-1)) {
    for (j in 1:(length(elevBreaks)-1)) {
      zone = data.frame(elevation=c(elevBreaks[j],elevBreaks[j+1],elevBreaks[j+1],elevBreaks[j],  elevBreaks[j]),
                         azimuth =c(azimBreaks[i],azimBreaks[i],  azimBreaks[i+1],azimBreaks[i+1],azimBreaks[i]))
      g = g + geom_polygon(data=zone, colour="gray",size=1,fill=NA) + 
              geom_text(x=mean(c(azimBreaks[i],azimBreaks[i+1])),
                        y=-1*mean(c(elevBreaks[j],elevBreaks[j+1])), # no reverse axis for text make negative instead
                        label=paste('z',i,'.',j,sep=''),color='grey')
      
    }
  }
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  g = g + geom_point(data=solarGeom,aes(x=azimuth,y=elevation,color=zone)) + scale_color_brewer(palette='YlGnBu') #color=hour(dates) )) color=zone))
  grid.arrange(g)
  #return(g)
}


# View of whole path, including below horizon 
# ggplot(solarGeom,aes(x=azimuth,y=elevation)) + 
#         #geom_point(aes(color=factor(aman))) + 
#         #scale_colour_manual(values=c('grey','black')) + 
#         geom_point(aes(color=factor(month))) + 
#         coord_polar() + 
#         scale_x_continuous(breaks=c(seq(0,330,by=30)),limits=c(0,360)) + 
#         scale_y_continuous(breaks=c(seq(-90,90,by=15))) + 
#         geom_line(data=horizon, colour="gray",size=2) + 
#         labs(title="Annual solar path for each of 24 hours of the day")
#          # + theme(panel.grid.minor.y = element_line(colour="white", size=0.5), panel.grid.major.y = element_line(colour="blue", size=0.5))

  
