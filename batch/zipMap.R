# Orig from MTabone
library('sp') # spatial data 
library('maps')
library('maptools')
library('plyr') # tools for splitting and combining data
library('fields') # curve, surface, and function fitting
library('classInt')
library('stringr') # contains padding function
data(countyMapEnv) # loads up county and state map data (?)
data(stateMapEnv)
 
# NOITE: the shape files used by this code must be in the current working directory.

# Merge outside data with area spatial data frame
append.spdf <- function(sp.df,sp.df.col,df.only,df.col){
  # sp.df     - spatial polygons data frame (SPDF) e.g. census shapefile data
  # sp.df.col - linking column of the spatial polygons data frame
  # df.only   - other data frame (e.g. PG&E CSV file)
  # df.col    - key column for the data frame

  x <- merge(sp.df@data,df.only,by.x = sp.df.col, by.y = df.col, all.x = TRUE, all.y = FALSE)
  m1 <- with(sp.df@data,get(sp.df.col))
  m2 <- with(x,get(sp.df.col))
  id <- match(m1,m2)
  x <- x[id,]
  return (x)
}

# PLOT: Simple function that plots a column of continuous values on a color scale
calZipPlot <- function (calzip,color.column,main=NULL,precis=0,colorMap=NULL,legend.cex=0.8,legend.title=NULL,legend.bg=NULL){
  # calzip is the california zipcodes SPDF
  # color.column is the column of the dataframe to be colorcoded
  # main is an optional plot title
  # precis sets the precision of the color interval boundaries
  if(is.null(colorMap)) { colorMap = brewer.pal(9,"Reds") } #default color pallette
  
  colorval <- with(calzip@data,get(color.column))
  
  if(is.character(colorval)) { colorval = factor(colorval) }
  if(is.factor(colorval))    { col = mapColors(colorval,colorMap) }
  else {
    quantiles <- classIntervals(colorval,n=7,style="quantile",dataPrecision=precis,na.rm=T)
    col <- findColours(quantiles, colorMap) 
    # plot(quantiles, pal = colorMap) 
  }
  
  # make a special color for NAN entries
  col[is.na(colorval)] <- 'grey80'
  col[grepl(c('XX'),calzip$NAME)] <- 'grey85' # XX means "large land areas such as national parks"
  col[grepl(c('HH'),calzip$NAME)] <- 'grey60'   # HH means body of water
  
  # this is the part that actually plots
  op = par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.1,1.1,1.1,1.1))
  plot(calzip,col=col,border=col)
  map("county",region="california",add=TRUE)
  if(is.factor(colorval)) {
    #colorval = factor(colorval[colorval != '']) # kill blank assignments
    legend(x=-117.52516,y=41.94532,fill=colorMap,legend=levels(colorval),
           cex=legend.cex,title=legend.title,y.intersp=1,bty='n',bg=legend.bg)
  }
  else {
    op = par(family='mono') # fixed width font
    legend(x=-117.52516,y=41.94532,fill=attr(col,"palette"),legend=names(attr(col,"table")),
           cex=legend.cex,title=legend.title,y.intersp=1,bty='n',bg=legend.bg)
    par(op)
  }

  
  if (is.null(main)){
    main = 'Map by zip code'
  }
  title(main)
  par(op)
}

getCalSPDF = function() {
  shpData = readShapePoly('zt06_d00.shp')
  proj4string(shpData) <- "+proj=longlat +datum=WGS84" # projection instructions
  return(shpData)
}
calMap = function(df,plotCol,zipCol=NULL,calZipData=NULL,main=NULL,colorMap=NULL,legend.cex=0.8,legend.title=NULL) {
  if(is.null(colorMap))   { colorMap   <- brewer.pal(7,"Reds") }
  if(is.null(calZipData)) { calZipData <- getCalSPDF()         }
  if(is.null(zipCol)) {
    zipCol = grep('^zip',colnames(df),value=T)[1]
  }
  calZipData@data<-append.spdf(calZipData,'NAME',df,zipCol)
  calZipPlot(calZipData,plotCol,main,colorMap=colorMap,legend.cex=legend.cex,legend.title=legend.title)
}

ggCalMap = function(df,zipCol,plotCol,calZipData=NULL,main=NULL) {
  if(is.null(calZipData)) {
    calZipData <- getCalSPDF()
  }
  dff = fortify(calZipData)
  print(colnames(dff))
  ca <- get_map('california', zoom = 6)
  caMap = ggmap(ca,legend = 'topleft')
  caMap + geom_polygon(aes(x = long, y = lat, group=group,fill=lag1,alpha = .4),data=dff)
}

# example usage: create data frame with zip codes and corresponding values
# pass the data frame into the calMap function with the name of the column to use
# for coloring zip code regions on the map. Color map and various formatting options  can be supplied as well.
demo=F
if(demo) {
  sample = data.frame(zip=93000:95000,vals=1:2001)
  calMap(sample,'vals',main='Sample plot',colorMap=rev(brewer.pal(9,"Blues")) )
}


#calMap(DATA_SOURCE$getZipCounts(),'count','zip5',main='Meter count by zip code',colorMap=brewer.pal(9,"Blues") )
#calMap(DATA_SOURCE$getZipData(),'cecclmzn','zip5',main='CEC climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )
#calMap(DATA_SOURCE$getZipData(),'climate' ,'zip5',main='PGE climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )
#ws = getWeatherSummary()
#calMap(ws,'tout','zip5',main='Mean annual temperature',colorMap=rev(brewer.pal(9,"RdBu")) )
#calMap(ws,'rain','zip5',main='Total annual rain (inches)',colorMap=brewer.pal(9,"Blues") )
#calMap(ws,'dp','zip5',main='Mean annual dew point (F)',colorMap=brewer.pal(9,"Purples") )
