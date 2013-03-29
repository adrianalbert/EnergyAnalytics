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
 
## Set the working directory
conf.basePath = '/Users/mehdyson/Documents/wharton/mapping'
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
}
setwd(conf.basePath)

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
calzipplot <- function (calzip,color.column,main=NULL,precis=0){
  # calzip is the california zipcodes SPDF
  # color.column is the column of the dataframe to be colorcoded
  # main is an optional plot title
  # precis sets the precision of the color interval boundaries
  
  colorval <- with(calzip@data,get(color.column))
  ## This section simply makes a vector of colors that correspond to each 
  ## row in the dataframe
  
  colorMap = brewer.pal(7,"Reds")
  col = mapColors(colorval,colorMap)
  idx = ceiling(colorval / max(colorval,na.rm=TRUE) * length(colorMap))
  quantiles <- classIntervals(colorval,n=7,style="quantile",dataPrecision=precis,na.rm=T)
  col <- findColours(quantiles, colorMap) # plot(quantiles, pal = colorMap)
  # make a special color for NAN entries
  col[is.na(colorval)] <- 'grey70'
  col[grepl(c('XX'),calzip$NAME)] <- 'grey80' # XX means "large land areas such as national parks"
  col[grepl(c('HH'),calzip$NAME)] <- 'blue'   # HH means body of water
  
  # this is the part that actually plots
  op = par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.1,1.1,1.1,1.1))
  plot(calzip, col = col, border = col)
  map("county",region="california",add=TRUE)
 
  legend("topright", fill = attr(col, "palette"),
        legend = names(attr(col, "table")),
        bty="n", cex = 0.8, y.intersp = 1)
  
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
calMap = function(df,plotCol,zipCol=NULL,calZipData=NULL,main=NULL) {
  if(is.null(calZipData)) {
    calZipData <- getCalSPDF()
  }
  if(is.null(zipCol)) {
    zipCol = grep('^zip',colnames(df),value=T)[1]
  }
  calZipData@data<-append.spdf(calZipData,'NAME',df,zipCol)
  calzipplot(calZipData,plotCol,main)
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

#calMap(db.getZipCounts(),'zip5','count',main='Meter count by zip code')

