# agg_OccupancyStatesAnalysis.r
#
# Occupancy states analysis on zipcode aggregates.
# 
# Adrian Albert
# Last modified: June 2013.
# -----------------------------------------------------------------------

# ________________________________________
# Initializations....

rm(list = ls())

library(maps);
library(maptools);
library(zipcode);
library(UScensus2010)
library(data.table)
library(scales)

setwd('~/Dropbox/OccupancyStates/')
source('code/utils/timing.r')
source('code/OccupancyStatesAnalysis.r')

# ________________________________________
# Set up analysis 

# populate OccupancyStatesAnalysis object
test      = new(Class='OccupancyStatesAnalysis', path = '~/Dropbox/OccupancyStates/fits/aggregate/')
test      = computeEffThermalResponse(test)

# test      = greedyThermalTargeting(test, method = 'portfolio')

# ______________
# Analysis plots

# set current directory for plots
setwd('~/Dropbox/OccupancyStates/plots/aggregate/')

# model size
png(filename = 'model_size_zip.png', width = 800, height = 400, res = 150)
print(plot(test, type = 'model-size'))
dev.off()

# performance plot
png(filename = 'performance_population_zip.png', width = 1000, height = 400, res = 110)
print(plot(test))
dev.off()

# heatmap of average response by temperature
png(filename = 'avg_response_temperature_heatmap.png', width = 1200, height = 600, res = 100)
print(plot(test, type = 'avg-response-temperature'))
dev.off()

# heatmap of average response by temperature
png(filename = 'thermal_duration_segmentation.png', width = 1200, height = 600, res = 100)
print(plot(test, type = 'thermal-duration-segmentation'))
dev.off()

png(filename = 'avg-response-query.png', width = 800, height = 400, res = 100)
print(plot(test, type = 'avg-response-query', query = c(1, 445)))
dev.off()

# ________________________________________
# Scale plots: performance vs group size

# load up aggregated data
load('~/Dropbox/OccupancyStates/data/zip_aggregate_120k_1yr.RData')
noUsers = sapply(aggrZip, function(l) unique(l$Count))
df      = data.frame(UID = 1:length(noUsers), noUsers = noUsers)
df      = merge(df, subset(test@OCCUP_STATS, select = c('MAPE.cv', 'entropy.uncr', 'UID')))
names(df)[c(3,4)] = c('MAPE (CV)', 'Entropy (IID)')
bins    = c(1, 5, 10, 15, 20, 25, 30, 50, 75, 100, 200, 300, 400, 500, 1000, 3000)
bin.str = paste(bins[1:(length(bins)-1)], bins[2:length(bins)], sep = '-')
df$Size = bin.str[findInterval(df$noUsers, bins, all.inside = T)]
df$Size = factor(df$Size, levels = bin.str)
df$UID  = NULL
df$noUsers = NULL
df      = melt(df, id.vars = c('Size'))

p = ggplot(df, aes(x = Size, y = value, fill = variable)) + 
  geom_boxplot()
p = p + facet_wrap(~variable, scales = 'free', ncol = 2)
p = p + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18, angle = 330),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.position  = 'none',
        axis.ticks = element_blank()) + 
  ggtitle("Neighborhood Size and Reconstruction Performance") + xlab('Number of Users') + ylab('Performance Metric')

png(filename = 'scale_performance_aggregate.png', width = 1400, height = 600, res = 100)
print(p)
dev.off()

# ________________________________________
# Zipcode data for maps

# set current directory for plots
setwd('~/Dropbox/OccupancyStates/plots/aggregate/')

# load up zipcode-climate zone correspondence
zones   = read.csv('~/Dropbox/OccupancyStates/data/metadata/zip_metadata.csv')
zips    = zones[,c('zip5', 'cecclmzn')]
names(zips) = c('ZIPCODE', 'ZONE')

# load up zipcode-ZCTA correspondence
zip.zta = read.csv('~/Dropbox/OccupancyStates/data/metadata/Zip_to_ZCTA_Crosswalk_2011_JSI.csv')
zip.zta = merge(zips, zip.zta[,c('ZIP', 'ZCTA_USE')], by.x = 'ZIPCODE', by.y = 'ZIP')

# load ZCTA lat/lon info
zip.gaz = read.table('~/Dropbox/OccupancyStates/data/metadata/Gaz_zcta_national.txt', 
                     header=TRUE, colClasses=c("character", rep("numeric", 8)))
# merge zipcode information
zcta.zip= merge(zip.zta, zip.gaz[,c('GEOID', 'INTPTLAT', 'INTPTLONG')], by.x=c("ZCTA_USE"), by.y=c("GEOID"))

# add in data to zipcodes
df.resp = subset(test@EFF_RESPONSE, TemperatureF %in% c(40,60,80,100))
df.resp = merge(df.resp, zcta.zip)
df.resp = aggregate(data = df.resp, cbind(Thermal.Component, Var.Thermal) ~ ZCTA_USE + INTPTLAT + INTPTLONG + TemperatureF, FUN = mean)
df.resp$ZCTA_USE = as.character(df.resp$ZCTA_USE)
df.resp = split(df.resp, df.resp$TemperatureF)

# load zipcode boundary shapes
shapes  = readShapePoly("~/Dropbox/OccupancyStates/data/metadata/tl_2010_06_zcta510.shp")
#convert the zip from a factor to a character
shapes@data$ZCTA_CHAR<-levels(shapes@data$ZCTA5CE10)[shapes@data$ZCTA5CE10]

# clean up data
shapes@data = subset(shapes@data, select = c('ZCTA5CE10', 'ZCTA_CHAR'))

# convert map to data frame for plotting
mapDf = fortify(shapes)

for (t in names(df.resp)) {
  
  df = df.resp[[t]]
  
  # merge map and data
  mapDT    = data.table(mapDf)
  respDT   = data.table(df)
  mapDT$id = levels(shapes@data$ZCTA5CE10)[as.numeric(mapDT$id)+1]
  names(respDT)[1] = 'id'
  mapDT    = merge(mapDT, respDT, by="id")
  mapDT    = mapDT[order(mapDT$order),]
  
  # limit map to the CA region for which we have data
  # ca.limits <- geocode(c("Redding, CA", "Bakersfield, CA", "Healdsburg, CA", "Santa Maria, CA", "Chico, CA", 'Sonora, CA', 'Oakhurst, CA'))
  # MapDT <- subset(mapDT, long > min(ca.limits$lon) & long < max(ca.limits$lon) & lat > min(ca.limits$lat) & lat < max(ca.limits$lat))
  bbox = c(min(df$INTPTLONG) - 1, max(df$INTPTLONG) +1, min(df$INTPTLAT) - 1, max(df$INTPTLAT) + 1)
  mapDT_sel <- subset(mapDT, long >= bbox[1] & long <= bbox[2] & 
                     lat >= bbox[3] & lat <= bbox[4])
  
  # ggplot mapping
  # over a GoogleMap (not working if not correctly projected)
  # map <- get_map(location = c(mean(df$INTPTLONG), mean(df$INTPTLAT)), 
  #                zoom=7, source  = 'stamen', maptype = 'watercolor')
  map <- get_map(location = c(mean(df$INTPTLONG), mean(df$INTPTLAT)), 
                 zoom=7,  maptype = 'roadmap')
  m0 <- ggmap(map)
  m1 <- m0 + geom_polygon(aes(x=long, y=lat, group=group, fill=Thermal.Component), 
                          data=mapDT_sel, alpha=1)
  m2 <- m1 + geom_path(aes(x=long, y=lat, group=group), data=mapDT_sel, color='gray')
  
  # add custom color gradient
  m2 <- m2 + scale_fill_gradient2(low="red", high="blue", limits=c(-0.025, 0.06))
  
#   # add points
#   m2 <- m2 + geom_point(data = df, aes(INTPTLONG, INTPTLAT), alpha = 0.8, size = 2) 
#   
#   # add text
#   m3 <- m2 + geom_text(aes(x=INTPTLONG, y=INTPTLAT, label=ZCTA_USE), data=df, col="yellow", cex=3)
  
  # format plot
  m2 = m2 + theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),
          legend.title     = element_text('Thermal\nResponse',size=18),
          legend.position  = c(0.82, 0.85),
          axis.ticks = element_blank()) + 
    ggtitle(paste("Distribution of Thermal Regimes T =", t))
  if (t != '100') m2 = m2 + theme(legend.position = 'none')

  png(filename = paste(t, 'response_agg_zipcode.png', sep = '_'), width = 1000, height = 1000, res = 130)
  print(m2)
  dev.off()  
}

# png(filename = 'response_agg_zipcode.png', width = 2000, height = 600, res = 100)
# multiplot(plotlist = plot_list, cols = 4)
# dev.off()
