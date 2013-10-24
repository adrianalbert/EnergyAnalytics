# sample10k_OccupancyStatesAnalysis.r
#
# Analysis on 10,000 users.
# 
# Adrian Albert
# Last modified: June 2013.
# -----------------------------------------------------------------------

rm(list = ls())

setwd('~/Dropbox/OccupancyStates/')
source('code/utils/timing.r')
source('code/OccupancyStatesAnalysis.r')

# load up zipcode-climate zone correspondence
zones   = read.csv('~/Dropbox/OccupancyStates/data/metadata/zip_metadata.csv')
zips    = zones[,c('zip5', 'cecclmzn')]
names(zips) = c('ZIPCODE', 'ZONE')

# populate OccupancyStatesAnalysis object
test      = new(Class='OccupancyStatesAnalysis', path = '~/Dropbox/OccupancyStates/fits/10k/')
test      = computeEffThermalResponse(test)
#test      = computeEffThermalDuration(test)
test      = computeEffTODResponse(test)

# ______________
# Analysis plots

setwd('~/Dropbox/OccupancyStates/plots/10k/')

# model size
png(filename = 'model_size.png', width = 800, height = 400, res = 150)
print(plot(test, type = 'model-size'))
dev.off()

# performance plot
png(filename = 'performance_population.png', width = 1000, height = 400, res = 100)
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

png(filename = 'avg-response-query_2users.png', width = 800, height = 400, res = 100)
print(plot(test, type = 'avg-response-query', query = c('2433956', '267098')))
dev.off()

png(filename = 'tod-response-query_2users.png', width = 1200, height = 400, res = 90)
print(plot(test, type = 'tod-response-query', query = c('2433956', '267098')))
dev.off()

png(filename = 'avg-duration-query.png', width = 800, height = 400, res = 100)
print(plot(test, type = 'avg-duration-query', query = c('2433956', '267098')))
dev.off()

png(filename = 'tod-response-distr.png', width = 800, height = 400, res = 100)
print(plot(test, type = 'tod-response-distr', query = zips))
dev.off()

