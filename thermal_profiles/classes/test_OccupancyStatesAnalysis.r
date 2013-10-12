# test_OccupancyStatesAnalysis.r
#
# Test for OccupancyStatesAnalysis.
# 
# Adrian Albert
# Last modified: September 2013.
# -----------------------------------------------------------------------

rm(list = ls())
options(error = recover)

source('../utils/timing.r')
source('classes/OccupancyStatesAnalysis.r')

# populate OccupancyStatesAnalysis object
test      = new(Class='OccupancyStatesAnalysis', 
                path = 'D:/Dropbox/OccupancyStates/fits/fits_sel/')

test      = computeEffThermalResponse(test)
test      = computeEffTODResponse(test)

test      = thermalSegmentation(test, type = 'temperature')

test      = greedyThermalTargeting(test, method = 'temperature')
test      = greedyThermalTargeting(test, method = 'random')

# ______________
# Analysis plots

setwd('~/Dropbox/OccupancyStates/plots/plots_sel/')

# model size
png(filename = 'model_size.png', width = 800, height = 400, res = 150)
print(plot(test, type = 'model-size'))
dev.off()

# performance plot
png(filename = 'performance_population.png', width = 1000, height = 400, res = 130)
print(plot(test))
dev.off()

# entropy plot
png(filename = 'entropy_population.png', width = 800, height = 400, res = 100)
print(plot(test, type = 'entropy'))
dev.off()

# predictability plot
png(filename = 'predictability_population.png', width = 800, height = 400, res = 100)
print(plot(test, type = 'predictability'))
dev.off()

#   # typical states by components (hard clustering)
#   png(filename = 'states_comp_types_kmeans.png', width = 1800, height = 600, res = 100)
#   print(plot(test, type = 'states_comp_types_hard'))
#   dev.off()
#   
#   # typical states by components (soft clustering)
#   png(filename = 'states_comp_types_cmeans.png', width = 1800, height = 600, res = 100)
#   print(plot(test, type = 'states_comp_types_soft'))
#   dev.off()

# classification of response per zipcode
png(filename = 'response_zipcode.png', width = 1200, height = 600, res = 100)
print(plot(test, type = 'response'))
dev.off()

# heatmap of average response by temperature
png(filename = 'avg_response_temperature_heatmap.png', width = 1200, height = 600, res = 100)
print(plot(test, type = 'avg-response-temperature'))
dev.off()

# heatmap of average response by temperature
png(filename = 'thermal_duration_segmentation.png', width = 1200, height = 600, res = 100)
print(plot(test, type = 'thermal-duration-segmentation'))
dev.off()

png(filename = 'targeting-temperature-cooling-margin.png', width = 1200, height = 400, res = 100)
print(plot(test, type = 'greedy-targeting', focus = 'cooling'))
dev.off()

png(filename = 'targeting-temperature-heating-margin.png', width = 1200, height = 400, res = 100)
print(plot(test, type = 'greedy-targeting', focus = 'heating'))
dev.off()

png(filename = 'avg-response-query.png', width = 800, height = 400, res = 100)
print(plot(test, type = 'avg-response-query', query = c(3284167, 3675267)))
dev.off()


