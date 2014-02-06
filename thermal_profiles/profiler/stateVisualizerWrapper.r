# stateVisualizerWrapper.r
#
# Wrapper function to implement visualizations for HMM decoding.
#
# Adrian Albert
# Last modified: December 2013.

# ---------------------------------------------------------
# Wrapper to perform analysis on a given time series data.
# ---------------------------------------------------------

options(error = recover)

source('classes/DataFormatter.r')
source('classes/StateDecoder.r')
source('classes/Interpreter.r')
source('classes/Visualizer.r')

stateVisualizerWrapper = function(decoder, interpreter, 
                                  plots_path = NULL, interval = NULL) {
  
  if (!is.null(plots_path)) if (!file.exists(plots_path))  dir.create(plots_path, recursive = TRUE)    
  
  # _______________________________
  # Produce analysis plots
  
  if (!is.null(plots_path)) {
        
    # visualize decoded & interpreted data
		visualizer     = new(Class = "Visualizer", decoder, interpreter,
		                    interval = interval, 
                         t.covar = 'TemperatureF',
                         r.covar = 'TemperatureD')
		visualizer.all = new(Class = "Visualizer", decoder, interpreter,
		                    interval = NULL, 
		                     t.covar = 'TemperatureF',
		                     r.covar = 'TemperatureD')

    # place plots in a directory
		dir.create(file.path(plots_path, visualizer@UID))
    plots_path = paste(plots_path, '/', visualizer@UID, '/', sep = '')
    
    # just the data
		png(paste(plots_path, visualizer@UID, '_raw.png', sep=''), width=1400, height=600, res = 150)
		plot(visualizer)
		dev.off()

		# HMM fit
    p1 = plot(visualizer, type='HMM-ts')
		png(paste(plots_path, visualizer@UID, '_HMM_fit.png', sep=''), width=1400, height=600, res = 150)
		print(p1)
		dev.off()

		# heatmap plots (entire data)
    png(paste(plots_path, visualizer.all@UID, '_HMM-consumption-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
		plot(visualizer.all, type = 'HMM-heatmap-obs')
		dev.off()

		# state breakdown heatmap
		png(paste(plots_path, visualizer.all@UID, '_HMM-state-heatmap2.png', sep=''), width = 2000, height = 1000, res = 200)
		print(plot(visualizer.all, type = 'HMM-state-heatmap-states'))
		dev.off()      

		# state breakdown by time of day
		png(paste(plots_path, visualizer.all@UID,'_HMM-state-breakdown.png', sep=''), width = 1400, height = 800, res = 150)
		print(plot(visualizer.all, type = 'HMM-state-breakdown'))
		dev.off()

		# HMM PACF
		png(paste(plots_path, visualizer.all@UID, '_HMM_pacf.png', sep=''), width=1000, height=600, res = 100)
		print(plot(visualizer.all, type='HMM-pacf'))
		dev.off()

		# HMM residuals
		png(paste(plots_path, visualizer.all@UID, '_HMM_res.png', sep=''), width=1200, height=500, res = 100)
		plot(visualizer.all, type = 'HMM-res')
		dev.off()

		# coefficient time series
		png(paste(plots_path, visualizer@UID, '_HMM-coefs-ts.png', sep=''), width = 1400, height = 800, res = 100)
		print(plot(visualizer, type = 'HMM-coefs-ts'))
		dev.off()      

		# temperature dependence
    png(paste(plots_path, visualizer.all@UID, '_HMM-dep-covar.png', sep = ''), width = 1000, height = 600, res = 150)
		print(plot(visualizer.all, type = 'HMM-dep-covar'))
		dev.off()      

		# separate states into their own panels
    png(paste(plots_path, visualizer.all@UID, 'HMM-dep-covar-sep.png', sep=''), width = 1000, height = 600, res = 150)
		print(plot(visualizer.all, type = 'HMM-dep-covar-sep'))
		dev.off()

		# state-dependent transition profiles (vs temperature)
    png(paste(plots_path, visualizer.all@UID,'HMM-trans-prob.png',sep=''), width = 1600, height = 800, res = 150)
		print(plot(visualizer.all, type = 'HMM-trans-prob'))
		dev.off()

    png(paste(plots_path, visualizer.all@UID,'HMM-stationary-prob.png',sep=''), width = 900, height = 500, res = 150)
		print(plot(visualizer.all, type = 'HMM-stationary-prob'))
		dev.off()

		# error forecasts
#     png(paste(plots_path, visualizer.all@UID,'HMM-err-horiz.png',sep=''), width = 1000, height = 400, res = 150)
# 		print(plot(visualizer.all, type = 'HMM-err-horiz'))
# 		dev.off()

		# aggregate seasonal contributions
    png(paste(plots_path, visualizer.all@UID,'HMM-aggregate-season.png',sep=''), width = 1600, height = 600, res = 150)
		print(plot(visualizer.all, type = 'HMM-aggregate-season'))
		dev.off()

		# HMM covariate contributions 
    png(paste(plots_path, visualizer@UID, '_HMM_contrib_ts.png', sep=''), width = 1400, height = 800, res = 150)
		plot(visualizer, type = 'HMM-contrib-ts')
		dev.off()
  }
  
  return(list(vis.all = visualizer.all, vis.int = visualizer))
}

