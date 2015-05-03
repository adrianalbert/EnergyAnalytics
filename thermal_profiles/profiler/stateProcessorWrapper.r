# stateProcessorWrapper.r
#
# Wrapper function to implement HMM decoding.
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

stateProcessorWrapper = function(cur_data, cur_covar, UID,  
                                 train.frac = 0.95, 
                                 tran.vars  = c('(Intercept)', 'TemperatureF'),
                                 resp.vars  = c('(Intercept)', 'TemperatureD'),                                 
                                 verbose    = F, dump_path = NULL, controls = NULL) {
  
  if (!is.null(dump_path)) if (!file.exists(dump_path))  dir.create(dump_path, recursive = TRUE)  

  timesteps = cur_covar$date
  cur_covar$date = NULL
  
  # ___________________________________
  # Construct DataFormatter object
  
  formatter = new(Class='DataFormatter', cur_data, UID)  
  formatter = addCovariates(formatter, timesteps, cur_covar)    
  good.data = extractFormattedData(formatter)    
    
  # ________________________________
  # Decode and interpret input data
    
  # initialize model
  decoder   = new(Class='StateDecoder', 
                  good.data$data, good.data$timestamps, good.data$UID,
                  train.frac = train.frac, 
                  tran.vars = tran.vars, 
                  resp.vars = resp.vars,
                  controls = controls)  
  # HMM analysis
  decoder   = learnStateDecoder(decoder, verbose = verbose)
  if (verbose) show(decoder)
  
  # perform interpretation and feature extraction
  interpreter  = new(Class = "Interpreter", decoder)
  
  # ______________________________
  # Save analysis results
  
  decode.dump  = dumpDecodedData(decoder, path = dump_path)
  interp.dump  = dumpInterpretedData(interpreter, path = dump_path)    
	
  return(list(decoder      = decoder, 
              interpreter  = interpreter, 
              decoded_data = decode.dump, 
              interp_data  = interp.dump))
}

