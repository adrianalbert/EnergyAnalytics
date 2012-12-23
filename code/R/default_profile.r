# Profile file for PGE data analysis
# default options

# default options for this run
NOBS_THRESH = 180   # discard users with less than a 1/2 year at a given premise
N_PROC      = 8     # number of cores to run process on 
DATA_ACCESS = 'sql' # flag to indicate type of media used for data
PLOTS       = FALSE
VERBOSE     = TRUE
USER        = 'adalbert'
PASSWORD    = 'adrian'
ZIP_START   = 201
ZIP_STOP    = 400
PROFILE     = FALSE
LOG_OUTPUT  = TRUE # should the output messages be diverted to disk?
SAVE_TO_SQL = TRUE

# paths to load data/save computation
plots_path  = '~/Dropbox/ControlPatterns/plots2/'
save_path   = '~/EnergyAnalytics/fits/'
logs_path   = '~/EnergyAnalytics/logs/'
prof_file   = '~/EnergyAnalytics/Rprof.out'
zips_file   = '~/EnergyAnalytics/data/metadata/zipcode_res.csv'

