# Profile file for PGE data analysis
# default options

# default options for this run
NOBS_THRESH = 180   # discard users with less than a 1/2 year at a given premise
N_CORES     = 1     # number of cores to run process on 
DATA_ACCESS = 'sql' # flag to indicate type of media used for data
PLOTS       = FALSE
VERBOSE     = TRUE
USER        = 'adalbert'
PASSWORD    = 'adrian'
ZIP_START   = 2
ZIP_STOP    = 2
PROFILE     = F

# paths to load data/save computation
plots_path  = '~/Dropbox/ControlPatterns/plots/'
save_path   = '~/EnergyAnalytics/fits/'
kwh_path    = ''
wthr_path   = ''
prof_file   = '~/EnergyAnalytics/Rprof.out'
zips_file   = '~/EnergyAnalytics/data/metadata/zipcode_res.csv'

