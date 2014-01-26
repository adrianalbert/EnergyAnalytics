# install_packages.r
#
# Install required packages for analysis on PGE data.
# 
# Adrian Albert
# Last modified: October 2013.

# define required pacakges
reqed.packages = c('RMySQL', 'zoo', 'depmixS4', 'ggplot2', 'utils', 
                   'multicore', 'parallel', 'segmented', 
                   'methods', 'lmtest', 'Amelia', 'imputation', 
                   'timeDate', 'lubridate', 
                   'dummies', 'linkcomm', 'useful', 'reshape', 
                   'grid', 'profr', 'R.utils', 
                   'proftools', 'data.table', 'igraph', 'useful', 
                   'RColorBrewer', 'reshape2', 'Rcpp', 'inline',
                   'lubridate', 'ggsubplot', 'ggmap', 'pls', 'zipcode')

# get already installed packages
avail.packages = installed.packages()

# what packages do we need?
need.packages  = setdiff(reqed.packages, avail.packages)

# install packages that are missing
install.packages(need.packages)
