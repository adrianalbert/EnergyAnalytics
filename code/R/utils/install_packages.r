# install_packages.r
#
# Install required packages for analysis on PGE data.
# 
# Adrian Albert
# Last modified: December 2012.

.libPaths('~/R/library') # use my local R library even from the comand line

# see here for RMySQL source install: http://stackoverflow.com/questions/4785933/adding-rmysql-package-to-r-fails
# define required pacakges
reqed.packages = c('RMySQL','zoo', 'depmixS4', 'ggplot2', 'utils', 'multicore', 'parallel',
                   'methods', 'lmtest', 'Amelia', 'imputation', 'timeDate', 'lubridate', 
                   'dummies', 'linkcomm', 'useful', 'reshape', 'grid', 'profr', 'R.utils', 
                   'proftools', 'data.table', 'igraph', 'useful', 'RColorBrewer', 'reshape2',
                   'solaR','mapproj','ggmap','gpclib','UScensus2010','colorspace','gridExtra',
                   'sp','gptk','gtools','hexbin','classInt','rgeos','cvTools')

# get already installed packages
avail.packages = installed.packages()

# what packages do we need?
need.packages  = setdiff(reqed.packages, avail.packages)

# install packages that are missing
install.packages(need.packages)

#library(UScensus2010)
#install.tract("windows") # or linux or osx

