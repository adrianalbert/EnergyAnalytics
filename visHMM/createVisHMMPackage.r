# createvisHMMPackage.r
#
# Create R package "ggHMM"
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------
# 
rm(list = ls())
files = c('acf_ggplot.r', 'plot_utils.r')
package.skeleton(name = 'ggHMM', code_files = files)