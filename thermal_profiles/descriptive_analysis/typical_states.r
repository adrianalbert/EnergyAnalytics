#!/usr/bin/Rscript

# typical_states.r
# 
# Cluster states and transition matrices
#
# - 
# Adrian Albert
# Last modified: March 2013.

# -----------------------------------------------
# Define constants, load libraries and functions
# -----------------------------------------------

rm(list=ls())

options(error = recover)

# load classes, utils
setwd('~/Dropbox/OccupancyStates/')
source('code/utils/timing.r')
source('code/utils/acf_ggplot.r')
source('code/Person.r')
source('code/personAnalysis.r')

# -----------------------------------------------
# 
# -----------------------------------------------

