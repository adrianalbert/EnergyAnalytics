# #########################################################################
# OccupancyStateAnalysis.r
# -------------------------
#
# OccupancyStatesAnalysis.r
#
# Class to encapsulate post-estimation analysis on occupancy states. 
# 
# Adrian Albert
# Last modified: September 2013.
# #########################################################################

library('methods')
library('timeDate')
library('zoo')
library('lubridate')
require('lmtest')
library('pls')
library('reshape')
library('ggplot2')
library('ggmap')
library('zipcode')
data('zipcode')

source('../utils/viz/plot_utils.r')
source('./clustering/kError.r')
source('./clustering/d_err.r')

# clean-up previous definitions of methods for class Person
removeClass('OccupancyStatesAnalysis')

# ________________________
# Class definition

setClass(
  Class = "OccupancyStatesAnalysis",
  representation = representation(
    PATH_DATA        = "character",           # path to where occupancy data is stored
    STATES_ATTR      = "data.frame",          # state attributes
    BENCHMARKS       = "data.frame",          # benchmarks 
    STATES_COMP      = "data.frame",          # variance components by state
    STATES_SEAS      = "data.frame",          # seasonal breakdown of states
    OCCUP_STATS      = "data.frame",          # stats on occupancy patterns
    SEGMENTATION     = 'list',                # clustering results for users
    TARGETING        = 'list',                # simple targeting results
    EFF_RESPONSE     = 'data.frame',          # effective thermal response vs temperature
    TOD_RESPONSE     = 'data.frame',          # effective thermal response vs temperature
    COVMAT           = 'list'                 # thermal covariance matrix of users
  )
)

# function to read in data from files
read_files = function(path, pattern) {
  cur_files  = list.files(path = path, pattern = pattern, full.names = T, recursive = T)
  res        = data.frame()
  for (f in cur_files) {
    if (file.info(f)$size < 10) next
    tmp = read.csv(f)
    if (nrow(res) == 0) res = tmp else res = rbind(res, tmp)
  }
  res <- res[,colSums(is.na(res))<nrow(res)]
  res$X = NULL
  
  return(list(data = res, noChunks = length(cur_files)))
}

# _____________________________________________________
# Constructor method for class OccupancyStatesAnalysis

setMethod(f = "initialize", 
          signature = "OccupancyStatesAnalysis",
          definition = function(.Object, 
                                path = './', 
                                states_attr = 'response',
                                states_comp = 'hmmcomps',
                                states_stat = 'hmmstats',
                                benchmark   = 'transitn',
                                states_seas = 'hmmseasn',
                                verbose = T) {
            
            if (verbose) cat('*** Reading in data ***\n')
            
            .Object@PATH_DATA = path
                                  
            # ____________________________
            # Format state attribute data
            
            res = read_files(path, states_attr)
            .Object@STATES_ATTR = res[['data']]
            if (verbose) cat(paste('State attributes:\tRead in', res[['noChunks']], 'chunks\n'))
            
            # ____________________________
            # Format state component data
            
            res = read_files(path, states_comp)
            .Object@STATES_COMP = res[['data']]
            if (verbose) cat(paste('State components:\tRead in', res[['noChunks']], 'chunks\n'))
            
            # ____________________________
            # Format HMM stats
            
            res = read_files(path, states_stat)
            .Object@OCCUP_STATS = res[['data']]
            if (verbose) cat(paste('Occupancy stats:\tRead in', res[['noChunks']], 'chunks\n'))
            
            # ____________________________
            # Format benchmark data
            
            res = read_files(path, benchmark)
            .Object@BENCHMARKS = res[['data']]
            if (verbose) cat(paste('Benchmarks:\tRead in', res[['noChunks']], 'chunks\n'))              
            # ____________________________
            # Format benchmark data
            
            res = read_files(path, states_seas)
            .Object@STATES_SEAS = res[['data']]
            if (verbose) cat(paste('Seasonal breakdown:\tRead in', res[['noChunks']], 'chunks\n'))              

            return(.Object)
          })

# ______________________________________________________________
# Compute effective thermal response for users (by temperature)

setGeneric(
  name = "computeEffThermalResponse",
  def = function(.Object, verbose = T){standardGeneric("computeEffThermalResponse")}
)
setMethod('computeEffThermalResponse',
          signature  = 'OccupancyStatesAnalysis',
          definition = function(.Object, verbose = T){
            
            if (verbose) cat('Computing effective thermal response...')
            
            # select data
            df.distr = subset(.Object@BENCHMARKS, variable == 'Distribution')
            df.resp  = .Object@STATES_ATTR[, c('UID', 'ZIPCODE', 'State', 'TemperatureF', 'Sigma', 'X.Intercept.')]
            
            # compute average response
            df.resp          = merge(df.distr, df.resp)
            df.resp$variable = NULL
            df.resp          = melt(df.resp, id.vars = c('UID', 'ZIPCODE', 'State', 'Sigma', 'TemperatureF', 'X.Intercept.'))
            df.resp$variable = as.numeric(gsub('X', '', df.resp$variable))
            df.resp$AvgTempF = df.resp$TemperatureF * df.resp$value
            df.resp$VarTempF = df.resp$TemperatureF^2 * df.resp$value
            df.resp$value    = NULL
            df.resp$TemperatureF = NULL
            names(df.resp)[c(5:8)] = c('Activity.Component', 'TemperatureF', 'Thermal.Component', 'Var.Thermal')
            
            # compute aggregated response              
            df.resp  = aggregate(data = df.resp, cbind(Thermal.Component, Var.Thermal) ~ UID + ZIPCODE + TemperatureF, FUN = sum)              
            df.resp$Var.Thermal = sqrt(df.resp$Var.Thermal - (df.resp$Thermal.Component)^2)
            df.resp$UID = as.factor(df.resp$UID)
            
            # store computation
            .Object@EFF_RESPONSE = df.resp
            
            return(.Object)
            
          }
)

# ______________________________________________________________________
# Compute effective thermal response for users by time of day & season

setGeneric(
  name = "computeEffTODResponse",
  def = function(.Object, verbose = T){standardGeneric("computeEffTODResponse")}
)
setMethod('computeEffTODResponse',
          signature  = 'OccupancyStatesAnalysis',
          definition = function(.Object, verbose = T){
            
            if (verbose) cat('Computing effective time-of-day response...')
            
            # select data
            df.distr = .Object@STATES_SEAS
            df.resp  = .Object@STATES_ATTR[, c('UID', 'ZIPCODE', 'State', 'TemperatureF', 'Sigma', 'X.Intercept.')]
            
            # compute average response
            df.resp          = merge(df.distr, df.resp)
            df.resp          = melt(df.resp, id.vars = c('UID', 'ZIPCODE', 'Season', 'State', 'Sigma', 'TemperatureF', 'X.Intercept.'))
            df.resp$variable = as.numeric(gsub('X', '', df.resp$variable))
            df.resp$AvgTempF = df.resp$TemperatureF * df.resp$value
            df.resp$VarTempF = df.resp$TemperatureF^2 * df.resp$value
            df.resp$value    = NULL
            df.resp$TemperatureF = NULL
            names(df.resp)[c(6:9)] = c('Activity.Component', 'Hour', 'Thermal.Component', 'Var.Thermal')
            
            # compute aggregated response              
            df.resp  = aggregate(data = df.resp, cbind(Thermal.Component, Var.Thermal) ~ UID + ZIPCODE + Hour + Season, FUN = sum)              
            df.resp$Var.Thermal = sqrt(df.resp$Var.Thermal - (df.resp$Thermal.Component)^2)
            df.resp$UID = as.factor(df.resp$UID)
            
            # store computation
            .Object@TOD_RESPONSE = df.resp
            
            return(.Object)
            
          }
)

# ______________________________
# Thermal segmentation for users

setGeneric(
  name = "thermalSegmentation",
  def = function(.Object, Kmin = 3, Kmax = 3, verbose = T, type = 'temperature')
  {standardGeneric("thermalSegmentation")}
)
setMethod('thermalSegmentation',
          signature  = 'OccupancyStatesAnalysis',
          definition = function(.Object, Kmin = 3, Kmax = 3, verbose = T, type = 'temperature'){
            
            if (verbose) cat('Thermal Segmentation:\n')
            
            # prepare data
            K = Kmin
            if (type == 'temperature') {
              # observations
              df = subset(.Object@EFF_RESPONSE, select = c('UID', 'TemperatureF', 'Thermal.Component'))              
              df = df[with(df, order(UID, TemperatureF)), ]
              uid= df[,'UID']
              df = matrix(as.matrix(df[,-1]), ncol = length(unique(df[,'TemperatureF'])))
              df = as.data.frame(df)
              X  = cbind(unique(uid), df)
              # errors
              df = subset(.Object@EFF_RESPONSE, select = c('UID', 'TemperatureF', 'Var.Thermal'))              
              df = df[with(df, order(UID, TemperatureF)), ]
              uid= df[,'UID']
              df = matrix(as.matrix(df[,-1]), ncol = length(unique(df[,'TemperatureF'])))
              df = as.data.frame(df)
              S  = cbind(unique(uid), df)              
            } else {
              X = .Object@TOD_RESPONSE
              # TODO: per-season analysis? or put both seasons together?
            }
            
            res = kError(X, S, K, iter = 10)
            return(res)
            
          }
)

# _____________________________________________
# Compute thermal correlations for all users

setGeneric(
  name = "computeThermalCorrelations",
  def = function(.Object, verbose = T){standardGeneric("computeThermalCorrelations")}
)
setMethod('computeThermalCorrelations',
          signature  = 'OccupancyStatesAnalysis',
          definition = function(.Object, verbose = T){
            
            if (verbose) cat('Computing effective thermal response...')
            
            # TODO!!
            
            # store computation
            .Object@COVMAT = df.resp
            
            return(.Object)
            
          }
)

# ________________________________________
# Simple targeting

setGeneric(
  name = "greedyThermalTargeting",
  def = function(.Object, method = 'temperature', verbose = T){standardGeneric("greedyThermalTargeting")}
)
setMethod('greedyThermalTargeting',
          signature  = 'OccupancyStatesAnalysis',
          definition = function(.Object, method = 'temperature', verbose = T){
              
            if (method == 'temperature') {
              
              if (verbose) cat('Greedy thermal targeting:\n')
              
              df.resp = subset(.Object@EFF_RESPONSE, TemperatureF %in% c(45, 60, 75, 95))              
              
              # classify response into low/medium/high
              df.resp$Type  = 'cooling'
              df.resp$Type[df.resp$Thermal.Component < 0] = 'heating'
              qnt      = quantile(abs(df.resp$Thermal.Component), probs = c(0,0.25,0.75, 1))
              type     = c('none', 'low', 'high')[findInterval(abs(df.resp$Thermal.Component), qnt, all.inside = T)]
              df.resp$Type[type == 'none'] = ''
              df.resp$Type  = paste(type, df.resp$Type, sep='-')     
              
              # compute cumulative sum of averted consumption
              df.resp$TemperatureF = as.factor(df.resp$TemperatureF)
              avert.c = lapply(levels(df.resp$TemperatureF), function(l) {
                cur = subset(df.resp, (TemperatureF == l) & (Thermal.Component > 0) & Type %in% c('high-cooling', 'low-cooling'))
                cur = cur[with(cur, order(-abs(Thermal.Component))), ]
                # cumsum(abs(cur$Thermal.Component))
                abs(cur$Thermal.Component)
              })
              names(avert.c) = levels(df.resp$TemperatureF)                            
              avert.h = lapply(levels(df.resp$TemperatureF), function(l) {
                cur = subset(df.resp, (TemperatureF == l) & (Thermal.Component < 0) & Type %in% c('high-heating', 'low-heating'))
                cur = cur[with(cur, order(-abs(Thermal.Component))), ]
                # cumsum(abs(cur$Thermal.Component))
                abs(cur$Thermal.Component)
              })
              names(avert.h) = levels(df.resp$TemperatureF)    
              
              .Object@TARGETING[['avert.temp.cool']]  = avert.c
              .Object@TARGETING[['avert.temp.heat']]  = avert.h            
            }                      

          if (method == 'random') {

            if (verbose) cat('Random thermal targeting:\n')
            
            df.resp = subset(.Object@EFF_RESPONSE, TemperatureF %in% c(45, 60, 75, 95))              
            
            # compute random performance by temperature level
            avert = lapply(unique(df.resp$TemperatureF), function(l) {    
              Kvec = seq(5,length(unique(df.resp$UID)), by=5)
              avert.r = lapply(Kvec, function(K) {                  
                  mu = c()             
                  for (iter in 1:20) {
                    # select users at random
                    uid.sel = sample(unique(df.resp$UID), K)                    
                    cur     = subset(df.resp, TemperatureF == l & UID %in% uid.sel)
                    mu[iter]= mean(cur$Thermal.Component)
                    # mu[iter]= sum(cur$Thermal.Component)
                  }
                  std = sd(mu)
                  mu  = mean(mu)
                  return(c(mu = mu, sd = std))
                })              
              
              df = as.data.frame(do.call('rbind', avert.r))
              df$No.Users = Kvec
              df$TemperatureF = l
              return(df)
            })
            avert = do.call('rbind', avert)  
            
            .Object@TARGETING[['random']] = avert
          }                           
          return(.Object)
        
        }          
)

# ____________________________
# Plot OccupancyStates object

setMethod('plot',
          signature  = 'OccupancyStatesAnalysis',
          definition = function(x, type = 'default', focus = 'cooling', query = NULL, center = NULL, verbose = T, ...){
            
            if (verbose) cat('*** Plot for OccupancyStates ***\n')
            
            if (type == 'default' | type == 'performance') {
              
              df.1 = subset(x@OCCUP_STATS, select = c('R2.cv', 'MAPE.cv'))
              names(df.1) = c('R2', 'MAPE')
              df.1 = melt(df.1)
              df.1$Type = 'Out.Of.Sample'
              df.2 = subset(x@OCCUP_STATS, select = c('R2', 'MAPE'))
              df.2 = melt(df.2)
              df.2$Type = 'In.Sample'
              df = rbind(df.1, df.2)
              plt = ggplot(df) + 
                geom_density(aes(x = value, colour = Type), size=2)
              plt = plt + facet_wrap( ~ variable, nrow=1, scales = 'free') + 
                theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18, angle = 330),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),   
                      legend.position  = c(0.8, 0.5),
                      axis.ticks = element_blank()) + 
                ylab('pdf') + ggtitle("Occupancy States: In-Sample Performance")       
              
              return(plt)
            }
            
            if (type == 'model-size') {
              tb = table(x@STATES_ATTR$UID)
              df = data.frame(UID = names(tb), No.States = as.factor(as.numeric(tb)))
              plt = ggplot(df) + 
                geom_bar(aes(x = No.States), size=2)
              plt = plt + 
                theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      axis.ticks = element_blank()) + 
                ylab('pdf') + ggtitle("Number of Occupancy States")       
              
              return(plt)
              
            }
            if (type == 'entropy') {
              df = subset(x@OCCUP_STATS, 
                          select = c('ZIPCODE', 'entropy.rand', 'entropy.uncr', 'entropy.lziv'))
              names(df) = c('ZIPCODE', 'Random', 'Uncorrelated', 'Lempel-Ziv')
              df = melt(df, id.vars = 'ZIPCODE')
              plt = ggplot(df) + 
                geom_density(aes(x = value, colour = variable), size=2)
              # plt = plt + facet_wrap(~ZIPCODE, ncol = 2)
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=15), 
                      axis.text.x      = element_text(size=15),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      axis.ticks = element_blank()) + 
                ylab('pdf') + ggtitle("Entropy Benchmarks")       
              
              return(plt)
              
            }
            
            if (type == 'predictability') {
              df = subset(x@OCCUP_STATS, 
                          select = c('ZIPCODE', 'predict.rand', 'predict.uncr', 'predict.gzip'))
              names(df) = c('ZIPCODE', 'Random', 'Uncorrelated', 'Lempel-Ziv')              
              df = melt(df, id.vars = 'ZIPCODE')
              plt = ggplot(df) + 
                geom_density(aes(x = value, colour = variable), size=2) + 
              # plt = plt + facet_wrap(~ZIPCODE, ncol = 2) + 
                theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=15), 
                      axis.text.x      = element_text(size=15),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      axis.ticks = element_blank()) + 
                ylab('pdf') + ggtitle("Predictability Bound Benchmarks")       
              
              return(plt)            
            }
            
            if (type == 'states_comp_types_hard' | type == 'states_comp_types_soft') {
              if(type == 'states_comp_types_hard') { 
                df = as.data.frame(x@COMPS_KMEAN$centers)
              } else 
                df = as.data.frame(x@COMPS_CMEAN$centers)          
              df$State = 1:nrow(df)
              df = melt(df, id.vars = 'State')
              plt = ggplot(df) + 
                geom_bar(aes(x = variable, y = value, fill = variable), stat = 'identity')
              plt = plt + facet_wrap(~State, ncol = 5) + 
                theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.text.y      = element_text(size=15), 
                      strip.text.x     = element_text(size=18),
                      axis.text.x      = element_text(size=15, angle=30),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.position  = "none",
                      axis.ticks = element_blank()) + 
                ylab('% Variance Explained') + ggtitle("Typical States (90% coverage)")       
              
              return(plt)              
            }
            
            # distribution of temperature response per zipcode
            if (type == 'response') {
              df = subset(x@STATES_ATTR, 
                          select = c('ZIPCODE', 'TemperatureF'))
              df$Type  = 'cooling'
              df$Type[df$TemperatureF < 0] = 'heating'
              qnt      = quantile(abs(df$TemperatureF), probs = c(0,0.25,0.75, 1))
              type     = c('none', 'low', 'high')[findInterval(abs(df$TemperatureF), qnt, all.inside = T)]
              df$Type[type == 'none'] = ''
              df$Type  = paste(type, df$Type, sep='-')
              plt = ggplot(subset(df, Type != 'none-')) + 
                geom_density(aes(x = TemperatureF, fill = Type))
              # plt = plt + facet_wrap(~ZIPCODE, ncol = 2)
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=15), 
                      axis.text.x      = element_text(size=15),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      axis.ticks = element_blank()) + 
                ylab('pdf') + ggtitle("Distribution of Thermal Regimes") + xlab('Thermal Response Rate')
              
              return(plt)            
            }    
            
            # density heatmap of average response by temperature
            if (type == 'avg-response-temperature') {

              # select data
              df.distr = subset(x@BENCHMARKS, variable == 'Distribution')
              df.resp  = x@STATES_ATTR[, c('UID', 'ZIPCODE', 'State', 'TemperatureF', 'Sigma', 'X.Intercept.')]
              temp.sel = paste('X', c(45, 60, 75, 95), sep = '')
              df.distr = df.distr[, c('UID', 'ZIPCODE', 'State', temp.sel)]
              
              # compute average response
              df.resp          = merge(df.distr, df.resp)
              df.resp$variable = NULL
              df.resp          = melt(df.resp, id.vars = c('UID', 'ZIPCODE', 'State', 'Sigma', 'TemperatureF', 'X.Intercept.'))
              df.resp$variable = as.numeric(gsub('X', '', df.resp$variable))
              df.resp$AvgTempF = df.resp$TemperatureF * df.resp$value
              df.resp$AvgSigma = ((df.resp$TemperatureF * df.resp$variable + df.resp$X.Intercept.)^2 + df.resp$Sigma^2) * df.resp$value
              df.resp$AvgFit   = (df.resp$TemperatureF * df.resp$variable + df.resp$X.Intercept.) * df.resp$value
              df.resp$value    = NULL
              df.resp$TemperatureF = NULL
              names(df.resp)[c(5:7)] = c('Activity.Component', 'TemperatureF', 'Thermal.Component')
              
              # compute variance of aggregated response              
              df.resp  = aggregate(data = df.resp, cbind(Thermal.Component, AvgSigma, AvgFit) ~ UID + ZIPCODE + TemperatureF, FUN = sum)              
              df.resp$AvgSigma = sqrt(df.resp$AvgSigma - (df.resp$AvgFit)^2)
            
              # construct plot
              plt = ggplot(df.resp, aes(x = Thermal.Component, y = AvgSigma)) + 
                stat_density2d(aes(fill = ..level..), geom="polygon")
              plt = plt + facet_wrap(~TemperatureF, ncol = 2)
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      legend.title     = element_text(size=18),
                      axis.ticks = element_blank()) + 
                ggtitle("Distribution of Thermal Regimes") + xlab('Effective Thermal Response Rate [kWh/F]') + ylab('Effective Standard Error [kWh]')
              
              return(plt)            
          }                        
          
          # segmentation by thermal and activity responses
          if (type == 'thermal-duration-segmentation') {
            
            # prepare data 
            df.distr = subset(x@BENCHMARKS, variable == 'Distribution')
            df.durat = subset(x@BENCHMARKS, variable == 'Duration')
            df.resp  = x@STATES_ATTR[, c('UID', 'ZIPCODE', 'State', 'TemperatureF')]
            temp.sel = paste('X', c(45, 60, 75, 95), sep = '')
            df.distr = df.distr[, c('UID', 'ZIPCODE', 'State', temp.sel)]
            df.durat = df.durat[, c('UID', 'ZIPCODE', 'State', temp.sel)]
            df.durat[,temp.sel] = apply(df.durat[,temp.sel], 2, function(x) {x[which(x>24)] = 24; x[which(x<1)] = 1; return(x)})
            
            # compute average response
            df.resp  = merge(df.distr, df.resp)
            df.resp$variable = NULL
            df.resp  = melt(df.resp, id.vars = c('UID', 'ZIPCODE', 'State', 'TemperatureF'))
            df.resp$TemperatureF = df.resp$TemperatureF * df.resp$value
            df.resp$value = NULL
            names(df.resp)[c(4,5)] = c('Thermal.Component', 'TemperatureF')
            df.resp$TemperatureF = gsub('X', 'T', df.resp$TemperatureF)
            
            # compute average duration
            df.durat  = cbind(df.distr[,1:3], df.distr[,temp.sel] * df.durat[,temp.sel])
            df.durat  = melt(df.durat, id.vars = c('UID', 'ZIPCODE', 'State'))
            names(df.durat)[c(4,5)] = c('TemperatureF','Duration')
            df.durat$TemperatureF = gsub('X', 'T', df.durat$TemperatureF)
            
            # compute user averages
            df = merge(df.durat, df.resp)
            df$UID = as.factor(df$UID)
            df$ZIPCODE = as.factor(df$ZIPCODE)
            df$TemperatureF = as.factor(df$TemperatureF)
            df$State = as.factor(df$State)
            df = aggregate(data = df, cbind(Thermal.Component, Duration) ~ UID + ZIPCODE + TemperatureF, FUN = sum)
            
            # classify states into low/medium/high by thermal response
            seg = subset(df, select = c('Thermal.Component', 'Duration'))
            seg$Thermal.Type  = 'cooling'
            seg$Thermal.Type[seg$Thermal.Component < 0] = 'heating'
            qnt.resp = quantile(abs(seg$Thermal.Component), probs = c(0,0.25,0.75, 1))
            type     = c('none', 'low', 'high')[findInterval(abs(seg$Thermal.Component), qnt.resp, all.inside = T)]
            seg$Thermal.Type[type == 'none'] = ''
            seg$Thermal.Type  = paste(type, seg$Thermal.Type, sep='-')  

            # segmentation by duration
            qnt.time = quantile(abs(seg$Duration), probs = c(0,0.25,0.75, 1))
            type     = c('short', 'medium', 'long')[findInterval(abs(seg$Duration), qnt.time, all.inside = T)]
            seg$Duration.Type = type            
            df = cbind(df, seg[,c('Thermal.Type', 'Duration.Type')])
            
            # construct plot
            plt = ggplot(df, aes(x = Thermal.Component, y = Duration, color = Thermal.Type, shape = Thermal.Type)) + 
              #stat_density2d(aes(fill = ..level..))
              geom_point(size=2)
            plt = plt + facet_wrap(~TemperatureF, ncol = 2)
            plt = plt + geom_vline(xintercept = qnt.resp[c(2,3)], color = 'blue', size = 1.5)
            plt = plt + geom_hline(yintercept = qnt.time[c(2,3)], color = 'red', size = 1.5)
            plt = plt + theme_bw() + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    strip.text.x     = element_text(size=18),
                    axis.text.y      = element_text(size=18), 
                    axis.text.x      = element_text(size=18),
                    axis.title.y     = element_text(size=18),
                    axis.title.x     = element_text(size=18),
                    plot.title       = element_text(size=20),            
                    legend.text      = element_text(size=18),
                    axis.ticks = element_blank()) + 
              ylab('Avg. Duration [hrs]') + ggtitle("Segmentation by Thermal Response and Duration") + xlab('Avg. Thermal Response [kWh/F]')
            
            return(plt)            
          } 
            
          # density heatmap of average response by temperature
          if (type == 'greedy-targeting') {

            # thermal targeting: cooling scenario
            df.c = list()
            for (v in names(x@TARGETING$avert.temp.cool)) {
              cur.c = x@TARGETING$avert.temp.cool[[v]]
              df.c[[v]] = data.frame(No.Users = 1:length(cur.c), Averted.kWh = cur.c, TemperatureF = v, Type = 'Cooling')
            }
            df.c = do.call('rbind', df.c)
            df.c$TemperatureF = as.numeric(gsub('T', '', df.c$TemperatureF))
            df.c$sd = 0
            
            # thermal targeting: heating scenario
            df.h = list()
            for (v in names(x@TARGETING$avert.temp.heat)) {
              cur.h = x@TARGETING$avert.temp.heat[[v]]
              df.h[[v]] = data.frame(No.Users = 1:length(cur.h), Averted.kWh = cur.h, TemperatureF = v, Type = 'Heating')
            }
            df.h = do.call('rbind', df.h)
            df.h$TemperatureF = as.numeric(gsub('T', '', df.h$TemperatureF))
            df.h$sd = 0
            
            # random targeting
            df2 = x@TARGETING$random
            df2$Type = 'Random'
            names(df2)[1] = 'Averted.kWh'
            df2 = df2[,names(df.h)]
            
            # put together dataframes
            df.c = rbind(df.c, subset(df2, No.Users <= max(df.c$No.Users)))
            df.c$TemperatureF = as.factor(paste('T =',df.c$TemperatureF))
#            df.c$TemperatureF = factor(df.c$TemperatureF, levels = levels(df.c$TemperatureF)[c(2,3,4,1)])            
            df.h = rbind(df.h, subset(df2, No.Users <= max(df.h$No.Users)))
            df.h$TemperatureF = as.factor(paste('T =',df.h$TemperatureF))
#            df.h$TemperatureF = factor(df.h$TemperatureF, levels = levels(df.h$TemperatureF)[c(2,3,4,1)])            
            
            # construct plot
            if (focus == 'cooling') mydf = df.c else mydf = df.h            
            plt1 = ggplot(mydf, aes(x = No.Users, y = Averted.kWh, color = Type)) + 
              geom_point(size = 3) + geom_line(size = 1.5) + facet_wrap(~TemperatureF, nrow = 1)
            plt1 = plt1 + geom_errorbar(aes(ymax = Averted.kWh + sd, ymin=Averted.kWh - sd))
            plt1 = plt1 + theme_bw() + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    strip.text.x     = element_text(size=18),
                    axis.text.y      = element_text(size=18), 
                    axis.text.x      = element_text(size=18, angle = 330),
                    axis.title.y     = element_text(size=18),
                    axis.title.x     = element_text(size=18),
                    plot.title       = element_text(size=20),            
                    legend.text      = element_text(size=18),
                    axis.ticks = element_blank()) + 
              ggtitle(paste("Marginal Benefit:", capitalize(focus), "Focus")) + xlab('Number of Users') + ylab('Averted Consumption [kWh]')
            
            plt = plt1 # multiplot(plt1, plt2)
            return(plt)            
          }                        
            
            # query a user ID for average temperature response
            if (type == 'avg-response-query') {
              
              df.resp = subset(x@EFF_RESPONSE, UID %in% query)
              
              # construct plot
              plt = ggplot(df.resp, aes(y = Thermal.Component, x = TemperatureF, color = UID)) + 
                geom_point(size = 2) + geom_line(size=1.5)
              plt = plt + geom_errorbar(aes(ymax = Thermal.Component + Var.Thermal, ymin=Thermal.Component - Var.Thermal))
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      legend.title     = element_text(size=18),
                      axis.ticks = element_blank()) + 
                ggtitle("Distribution of Thermal Regimes") + ylab('Effective Thermal Response [kWh/F]') + xlab('Temperature [F]')
              
              return(plt)            
            }
            
            # query a user ID for average temperature response
            if (type == 'avg-duration-query') {
              
              df.resp = subset(x@EFF_DURATION, UID %in% query)
              
              # construct plot
              plt = ggplot(df.resp, aes(y = Thermal.Duration, x = TemperatureF, color = UID)) + 
                geom_point(size = 2) + geom_line(size=1.5)
              plt = plt + geom_errorbar(aes(ymax = Thermal.Duration + Var.Thermal, ymin=Thermal.Duration - Var.Thermal))
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      legend.title     = element_text(size=18),
                      axis.ticks = element_blank()) + 
                ggtitle("Distribution of Thermal Regimes") + ylab('Effective Thermal Duration [hrs]') + xlab('Temperature [F]')
              
              return(plt)            
            }
            
          # query a user ID for average temperature response; plot on map
          if (type == 'avg-response-map') {
            
            df.resp = subset(x@EFF_RESPONSE, TemperatureF %in% c(40, 60, 80, 100))
            df.resp$Var.Thermal = NULL
            
            # get general location for zipcodes                      
            myzips = subset(zipcode, zip %in% df.resp$ZIPCODE, select = c('zip', 'latitude', 'longitude'))
            names(myzips) = c('ZIPCODE', 'lat', 'lon')
            df = merge(df.resp, myzips)
            
            map = get_googlemap(center  = c(mean(df$lon), mean(df$lat)), 
                                zoom    = 8,
                                source  = 'stamen',
                                maptype = 'roadmap')
            p   = ggmap(map, extent = 'device')
            p   = p + stat_density2d(data = df, 
                                     aes(x = lon, y = lat, 
                                         size = ..level..,  alpha = ..level..),
                                     size = 2, bins = 10, geom = 'polygon')  
            p   = p + geom_point(data = df, aes(x = lon, y = lat, color = Thermal.Component), size = 2)
            p   = p + scale_color_gradient(low = "blue", high = "red")
            p   = p + facet_wrap(~TemperatureF, nrow = 1)
            # format plot          
            p = p + theme_bw() + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    strip.text.x     = element_text(size=18),
                    axis.text.y      = element_text(size=18), 
                    axis.text.x      = element_text(size=18),
                    axis.title.y     = element_text(size=18),
                    axis.title.x     = element_text(size=18),
                    plot.title       = element_text(size=20),            
                    legend.text      = element_text(size=18),
                    legend.title     = element_text(size=18),
                    axis.ticks = element_blank()) + 
              ggtitle("Distribution of Thermal Regimes")
            
            return(p)            
          }                        
            
          # query a user ID for average temperature response
          if (type == 'tod-response-query') {
            
            df.resp = subset(x@TOD_RESPONSE, UID %in% query)
            
            # construct plot
            plt = ggplot(df.resp, aes(y = Thermal.Component, x = Hour, color = Season)) + 
              geom_point(size = 2) + geom_line(size=1.5)
            plt = plt + geom_errorbar(aes(ymax = Thermal.Component + Var.Thermal, ymin=Thermal.Component - Var.Thermal))
            plt = plt + facet_wrap(~UID, ncol = 2)
            plt = plt + theme_bw() + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    strip.text.x     = element_text(size=18),
                    axis.text.y      = element_text(size=18), 
                    axis.text.x      = element_text(size=18),
                    axis.title.y     = element_text(size=18),
                    axis.title.x     = element_text(size=18),
                    plot.title       = element_text(size=20),            
                    legend.text      = element_text(size=18),
                    legend.title     = element_text(size=18),
                    axis.ticks = element_blank()) + 
              ggtitle("Distribution of Thermal Regimes") + ylab('Effective Thermal Response [kWh/F]') + xlab('Time of Day')
            
            return(plt)            
          }
          
          # query a user ID for average temperature response
          if (type == 'tod-response-distr') {
            
            # classify response in high/medium/low
            df.resp = x@TOD_RESPONSE
            df.resp = merge(df.resp, query)
            df.resp$ZONE = as.character(df.resp$ZONE)
                        
            df.resp$Type  = 'cooling'
            df.resp$Type[df.resp$Thermal.Component < 0] = 'heating'
            qnt      = quantile(abs(df.resp$Thermal.Component), probs = c(0,0.25,0.75, 1))
            type     = c('none', 'low', 'high')[findInterval(abs(df.resp$Thermal.Component), qnt, all.inside = T)]
            df.resp$Type[type == 'none'] = ''
            df.resp$Type  = paste(type, df.resp$Type, sep='-')     
                    
            # compute stats
            tab = t(table(df.resp$Type, df.resp$ZONE))
            tab = tab / rowSums(tab)
            df  = as.data.frame(tab)
            names(df) = c('Zone', 'Type', 'Breakdown')
            
            # construct plot
            plt = ggplot(df, aes(x=Zone, fill = Type, y = Breakdown)) +
              geom_bar(stat = 'identity')
            plt = plt + theme_bw() + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    strip.text.x     = element_text(size=18),
                    axis.text.y      = element_text(size=18), 
                    axis.text.x      = element_text(size=18),
                    axis.title.y     = element_text(size=18),
                    axis.title.x     = element_text(size=18),
                    plot.title       = element_text(size=20),            
                    legend.text      = element_text(size=18),
                    legend.title     = element_text(size=18),
                    axis.ticks = element_blank()) + 
              ggtitle("Distribution of Thermal Regimes by Zone") 
            
            return(plt)            
          }                    
})
            
  