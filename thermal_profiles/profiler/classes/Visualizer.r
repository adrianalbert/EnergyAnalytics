# Visualizer.r
#
# Plot capabilities for non-homogenous HMMs.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------


source('../../visHMM/plot_utils.r')
source('../../visHMM/acf_ggplot.r')

# ________________________
# Class definition

setClass(
  Class = "Visualizer",
  representation  = representation(
    UID           = "character",         # unique person ID
    timestamps    = "POSIXct",
    obs           = "numeric",
    t.covar       = "data.frame",    
    r.covar       = "data.frame",    
    HMM           = "list"    
  )
)

# ________________________
# Initializer 

setMethod(f = "initialize", 
          signature = "Visualizer",
          definition = function(.Object, decoder, interpreter, interval = c(), verbose = T, 
                                t.covar = NULL, r.covar = NULL) {
            
            if (verbose) {
              cat(paste('*** Initializing Visualizer (', .Object@UID, ') ***\n', sep=''))
              t0 = proc.time()
            }
            
            .Object@UID         = decoder@UID
            .Object@timestamps  = as.POSIXct(decoder@data.train$timestamps)
            .Object@obs         = decoder@data.train$obs
            HMM         = list()    
            HMM$types   = interpreter@regime.type
            HMM$means   = decoder@HMM$fit
            
            if (!is.null(r.covar)){            
              .Object@r.covar   = as.data.frame(decoder@data.train[,r.covar])
              colnames(.Object@r.covar) = r.covar
              HMM$r.covar   = t(decoder@HMM$response$means[r.covar,])[decoder@HMM$states[,1]]
            } else {
              .Object@r.covar = NULL
              HMM$r.covar = NULL
            }            
            if (!is.null(t.covar)){            
              .Object@t.covar   = as.data.frame(decoder@data.train[,t.covar])
              colnames(.Object@t.covar) = t.covar
            } else {
              .Object@t.covar = NULL
            }            
            
            HMM$states  = interpreter@regime.type[decoder@HMM$states[,1]]
            HMM$sigma   = decoder@HMM$response$stdev[HMM$states]
            HMM$residual= decoder@HMM$residual      
            HMM$probProfile = interpreter@benchmarks$probProfile
            HMM$steadyDistr = interpreter@benchmarks$steadyDistr
            HMM$performance = decoder@performance
            HMM$contributions = interpreter@contributions$agg
            HMM$contrib.ts  = interpreter@contributions$ts
            HMM$fit         = decoder@HMM$fit
            
            if (!is.null(interval)) {
              idx_ok      = which(.Object@timestamps >= as.POSIXct(interval[1]) & 
                                  .Object@timestamps <= as.POSIXct(interval[2]))
              .Object@timestamps  = .Object@timestamps[idx_ok]
              .Object@obs = .Object@obs[idx_ok]
              HMM$means   = HMM$means[idx_ok]
              HMM$sigma   = HMM$sigma[idx_ok]
              HMM$states  = HMM$states[idx_ok]
              HMM$residual= HMM$residual[idx_ok]
              HMM$r.covar   = HMM$r.covar[idx_ok]
              HMM$contrib.ts  = HMM$contrib.ts[idx_ok,]
              HMM$fit         = HMM$fit[idx_ok]
              .Object@t.covar= as.data.frame(.Object@t.covar[idx_ok,])
              colnames(.Object@t.covar) = t.covar
              .Object@r.covar= as.data.frame(.Object@r.covar[idx_ok,])
              colnames(.Object@r.covar) = r.covar
            }  
            .Object@HMM = HMM
            
            return(.Object)
            
          })
# _______________________
# Plot Visualizer object

setMethod('plot',
          signature  = 'Visualizer',
          definition = function(x, type = 'default', verbose = T){
                        
            if (type == 'default') {
              plot(x@timestamps, x@obs, 
                   main = paste('(', x@UID,')'),
                   xlab = 'Time', ylab = 'kWh', type = 'l', lwd = 2)
            }
                                    
            if (type == 'HMM-ts') {
              title = paste("Zoom-In: Data and Decoded States (", x@UID,')',sep='')                       
              plt = plot_hmm_ts(x@HMM$means, x@HMM$sigma, x@HMM$states, x@timestamps, x@obs, 
                                y.lab = 'kWh', title = title)
              return(plt)
            }
            
            # plot time series of state parameters
            if (type == 'HMM-coefs-ts') {
              
              r.covar = colnames(x@r.covar)
              covar_state = as.data.frame(as.matrix(x@HMM$r.covar))
              colnames(covar_state) = r.covar
              
              title = paste("Zoom-In: States Parameters (", x@UID,')',sep='')                       
              plt = plot_hmm_coefs_ts(covar_state, x@HMM$states, x@obs, x@timestamps,  
                                      title = title)
              return(plt)
            }
            
            # plot ToD/Dow breakdown of occupancy states
            if (type == 'HMM-state-breakdown') {
              title = paste("States breakdown (", x@UID,')',sep='')                       
              plt   = plot_state_breakdown(as.factor(x@HMM$states), x@timestamps, 
                                           state.type =  NULL, title = title)
              return(plt)
            }
            
            # plot heatmap of quantity of interest
            if (type == 'HMM-heatmap-obs') {
              myMat = -x@obs 
              title = paste("Heatmap (", x@UID,')',sep='')                       
              plot_state_heatmap(myMat, x@timestamps, title = title)
            }            
            
            # plot heatmap of occupancy states (ggplot)
            if (type == 'HMM-state-heatmap-states') {
              myMat = as.factor(x@HMM$states) 
              title = paste("States heatmap (", x@UID,')',sep='')                       
              plt = plot_state_heatmap2(myMat, x@timestamps, title = title)
              return(plt)
            }
            
#             if (type == 'HMM-MC-cov') {
#               title             = paste("State Space Diagram (", x@UID, ')',sep='')
#               contrib           = x@HMM$components$stat
#               contrib$mu        = as.numeric(x@HMM$response$means[1,])
#               contrib$sigma2    = as.numeric(x@HMM$response$stdev)
#               plt               = plot_HMM_MC_cov(x@HMM$P.trans, x@HMM$state.duration, contrib, title=title)
#               return(plt)
#             }
            
            if (type =='HMM-pacf') {
              
              plt = acf_ggplot(na.omit(x@obs), na.omit(x@HMM$means), 
                               title = 'PACF: Empirical vs Model', PACF = T)
              return(plt)
            }
            
            if (type == 'HMM-res') { 
              plt = plot_HMM_res(x@HMM$residual, x@HMM$states, x@HMM$response$stdev)
              plt
            }             
            
            if (type == 'HMM-dep-covar') {
              t.covar = colnames(x@t.covar)
              title = paste(paste(t.covar, " dependence",sep=''))
              states = x@HMM$states
#               states = paste(states, ': slope=', 
#                              round(x@HMM$response$means[covar,states]*100, digits=3),'(/100)',sep='')
              
              p     = plot_dep_covar(as.numeric(x@t.covar[,1]),
                                     x@obs, 
                                     states, 
                                     title = title, x.lab = t.covar, y.lab = 'kWh')  	
              return(p)
            }            
            
            if (type == 'HMM-dep-covar-sep') {
              t.covar = colnames(x@t.covar)
              title = paste(paste(t.covar, " dependence (", x@UID,')',sep=''))
              states = x@HMM$states
              p     = plot_dep_covar(as.numeric(x@t.covar[,1]), x@obs, states, 
                                     title = title, x.lab = t.covar, y.lab = 'kWh', separate=T)    
              return(p)
            }            
            
            if (type == 'HMM-trans-prob') {
              t.covar = colnames(x@t.covar)
              title = paste("Transition probabilities" ,sep='')  
              x_var = as.data.frame(1:nrow(x@HMM$probProfile[[1]]))
              names(x_var) = t.covar
              p     = plot_tran_covar(x_var, x@HMM$probProfile, 
                                      title = title, x.lab = t.covar, y.lab = 'P', markers = NULL)  	
              return(p)
            }            
            
            if (type == 'HMM-stationary-prob') {
              t.covar = colnames(x@t.covar)
              title = paste("Stationary distribution",sep='')
              x_var = as.data.frame(1:nrow(x@HMM$probProfile[[1]]))
              names(x_var) = t.covar
              dep   = x@HMM$steadyDistr
              dep   = list(as.matrix(dep))
              type  = x@HMM$types

              # plot 
              p     = plot_tran_covar(x_var, dep, type = type,
                                      title = title, x.lab = t.covar, y.lab = 'P')    
              return(p)
            }            
            
            if (type == 'HMM-err-horiz') {
              title = paste("Out-of-sample decoding performance (", x@UID,')',sep='')
              perf  = x@HMM$performance
              x_var = rownames(x@HMM$performance)
              perf  = data.frame(perf)
              perf$Horizon = rownames(perf)
              df    = melt(perf, id.vars = 'Horizon')
              df$Horizon = as.numeric(df$Horizon)
              # plot 
              plt = ggplot(df, aes(x = Horizon, y = value))
              plt = plt + geom_point(aes(color = variable, shape = variable), show_guide = FALSE, size=4)
              plt = plt + geom_line(aes(color = variable), show_guide = FALSE, size=2.5)
              plt = plt + facet_wrap(~variable, ncol = 2, scales = 'free')
              
              # formatting
              plt = plt + 
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
                      legend.text      = element_text(size=15),
                      legend.title     = element_text(size=15),
                      axis.ticks       = element_blank() ) + 
                ylab('Performance') + xlab('Time Horizon [hours]') + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( title )
              return(plt)
            }            
                        
            if (type == 'HMM-aggregate-season') {
              title = paste("Seasonal analysis (", x@UID,')',sep='')              
              df.mlt = melt(x@HMM$contributions, id.vars = c('Day', 'Hour', 'Season'))
              t.covar = colnames(x@t.covar)
              
              # plot 
              plt = ggplot(df.mlt, aes(x = Hour, y = value, color = Season, shape = Day))
              plt = plt + geom_point(size = 4) + geom_line()
              plt = plt + facet_wrap(~variable, scales = 'free')
              
              # formatting
              plt = plt + 
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
                      legend.text      = element_text(size=15),
                      legend.title     = element_text(size=15),
                      axis.ticks       = element_blank() ) + 
                ylab('kWh') + xlab('Hour of Day') + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( title )
              return(plt)
            }          
            
            if (type == 'HMM-contrib-ts'){
              title = paste("Zoom-In: Contributions (", x@UID,')',sep='')
              df = x@HMM$contrib.ts
              df$obs = x@obs
              df$fit = x@HMM$fit
              
              r.covar = colnames(x@r.covar)
              covar_state = as.data.frame(as.matrix(x@HMM$r.covar))
              colnames(covar_state) = r.covar 
              
              plot_components_ts(df, x@timestamps,
                                 title = title, states = x@HMM$states, covars = covar_state)              
            }
            
          })
