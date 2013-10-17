# Visualizer.r
#
# Plot capabilities for non-homogenous HMMs.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------


# source('../utils/viz/plot_utils.r')
source('../utils/viz/acf_ggplot.r')

# ________________________
# Class definition

setClass(
  Class = "Visualizer",
  representation  = representation(
    UID           = "character",         # unique person ID
    timestamps    = "POSIXct",
    obs           = "numeric",
    covar         = "numeric",    
    HMM           = "list"    
  )
)

# ________________________
# Initializer 

setMethod(f = "initialize", 
          signature = "Visualizer",
          definition = function(.Object, decoder, interpreter, interval = c(), verbose = T, covar = NULL) {
            
            if (verbose) {
              cat(paste('*** Initializing Visualizer (', .Object@UID, ') ***\n', sep=''))
              t0 = proc.time()
            }
            
            .Object@UID         = decoder@UID
            .Object@timestamps  = as.POSIXct(decoder@data.train$timestamps)
            .Object@obs         = decoder@data.train$obs
            .Object@covar       = decoder@data.train[,covar]
            
            HMM         = list()            
            HMM$means   = decoder@HMM$fit
            HMM$states  = decoder@HMM$states[,1]
            HMM$sigma   = decoder@HMM$response$stdev[HMM$states]
            HMM$residual= decoder@HMM$residual      
                        
            if (!is.null(interval)) {
              idx_ok      = which(.Object@timestamps >= as.POSIXct(interval[1]) & 
                                  .Object@timestamps <= as.POSIXct(interval[2]))
              .Object@timestamps  = .Object@timestamps[idx_ok]
              .Object@obs = .Object@obs[idx_ok]
              HMM$means   = HMM$means[idx_ok]
              HMM$sigma   = HMM$sigma[idx_ok]
              HMM$states  = HMM$states[idx_ok]
              HMM$residual= HMM$residual[idx_ok]
              .Object@covar= .Object@covar[idx_ok]
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
              title = paste("Zoom-In: States and Observed Emissions (", x@UID,')',sep='')                       
              plt = plot_hmm_ts(x@HMM$means, x@HMM$sigma, x@HMM$states, x@timestamps, x@obs, 
                                y.lab = 'kWh', title = title)
              return(plt)
            }
            
            # plot time series of state parameters
            if (type == 'HMM-coefs-ts') {
              covar_state = as.data.frame(t(as.matrix(x@HMM$response$means[covar,])))
              names(covar_state) = covar
              title = paste("Zoom-In: States Parameters (", x@UID,')',sep='')                       
              plt = plot_hmm_coefs_ts(covar_state, x@HMM$states, x@obs, x@timestamps,  
                                      title = title)
              return(plt)
            }
            
            # plot ToD/Dow breakdown of occupancy states
            if (type == 'HMM-state-breakdown') {
              title = paste(covar, " states breakdown (", x@UID,')',sep='')                       
              plt   = plot_state_breakdown(as.factor(states), timestamps, 
                                           state.type =  NULL, #x@HMM$statMetrics$StateType,
                                           title = title)
              return(plt)
            }
            
            # plot heatmap of quantity of interest
            if (type == 'HMM-heatmap-obs') {
              myMat = -x@obs 
              title = paste("Heatmap (", x@UID,')',sep='')                       
              plot_state_heatmap(myMat, x@timestamps, title = title)
            }            
            
            # plot heatmap of occupancy states (ggplot)
            if (type == 'HMM-state-heatmap2') {
              if (covar == 'States') myMat = as.factor(states) else 
                if (covar == 'kWh') myMat = -kWh else {
                  myMat = -t(x@HMM$response$means[covar,states])
                  colnames(myMat) = NULL
                }
              title = paste(covar, " heatmap (", x@UID,')',sep='')                       
              plt = plot_state_heatmap2(myMat, timestamps, title = title)
              return(plt)
            }
            
            if (type == 'HMM-MC-cov') {
              title             = paste("State Space Diagram (", x@UID, ')',sep='')
              contrib           = x@HMM$components$stat
              contrib$mu        = as.numeric(x@HMM$response$means[1,])
              contrib$sigma2    = as.numeric(x@HMM$response$stdev)
              plt               = plot_HMM_MC_cov(x@HMM$P.trans, x@HMM$state.duration, contrib, title=title)
              return(plt)
            }
            
            if (type =='HMM-pacf') {
              
              plt = acf_ggplot(na.omit(kWh), na.omit(hmm.means), 
                               title = 'PACF: Empirical vs Model', PACF = PACF)
              return(plt)
            }
            
            if (type == 'HMM-res') { 
              plt = plot_HMM_res(hmm.residual, states, x@HMM$response$stdev, sw_test = x@HMM$sw_test)
              plt
            }             
            
            if (type == 'HMM-contrib-ts'){
              title = paste("Zoom-In: HMM Covariate Contributions (", x@UID,')',sep='')
              covar_state = as.data.frame(t(as.matrix(x@HMM$response$means[covar,])))
              names(covar_state) = covar              
              plot_components_ts(comps.hmm, timestamps, 
                                 title = title, states = states, 
                                 covars = covar_state)              
            }
            
            if (type == 'HMM-dep-covar') {
              title = paste(covar, paste("dependence (", x@UID,')',sep=''))
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H')
              states = x@HMM$states[,1]
              #              states = paste(states, x@HMM$statMetrics$StateType[states],sep=':')
              #               states = paste(states, ': slope=', 
              #                              round(x@HMM$response$means[covar,states]*100, digits=3),'+/-',
              #                              round(x@HMM$response$error[covar,states]*100, digits=3),'(/100)',
              #                              sep='')
              states = paste(states, ': slope=', 
                             round(x@HMM$response$means[covar,states]*100, digits=3),'(/100)',sep='')
              p     = plot_dep_covar(x@data.train[,covar], x@data.train$kWh, states, 
                                     title = title, x.lab = covar, y.lab = 'kWh', highlight = highlight,
                                     markers = marks, separate = separate)  	
              return(p)
            }            
            
            if (type == 'HMM-trans-prob') {
              title = paste("Transition probabilities (", x@UID,')',sep='')
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H') 
              
              # compute transition probabilities over a predefined grid
              x_var = x@HMM$model@transition[[1]]@x
              x_var = cbind(rep(1,100), seq(min(x_var[,covar]), max(x_var[,covar]), length.out=100))
              colnames(x_var) = colnames(x@HMM$model@transition[[1]]@x)
              dep = lapply(x@HMM$model@transition, function(l) {
                y = x_var %*% l@parameters$coefficients 
                y = exp(y) / rowSums(exp(y))
                return(y)
              })
              
              # plot 
              x_var = as.data.frame(x_var[,covar])
              names(x_var) = covar
              p     = plot_tran_covar(x_var, dep, 
                                      title = title, x.lab = covar, y.lab = 'P', markers = marks)  	
              return(p)
            }            
            
            if (type == 'HMM-stationary-prob') {
              title = paste("Stationary transition probabilities (", x@UID,')',sep='')
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H') 
              x_var = data.frame(x@HMM$statMetrics$distr.covar[,covar])
              names(x_var) = covar
              dep   = x@HMM$statMetrics$distr.covar[,-which(names(x@HMM$statMetrics$distr.covar) == covar)]
              dep   = list(as.matrix(dep))
              
              # plot 
              p     = plot_tran_covar(x_var, dep,
                                      title = title, x.lab = covar, y.lab = 'P', markers = marks)    
              return(p)
            }            
            
            if (type == 'HMM-err-horiz') {
              title = paste("Out-of-sample decoding performance (", x@UID,')',sep='')
              perf  = x@HMM$PredictStats
              x_var = rownames(x@HMM$PredictStats)
              perf  = data.frame(perf)
              perf$Horizon = rownames(perf)
              df    = melt(perf, id.vars = 'Horizon')
              df$Horizon = as.numeric(df$Horizon)
              # plot 
              plt = ggplot(df, aes(x = Horizon, y = value))
              plt = plt + geom_point(aes(color = variable, shape = variable), show_guide = FALSE, size=4)
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
            
            if (type == 'HMM-time-temp') {
              title = paste("State duration (", x@UID,')',sep='')
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H') 
              tstart= round(marks[1])
              tstop = round(marks[2])              
              dep   = as.data.frame(x@HMM$statMetrics$time.temp[tstart:tstop,])
              names(dep) = paste('State', 1:ncol(dep))
              dep[,covar] = tstart:tstop
              df.mlt = melt(dep, id.vars = covar)
              
              # plot 
              plt = ggplot(df.mlt, aes_string(x = covar, y = 'value'))
              plt = plt + geom_point(aes(color = variable, shape = variable), size = 4)
              
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
                ylab('Duration') + xlab(covar) + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( title )
              return(plt)
            }              
            
            if (type == 'HMM-aggregate-season') {
              title = paste("Seasonal analysis (", x@UID,')',sep='')              
              df.mlt = melt(x@HMM$components$agg.tod, id.vars = c('Day', 'Hour', 'Season'))
              
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
                ylab('kWh') + xlab(covar) + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( title )
              return(plt)
            }                                  
          })
