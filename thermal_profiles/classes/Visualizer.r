# Plotter.r
#
# Interpret output from HMM decoder.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------

# _______________________
# Plot Visualizer object

setMethod('plot',
          signature  = 'Visualizer',
          definition = function(x, type = 'default', interval = NULL, verbose = T, 
                                covar = NULL, PACF = T, highlight = T, separate = F, ...){
            
            if (verbose)
              cat(paste('*** ', type, ' Plot for Visualizer (', x@UID, ') ***\n', sep=''))
            
            timestamps = as.POSIXct(x@data.train$timestamps)
            ols.resids = x@OLS[['fit.summary']]$residuals
            ols.resids = naresid(x@OLS$fit.summary$na.action, ols.resids)
            kWh= x@data.train$kWh
            
            ols.means  = x@OLS[['fitted.vals']]
            hmm.means  = x@HMM$fit
            states     = x@HMM$states[,1]
            hmm.sigma  = x@HMM$response$stdev[states]
            hmm.residual=x@HMM$residual      
            wthr_data   =x@weather
            comps.hmm   =x@HMM$components$ts
            if (!is.null(interval)) {
              idx_ok      = which(x@timestamps >= as.POSIXct(interval[1]) & 
                                    x@timestamps <= as.POSIXct(interval[2]))
              timestamps  = timestamps[idx_ok]
              ols.resids  = ols.resids[idx_ok]
              kWh = kWh[idx_ok]
              ols.means   = ols.means[idx_ok]
              hmm.means   = hmm.means[idx_ok]
              hmm.sigma   = hmm.sigma[idx_ok]
              states      = states[idx_ok]
              hmm.residual= hmm.residual[idx_ok]
              wthr_data   = wthr_data[idx_ok,]
              comps.hmm   = comps.hmm[idx_ok,]
            }  
            
            if (type == 'default') {
              plot(timestamps, kWh, 
                   main = paste('kWh Visualizer-SPID (', x@UID,')'),
                   xlab = 'Time', ylab = 'kWh', type = 'l', lwd = 2)
            }
            
            if (type == 'weather') {
              if (is.null(wthr_data) | nrow(wthr_data)==0) {
                cat('No weather data!\n')
              } else {
                wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                                'HourlyPrecip', 'SolarRadiation')
                wthr = wthr_data[,names(wthr_data) %in% wthr_names]
                wthr$kWh = kWh
                wthr = zoo(wthr, order.by = timestamps)  
                plot(wthr, 
                     main = paste('kWh Visualizer-SPID (', x@UID,')'),
                     xlab = 'Time', type = 'l', lwd = 2)                
              }
            }
            
            if (type == 'OLS-res') {
              layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))            
              plot(timestamps, ols.resids, 
                   main = 'OLS Residuals', ylab = 'Residuals', xlab = 'Time', lwd = 3, type='l')
              plot(density(na.omit(ols.resids)), 
                   main = 'OLS Residuals: Density', ylab = 'pdf', xlab = 'Residuals', lwd = 3, type='l')
              acf(na.omit(ols.resids), 
                  main = 'OLS Residuals', ylab = 'Residuals', xlab = 'Time', lwd = 3)
            }    
            
            if (type == 'OLS-fit') {
              yrange = range(c(ols.means, kWh))
              plot(timestamps, ols.means, 
                   main = 'OLS Fitted Values', ylab = 'Values', xlab = 'Time', 
                   lwd = 3, type='l', col='red', ylim = yrange)
              points(timestamps, kWh, type='b', pch=20, col = 'black')
              legend('topleft', c('OLS Fit', 'Observations'), col=c('red', 'black'))
            }    
            
            if (type == 'HMM-ts') {
              title = paste("Zoom-In: States and Observed Emissions (", x@UID,')',sep='')                       
              plt = plot_hmm_ts(hmm.means, hmm.sigma, states, timestamps, kWh, 
                                y.lab = 'kWh', title = title)
              return(plt)
            }
            
            # plot time series of state parameters
            if (type == 'HMM-coefs-ts') {
              covar_state = as.data.frame(t(as.matrix(x@HMM$response$means[covar,])))
              names(covar_state) = covar
              title = paste("Zoom-In: States Parameters (", x@UID,')',sep='')                       
              plt = plot_hmm_coefs_ts(covar_state, states, kWh, timestamps,  
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
            
            # plot heatmap of occupancy states
            if (type == 'HMM-state-heatmap') {
              if (covar == 'States') myMat = as.factor(states) else 
                if (covar == 'kWh') myMat = -kWh else {
                  myMat = -t(x@HMM$response$means[covar,states])
                  colnames(myMat) = NULL
                }
              title = paste(covar, " heatmap (", x@UID,')',sep='')                       
              plot_state_heatmap(myMat, timestamps, title = title)
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
            
            if (type == 'HMM-MC2') {
              title = paste("HMM MC structure (", x@UID, ')',sep='') 
              # organize MC information
              plt = plot_HMM_MC_NET(as.numeric(x@HMM$response$means[1,]), as.numeric(x@HMM$response$stdev), x@HMM$transition)
              plot(plt, main = title) 
            }
            
            if (type == 'HMM-contrib-ts'){
              title = paste("Zoom-In: HMM Covariate Contributions (", x@UID,')',sep='')
              covar_state = as.data.frame(t(as.matrix(x@HMM$response$means[covar,])))
              names(covar_state) = covar              
              plot_components_ts(comps.hmm, timestamps, 
                                 title = title, states = states, 
                                 covars = covar_state)              
            }
            
            if (type == 'HMM-contrib-tot') {
              df = x@HMM$components$stat
              title = paste("HMM Covariate Contributions (", x@UID,')',sep='')
              p     = plot_tornado(df, title = title, nrow = 2)		
              return(p)
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

# _________________________________________
# Method to dump Visualizer model data to file
setGeneric(
  name = "dumpComputationToFile",
  def = function(.Object, verbose = T, path = NULL){standardGeneric("dumpComputationToFile")}
)
setMethod('dumpComputationToFile',
          signature  = 'Visualizer',
          definition = function(.Object, verbose = T, path = NULL){
            if (verbose) 
              cat(paste('*** Dumping model data to file (', .Object@UID, ')***\n', sep=''))
            
            # ______________________________
            # Format analysis results: OLS
            
            df     = as.data.frame(t(coef(.Object@OLS$fit.summary)[,1]))
            df$R2  = .Object@OLS$fit.summary$adj.r.squared
            df$GR2 = .Object@OLS$GR2
            df$UID = .Object@UID
            df$ZIPCODE= .Object@ZIPCODE
            df$dwTest = .Object@OLS$dw.test
            df$bpTest = .Object@OLS$bp.test 
            df.ols = df
            
            # ______________________________
            # Format analysis results: HMM
            
            # stats
            df        = as.data.frame(.Object@HMM$entropy)
            names(df) = paste('entropy',names(df), sep='.')          
            df$UID    = .Object@UID
            df$ZIPCODE= .Object@ZIPCODE
            df$GR2    = .Object@HMM$GR2
            df$R2     = .Object@HMM$R2
            df$MAPE   = .Object@HMM$MAPE
            df$MAPE.cv= .Object@HMM$cv.metrics[as.character(.Object@HMM$nStates),'MAPE.cv']
            df$R2.cv  = .Object@HMM$cv.metrics[as.character(.Object@HMM$nStates),'R2.cv']
            tmp       = as.data.frame(.Object@HMM$predict)
            names(tmp)= paste('predict',names(tmp), sep='.')
            df        = cbind(df, tmp)
            df        = cbind(df, data.frame(t(.Object@breakpoint$psi)))
            df.stat   = df
            
            # response
            
            df        = as.data.frame(t(.Object@HMM$response$means))
            df$Sigma  = .Object@HMM$response$stdev
            df$UID    = .Object@UID
            df$ZIPCODE= .Object@ZIPCODE
            df$State  = 1:.Object@HMM$nStates
            df.resp   = df
            
            # benchmarks: temperature-dependent state distribution
            x_var     = .Object@HMM$statMetrics$distr.covar[,1]
            df        = as.data.frame(t(.Object@HMM$statMetrics$distr.covar[,-1]))
            names(df) = x_var
            df.distr  = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), variable = 'Distribution', df)
            
            # benchmarks: temperature-dependent state duration
            df        = as.data.frame(t(.Object@HMM$statMetrics$time.temp))
            names(df) = x_var
            df.time   = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), variable = 'Duration', df)
            
            df.tran   = rbind(df.distr, df.time)
            rownames(df.tran) = NULL
            
            # benchmarks: state breakdown by time-of-day
            df = data.frame(matrix(.Object@HMM$statMetrics$breakdown$Summer, ncol = 24, byrow=T))
            names(df) = 0:23
            df.summer = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), Season = 'Summer', df)
            
            df = data.frame(matrix(.Object@HMM$statMetrics$breakdown$Winter, ncol = 24, byrow=T))
            names(df) = 0:23
            df.winter = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), Season = 'Winter', df)
            
            df.season = rbind(df.summer, df.winter)
            
            # ________________________
            # Component contributions
            
            df.comp  = .Object@HMM$components$agg.tod
            df.comp  = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, df.comp)
            
            # ______________________________
            # Dump to file
            
            if (!is.null(path)) {
              write.csv(df.ols,    file = paste(path,.Object@UID, '_ols.csv', sep=''), quote = F, row.names = F)
              write.csv(df.resp,   file = paste(path,.Object@UID, '_hmm_resp.csv', sep=''), quote = F, row.names = F)
              write.csv(df.tran,   file = paste(path,.Object@UID, '_hmm_tran.csv', sep=''), quote = F, row.names = F)
              write.csv(df.season, file = paste(path,.Object@UID, '_hmm_seas.csv', sep=''), quote = F, row.names = F)
              write.csv(df.comp,   file = paste(path,.Object@UID, '_hmm_comp.csv', sep=''), quote = F, row.names = F)
            }
            
            return(list(OLS = df.ols, 
                        HMM.response = df.resp, HMM.transition = df.tran, HMM.stats = df.stat,
                        HMM.season = df.season, HMM.comp = df.comp))
          })

# _________________________________________
# Method to dump Visualizer model data to file
setGeneric(
  name = "dumpComputation",
  def = function(.Object, verbose = T, path = NULL){standardGeneric("dumpComputation")}
)
setMethod('dumpComputation',
          signature  = 'Visualizer',
          definition = function(.Object, verbose = T, path = NULL){
            if (verbose) 
              cat(paste('*** Dumping model data to RData file (', .Object@UID, ')***\n', sep=''))
            
            # trim down object for saving to disc
            tmp         = list()
            tmp$UID     = .Object@UID
            tmp$ZIP     = .Object@ZIPCODE
            tmp$OLS     = .Object@OLS[c('heterosked', 'serial.corr', 'R2', 'GR2')]
            tmp$breakpoint = .Object@breakpoint
            tmp$HMM     = .Object@HMM[c('nStates', 'cv.metrics', 'PredictStats',
                                        'response', 'BIC', 'sw_test', 'MAPE',
                                        'GR2', 'R2', 'state.duration', 'time.temp',
                                        'P.trans', 'statMetrics', 'fit')]
            tmp$HMM$states = data.frame(states = .Object@HMM$states[,1],
                                        DateTime = .Object@data.train$timestamps)
            tmp$HMM$entropy = unlist(.Object@HMM$entropy)
            tmp$HMM$predict = unlist(.Object@HMM$predict)
            tmp$HMM$comps   = .Object@HMM$components
            
            # save to disc if requested
            if (!is.null(path)) 
              save(file = paste(path, .Object@UID, '.RData', sep=''), list = 'tmp')
            
            return(tmp)
          })

