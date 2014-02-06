# #########################################################################
# Scheduler.r
# -------------------------
#
# Interface to implement scheduling problems for thermally-sensitive users.
# 
# Adrian Albert
# Last modified: February 2014.
# #########################################################################

library('timeDate')
library('lubridate')
library('reshape')
library('CVXfromR')
library('ggplot2')

# clean-up previous definitions of methods for class Person
removeClass('Scheduler')

# ________________________
# Class definition

setClass(
  Class = "Scheduler",
  representation = representation(
    INPUTS       = "list",             # problem setup: exogenous inputs
    SETUP        = "list",             # variables for problem setup
    OUTPUT       = "list",             # problem output
    MODELS       = "list",             # model parameters for all users in the problem
    SOLVER       = "list"              # options for optimization solver
  )
)

# _____________________________________________________
# Constructor method for class Scheduler

setMethod(f = "initialize", 
          signature = "Scheduler",
          definition = function(.Object, importer, profile=NULL, target = NULL, setup = list(), verbose = T) {
            
            if (verbose) cat('*** Initializing Scheduler ***\n')
            
            # inputs
            .Object@INPUTS$HORIZON    = nrow(profile)
            .Object@INPUTS$DT         = setup$DT
            .Object@INPUTS$BUDGET     = setup$budget            
            .Object@INPUTS$PROFILE    = profile            
            .Object@INPUTS$TARGET     = target       
            
            # models
            .Object@MODELS$RESPONSE   = importer@STATES_PAR            
            
            # initialize states over the planning horizon
            .Object@MODELS$STATES     = propagateDistribution(.Object, importer@TRANSITION)
            
            # setup problem
            .Object@SETUP             = setupProblem(.Object)
            
            # solver parameters
            .Object@SOLVER$cvx.setup.dir = setup$cvx.setup.dir          
            
            return(.Object)
          })

# __________________________________________________________
# Method to set up & solve deterministic scheduling problem

# Puts the problem into CVX and calls matlab. 

setGeneric(
  name = "setupProblem",
  def = function(.Object, verbose = T){standardGeneric("setupProblem")}
)
setMethod('setupProblem',
          signature  = 'Scheduler',
          definition = function(.Object, verbose = T) {
            
            S = .Object@MODELS$STATES[,-1]
            
            # regimes
            R = unstack(S, State ~ Time)
            R = cbind(UID = unique(.Object@MODELS$STATES$UID), R)
            
            # response part
            A = subset(.Object@MODELS$RESPONSE, select = c('UID', 'State', 'TemperatureF'))
            A = merge(S, A, by = c('UID', 'State'))
            A = A[with(A, order(UID, Time)), ]
            A = unstack(A, TemperatureF ~ Time)
            # A = cbind(UID = unique(.Object@STATES$UID), A)
            A = as.matrix(A)
            
            # baseload part
            B = subset(.Object@MODELS$RESPONSE, select = c('UID', 'State', 'X.Intercept.'))
            names(B)[3] = 'Baseload'
            B = merge(S, B, by = c('UID', 'State'))
            B = B[with(B, order(UID, Time)), ]
            B = unstack(B, Baseload ~ Time)
            # B = cbind(UID = unique(.Object@STATES$UID), B)
            B = as.matrix(B)              
            b = colSums(B)
            
            # some constants
            p = nrow(.Object@INPUTS$TARGET)
            N = nrow(B)
            
            # temperature profiles for all users
            # this can be easily replaced with different temperature profiles for each user
            TF = .Object@INPUTS$PROFILE[,1] 
            TF = matrix(rep(TF, N), nrow = N, byrow=T)
            
            # save setup
            setup = list()                      
            setup$p   = p
            setup$N   = N
            setup$g   = .Object@INPUTS$TARGET[,1]
            setup$A   = A
            setup$B   = B
            setup$b   = b
            setup$UID = unique(.Object@MODELS$STATES$UID)
            setup$TF  = TF
            setup$DT  = .Object@INPUTS$DT
            setup$BUDGET = .Object@INPUTS$BUDGET
            return(setup)
          })


# __________________________________________________________
# Method to set up & solve deterministic scheduling problem

# Puts the problem into CVX and calls matlab. 

setGeneric(
  name = "solveDeterministic",
  def = function(.Object, verbose = T){standardGeneric("solveDeterministic")}
)
setMethod('solveDeterministic',
          signature  = 'Scheduler',
          definition = function(.Object, verbose = T) {
                        
            if (!is.null(.Object@SOLVER$cvx.setup.dir)) {
              
              A.T   = t(.Object@SETUP$A) %*% .Object@SETUP$TF
              A.hat = .Object@SETUP$A * .Object@SETUP$DT # abs(.Object@SETUP$A * .Object@SETUP$DT)

              # form CVX problem              
              cvxcode <- paste("variables U(N,p)",
                               "minimize( norm(g - b - trace(A.T - A.hat' * U), 2) )",
                               "subject to",
                               "abs(U) * ones(p,1) <= BUDGET",
                               "-1 <= U <= 1",
                               sep=";")
              const.vars = list(p      = .Object@SETUP$p, 
                                N      = .Object@SETUP$N, 
                                A.T    = A.T, 
                                A.hat  = A.hat, 
                                g      = .Object@SETUP$g, 
                                b      = .Object@SETUP$b, 
                                BUDGET = .Object@SETUP$BUDGET)
              opt <- CallCVX(cvxcode, 
                             const.vars=const.vars,
                             opt.var.names="U", setup.dir=.Object@SOLVER$cvx.setup.dir)
            }
            
            .Object@OUTPUT$deterministic$U = cbind(UID = .Object@SETUP$UID, as.data.frame(opt$U))
            .Object@OUTPUT$deterministic$L = .Object@SETUP$A * .Object@SETUP$TF
            .Object@OUTPUT$deterministic$C = A.hat * opt$U
            return(.Object)
          })


# _________________________
# Method to plot schedules 

setMethod('plot',
          signature  = 'Scheduler',
          definition = function(x, type = 'heatmap', selected = NULL){
            
            U = x@OUTPUT$deterministic$U
            if (!is.null(selected)) U = subset(U, UID %in% selected)            
            uids  = as.character(U$UID)
            U.mat = as.matrix(U[,-1])
#             idx   = which(rowSums(abs(U.mat)) < 1e-3)
#             U     = U[-idx,]
#             U.mat = U.mat[-idx,]
            
            # heatmap plot for lots of users
            if (type == 'heatmap') {           
              
              hmcol<-brewer.pal(50,'YlOrBr')
              heatmap(U.mat, Rowv = NA, Colv=NA, 
                      labRow = uids, 
                      labCol = paste('Hour', 1:x@SETUP$p), 
                      cexRow = 1.5, cexCol = 1.3, col = hmcol, main = 'Opportunistic Schedules')                          
            }
            
            # selected profiles (good for up to 10 users)
            if (type == 'profile') {   
              
              df = melt(U, id.vars = 'UID')
              df$variable = as.numeric(gsub('V', '', as.character(df$variable)))
              names(df) = c('UID', 'Hour', 'Action')
              df$UID = as.factor(df$UID)
              df$Action = df$Action * x@SETUP$DT
              
              if (!is.null(names(selected))) {
                main = paste(names(selected), collapse = ' and ')
                levels(df$UID) = names(selected)
              }
              
              # construct plot
              plt = ggplot(df, aes(y = Action, x = Hour, color = UID)) + 
                geom_point(size = 2.5) + geom_line(size=1.5)
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
                ggtitle(paste("Selected Schedules")) + ylab('Action [deg F]')
              
              return(plt)                                      
            }
            
            # density plot
            if (type == 'density') {   
              
              df = data.frame(UID = U$UID, Action = rowSums(U.mat))              
              df$Action = df$Action * x@SETUP$DT
              
              # construct plot
              plt = ggplot(df, aes(x = Action)) + 
                geom_histogram(size = 2) 
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
                ggtitle(paste("Distribution of Total Daily Effort:", nrow(df), "Users")) + xlab('Effort [deg F]')
              
              return(plt)                                        
            }
            
            # density plot (hourly)
            if (type == 'density-hourly') {   
              
              df = melt(U, id.vars = 'UID')
              df$variable = as.numeric(gsub('V', '', as.character(df$variable)))
              names(df) = c('UID', 'Hour', 'Action')
              df$UID = as.factor(df$UID)
              df$Action = df$Action * x@SETUP$DT
              
              # construct plot
              plt = ggplot(df, aes(x = Action)) + 
                geom_histogram(size = 2) +
                facet_wrap(~Hour, ncol = 6)
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
                ggtitle(paste("Distribution of Hourly Effort:", nrow(U), "Users")) + xlab('Effort [deg F]')
              
              return(plt)                                        
            }
            
            # input profiles
            if (type == 'inputs') {   
              
              df = cbind(x@INPUTS$PROFILE, Generation.kWh = x@SETUP$g, Time = 1:x@INPUTS$HORIZON)
              df = melt(df, id.vars = 'Time')
              
              # construct plot
              plt = ggplot(df, aes(x = Time, y = value)) + 
                geom_point(size = 2.5) + geom_line(size=1.5) + 
                facet_wrap(~variable, ncol = 2, scales = 'free')
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
                ggtitle('Input Profiles') + xlab('Hour of Day')
              
              return(plt)                                        
            }
            
            # input profiles
            if (type == 'aggregates') {   
              
              cntrl = colSums(x@OUTPUT$deterministic$C)
              thrml = colSums(x@OUTPUT$deterministic$L)
              df = data.frame(Time       = 1:x@INPUTS$HORIZON, 
                              Generation = x@SETUP$g, 
                              Baseload   = x@SETUP$b + thrml)
              df = cbind(df, df$Baseload + cntrl)
              names(df)[4] = 'Baseload+Controlled'
              df = melt(df, id.vars = 'Time')
              
              # construct plot
              plt = ggplot(df, aes(x = Time, y = value, color = variable, shape = variable)) + 
                geom_point(size = 2.5) + geom_line(size=1.5)
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
                ggtitle('Aggregate Profiles') + xlab('Hour of Day')
              
              return(plt)            
              
              
            }
          })
            
                  