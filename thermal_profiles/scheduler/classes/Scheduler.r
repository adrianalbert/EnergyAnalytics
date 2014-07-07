# #########################################################################
# Scheduler.r
# -------------------------
#
# Interface to implement scheduling problems for thermally-sensitive users.
# 
# Adrian Albert
# Last modified: July 2014.
# #########################################################################

library('timeDate')
library('lubridate')
library('reshape')
library('CVXfromR')
library('ggplot2')
library('nloptr')
library('Matrix')
library('quadprog')
library('bigmemory')
library("bigalgebra")
library('Matrix')
library('Rmosek')

# clean-up previous definitions of methods for class Person
removeClass('Scheduler')

# ________________________
# Class definition

setClass(
  Class = "Scheduler",
  representation = representation(
    SETUP        = "list",             # variables for problem setup
    OUTPUT       = "list",             # problem output
    SOLVER       = "list"              # options for optimization solver
  )
)

# _____________________________________________________
# Constructor method for class Scheduler

setMethod(f = "initialize", 
          signature = "Scheduler",
          definition = function(.Object, inputs, setup = list(), verbose = T) {
            
            if (verbose) cat('*** Initializing Scheduler ***\n')
            
            # parse inputs into form required by optimization
            .Object@SETUP$tau        = nrow(setup$goal)
            .Object@SETUP$N          = length(inputs) 
            .Object@SETUP$DT         = setup$DT
            .Object@SETUP$g          = setup$goal$Goal   
            .Object@SETUP$q          = setup$tou_rates               
                                    
            # response rates and baselines & scale rsponse by DT
            .Object@SETUP$Abar       = do.call('rbind', lapply(inputs, function(l) as.numeric(l$a$mu[,1]  * .Object@SETUP$DT)))            
            rownames(.Object@SETUP$Abar)= names(inputs)
            .Object@SETUP$W          = lapply(inputs, function(l) l$a$covmat * .Object@SETUP$DT^2)
            
            return(.Object)
          })

# _________________________________
# Method to print scheduler object

setMethod('show',
          signature  = 'Scheduler',
          definition = function(object){
            cat('*** Scheduler Object ***\n')
            
            # basic info
            cat(sprintf('No. Users:   %d\n', length(unique(object@N))))            
            
            cat('*** END: Scheduler Object ***\n')
          })

# __________________________________________________________
# Method to set up & solve full scheduling problem

setGeneric(
  name = "solveSchedules",
  def = function(.Object, options = NULL, verbose = T){standardGeneric("solveSchedules")}
)
setMethod('solveSchedules',
          signature  = 'Scheduler',
          definition = function(.Object, options = NULL, verbose = T) {
            
            # access model parameters and data
            N = .Object@SETUP$N
            Abar = .Object@SETUP$Abar
            tau  = .Object@SETUP$tau
            W    = .Object@SETUP$W
            w    = t(sapply(W, function(w) diag(w)))
            g    = .Object@SETUP$g
            q    = .Object@SETUP$q
            Q    = diag(q)
            if (length(options$budget) == 1) beta = rep(options$budget, N) else beta = options$budget
            gamma= options$gamma
            npars = N*tau
            if (is.null(options$Nr)) Nr = rep(1,N) else Nr = options$Nr
            if (!is.null(options$presaved)) presaved = options$presaved else presaved = FALSE
              
            # do not consider negative rates
            idx0 = which(Abar < 0)
            if (length(idx0)>0) Abar[idx0] = 0
            
            # form optimization variables            
            if (presaved & file.exists('presaved.RData')) {
              cat('Loading presaved objects...\n')
              load('presaved.RData')
            } else {
              cat('Creating QP objects...\n')
              for (i in 1:N){
                cat(paste('i=', i, ','))
                t0 = proc.time()
                for (j in 1:N){
                  if (i == j) {
                    tmp = Q %*% diag(w[j,] + Abar[j,]^2)
                  } else {
                    tmp = diag(Abar[i,]) %*% Q %*% diag(Abar[j,]) 
                  }
                  if (j == 1) row = Matrix(tmp) else row = cBind(row, tmp)
                }
                dt = proc.time() - t0
                print(dt)
                if (i == 1) H = Matrix(row) else H = rBind(H, row)
              }
              H = (H + t(H)) / 2
              Amean = Matrix(diag(Abar[1,])); for (i in 2:N) Amean = cBind(Amean, diag(Abar[i,]));            
              dvec  = -as.numeric(t(g) %*% Q %*% Amean) 
              Amat = kronecker(Matrix(diag(N)), Matrix(t(rep(1,tau))))           # (u_i)^T 1 < beta
              save(file = 'presaved.RData', list = c('H', 'Amat', 'Amean', 'dvec'))
            }
            cat('... done!\n')
            
            # form constraints
            lbnd = rep(0, npars)
            ubnd = rep(Nr, each = tau)
            if (length(beta)==1) bvec = rep(beta, N) else bvec = beta
            bvec = bvec * Nr
            
            # Choleski decomposition of H
            R = chol(2*H); R1 = solve(R)

            # solve QP and return results
            cat('Solving QP...')
            
            # Construct and solve problem 
            prob <- mosek_qptoprob(R, 2*dvec, Amat, bvec, NA, NA, lbnd, ubnd);
            fit  <- mosek(prob);
            
            # retrieve variables
            status <- fit$sol$itr$solsta
            coefs  <- fit$sol$itr$xx
            u      <- coefs[1:(npars)]# - coefs[(npars+1):(2*npars)]
            u      <- matrix(u, nrow = N, byrow = TRUE)
            resid  <- coefs[(2*(npars) + 1):(2*(npars) + npars)]            
            val = as.numeric(prob$c %*% coefs + t(g) %*% Q %*% g)
            
            cat('done!\n')
                        
            # compute aggregate profiles
            Delta.bar = colSums(Abar * u)            
            Delta.Var = matrix(0, nrow = tau, ncol = tau)
            for (i in 1:length(W)) {
              Delta.Var = Delta.Var + diag(u[i,]) %*% W[[i]] %*% diag(u[i,])
            }
            
            # normalize effort per each user
            nr= round(rowSums(u / beta))          
            u = u / matrix(rep(nr, each = ncol(u)), nrow = nrow(u), byrow=T)   
            i = which(!is.finite(u) | is.nan(u))
            if (length(i)>0) u[i] = 0
            
            # pass solution further
            .Object@OUTPUT$EC        = val
            .Object@OUTPUT$U         = u
            .Object@OUTPUT$Delta.bar = Delta.bar
            .Object@OUTPUT$Delta.Var = Delta.Var
            .Object@OUTPUT$nr        = nr
            .Object@OUTPUT$Nr        = Nr
            
            rm(list = c('H', 'R', 'R1', 'Amat')); gc()
            
            return(.Object)
          })

# _________________________
# Method to plot schedules 

setMethod('plot',
          signature  = 'Scheduler',
          definition = function(x, type = 'effort-profiles', selected = NULL, compare = NULL){
            
            if (is.null(selected)) selected = 1:x@SETUP$N
            if (length(selected) == 1) selected = sample(x@SETUP$N, selected) 
            
            if (is.character(selected)) selected = which(rownames(x@SETUP$Abar) %in% selected)
            
            U = as.data.frame(x@OUTPUT$U[selected,])
            U$name = rownames(x@SETUP$Abar[selected,])
            nr= x@OUTPUT$nr[selected]
            Nr= x@OUTPUT$Nr[selected]
            
            # effort profiles
            if (type == 'effort-profiles') {   
              
              U$name = sprintf("%02s", U$name)
              
              df = melt(U, id.vars = 'name')              
              df$variable = as.numeric(df$variable)
              
              # construct plot
              plt = ggplot(df, aes(y = value, x = variable)) + 
                geom_point(size = 2.5) + geom_line(size=1.5)
              plt = plt + facet_wrap(~name, ncol = 4)
              plt = plt + scale_x_continuous(limits = c(1, 24))
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=22),
                      axis.text.y      = element_text(size=20), 
                      axis.text.x      = element_text(size=20),
                      axis.title.y     = element_text(size=20),
                      axis.title.x     = element_text(size=20),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      legend.position  = 'none',
                      axis.ticks = element_blank()) + 
                ggtitle(paste("Effort Schedules")) + ylab('Requested Effort [deg F]') + xlab('Hour of Day')
              
              return(plt)                                      
            }            
            
            # goal match
            if (type == 'goal-match') {   
              
              df = data.frame(value = x@OUTPUT$Delta.bar, Variance = sqrt(abs(diag(x@OUTPUT$Delta.Var))))
              dh = data.frame(value = x@SETUP$g, Variance = 0)              
              df$Hour = 1:nrow(df)
              df$variable = 'Delta'
              dh$variable = 'Goal'
              dh$Hour = 1:nrow(dh)
              df = rbind(df, dh)
              if (!is.null(compare)) {
                dc = data.frame(value = compare$mu, Variance = sqrt(abs(diag(compare$var))))
                dc$variable = 'Actual Delta'
                dc$Hour = 1:nrow(dc)
                df = rbind(df, dc)
              }
              
              # construct plot
              plt = ggplot(df, aes(y = value, x = Hour, color = variable)) + 
                geom_point(size = 2.5) + geom_line(size=1.5)
              plt = plt + geom_errorbar(aes(ymax = value + Variance, ymin = value - Variance))        
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
                      legend.title      = element_text(size=18),
                      legend.position  = c(0.25, 0.7),
                      axis.ticks = element_blank()) + 
                ggtitle(paste("Matching the reductions goal profile")) + ylab('kWh') + xlab('Hour of Day')
              
              return(plt)                                      
            }
            
            # goal match
            if (type == 'number-selected') {   
              
              df = data.frame(Available = Nr, Selected = nr, Consumer = sprintf('%02d', 1:length(Nr)))
              df = melt(df, id.vars = 'Consumer')
              
              # construct plot
              plt = ggplot(df, aes(y = value, x = Consumer, color = variable, fill = variable)) + 
                geom_bar(stat = 'identity', position = 'dodge', width = 0.8)
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18, angle = -20),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      legend.title      = element_text(size=18),
                      legend.position  = c(0.65, 0.8),
                      axis.ticks = element_blank()) + 
                ggtitle(paste("Selected Consumers")) + ylab('No. Consumers') + xlab('Consumer')
              
              return(plt)                                      
            }
            
            # effort distribution
            if (type == 'effort-distribution') {   
              
              id = which(rowSums(U[,-ncol(U)])>0)
              df = melt(U[id,-ncol(U)])
              
              # construct plot
              plt = ggplot(df, aes(y = value, x = variable)) + 
                geom_boxplot(size = 1.5)
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18, angle = -20),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),            
                      legend.text      = element_text(size=18),
                      legend.title      = element_text(size=18),
                      legend.position  = c(0.65, 0.8),
                      axis.ticks = element_blank()) + 
                ggtitle(paste("Effort distribution")) + ylab('Effort') + xlab('Hour of day')
              
              return(plt)                                      
            }
          })
