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
library('nloptr')
library('Matrix')

source('../clustering/kError.r')
# source('./init_objective.r')

# clean-up previous definitions of methods for class Person
removeClass('Scheduler')

# ________________________
# Class definition

setClass(
  Class = "Scheduler",
  representation = representation(
    SETUP        = "list",             # variables for problem setup
    SEGMENTS     = "list",             # segments
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
            .Object@SETUP$O          = setup$cov.mat
            
            # response rates and baselines & scale rsponse by DT
            .Object@SETUP$Abar       = do.call('rbind', lapply(inputs, function(l) as.numeric(l$a$mu[,1]  * .Object@SETUP$DT)))            
            rownames(.Object@SETUP$Abar)= names(inputs)
            .Object@SETUP$Bbar       = do.call('rbind', lapply(inputs, function(l) as.numeric(l$b$mu[,1])))
            rownames(.Object@SETUP$Bbar)= names(inputs)
            .Object@SETUP$W          = lapply(inputs, function(l) l$a$covmat * .Object@SETUP$DT^2)
            .Object@SETUP$V          = lapply(inputs, function(l) l$b$covmat)
            .Object@SETUP$sigma      = do.call('rbind', lapply(inputs, function(l) as.numeric(l$sigma)))
            rownames(.Object@SETUP$sigma)= names(inputs)
            
            # solver parameters
            .Object@SOLVER$cvx.setup.dir = setup$cvx.setup.dir    

            # if clustering is requested
            if (!is.null(setup$noClusters)) {
              cat('Segmentation request detected...\n')
              .Object@SETUP$K = setup$noCluster
              .Object@SOLVER$max.iter = setup$max.iter
            }
            
            return(.Object)
          })

# _________________________
# Method to segment users

setGeneric(
  name = "segmentProfiles",
  def = function(.Object, verbose = T){standardGeneric("segmentProfiles")}
)
setMethod('segmentProfiles',
          signature  = 'Scheduler',
          definition = function(.Object, verbose = T) {
            
            if (is.null(.Object@SETUP$K)) return(0) 
              
            if (length(.Object@SETUP$K) == 1) {
              Kmin = .Object@SETUP$K
              Kmax = .Object@SETUP$K
            } else {
              Kmin = .Object@SETUP$K[1]
              Kmax = .Object@SETUP$K[2]
            }
              
            W = lapply(.Object@SETUP$W, function(w) abs(w))
            if (Kmin <= Kmax) {
              # compute K-means solutions in parallel for model (K) selection
              res = lapply(Kmin:Kmax, #mc.cores = 5, 
                           function(k){
                fit = kError(.Object@SETUP$Abar, W, k, iter = .Object@SOLVER$max.iter)
                return(fit)              
              })                            
              .Object@SEGMENTS = res
            } else return(0)
                        
            return(.Object)
            
          })

# # __________________________________________________________
# # Method to set up & solve deterministic scheduling problem
# 
# # Puts the problem into CVX and calls matlab. 
# 
# setGeneric(
#   name = "solveSchedules",
#   def = function(.Object, verbose = T){standardGeneric("solveSchedules")}
# )
# setMethod('solveSchedules',
#           signature  = 'Scheduler',
#           definition = function(.Object, verbose = T) {
#             
#             if (is.null(.Object@SOLVER$cvx.setup.dir)) return(0)
#             
#             # form objective and constraints
#             Abar = .Object@SETUP$Abar
#             N    = .Object@SETUP$N
#             tau  = .Object@SETUP$tau
#             w    = sapply(.Object@SETUP$W, function(w) diag(w))
#             W    = array(dim = c(N, tau, tau))
#             for (i in 1:N) W[i,,] = .Object@SETUP$W[[i]]
#             g    = .Object@SETUP$g
#             q    = .Object@SETUP$q
#             Q    = diag(q)
#             beta = rep(.Object@SETUP$beta, N)
#             O    = list()
#             for (t in 1:tau) O[[1+length(O)]] = matrix(5e-5*(-1+2*runif(N^2)), nrow=N, ncol=N)
#             #             O    = list()
#             #             for (t in 1:tau) O[[1+length(O)]] = matrix(0, nrow=N, ncol=N)
#             
#             print(bla)
#             
#             f0 = f_reg(matrix(0, nrow = N, ncol = tau), g, q, Abar, O, w)   
#             # objective function
#             f_U = function(u, beta) {
#               U = matrix(u, ncol = tau, byrow = T)
#               #f_reg(U, g, q, Abar, O, w) 
#               f_obj(U, g, q, Abar, O, w)          
#             }
#             # gradient
#             g_U = function(u, beta) {
#               U  = matrix(u, ncol = tau, byrow = T)
#               gU = f_grad(U, g, q, Abar, O, w)            
#               #gU = f_grad(U, g, q, Abar, w)            
#               gU = as.numeric(t(gU))
#             }
#             
#             # bound constraints
#             lb = rep(0, N*tau) 
#             ub = rep(1, N*tau)
#             # inequality constraints
#             constr_U = function(u, beta) {
#               U = matrix(abs(u), ncol = tau, byrow=T)
#               return(as.numeric(U %*% as.matrix(rep(1, tau)) - beta))
#             }
#             # inequality constraints gradient
#             constr_U_grad = function(u, beta) {
#               idx.1 = matrix(1:(N*tau), ncol = tau, byrow=T)
#               J = matrix(0, ncol = N*tau, nrow = N)
#               for (i in 1:N) J[i,idx.1[i,]] = 1
#               return(J)
#             }
#             
#             # initial guesses
#             U0   = 0.1 * matrix(runif(N*tau), ncol = tau)
#             u0  = c(t(U0))
#             z   = f_U(u0, beta)
#             cz  = constr_U(u0, beta)
#             czz = constr_U_grad(u0, beta)
#             
#             # run nonlinear optimizer
#             # 
#             opts <- list(algorithm="NLOPT_LD_MMA",
#                          xtol_rel=1.0e-6,
#                          print_level = 1,
#                          maxeval = 50,
#                          check_derivatives = F)
#             res <- nloptr( x0          = u0,
#                            eval_f      = f_U,
#                            eval_grad_f = g_U,
#                            lb          = lb,
#                            ub          = ub,
#                            eval_g_ineq = constr_U,
#                            eval_jac_g_ineq = constr_U_grad,
#                            opts        = opts,
#                            beta        = beta)
#             U = matrix(res$solution, ncol = tau, byrow = T)
#             
#             #             .Object@OUTPUT$original$U = cbind(UID = .Object@SETUP$UID, as.data.frame(opt$U))
#             #             .Object@OUTPUT$original$L = .Object@SETUP$A * .Object@SETUP$TF
#             #             .Object@OUTPUT$original$C = A.hat * opt$U
#             return(.Object)
#           })

# __________________________________________________________
# Simplified scheduling using clusters

setGeneric(
  name = "solveSchedulesClusters",
  def = function(.Object, options = NULL, verbose = T){standardGeneric("solveSchedulesClusters")}
)
setMethod('solveSchedulesClusters',
          signature  = 'Scheduler',
          definition = function(.Object, options = NULL, verbose = T) {
            
            if (is.null(.Object@SOLVER$cvx.setup.dir)) return(0)
            
            # form objective and constraints
            Abar = .Object@SEGMENTS[[1]]$centers
            N    = nrow(Abar)
            Nr   = table(.Object@SEGMENTS[[1]]$assignment)
            tau  = ncol(Abar)
            w    = t(sapply(.Object@SEGMENTS[[1]]$error, function(w) diag(w)))
            W    = .Object@SEGMENTS[[1]]$error
            g    = .Object@SETUP$g
            q    = .Object@SETUP$q
            Q    = diag(q)
            .Object@SETUP$beta  = options$budget
            .Object@SETUP$gamma = options$gamma
            
            # simulation parameters
            beta = Nr * options$budget
            gamma= options$gamma
                      
            # expected cost
            compute_EC = function(Abar, W, U) {
              Delta.bar = colSums(Abar * U)
              Delta.Var = matrix(0, nrow = tau, ncol = tau)
              for (i in 1:length(W)) {
                Delta.Var = Delta.Var + diag(U[i,]) %*% W[[i]] %*% diag(U[i,])
              }
              EC = sum(diag(Q %*% Delta.Var)) + t(Delta.bar - g) %*% Q %*% (Delta.bar - g)
              # add L1 penalty
              EC = EC + gamma * sum(abs(U))
              return(as.numeric(EC))
            }
            
            # objective function
            f_U = function(u, beta) {
              U = matrix(u, ncol = tau, byrow = T)
              compute_EC(Abar, W, U) 
            }
                        
            # bound constraints
            lb = rep(0, N*tau) 
            ub = rep(as.numeric(Nr), each = tau)
            # inequality constraints
            constr_U = function(u, beta) {
              U = matrix(abs(u), ncol = tau, byrow=T)
              return(as.numeric(t(U %*% as.matrix(rep(1, tau))) - t(as.matrix(beta))))
            }
            
            # print(bla)
            
            # initial guesses
            U0   = min(Nr) * options$budget *  matrix(0.1, ncol = tau, nrow = N)                      
            u0  = c(t(U0))
            u0[sample(1:length(u0), N)] = 0
            constr_U(u0, beta)
            constr_U(u0, beta)
            
            # run nonlinear optimizer
            # NLOPT_LN_BOBYQA
            opts <- list(algorithm="NLOPT_LN_COBYLA",
                         xtol_rel=1.0e-6,
                         print_level = 0,
                         maxeval = 500,
                         check_derivatives = F)
            res <- nloptr( x0          = u0,
                           eval_f      = f_U,
                           lb          = lb,
                           ub          = ub,
                           eval_g_ineq = constr_U,
                           opts        = opts,
                           beta        = beta)
            
            # extract and form solution
            U = matrix(res$solution, ncol = tau, byrow = T)
            nr= apply(U, 1, function(x) trunc(sum(x) / .Object@SETUP$beta))
            
            # compute solution profiles
            Delta.bar = colSums(Abar * U)            
            Delta.Var = matrix(0, nrow = tau, ncol = tau)
            for (i in 1:length(W)) {
              Delta.Var = Delta.Var + diag(U[i,]) %*% W[[i]] %*% diag(U[i,])
            }
            
            # normalize effort per class
            U = U / matrix(rep(nr, each = ncol(U)), nrow = nrow(U), byrow=T)                                           
            
            # pass solution further
            .Object@OUTPUT$segments$EC = res$objective
            .Object@OUTPUT$segments$U  = U
            .Object@OUTPUT$segments$nr = nr
            .Object@OUTPUT$segments$Delta.bar = Delta.bar
            .Object@OUTPUT$segments$Delta.Var = Delta.Var
            return(.Object)
          })

# _________________________
# Method to plot schedules 

setMethod('plot',
          signature  = 'Scheduler',
          definition = function(x, type = 'heatmap', selected = NULL){
            
            if (type == 'cluster-centers' & length(x@SEGMENTS) > 0) {
              
              if (is.null(selected)) selected = 1
              res = x@SEGMENTS[[selected]]
              
              # extract error diagonals
              sd = as.data.frame(t(sapply(res$errors, function(x) sqrt(abs(diag(x))))))
              names(sd) = 1:ncol(sd)
              sd$Segment = 1:nrow(sd)
              sd = melt(sd, id.vars = 'Segment')
              names(sd)[3] = 'sd'

              # extract centers
              df = as.data.frame(res$centers)                            
              names(df) = 1:ncol(df)
              df$Segment = 1:nrow(df)
              df = melt(df, id.vars = 'Segment')
              names(df)[3] = 'mu'
              
              # form plotting data
              dfp = merge(df, sd, by = c('Segment', 'variable'))
              dfp$variable = as.numeric(dfp$variable)
              dfp$Segment = as.factor(dfp$Segment)
                            
              # percentages of membership              
              tab = table(res$assignment)
              tab = round(tab / sum(tab), digits = 4)
              tab = tab * 100
              dfp$Segment = paste(as.numeric(dfp$Segment), ': ', tab[as.numeric(dfp$Segment)], '%', sep = '')
              
              plt = ggplot(dfp, aes(y = mu, x = variable)) + 
                geom_point(size = 2.5) + geom_line(size=1.5)
              plt = plt + geom_errorbar(aes(ymax = mu + sd, ymin = mu - sd))        
              plt = plt + facet_wrap(~Segment, ncol = 3, scales = 'free')
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
                      legend.position  = "none",
                      axis.ticks = element_blank()) + 
                ggtitle(paste("Segments")) + ylab('Response [kWh/F]') + xlab('Time of day')
              
              return(plt)
            }
            
            # effort profiles
            if (type == 'effort-profiles') {   
              
              if (is.null(selected)) selected = 1
              seg = x@SEGMENTS[[selected]]
              eff = x@OUTPUT$segments
              
              U = as.data.frame(eff$U)
              nr= eff$nr
#               U = U / matrix(rep(nr, each = ncol(U)), nrow = nrow(U), byrow=T)                                           
              names(U)  = 1:ncol(U)
              U$Segment = paste(1:nrow(U), '(Selected: ', nr[1:nrow(U)], ')', sep = '')
              df = melt(U, id.vars = 'Segment')              
              df$variable = as.numeric(df$variable)
              
              # construct plot
              plt = ggplot(df, aes(y = value, x = variable, color = Segment)) + 
                geom_point(size = 2.5) + geom_line(size=1.5)
              plt = plt + facet_wrap(~Segment, ncol = 3)
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
                      legend.position  = 'none',
                      axis.ticks = element_blank()) + 
                ggtitle(paste("Typical Schedules")) + ylab('Requested Effort [deg F]') + xlab('Hour of Day')
              
              return(plt)                                      
            }            
            
            # goal match
            if (type == 'goal-match') {   
              
              if (is.null(selected)) selected = 1
              seg = x@SEGMENTS[[selected]]
              eff = x@OUTPUT$segments
              
              df = data.frame(value = eff$Delta.bar, Variance = sqrt(abs(diag(eff$Delta.Var))))
              dh = data.frame(value = x@SETUP$g, Variance = 0)
              df$Hour = 1:nrow(df)
              df$variable = 'Delta'
              dh$variable = 'Goal'
              dh$Hour = 1:nrow(dh)
              df = rbind(df, dh)
              
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
                      legend.position  = c(0.15, 0.8),
                      axis.ticks = element_blank()) + 
                ggtitle(paste("Matching the reductions goal profile")) + ylab('kWh') + xlab('Hour of Day')
              
              return(plt)                                      
            }
          })
            
                  
# # form CVX problem              
# cvxcode <- paste("variables U(N,tau)",
#                  "minimize(  )",
#                  "subject to",
#                  "abs(U) * ones(tau,1) <= beta",
#                  "-1 <= U <= 1",
#                  sep=";")
# const.vars = list(p      = .Object@SETUP$p, 
#                   N      = .Object@SETUP$N, 
#                   A.T    = A.T, 
#                   A.hat  = A.hat, 
#                   g      = .Object@SETUP$g, 
#                   b      = .Object@SETUP$b, 
#                   BUDGET = .Object@SETUP$BUDGET)
# opt <- CallCVX(cvxcode, 
#                const.vars=const.vars,
#                opt.var.names="U", setup.dir=.Object@SOLVER$cvx.setup.dir)

# define objective function
# idx = expand.grid(i = 1:N, j = 1:N)
# f_obj = function(U) {
#   obj = mclapply(1:tau, mc.cores = 5,
#                  function(t) {
#                    res = sapply(1:length(idx), function(k){
#                      i = idx[k,'i']; j = idx[k,'j'];
#                      if (i == j) v = w[i,t] else v = 0#omega[[i]][[j]][t,t]
#                      r = U[i,t]*U[j,t]*( v + (Abar[i,t] - g[t])^2) 
#                      return(r)
#                    })
#                    return(sum(unlist(res)))
#                  })
#   return(unlist(obj))              
# }
# 
# # define gradient function
