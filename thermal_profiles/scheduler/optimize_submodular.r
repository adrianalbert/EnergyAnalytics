# #########################################################################
# optimize_submodular.r
# -------------------------
#
# Implements greedy algorithm to optimize submodular non-monotone function.
# 
# Adrian Albert
# Last modified: July 2014.
# #########################################################################

# Omega:   set of elements (list)
# U_list:  set of sets of schedules for elements in Omega (list)
# f(A, U): function to optimize with parameters
#             A: set of elements
#             U: set of actions for each element in A
optimize_submodular = function(f, Omega, U_list, eps = 0.001) {
  
  # define marginal payoff function
  rho = function(A, U, e, u) {
    A_new = c(A, e); U_new = c(U, u)
    res = f(A_new, U_new) - f(A, U)
    return(res)
  }
  
  Y = U_list
  
  # initialize sets
  A = list();
  U = list()
  d = as.data.frame(t(sapply(names(Omega), function(e) {
    Ue = U_list[[e]]
    r  = sapply(1:length(Ue), function(j) rho(A, U, Omega[e], Ue[j]))
    rs = which.max(r)
    return(c(rho = r[rs], u.idx = rs))
  })))
  
  e         = which.max(d[,'rho']); 
  A         = Omega[e]
  Ue        = U_list[[e]]
  U         = Ue[d[e,'u.idx']]
  Omega_cur = Omega[-e]
  Y_cur     = Y[-e]
  
  # main loop
  iter = 0
  repeat {
    
    n = length(A)
    terminate = FALSE
    iter = iter + 1
    
    cat(paste('Iteration', iter,'...\n'))
    
    # step forward
    cat(paste("Stepping forward...|A|=", length(A), '\n'))
    
    # find element in Omega' for which f(c(A, ae), c(U, u)) > (1 + eps/n^2) * f(A, U)
    ide = sapply(names(Omega_cur), function(e) {
      ae = Omega_cur[e]
      Ue = Y_cur[[e]]
      id = sum(sapply(1:length(Ue), function(j) f(c(A, ae), c(U, Ue[j])) > (1 + eps/n^2) * f(A, U)))
      return(id)
    })
    ide = names(Omega_cur)[which(ide > 0)]
    
    if (length(ide) > 0) {
      d = as.data.frame(t(sapply(names(Omega_cur[ide]), function(e) {
        Ue = Y_cur[[e]]
        r  = sapply(1:length(Ue), function(j) rho(A, U, Omega_cur[e], Ue[j]))
        rs = which.max(r)
        return(c(rho = r[rs], u.idx = rs))
      })))       
      es = which.max(d[,'rho'])
      e  = ide[es]; 
      ae = Omega_cur[e]; 
      Ue = Y_cur[[e]]
      u  = Ue[d[e,'u.idx']]
      
      A = c(A, ae)
      U = c(U, u)
      Y_cur = Y_cur[-es]   
      Omega_cur = Omega_cur[-es]
    } else {    
      
      # step backward 
      cat(paste("Stepping backward...|A|=", length(A), '\n'))

      # find element in A for which f(A[-e], U[-e]) > (1 + eps/n^2) * f(A, U))
      idf = which(sapply(1:length(A), function(e) f(A[-e], U[-e]) > (1 + eps/n^2) * f(A, U)))
      if (length(idf)>0) {
        d = sapply(1:length(A), function(e) -rho(A[-e], U[-e], Omega_cur[e], U[e]))
        e = which.max(d)
        e_name = names(A)[e]
        A = A[-e]
        U = U[-e]
        Y_cur = c(Y_cur, U_list[e_name]); names(Y_cur)[length(Y_cur)] = e_name
        Omega_cur = c(Omega_cur, Omega[e_name]); names(Omega_cur)[length(Omega_cur)] = e_name
      } else terminate = TRUE 
    }
    
    # terminate condition
    if (terminate | length(Y_cur) == 0) break
    
  }
  return(list(A = A, U = U, val = f(A,U)))  
}

test = FALSE

if (test) {
    
  rm(list = setdiff(ls(), 'optimize_submodular'))
  # options(error = recover)
  
  # generate fake data
  g = runif(24)
  N = 100;   
  zeta = sample(2:6, N, replace = T)
  usr_names = paste('Name', 1:N)
  Omega = lapply(usr_names, function(s) runif(24))
  names(Omega) = usr_names
  U_list= lapply(zeta, function(z) {
    eta   = sample(1:23, z, replace = T)
    gamma = sample(1:5, z, replace = T)
    U     = lapply(1:z, function(s) {
      u = rep(0, 24)
      u[eta[s]:min(eta[s] + gamma[s], 24)] = 1
      return(u)
    })
    return(U)
  })
  names(U_list) = usr_names
  
  # define submodular function
  f = function(A, U) {
    if (length(A) == 0) return(sum(g^2))
    D = sapply(1:length(A), function(i) A[[i]] * U[[i]])
    C = sum(g^2) - sum((colSums(D) - g)^2)
    return(C)
  }
  
  # optimize submodular function
  Rprof()
  res = optimize_submodular(f, Omega, U_list)
  Rprof(NULL)
  tmp = summaryRprof()
}