# 
# Maarten Speekenbrink, 23-3-2008
# 

viterbi_states <-
function(object, new.data.x, new.data.y, na.allow=TRUE) {
	# returns the most likely state sequence
	nt <- nrow(new.data.x)
	lt <- length(object@ntimes)
	et <- nrow(new.data.x)
	bt <- c(1,et[-lt]+1)		
	ns <- object@nstates
  
  # extract distribution over last state in training set
  pi.0 = object@posterior[nrow(object@posterior),]
  pi.0 = t(t(as.numeric(pi.0[-1])))
	
	for (s in 1:ns) {
    select_vars = colnames(object@response[[s]][[1]]@x)    
    if ('(Intercept)' %in% select_vars) {
      tmp = cbind(X = 1, new.data.x[,setdiff(select_vars, '(Intercept)')])
      colnames(tmp)[1] = '(Intercept)'
    } else tmp = new.data.x[,setdiff(select_vars, '(Intercept)')]
    
    object@response[[s]][[1]]@x = as.matrix(tmp)
	  object@response[[s]][[1]]@y = as.matrix(new.data.y)
	}
	
  # compute density of new data under given model in each state
	nr <- nresp(object)	
	B  <- array(,c(nt,nr,ns))
	R  <- array(,c(nt,nr,ns))
	for(i in 1:ns) {
	  for(j in 1:nr) {
	    B[,j,i] <- dens(object@response[[i]][[j]]) # remove this response as an argument from the call to setpars
	    R[,j,i] <- predict(object@response[[i]][[j]])
	  }
	}
	  
	delta <- psi <- matrix(nrow=nt,ncol=ns)
	state <- vector(length=nt)
	
	prior <- object@init
	
	if(max(ntimes(object)>1)) A <- object@trDens
# 	B <- object@dens
	if(na.allow) B <- replace(B,is.na(B),1)
	B <- apply(B,c(1,3),prod)
	
	for(case in 1:lt) {
		# initialization
		delta[bt[case],] <- prior[case,]*B[bt[case],]
		delta[bt[case],] <- delta[bt[case],]/(sum(delta[bt[case],]))
		psi[bt[case],] <- 0
		# recursion
		if(object@ntimes[case]>1) {
			for(tt in ((bt[case]+1):et[case])) {
				for(j in 1:ns) {
					if(!object@stationary) {
						delta[tt,j] <- max(delta[tt-1,]*(A[tt,j,]))*B[tt,j]
						k <- which.max(delta[tt-1,]*A[tt,j,])
					} else {
						delta[tt,j] <- max(delta[tt-1,]*(A[1,j,]))*B[tt,j]
						k <- which.max(delta[tt-1,]*A[1,j,])
					}
					if(length(k) == 0) k <- 0 # what's this doing here??? can this ever occur? FIX ME
					psi[tt,j] <- k
				}
				delta[tt,] <- delta[tt,]/(sum(delta[tt,]))

			}
		}
		
		# trace maximum likely state
		state[et[case]] <- which.max(delta[et[case],])
		
		# this doesn't need a for loop does it???? FIX ME
		if(object@ntimes[case]>1) {
			for(i in (et[case]-1):bt[case]) {
				state[i] <- psi[i+1,state[i+1]]
			}
		}
	}
	
  # compute predicted fit using Viterbi-decoded sequence
	resp_fit     = sapply(1:length(state), function(j) R[j,,state[j]])
	resp_fit_avg = sapply(1:length(state), function(j) sum(delta[j,] * R[j,,]))
	
  # compute predicted fit using Markov-propagated state distribution
	P0 = t(A[1,,])
  P  = P0
  pi.cur = pi.0
  s      = which.max(pi.cur)
  pred = rep(0, nrow(new.data.x))
  for (t in 1:nrow(new.data.x)){
    #     pi.cur = t(t(pi.cur) %*% P)
    pi.cur = P[s,]
    pred[t]= R[t,,s]
    s      = which.max(pi.cur)    
    P      = P %*% P0
  }
  
	delta <- data.frame(state,delta, fit = resp_fit_avg, fit.max = resp_fit) 	
	return(delta)
}

