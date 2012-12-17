              if (constrMC) {
                # define constraints on transition matrix: telescope model
                A = matrix(0, nrow = K, ncol = K)
                for (i in 1:K) {
                  if (i<K) A[i,i+1] = 1
                  A[i,i]   = 1
                  A[i,1]   = 1
                }            
                
                # set constraints on transition matrix in depmix model & fit constrained model
                setpars(mod_unc, value = 1:npar(mod_unc))
                pars <- c(unlist(getpars(fm_unc)))
                conpat = rep(1,npar(mod_unc))
                for (i in 1:K)
                  for (j in 1:K)
                    if (A[i,j] == 0) {
                      idx_constr = K * (i-1) * (length(lifestyl_vars)+1) + 
                                   (1:(length(lifestyl_vars)+1)) * K + j
                      pars[idx_constr] = 0
                      conpat[idx_constr] = 0
                    }                      
                fm <- setpars(mod_unc, pars)
                set.seed(1)
                fm <- fit(fm, equal = conpat, verbose=T)                
              } else 
                
                
                
            # ______________________________________________
            # Clean up data: remove some obvious errors
            
            # remove very unfrequent Proxy IDs
            tab_psa_id = table(raw_data$PSA_ID) / nrow(raw_data)
            infreq_ids = which(tab_psa_id < 0.05)
            if (length(infreq_ids) > 0) {
              idx_rm   = which(raw_data$PSA_ID %in% names(infreq_ids))
              raw_data = raw_data[-idx_rm,]
            }
