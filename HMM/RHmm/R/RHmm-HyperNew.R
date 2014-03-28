 ###############################################################
 #### RHmm package                             
 ####                                                         
 #### File: RHmm-HyperNew.R 
 ####                                                         
 #### Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr>
 #### Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ####                                                         
 ###############################################################

tolMin <- .Machine$double.eps*100
# Contraintes sur les paramètres
eProba <- 1
eVar <- 2
eCor <- 3
eNone <- 0

sumList <- function(List, n)
{
    if (is.list(List))
    {   Res <- 0
        for (i in 1:n)
            Res <- Res + List[[i]]
    }
    else
        Res <- List
    return(Res)
}

GetNParam<-function(object) UseMethod("GetNParam")
GetNAllParam<-function(object) UseMethod("GetNAllParam")
GetVectParam <- function(object) UseMethod("GetVectParam")
GetVectAllParam <- function(object) UseMethod("GetVectAllParam")
SetVectParam <- function(object, Vect) UseMethod("SetVectParam")
SetVectAllParam <- function(object, Vect) UseMethod("SetVectAllParam")
GetNConstraint <- function(object) UseMethod("GetNConstraint")
GradConstraint <- function(object) UseMethod("GradConstraint")

GetNParam.default<-function(object)
{
    return(list(nParam=0,paramConstr=NULL))
}

GetNAllParam.default<-function(object)
{
    return(list(nParam=0, paramConstr=NULL))
}

GetVectParam.default<-function(object)
{
    return(NULL)
}

GetVectAllParam.default<-function(object)
{
    return(NULL)
}

SetVectParam.default<-function(object, Vect)
{
    return(NULL)
}


SetVectAllParam.default<-function(object, Vect)
{
    return(NULL)
}

GetNConstraint.default<-function(object)
{
    return(0)
}

GradConstraint.default<-function(object)
{
    return(NULL)
}


GetNParam.HMMFitClass<-function(object)
{
    return(GetNParam(object$HMM))
}

GetNAllParam.HMMFitClass<-function(object)
{
    return(GetNAllParam(object$HMM))
}

GetVectParam.HMMFitClass<-function(object)
{
    return(GetVectParam(object$HMM))
}

GetVectAllParam.HMMFitClass<-function(object)
{
    return(GetVectAllParam(object$HMM))
}

SetVectParam.HMMFitClass<-function(object, Vect)
{
    return(SetVectParam(object=object$HMM, Vect=Vect))
}

SetVectAllParam.HMMFitClass<-function(object, Vect)
{
    return(SetVectAllParam(object=object$HMM, Vect=Vect))
}

GetNConstraint.HMMFitClass<-function(object)
{
    return(GetNConstraint(object$HMM))
}

GradConstraint.HMMFitClass<-function(object)
{
    return(GradConstraint(object$HMM))
}


GetNParam.HMMClass<-function(object)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nStates <- object$distribution$nStates
    nOtherParam <- (nStates - 1)*(nStates + 1)
    Res1<-GetNParam(object$distribution)
    nParam <- Res1$nParam * nStates + nOtherParam
    paramConstr <- c(rep(eProba, nOtherParam), rep(Res1$paramConstr, nStates))
    Res <- list(nParam=nParam, paramConstr=paramConstr)
    return(Res)
}


GetNAllParam.HMMClass<-function(object)
{   nStates <- object$distribution$nStates
    Res1<-GetNAllParam(object$distribution)
    nOtherParam <- nStates*(nStates+1)
    nParam <- Res1$nParam * nStates + nOtherParam
    paramConstr <- c(rep(eProba, nOtherParam), rep(Res1$paramConstr, nStates))
    Res <- list(nParam=nParam, paramConstr=paramConstr)
    return(Res)
}

GetVectParam.HMMClass<-function(object)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nParam <- GetNParam(object)$nParam
    nStates <- object$distribution$nStates
    Res <- rep(0, nParam)
    #initProb
    for (i in 1:(nStates-1))
    {   Res[i] <- object$initProb[i]
    }
    k <- nStates
    for (i in 1:nStates)
    {   for (j in 1:(nStates -1))
        {   Res[k] <- object$transMat[i,j]
            k <- k + 1
        }
    }
    Res[k:nParam] <- GetVectParam(object$distribution)
    return(Res)
}

GetVectAllParam.HMMClass<-function(object)
{
    nStates <- object$distribution$nStates
    Res <- object$initProb

    for (i in 1:nStates)
        Res <- c(Res, object$transMat[i,])
    Res1 <- GetVectAllParam(object$distribution)
    return(c(Res, Res1))
}

SetVectParam.HMMClass<-function(object, Vect)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nStates <- object$distribution$nStates

    #initProb
    initProb <- rep(0, nStates)
    initProb[1:(nStates-1)] <- Vect[1:(nStates-1)]
    initProb[nStates] <- 1.0 - sum(initProb)

    #transMat
    transMat <- matrix(0, nrow=nStates, ncol=nStates)
    k <- nStates
    for (i in 1:nStates)
    {   transMat[i,1:(nStates-1)] <- Vect[k:(k+nStates-2)]
        k <- k + nStates-1
        transMat[i,nStates] <- 1.0 - sum(transMat[i,])
    }

    ResDistr <- SetVectParam(object=object$distribution, Vect=Vect[k:length(Vect)])
    HMM<-HMMSet(initProb=initProb, transMat=transMat, ResDistr)
    return(HMM)
}

SetVectAllParam.HMMClass<-function(object, Vect)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nStates <- object$distribution$nStates

    #initProb
    initProb <- Vect[1:nStates]

    #transMat
    transMat <- matrix(0, nrow=nStates, ncol=nStates)
    k <- nStates+1
    for (i in 1:nStates)
    {   transMat[i,1:nStates] <- Vect[k:(k+nStates-1)]
        k <- k + nStates
    }

    ResDistr <- SetVectAllParam(object=object$distribution, Vect=Vect[k:length(Vect)])
    HMM<-HMMSet(initProb=initProb, transMat=transMat, ResDistr)
    return(HMM)
}

GetNConstraint.HMMClass<-function(object)
{
    nStates <- object$distribution$nStates
    nConstr <- 1 + nStates + GetNConstraint(object=object$distribution)
    return(nConstr)
}

GradConstraint.HMMClass<-function(object)
{
    LParam <- GetNAllParam(object)
    nConstr <- GetNConstraint(object)
    Res <- matrix(0, nrow=nConstr, ncol=LParam$nParam)
    nStates <- object$distribution$nStates
    Grad1 <- rep(1, nStates)
#   probInit
    Res[1,1:nStates] <- Grad1
#   transMat
    for (i in 1:nStates)
        Res[i+1, (i*nStates+1):((i+1)*nStates)] <- Grad1
#   Distribution
    Aux <- GradConstraint(object$distribution)
    if (! is.null(Aux))
    {    Res[(nStates+2):nConstr,(nStates*(nStates+1) + 1):LParam$nParam] <- Aux
    }
    return(Res)
}


GetNParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetNParam", object=object))

}

GetNAllParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetNAllParam", object=object))

}

GetVectParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetVectParam", object=object))

}

GetVectAllParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetVectAllParam", object=object))

}

SetVectParam.distributionClass<-function(object, Vect)
{
    return(NextMethod(generic="SetVectParam", object=object, Vect=Vect))

}

SetVectAllParam.distributionClass<-function(object, Vect)
{
    return(NextMethod(generic="SetVectAllParam", object=object, Vect=Vect))

}

GetNConstraint.distributionClass<-function(object)
{
    return(NextMethod(generic="GetNConstraint", object=object))
}

GradConstraint.distributionClass<-function(object)
{
    return(NextMethod(generic="GradConstraint", object=object))
}


GetNParam.univariateNormalClass<-function(object)
{
    return(list(nParam=2, paramConstr=c(0,eVar)))
}

GetNAllParam.univariateNormalClass<-function(object)
{
    return(list(nParam=2, paramConstr=c(0,eVar)))
}

GetVectParam.univariateNormalClass<-function(object)
{   nStates <- object$nStates
    nParam <- nStates * 2
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   Res[k] <- object$mean[i]
        k<-k+1
        Res[k] <- object$var[i]
        k<-k+1
    }
    return(Res)
}

GetVectAllParam.univariateNormalClass<-function(object)
{   return(GetVectParam.univariateNormalClass(object))
}

SetVectParam.univariateNormalClass<-function(object, Vect)
{   nStates <- object$nStates
    mean<-var<-rep(0, nStates)
    k <- 1
    for (i in 1:nStates)
    {   mean[i] <- Vect[k]
        k<-k+1
        var[i] <- Vect[k]
        k<-k+1
    }
    Res<-distributionSet(dis='NORMAL', mean=mean, var=var, verif=FALSE)
    return(Res)
}

SetVectAllParam.univariateNormalClass<-function(object, Vect)
{   return(SetVectParam.univariateNormalClass(object, Vect))
}

GetNConstraint.univariateNormalClass<-function(object)
{
    return(0)
}

GradConstraint.univariateNormalClass<-function(object)
{
    return(NULL)
}

GetNParam.multivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nParam <- as.integer(dimObs + dimObs*(dimObs+1)/2)
    paramConstr <- c(rep(eNone,dimObs), rep(eVar, dimObs), rep(eCor, as.integer(dimObs*(dimObs-1)/2)))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetNAllParam.multivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nParam <- as.integer(dimObs + dimObs*(dimObs+1)/2)
    paramConstr <- c(rep(eNone,dimObs), rep(eVar, dimObs), rep(eCor, as.integer(dimObs*(dimObs-1)/2)))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetVectParam.multivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nStates <- object$nStates
    nParam <- as.integer(dimObs + dimObs*(dimObs+1)/2)*nStates
    Res <- rep(0, nParam)

    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+dimObs-1)] <- object$mean[[i]]
        k <- k + dimObs
#        Res[k:(k+dimObs-1)] <- diag(object$cov[[i]])
#        k <- k + dimObs
#        aa<-diag(1/sqrt(diag(object$cov[[i]])))
#        matCor <- aa%*%object$cov[[i]]%*%aa
        for (n in 1:dimObs)
        {   for (m in n:dimObs)
            {   Res[k] <- object$cov[[i]][m,n]
                k <- k + 1
            }
        }
    }
    return(Res)
}

GetVectAllParam.multivariateNormalClass<-function(object)
{
    return(GetVectParam.multivariateNormalClass(object))
}

SetVectParam.multivariateNormalClass<-function(object, Vect)
{
    dimObs <- object$dimObs
    nStates <- object$nStates
    k <- 1
    distr <- object
    for (i in 1:nStates)
    {   distr$mean[[i]] <- Vect[k:(k+dimObs-1)]
        k <- k + dimObs
#        aa <- sqrt(diag(Vect[k:(k+dimObs-1)]))
#        k <- k + dimObs
#        matCor <- diag(rep(1, dimObs))
        for (n in 1:dimObs)
        {   for (m in n:dimObs)
            {  distr$cov[[i]][m,n] <- distr$cov[[i]][n,m] <- Vect[k]
                k <- k + 1
            }
        }
    }
    return(distr)
}

SetVectAllParam.multivariateNormalClass<-function(object, Vect)
{
    return(SetVectParam.multivariateNormalClass(object, Vect))
}

GetNConstraint.multivariateNormalClass<-function(object)
{
    return(0)
}

GradConstraint.multivariateNormalClass<-function(object)
{
    return(NULL)
}

GetNParam.mixtureUnivariateNormalClass<-function(object)
{
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nProba <- nMixt - 1
    nParam <- nMean + nVar + nProba
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eProba, nProba))
    return(list(nParam=nParam, paramConstr=paramConstr))
}


GetNAllParam.mixtureUnivariateNormalClass<-function(object)
{
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nParam <- nMean + nVar + nMixt
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eProba, nMixt))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetVectParam.mixtureUnivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nProba <- nMixt - 1
    nParam <- (nMean + nVar + nProba)*nStates
    Res <- rep(0, nParam)

    k <- 1
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Res[k] <- object$mean[[i]][j]
            k <- k + 1
            Res[k] <- object$var[[i]][j]
            k <- k + 1
            if (j < nMixt)
            {   Res[k] <- object$prop[[i]][j]
                k <- k + 1
            }
        }
    }
    return(Res)
 }

GetVectAllParam.mixtureUnivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nParam <- (nMean + nVar + nMixt)*nStates
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Res[k] <- object$mean[[i]][j]
            Res[k+1] <- object$var[[i]][j]
            Res[k+2] <- object$prop[[i]][j]
            k <- k + 3
        }
    }
    return(Res)
 }

SetVectParam.mixtureUnivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nProba <- nMixt - 1
    distr <- object

    k <- 1
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   distr$mean[[i]][j] <- Vect[k]
            k <- k + 1
            distr$var[[i]][j] <- Vect[k]
            k <- k + 1
            if (j < nMixt)
            {   distr$prop[[i]][j] <- Vect[k]
                k <- k + 1
            }
            else
            {   distr$prop[[i]][nMixt] <- 1.0 -sum(distr$prop[[i]][1:(nMixt-1)])
            }
        }
    }
    return(distr)
 }

SetVectAllParam.mixtureUnivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    distr <- object
    k <- 1
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   distr$mean[[i]][j] <- Vect[k]
            distr$var[[i]][j] <- Vect[k+1]
            distr$prop[[i]][j] <- Vect[k+2]
            k <- k + 3
        }
    }
    return(distr)
 }

GetNConstraint.mixtureUnivariateNormalClass<-function(object)
{
    return(object$nStates)
}

GradConstraint.mixtureUnivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    LParam <- GetNAllParam(object)
    Res <- matrix(0, nrow=nStates, ncol=LParam$nParam*nStates)
    k <- 2*nMixt
    Grad1 <- rep(1, nMixt)
    for ( i in 1:nStates)
    {    Res[i,(k+1):(k+nMixt)] <- Grad1
        k <- k + 3*nMixt
    }
    return(Res)
}


GetNParam.mixtureMultivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nStates <- object$nStates
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    nProba <- (nMixt - 1)
    nParam <- nMean + nVar + nCor + nProba
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eCor, nCor), rep(eProba, nProba))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetNAllParam.mixtureMultivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nStates <- object$nStates
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    nParam <- nMean + nVar + nCor + nMixt
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eCor, nCor), rep(eProba, nMixt))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetVectParam.mixtureMultivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    nProba <- nMixt - 1
    nParam <- (nMean + nVar + nCor + nProba)*nStates
    Res <- rep(0, nParam)
    k <- 1

    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Res[k:(k+dimObs-1)] <- object$mean[[i]][[j]]
            k <- k + dimObs
            aa <- object$cov[[i]][[j]]
            for (n in 1:dimObs)
            {   for (m in n:dimObs)
                {   Res[k] <- aa[m,n]
                    k <- k + 1
                }
            }
            if ( j < nMixt)
            {   Res[k] <- object$proportion[[i]][j]
                k <- k + 1
            }
        }
    }
    return(Res)
}

GetVectAllParam.mixtureMultivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
   nParam <- (nMean + nVar + nCor + nMixt)*nStates
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Res[k:(k+dimObs-1)] <- object$mean[[i]][[j]]
            k <- k + dimObs
            aa <- object$cov[[i]][[j]]
            for (n in 1:dimObs)
            {   for (m in n:dimObs)
                {   Res[k] <- aa[m,n]
                    k <- k + 1
                }
            }
            Res[k] <- object$proportion[[i]][j]
            k <- k + 1
         }
    }
   return(Res)
}

SetVectParam.mixtureMultivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    k <- 1
    distr <- object
    aa <- matrix(0, nrow=dimObs, ncol=dimObs)
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   distr$mean[[i]][[j]] <- Vect[k:(k+dimObs-1)]
            k <- k + dimObs
            for (n in 1:dimObs)
            {   for (m in n:dimObs)
                {   aa[m,n] <- aa[n,m] <- Vect[k]
                    k <- k + 1
                }
            }
            distr$cov[[i]][[j]] <- aa
            if (j < nMixt)
            {   distr$proportion[[i]][j] <- Vect[k]
                k <- k + 1
            }
            else
            {   distr$proportion[[i]][nMixt] <- 1.0 - sum(distr$proportion[[i]][1:(nMixt-1)])
            }
        }
    }
    return(distr)
}

SetVectAllParam.mixtureMultivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    k <- 1
    distr <- object
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   distr$mean[[i]][[j]] <- Vect[k:(k+dimObs-1)]
            k <- k + dimObs
            for (n in 1:dimObs)
            {   for (m in n:dimObs)
                {   distr$cov[[i]][[j]][m,n] <- distr$cov[[i]][[j]][n,m] <- Vect[k]
                    k <- k + 1
                }
            }
           distr$proportion[[i]][j] <- Vect[k]
           k <- k + 1
        }
    }
    return(distr)
}

GetNConstraint.mixtureMultivariateNormalClass<-function(object)
{
    return(object$nStates)
}

GradConstraint.mixtureMultivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    LParam <- GetNAllParam(object)
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    k <- nOtherParam <- nMean + nVar + nCor
    Res <- matrix(0, nrow=nStates, ncol=LParam$nParam*nStates)
    Grad1 <- rep(1, nMixt)
    for ( i in 1:nStates)
    {   Res[i,(k+1):(k+nMixt)] <- Grad1
        k <- k + nMixt + nOtherParam
    }
    return(Res)
}

GetNParam.discreteClass <- function(object)
{
    nParam <- object$nLevels - 1
    return(list(nParam=nParam, paramConstr=rep(eProba, nParam)))
}

GetNAllParam.discreteClass <- function(object)
{
    nParam <- object$nLevels
    return(list(nParam=nParam, paramConstr=rep(eProba, nParam)))
}

GetVectParam.discreteClass <- function(object)
{   nStates <- object$nStates
    nLevels <- object$nLevels
    nParam <- (nLevels - 1)*nStates
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+nLevels-2)] <- as.vector(object$proba[[i]][1:(nLevels-1)])
        k <- k+nLevels-1
    }
    return(Res)
}

GetVectAllParam.discreteClass <- function(object)
{   nStates <- object$nStates
    nLevels <- object$nLevels
    nParam <- nLevels*nStates
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+nLevels-1)] <- as.vector(object$proba[[i]])
        k <- k+nLevels
    }
    return(Res)
}

SetVectParam.discreteClass <- function(object, Vect)
{
    nStates <- object$nStates
    nLevels <- object$nLevels
    distr <- object
    noms <- names(object$proba[[1]])
    k <- 1
    for (i in 1:nStates)
    {   distr$proba[[i]][1:(nLevels-1)] <-  Vect[k:(k+nLevels-2)]
        distr$proba[[i]][nLevels] <- 1.0 -sum(distr$proba[[i]][1:(nLevels-1)])
        names(distr$proba[[i]]) <- noms
        k <- k+nLevels-1
    }
    return(distr)
}

SetVectAllParam.discreteClass <- function(object, Vect)
{
    nStates <- object$nStates
    nLevels <- object$nLevels
    distr <- object
    noms <- names(object$proba[[1]])
    k <- 1
    for (i in 1:nStates)
    {   distr$proba[[i]] <-  Vect[k:(k+nLevels-1)]
        names(distr$proba[[i]]) <- noms
        k <- k+nLevels
    }
    return(distr)
}

GetNConstraint.discreteClass<-function(object)
{
    return(object$nStates)
}

computeScoreAndInformation <- function(HMM, obs)
{
    if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM

    if (is.data.frame(obs))
    {   dimObs <- dim(obs)[2]
        if (dimObs == 1)
            obs <- obs[,1]
        else
            obs <- as.matrix(obs[,1:dimObs])
    }

    maListe <- TransfListe(HMM$distribution, obs)
    HMM <- setStorageMode(HMM)

    Res <- .Call("RScoreAndInformation", HMM, maListe$Zt)
    score <- Res[[1]]
    information <- Res[[2]]
    noms <- NomsIndepParamHMM(HMM)
    names(score) <- noms
    row.names(information) <- noms
    colnames(information) <- noms
    return(list(score=score, information=information))
}

asymptoticCov <- function(HMM, obs)
{
    if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM

    if (is.data.frame(obs))
    {   dimObs <- dim(obs)[2]
        if (dimObs == 1)
            obs <- obs[,1]
        else
            obs <- as.matrix(obs[,1:dimObs])
    }

    maListe <- TransfListe(HMM$distribution, obs)
    HMM <- setStorageMode(HMM)

    Res <- .Call("RComputeCov", HMM, maListe$Zt)
    colnames(Res) <- rownames(Res) <- NomsParamHMM(HMM)
    return(Res)
}


asymptoticIterSimCovMat <- function(HMM, obs, nSimul, verbose=FALSE, oldCovMat=NULL)
{
    if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM
    nParam <- GetNAllParam(HMM)$nParam
    if (is.null(oldCovMat))
    {   esp2Mat <- matCov <- matrix(0, nParam, nParam)
        espVect <- rep(0, nParam)
        oldNSimul <- 0
    }
    else
    {   matCov <- oldCovMat$cov
        esp2Mat <- oldCovMat$esp2Mat
        espVect <- oldCovMat$mean
        oldNSimul <- oldCovMat$nSimul
    }

    if (is.list(obs))
    {   nList <- length(obs)
        nObs <- rep(0, nList)
        for (j in 1:nList)
        {   nObs[j] <- length(obs[[j]])
        }
    }
    else
    {   nList <- 1
        nObs <- c(length(obs))
    }

    i <- oldNSimul
    while (i < oldNSimul+nSimul)
    {   simObs <- NULL
        for (j in 1:nList)
        {   simObs <- c(simObs, HMMSim(nObs[j], HMM=HMM)$obs)
        }
        if (HMM$distribution$dis=="NORMAL")
        {   Res <- HMMFit(simObs, dis="NORMAL", nStates=HMM$distribution$nStates, asymptCov=FALSE, control=list(init="USER", initPoint=HMM))
            Teta <- GetVectAllParam(Res$HMM)
            if (!any(is.na(Teta)) && !any(is.infinite(Teta)))
            {   espVect <- (i*espVect + Teta)/(i+1)
                esp2Mat <- (i*esp2Mat + Teta%*%t(Teta))/(i+1)

                if (verbose)
                {   cat(sprintf("iteration %d\n", i+1))
                }
                i <- i + 1
            }
        }
    }
    matCov <- esp2Mat - espVect%*%t(espVect)
    colnames(matCov) <- rownames(matCov) <- NomsParamHMM(HMM)
    return(list(cov=matCov, esp2Mat=esp2Mat, mean = espVect, nSimul=oldNSimul+nSimul))
}



setAsymptoticCovMat<-function(HMMFit, asymptCovMat)
{
    if ( ( class(HMMFit) != "HMMFitClass" ) )
        stop("class(HMMFit) must be 'HMMFitClass'\n")
    if (! is_numeric_matrix(asymptCovMat) )
        stop("asymptCovMat must be a matrix\n")
 #   if (! is_positive_definite(asymptCovMat) )
 #       stop("asymptCovMat must be a definite positive matrix\n")
    colnames(asymptCovMat) <- rownames(asymptCovMat) <- NomsParamHMM(HMMFit)
    HMMFit$asymptCov <- asymptCovMat
    return(HMMFit)
}


NomsParamHMM <- function(object) UseMethod("NomsParamHMM")
NomsIndepParamHMM <- function(object) UseMethod("NomsIndepParamHMM")

NomsParamHMM.default <- function(object)
{   return(NULL)
}

NomsIndepParamHMM.default <- function(object)
{   return(NULL)
}

NomsParamHMM.distributionClass <- function(object)
{
    return(NextMethod(generic="NomsParamHMM", object=object))
}

NomsIndepParamHMM.distributionClass <- function(object)
{
    return(NextMethod(generic="NomsIndepParamHMM", object=object))
}

NomsParamHMM.discreteClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    nLevels <- object$nLevels
    namesProba <- names(object$proba[[1]])
    if (is.null(namesProba))
    {   for (j in 1:nLevels)
        {   Aux <- sprintf("p[%d]", j)
            namesProba <- c(Res, Aux)
        }
    }
    for (i in 1:nStates)
    {   for (j in 1:nLevels)
        {   Aux <- sprintf("State[%d]-%s", i, namesProba[j])
            Res <- c(Res, Aux)
        }
    }
    return(Res)
}

NomsIndepParamHMM.discreteClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    nLevels <- object$nLevels
    namesProba <- names(object$proba[[1]])
    if (is.null(namesProba))
    {   for (j in 1:(nLevels-1))
        {   Aux <- sprintf("p[%d]", j)
            namesProba <- c(Res, Aux)
        }
    }
    for (i in 1:nStates)
    {   for (j in 1:(nLevels-1))
        {   Aux <- sprintf("State[%d]-%s", i, namesProba[j])
            Res <- c(Res, Aux)
        }
    }
    return(Res)
}

NomsParamHMM.univariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    for (i in 1:nStates)
    {   Aux <- sprintf("mean[%d]", i)
        Res <- c(Res, Aux)
        Aux <- sprintf("var[%d]", i)
        Res <- c(Res, Aux)
    }
    return(Res)

}

NomsIndepParamHMM.univariateNormalClass <- function(object)
{
    return(NomsParamHMM(object))
}

NomsParamHMM.mixtureUnivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    nMixt <- object$nMixt
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Aux <- sprintf("state[%d]-mean[%d]", i, j)
            Res <- c(Res, Aux)
            Aux <- sprintf("state[%d]-var[%d]", i, j)
            Res <- c(Res, Aux)
            Aux <- sprintf("state[%d]-prop[%d]", i, j)
            Res <- c(Res, Aux)
        }
    }
    return(Res)
}

NomsIndepParamHMM.mixtureUnivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    nMixt <- object$nMixt
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Aux <- sprintf("state[%d]-mean[%d]", i, j)
            Res <- c(Res, Aux)
            Aux <- sprintf("state[%d]-var[%d]", i, j)
            Res <- c(Res, Aux)
            if (j < nMixt)
            {   Aux <- sprintf("state[%d]-prop[%d]", i, j)
                Res <- c(Res, Aux)
            }
        }
    }
    return(Res)
}

NomsParamHMM.multivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    dimObs <- object$dimObs
    for (i in 1:nStates)
    {   for (j in 1:dimObs)
        {   Aux <- sprintf("state[%d]-mean[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:dimObs)
        {   for (k in j:dimObs)
            {   if (j == k)
                {    Aux <- sprintf("state[%d]-var[%d,%d]", i, j, j)
                }
                else
                {   Aux <- sprintf("state[%d]-cov[%d,%d]", i, j, k)
                }
                Res <- c(Res, Aux)
            }
        }
    }
    return(Res)
}

NomsIndepParamHMM.multivariateNormalClass <- function(object)
{
    return(NomsParamHMM(object))
}

NomsParamHMM.mixtureMultivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    for (i in 1:nStates)
    {   for (m in 1:nMixt)
        {   for (j in 1:dimObs)
            {   Aux <- sprintf("state[%d]-Mixt[%d]-mean[%d]", i, m, j)
                Res <- c(Res, Aux)
            }
            for (j in 1:dimObs)
            {   for (k in j:dimObs)
                {   if (j == k)
                    {   Aux <- sprintf("state[%d]-Mixt[%d]-var[%d,%d]", i, m, j, j)
                    }
                    else
                    {   Aux <- sprintf("state[%d]-Mixt[%d]-Cov[%d,%d]", i, m, j, k)
                    }
                    Res <- c(Res, Aux)
                }
            }
            Aux <- sprintf("state[%d]-prop[%d]", i, m)
            Res <- c(Res, Aux)
        }
    }
    return(Res)
}

NomsIndepParamHMM.mixtureMultivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    for (i in 1:nStates)
    {   for (m in 1:nMixt)
        {   for (j in 1:dimObs)
            {   Aux <- sprintf("state[%d]-Mixt[%d]-mean[%d]", i, m, j)
                Res <- c(Res, Aux)
            }
            for (j in 1:dimObs)
            {   for (k in j:dimObs)
                {   if (j == k)
                    {   Aux <- sprintf("state[%d]-Mixt[%d]-var[%d,%d]", i, m, j, j)
                    }
                    else
                    {   Aux <- sprintf("state[%d]-Mixt[%d]-Cov[%d,%d]", i, m, j, k)
                    }
                    Res <- c(Res, Aux)
                }
            }
            if (m < nMixt)
            {   Aux <- sprintf("state[%d]-prop[%d]", i, m)
                Res <- c(Res, Aux)
            }
        }
    }
    return(Res)
}

NomsParamHMM.HMMClass <- function(object)
{
    nStates <- object$distribution$nStates
    Noms <- NULL
    for (i in 1:nStates)
    {   Aux <- sprintf("Pi[%d]", i)
        Noms <- c(Noms, Aux)
    }

    for (i in 1:nStates)
        for (j in 1:nStates)
        {   Aux <- sprintf("transMat[%d,%d]", i, j)
            Noms <- c(Noms, Aux)
        }
    Aux <- NomsParamHMM(object$distribution)
    Noms <- c(Noms, Aux)
    return(Noms)
}

NomsIndepParamHMM.HMMClass <- function(object)
{
    nStates <- object$distribution$nStates
    Noms <- NULL
    for (i in 1:(nStates-1))
    {   Aux <- sprintf("Pi[%d]", i)
        Noms <- c(Noms, Aux)
    }

    for (i in 1:nStates)
        for (j in 1:(nStates-1))
        {   Aux <- sprintf("transMat[%d,%d]", i, j)
            Noms <- c(Noms, Aux)
        }
    Aux <- NomsIndepParamHMM(object$distribution)
    Noms <- c(Noms, Aux)
    return(Noms)
}

NomsParamHMM.HMMFitClass <- function(object)
{
    return(NomsParamHMM(object$HMM))
}

NomsIndepParamHMM.HMMFitClass <- function(object)
{
    return(NomsIndepParamHMM(object$HMM))
}

summary.HMMFitClass <-function (object, ...)
{
    ans = NULL
    ans$call = object$call
    ans$nIter <- object$nIter
    ans$relVariation <- object$relVariation
    y <- object$HMM$distribution
    if (y$dis == "NORMAL")
    {    if (y$dimObs==1)
            nomloi <- "univariate gaussian"
        else
            nomloi <- sprintf("%d-d gaussian", y$dimObs)
    }
    if (y$dis == "DISCRETE")
        nomloi <- "discrete"
    if (y$dis == "MIXTURE")
    {   if (y$dimObs == 1)
            nomloi <- sprintf("mixture of %d gaussian", y$nMixt)
        else
            nomloi <- sprintf("mixture of %d %d-d gaussian", y$nMixt, y$dimObs)
    }
    Model <- sprintf("%d states HMM with %s distribution", y$nStates, nomloi)
    ans$model <- Model
    ans$LLH = object$LLH
    ans$BIC <- object$BIC
    ans$AIC <- object$AIC
    if (is.null(object$asymptCov))
    {   cat(sprintf("Computing the asymptotic covariance matrix of estimates\n"))
        object$asymptCov <- asymptoticCov(object, object$obs)
    }

    nAllParam <- GetNAllParam(object)
    Value <- GetVectAllParam(object)
    asymptVar <- diag(object$asymptCov)
    se.coef <- sqrt(asymptVar)
    tval = Value/se.coef
    prob = 2 * (1 - pnorm(abs(tval)))

    Noms <- NomsParamHMM(object$HMM)
    ans$coef = cbind(Value, se.coef, tval, prob)
    dimnames(ans$coef) = list(Noms, c(" Estimate",
        " Std. Error", " t value", "Pr(>|t|)"))
    class(ans) <- 'summary.HMMFitClass'
    return(ans)
}

print.summary.HMMFitClass <- function(x, ...)
{  # Description:
    #   Print summary method for an x of class "HMMFitClass".

    # FUNCTION:

    # Call and Model:
    cat("\nCall:", sep="\n")
    cat("----", sep="\n")
    cat(deparse(x$call), "\n", sep = "")

    cat("\nModel:", sep="\n")
    cat("------", sep="\n")
    cat(x$model, "\n", sep = "")
    cat("\nBaum-Welch algorithm status:", sep="\n")
    cat("----------------------------", sep="\n")
    cat(sprintf("Number of iterations : %d\n", x$nIter), sep="")
    cat(sprintf("Last relative variation of LLH function: %f\n", x$relVariation), sep="")

    # Coef Matrix:
    cat("\nCoefficient(s):\n")
    cat("---------------\n")
    signif.stars = getOption("show.signif.stars")
    digits = max(4, getOption("digits") - 4)
    printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)

     cat("\nLog Likelihood: ",
        format(round(x$LLH, 2)), "\n")

    cat("BIC Criterion: ",
        format(round(x$BIC, 2)), "\n")

    cat("AIC Criterion: ",
        format(round(x$AIC, 2)), "\n")

    # Return Value:
#    cat("\n")
    invisible()
}


HMMRank <- function(object) UseMethod("HMMRank")
HMMRank.HMMFitClass <- function(object)
{
    return(HMMRank(object$HMM))
}

HMMRank.HMMClass <- function(object)
{
    Mat <- cbind(rank(object$initProb), HMMGetRank(object$distribution))
    Res <- matrixRank(Mat)
    return(Res)
}

HMMGetRank <- function(object) UseMethod("HMMGetRank")
HMMGetRank.distributionClass <- function(object)
{
    return(NextMethod(generic="HMMGetRank", object=object))
}

HMMGetRank.univariateNormalClass <- function(object)
{

    return(cbind(rank(object$mean), rank(object$var)))

}

HMMGetRank.multivariateNormalClass <- function(object)
{   mat <- matrix(0, ncol=2*object$dimObs, nrow=object$nStates)
    for ( i in 1:object$nStates)
    {   for (j in 1:object$dimObs)
        {   mat[i,j] <- object$mean[[i]][j]
            mat[i, j+object$dimObs] <- object$cov[[i]][j][j]
        }
    }
    for (j in 1:(2*object$dimObs))
    {   mat[,j] <- rank(mat[,j])
    }
    return(mat)
}

HMMGetRank.mixtureUnivariateNormalClass<-function(object)
{
    mat <- matrix(0, ncol = 3, nrow=object$nStates)
    for (i in 1:object$nStates)
    {   mat[i,1] <- min(object$proportion[[i]])
        mat[i,2] <- as.double(t(object$proportion[[i]])%*%object$mean[[i]])
        mat[i,3] <- as.double(t(object$proportion[[i]])%*%object$var[[i]])
    }
    for(i in 1:3)
        mat[,i] <- rank(mat[,i])
    return(mat)
}

HMMGetRank.mixtureMultivariateNormalClass<-function(object)
{
    mat <- matrix(0, ncol = 1 + 2*object$dimObs, nrow=object$nStates)
    for (i in 1:object$nStates)
    {   mat[i,1] <- min(object$proportion[[i]])
        esp <- NULL
        vari <- NULL
        for (j in 1:object$nMixt)
        {   esp <- rbind(esp, object$mean[[i]][[j]])
            vari <- rbind(vari, diag(object$cov[[i]][[j]]))
        }

        for (j in 1:object$dimObs)
        {   mat[i,1 + j] <- as.double(t(object$proportion[[i]])%*%esp[,j])
            mat[i,1 + 2 * j] <- as.double(t(object$proportion[[i]])%*%vari[,j])
        }
    }

    for (j in 1:(1+2*object$dimObs))
    {   mat[,j] <- rank(mat[,j])
    }

    return(mat)
}

HMMGetRank.discreteClass<-function(object)
{
    Res <- matrix(0, nrow=object$nStates, ncol = object$nLevels)
    for ( i in 1:object$nStates )
    {   for (j in 1:object$nLevels)
        {   Res[i,j] <- object$proba[[i]][j]
        }
    }
    for (j in 1:object$nLevels )
    {   Res[,j] <- rank(Res[,j])
    }
    return(Res)
}

HMMSort <- function(object, r=NULL) UseMethod("HMMSort")
HMMSort.HMMFitClass <- function(object, r=NULL)
{
    Res <- object
    object$HMM <- HMMSort(object$HMM, r)
    return(Res)
}

HMMSort.HMMClass <- function(object, r=NULL)
{
    if (is.null(r))
        r <- HMMRank(object)
    initProb <- object$initProb
    initProb[r] <- object$initProb
    transMat <- object$transMat
    transMat[r,r] <- object$transMat
    distribution <- HMMSort(object$distribution, r=r)
    return(HMMSet(initProb=initProb, transMat = transMat, distribution=distribution))
}

HMMSort.distributionClass <- function(object, r)
{
    return(NextMethod(generic="HMMSort", object=object, r=r))
}

HMMSort.univariateNormalClass <- function(object, r)
{
    moy <- object$mean
    vari <- object$var
    moy[r] <- moy
    vari[r] <- vari
    return(distributionSet("NORMAL",mean=moy, var=vari, verif=FALSE))
}


HMMSort.multivariateNormalClass <- function(object, r)
{
    moy <- object$mean
    cova <- object$cov
    moy[r] <- moy
    cova[r] <- cova
    return(distributionSet("NORMAL",mean=moy, cov=cova, verif=FALSE))
}

HMMSort.mixtureUnivariateNormalClass <- function(object, r)
{
    prop <- object$proportion
    moy <- object$mean
    vari <- object$var

    for (i in 1:object$nStates)
    {
        mat <- cbind(prop[[i]],moy[[i]], vari[[i]])
        for (j in 1:dim(mat)[2])
        {    mat[,j] <- rank(mat[,j])
        }
        r1 <- matrixRank(mat)
        prop[[i]][r1] <- prop[[i]]
        moy[[i]][r1] <- moy[[i]]
        vari[[i]][r1] <- vari[[i]]
    }
    prop[r] <- prop
    moy[r] <- moy
    vari[r] <- vari

    return(distributionSet("MIXTURE", proportion=prop, mean=moy, var=vari, verif=FALSE))
}

HMMSort.mixtureMultivariateNormalClass <- function(object, r)
{
    prop <- object$proportion
    moy <- object$mean
    cova <- object$cov
    for (i in 1:object$nStates)
    {
        mat <- NULL
        for (j in 1:object$nMixt)
        {
            mat <- rbind(mat,c(prop[[i]][j],moy[[i]][[j]], diag(cova[[i]][[j]])))
        }
        for (j in 1:dim(mat)[2])
        {   mat[,j] <- rank(mat[,j])
        }
        r1 <- matrixRank(mat)
        prop[[i]][r1] <- prop[[i]]
        moy[[i]][r1] <- moy[[i]]
        cova[[i]][r1] <- cova[[i]]
    }
    prop[r] <- prop
    moy[r] <- moy
    cova[r] <- cova

    return(distributionSet(dis="MIXTURE", proportion=prop, mean=moy, cov=cova, verif=FALSE))
}


HMMSort.discreteClass <- function(object, r)
{
    proba <- object$proba
    name <- names(proba[[1]])
    proba[r] <- proba
    if ( name[1] == 'p 1')
    {   name = NULL
    }
    return(distributionSet(dis="DISCRETE", proba=proba, labels=name, verif=FALSE))
}


Exaeco <- function(x)
{   z <- sort(x)
    n <- length(z)
    return(any(z[2:n]==z[1:(n-1)]))
}


matrixRank <- function(theMat)
{
    if (is.null(theMat))
        return(NULL)
    dimens <- dim(theMat)
    nStates <- dimens[1]
    nCol <- dimens[2]
    Res <- rank(theMat[1:nStates,1])
    isExaeco <- Exaeco(Res)
    i <- 1
    while (isExaeco)
    {   ind <- c(i)
        for (j in (i+1):nStates)
        {   if (Res[j]==Res[i])
            {   ind <- c(ind, j)
            }
        }
        if (length(ind) > 1)
        {   Mat1 <- as.matrix(theMat[ind,2:nCol])
            Res1 <- matrixRank(Mat1)
            Somme <- sum(Res[ind])
            n <- length(ind)
            k <- as.integer((Somme - n*(n-1)/2)/n)
            Res[ind] <- k-1+Res1
            isExaeco <- Exaeco(Res)
        }
        i <- i + 1
    }
    return(Res)
}
