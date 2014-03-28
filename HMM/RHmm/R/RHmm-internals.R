 ###############################################################
 #### RHmm package                             
 ####                                                         
 #### File: RHmm-internals.R 
 ####                                                         
 #### Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr>
 #### Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ####                                                         
 ###############################################################

setStorageMode <- function(object) UseMethod("setStorageMode")

setStorageMode.paramHMM <- function(object)
{   x <- object
    storage.mode(x$nStates) <- "integer"
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$nMixt) <- "integer"
    storage.mode(x$nLevels) <- "integer"
    storage.mode(x) <- "list"
    class(x) <- "paramHMM"
    return(x)
}

setStorageMode.paramAlgoBW <- function(object)
{   x <- object
    storage.mode(x$iter) <- "integer"
    storage.mode(x$verbose) <- "integer"
    storage.mode(x$nInit) <- "integer"
    storage.mode(x$nIterInit) <- "integer"
    if (!is.null(x$initPoint))
       x$initPoint <- setStorageMode(x$initPoint)
    storage.mode(x) <- "list"
    class(x) <- "paramAlgoBW"
    return(x)
}

setStorageMode.HMMClass <- function(object)
{   x <- object
    storage.mode(x$initProb) <- "double"

    if (is.list(x$transMat)) lapply(x$transMat,function(x){storage.mode(x) <- "double"})
    else storage.mode(x$transMat) <- "double"
    x$distribution <- setStorageMode(object$distribution)
    class(x$distribution) <- class(object$distribution)
    storage.mode(x) <- "list"
    class(x) <- "HMMClass"
    return(x)
}

setStorageMode.distributionClass <- function(object)
{   x<-NextMethod("setStorageMode", object)
    storage.mode(x$nStates) <- "integer"
    storage.mode(x$dis) <-"character"
    class(x) <- c("distributionClass", class(x))
    return(x)
}

setStorageMode.univariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$mean) <- "double"
    storage.mode(x$var) <- "double"
    storage.mode(x) <- "list"
    class(x) <- "univariateNormalClass"
    return(x)
}

setStorageMode.multivariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$mean) <- "list"
    storage.mode(x$cov) <- "list"
    storage.mode(x) <- "list"
    class(x) <- "multivariateNormalClass"
    return(x)
}

setStorageMode.mixtureUnivariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$nMixt) <- "integer"
    storage.mode(x$mean) <- "list"
    storage.mode(x$var) <- "list"
    storage.mode(x$proportion) <- "list"
    storage.mode(x) <- "list"
    class(x) <- "mixtureUnivariateNormalClass"
    return(x)
}

setStorageMode.mixtureMultivariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$nMixt) <- "integer"
    storage.mode(x$mean) <- "list"
    storage.mode(x$cov) <- "list"
    storage.mode(x$proportion) <- "list"
    storage.mode(x) <- "list"
    class(x) <- "mixtureMultivariateNormalClass"
    return(x)
}

setStorageMode.discreteClass <- function(object)
{  x <- object
   storage.mode(x$nLevels) <- "integer"
   storage.mode(x$proba) <- "list"
   storage.mode(x$dimObs) <- "integer"
   storage.mode(x) <- "list"
   class(x) <- "discreteClass"
   return(x)
}

GetAllLevels <- function(obs)
{
    if (is.list(obs))
    {   if (is.factor(obs[[1]]))
        {   Aux1 <- unlist(lapply(obs, levels))
            Aux2 <- sort(Aux1)
            lAux <- length(Aux2)
            Aux3 <- c(TRUE,(Aux2[2:lAux]!=Aux2[1:(lAux-1)]))
            return(Aux2[Aux3])
        }
        else
        {   Y <- obs
            lY <- length(Y)
            for ( i in 1:lY)
                Y[[i]] <- as.factor(Y[[i]])
            return(GetAllLevels(Y))
        }
    }
    else
    {   if (is.factor(obs))
            return(levels(obs))
        else
            return(levels(as.factor(obs)))
    }
}

TransformeListe <- function(paramHMM, obs)
{
    if (paramHMM$dis=='DISCRETE')
    {   if (is.null(paramHMM$Levels))
            labels <- GetAllLevels(obs)
        else
            labels <- paramHMM$Levels
        nLevels <- length(labels)
        if (is.list(obs))
        {   Z <- obs
            lZ <- length(Z)
            Zt <- rep(list(0), lZ)
            for (i in 1:lZ)
            {   Aux1 <- as.character(Z[[i]])
                Aux2 <- factor(Aux1, levels=labels)
                Zt[[i]] <- as.double(Aux2)-1.0
            }
        }
        else
        {   Z <- as.character(obs)
            Z <- factor(Z, levels=labels)
            Zt <- list(as.double(Z) - 1.0)
        }
        return(list(Zt=Zt, nLevels=nLevels, labels=labels))
    }
    else
    {    if (is.list(obs))
         {   return(list(Zt=obs, nLevels=0, labels=NULL))
         }
         else
         {  return(list(Zt=list(obs), nLevels=0, labels=NULL))
         }
    }
}

TransfListe <- function(distribution, obs)
{
    if (distribution$dis == 'DISCRETE')
    {   Levels <- names(distribution$proba[[1]])
    }
    else
        Levels <- NULL
    paramHMM <- list(dis=distribution$dis, Levels=Levels)
    return(TransformeListe(paramHMM, obs))
}





########## end of RHmm_internals
