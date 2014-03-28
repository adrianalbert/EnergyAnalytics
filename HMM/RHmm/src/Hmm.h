/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: Hmm.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _HMM_H_
#define _HMM_H_
#pragma once

#ifdef _DEBOGAGE_
#include "Debogage.h"
        #define MESS_CREAT(p){Rprintf("Creation de %s\n", (p)) ;}
        #define MESS_DESTR(p){Rprintf("Destruction de %s\n", (p)) ;}
        #define MESS(p) {Rprintf("%s\n", (p)) ;}
#else
        #define MESS_CREAT(p)
        #define MESS_DESTR(p)
        #define MESS(p)
#endif //_DEBOGAGE

#ifdef _RDLL_
        #include <R.h>
        #include <Rinternals.h>
        #include <Rmath.h>
#else
        #include "Kmeans.h"
        #include "REquivalents.h"
        #include "cRandomGenerator.h"
#endif //_RDLL_


#ifndef SQRT_TWO_PI
        #define SQRT_TWO_PI     2.5066282746310002
#endif // SQRT_TWO_PI

#ifndef MAX
        #define MAX(p,q) ((p) > (q) ? (p) : (q))
#endif //MAX

typedef unsigned int uint ;
typedef enum distrDefinitionEnum
{       eUnknownDistr = -1,
        eNormalDistr = 0,
        eMultiNormalDistr = 1,
        eMixtUniNormalDistr = 2,
        eMixtMultiNormalDistr = 3,
        eDiscreteDistr=4
}distrDefinitionEnum;

typedef enum initEnum
{       eRandom = 0,
        eKMeans = 1,
        eUser = 2
}initEnum;
#endif //_HMM_H_
