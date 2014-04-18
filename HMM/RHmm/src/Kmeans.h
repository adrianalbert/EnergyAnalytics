/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: Kmeans.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _KMEANS_H_
#define _KMEANS_H_
#pragma once

#ifndef _RDLL_

#include "OTMathUtil.h"
#include "REquivalents.h"
#ifndef uint
        typedef unsigned int uint ;
#endif //uint


void mkmeans(cDVector& theYt, uint theNClass, int* theSeq) ;
void mkmeans(cDVector& theYt, uint theNClass, uint theDimObs, int* theSeq) ;

#endif //_RDLL_

#endif //_KMEANS_H_
