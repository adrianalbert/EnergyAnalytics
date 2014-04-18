/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: SamplesUtil.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _SAMPLESUTIL_H_
#define _SAMPLESUTIL_H_
#pragma once
#include "OTMathUtil.h"

void flatSamples(cDVector* theInVect, uint theNSample, uint theDimObs, uint theNObsAllSamples, cDVector& theOutVect) ;
void listSamples(cDVector& theInVect, uint theNSample, uint theDimObs, uint* theNObsSample, cDVector* theOutVect) ;

#endif //_SAMPLESUTIL_H_
