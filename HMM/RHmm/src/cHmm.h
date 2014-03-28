/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cHmm.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CHMM_H_
#define _CHMM_H_
#pragma once

#include "Hmm.h"
#include "cInParam.h"
#include "cCyclicVector.h"

class cDistribution ;
class cHmm
{       public :
                distrDefinitionEnum mDistrType  ;
                cDVector mInitProba ;
                cCyclicVector<cDMatrix> mTransMatVector ;
                cDistribution*   mDistrParam ;
        public :
                cHmm(distrDefinitionEnum theDistrType, uint theNClass , uint theDimObs=1 , uint theNMixt=0, uint theNProba=0) ;
                cHmm(const cInParam &theInParam) ;
                virtual ~cHmm() ;
                cHmm & operator =(cHmm& theSrc) ;
                void CopyHmm(cHmm& theSrc) { *this = theSrc ; } ;
                void Print() ;
                uint GetNParam(void) ;
				uint GetNFreeParam(void) ;
                void SetParam(cDVector& theParam) ;
                void GetParam(cDVector& theParam) ;            
} ;
#endif //_CHMM_H_
