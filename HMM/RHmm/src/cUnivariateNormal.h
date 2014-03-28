/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cUnivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CUNIVARIATENORMAL_H_
#define _CUNIVARIATENORMAL_H_
#pragma once
#include "cDistribution.h"

class cUnivariateNormal : public cDistribution
{       public :
                cDVector       mMean   ;
                cDVector       mVar    ;
        public :
                cUnivariateNormal(uint theNClass = 0) ;
				cUnivariateNormal(cDistribution& theSrc) ;
                virtual ~cUnivariateNormal() ;
                void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba) ;
                void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL) ;
                void ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess) ;
				void ComputeCov(cDMatrix& theCov) ;
				cDVector GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd) ;
				void InitParameters(cBaumWelchInParam &theInParam) ;
                void Print() ;
                void GetParam(uint theDeb, cDVector& theParam) ;
                void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetNParam(void){ return 2 ; } ;
                uint GetNFreeParam(void){ return 2 ; } ;
                void CopyDistr(cDistribution* theSrc) ;
} ;
#endif //_CUNIVARIATENORMAL_H_
