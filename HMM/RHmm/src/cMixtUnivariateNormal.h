/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cMixtUnivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CMIXTUNIVARIATENORMAL_H_
#define _CMIXTUNIVARIATENORMAL_H_
#pragma once

#include "cDistribution.h"

class cMixtUnivariateNormal : public cDistribution
{       private :
                uint    mvNClass        ;
                uint    mvNMixt         ;
        public :
                cDVector*      mMean   ;
                cDVector*      mVar    ;
                cDVector*      mp              ;
        public :
                cMixtUnivariateNormal(uint theNClass = 0, uint theNMixt = 1) ;
				cMixtUnivariateNormal(cDistribution& theSrc) ;
                virtual ~cMixtUnivariateNormal() ;
                void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)  ;
                void ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess) ;
				void ComputeCov(cDMatrix& theCov) ;
				cDVector GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd) ;
				void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL)  ;
                void InitParameters(cBaumWelchInParam &theInParam) ;
                void Print() ;
                void CopyDistr(cDistribution *theSrc) ;
                void GetParam(uint theDeb, cDVector& theParam) ;
                void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetNParam(void){ return mvNMixt * 3 ; }
                uint GetNFreeParam(void){ return mvNMixt * 3 - 1 ; }
} ;
#endif //_CMIXTUNIVARIATENORMAL_H_
