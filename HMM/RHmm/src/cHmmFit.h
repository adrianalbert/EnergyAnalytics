/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cHmmFit.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CHMMFIT_H_
#define _CHMMFIT_H_
#pragma once

#include "cBaumWelchInParam.h"
#include "cBaumWelch.h"
#include "cDistribution.h"
#include "cHmm.h"

class cHmmFit : public cBaumWelch, public cHmm
{       public :        
                double  mBic    ;
				double	mAic	;
                uint    mNIter  ;
                double  mTol    ;
                double  mLLH    ;
        public :
                cHmmFit(distrDefinitionEnum theDistrType, uint theNClass, uint theDimObs=1, uint theNMixt=0, uint theNProba=0, uint theNSample=1, uint* myT=NULL) ;
                cHmmFit(cInParam& theInParam) ;
                virtual ~cHmmFit() ;
                void BaumWelchAlgo(cBaumWelchInParam& theInParam) ;
                void BaumWelchAlgoInit(cBaumWelchInParam& theInParam) ;
                cHmmFit & operator = (cHmmFit& theSrc) ;
                double ComputeLLH(cBaumWelchInParam& theInParam, cDMatrix* theProbaCond) ;
                void ComputeFunction(cBaumWelchInParam& theInParam, cDVector& theValFunct, cDVector& theh, cDMatrix* theProbaCond, double theDelta=1e-3) ;
                void ComputeFunction(cBaumWelchInParam& theInParam, cDMatrix& theValFunct, cDVector& theh, cDMatrix* theProbaCond, double theDelta=1e-3) ;
                void ComputeGradient(cBaumWelchInParam& theInParam, cDVector& theGrad, double theDelta=1e-3) ;
                void ComputeHessian(cBaumWelchInParam& theInPram, cDMatrix& theHess, double theDelta=1e-3) ;
} ;
#endif //_CHMMFIT_H_
