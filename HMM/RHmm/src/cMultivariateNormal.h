/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cMultivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CMULTIVARIATENORMAL_H_
#define _CMULTIVARIATENORMAL_H_
#pragma once

#include "cDistribution.h"
#include "SamplesUtil.h"

class cMultivariateNormal : public cDistribution
{       private :
                uint    mvNClass        ;
        public :
                cDVector*      mMean   ;
                cDMatrix*      mCov    ;
        public :
                cMultivariateNormal(uint theNClass = 0, uint theDimObs = 1) ;
                cMultivariateNormal(cDistribution& theSrc) ;
				virtual ~cMultivariateNormal() ;
                void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba) ;
                void ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess) ;
				void ComputeCov(cDMatrix& theCov) ;
				cDVector GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd) ;
				void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL) ;
                void InitParameters(cBaumWelchInParam &theInParam) ;
                void Print() ;
                void GetParam(uint theDeb, cDVector& theParam) ;
                void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetDimObs() ;
                void CopyDistr(cDistribution* theSrc) ;
                uint GetNParam(void){ return mMean[0].mSize + mMean[0].mSize * (mMean[0].mSize + 1)/2 ; }
                uint GetNFreeParam(void){ return mMean[0].mSize + mMean[0].mSize * (mMean[0].mSize + 1)/2 ; }

} ;
#endif //_CMULTIVARIATENORMAL_H_
