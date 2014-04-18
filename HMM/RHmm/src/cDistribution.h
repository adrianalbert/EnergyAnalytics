/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDistribution.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CDISTRIBUTION_H_
#define _CDISTRIBUTION_H_
#pragma once

#include "cBaumWelchInParam.h"
#include "cBaumWelch.h"

class cDistribution  
{       public :
                virtual void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)=0 ;
                virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL)=0 ;
                virtual void InitParameters(cBaumWelchInParam& theInParam)=0 ;
                virtual void Print() = 0 ;
                virtual void CopyDistr(cDistribution* theSrc)=0 ;
                virtual uint GetNParam(void)=0;
                virtual uint GetNFreeParam(void)=0 ;
                virtual void GetParam(uint theDeb, cDVector& theParam)=0 ;
				virtual void SetParam(uint theDeb, cDVector& theParam)=0 ;
				virtual ~cDistribution(){};
				virtual void ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess) = 0 ;
				virtual void ComputeCov(cDMatrix& theCov) = 0 ;
				virtual cDVector GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd)=0;
#ifndef _RDLL_
                void KMeans(cDVector& theYt, uint theNClass, int* theSeq) {
                                mkmeans(theYt, theNClass, theSeq) ;
                } ;
                void KMeans(cDVector& theYt, uint theNClass, uint theDimObs, int* theSeq) {
                                mkmeans(theYt, theNClass, theDimObs, theSeq) ;
                } ;
#endif // _RDLL_

} ;


#endif //_CDISTRIBUTION_H_
