/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDiscrete.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CDISCRETE_H_
#define _CDISCRETE_H_
#pragma once

#include "cDistribution.h"
class cDiscrete : public cDistribution 
{       private :
                uint            mvNClass        ;
        public :
                cCyclicVector<cDMatrix>        mProbaMatVector;
        public :
                cDiscrete(uint theNClass, uint theNProba) ;
				cDiscrete(cDistribution& theSrc) ;
                virtual ~cDiscrete() ;
                void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)  ;
                void ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess) ;
				void ComputeCov(cDMatrix& theCov) ;
				cDVector GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd) ;
				void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL)  ;
                void InitParameters(cBaumWelchInParam &theInParam) ;
                void Print() ;
                uint GetNProba(void) ;
                void GetParam(uint theDeb, cDVector& theParam) ;
                void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetNParam(void){ return GetNProba() ; } ;
				uint GetNFreeParam(void) { return  GetNProba() - 1 ; } ;
                void CopyDistr(cDistribution *theSrc) ;
} ;

extern cDMatrix AddOneVariable(cDMatrix& theCovCour, cDVector& theVect) ;

#endif //_CDISCRETE_H_
