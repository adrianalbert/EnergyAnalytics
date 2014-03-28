/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDerivative.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CDERIVATIVE_H_
#define _CDERIVATIVE_H_
#pragma once

#include "cDistribution.h"
#include "cHmm.h"


class cDerivative /*: public cBaumWelch*/
{
	private :
		uint	mvNSample ;
		uint*	mvT ;
		uint	mvNFreeParam ;
		uint	mvNClass ;
	public :    
		cDVector***	mPsi ;
		cDMatrix***	mOmega ;
		cDVector* mScore ;
		cDMatrix* mInformation ;

	public :
		cDerivative(uint theNSample, uint* theT, uint theNClass, uint theNFreeParam)  ;
		cDerivative(const cInParam &theInParam, uint theNFreeParam) ;
		virtual ~cDerivative() ;
		void ComputeDerivative(cHmm& theHmm, cInParam& theParam) ;
		void ComputeScoreAndInformation(cDVector& theScore, cDMatrix& theInformation) ;
		void ComputeCov(cHmm& theHmm, cDMatrix& theCov) ;
} ;

#endif //_CDERIVATIVE_H_
