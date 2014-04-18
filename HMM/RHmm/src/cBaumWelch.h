/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cBaumWelch.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CBAUMWELCH_H_
#define _CBAUMWELCH_H_
#pragma once

#include "cInParam.h"
#include "cHmm.h"

class cBaumWelch
{       protected :
                        uint	mtNSample ;
                        uint*   mtT ;
        public :                                                
                        cDMatrix*   mAlpha ;
                        cDMatrix*   mBeta ;
                        cDVector*	mRho ;
                        cDMatrix*   mGamma ;
                        cDMatrix**  mXsi ;
                        cDMatrix*   mSumXsi ;
						cDMatrix*	mDelta ;
                        cDVector    mLogVrais ;
        public :
                cBaumWelch(uint theNSample, uint* theT, uint theNClass) ;
                cBaumWelch(const cInParam &theInParam) ;
                void ForwardBackward(cDMatrix* theCondProba, cHmm& theHMM) ;
				void OutForwardBackward(cDMatrix* theCondProba, cHmm& theHMM, bool theLogData=true) ;
                uint GetSampleSize(uint theN){ return mtT[theN] ;}
                virtual ~cBaumWelch() ;
} ;


#endif //_CBAUMWELCH_H_
