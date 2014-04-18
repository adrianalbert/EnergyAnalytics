/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cLogBaumWelch.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CLOGBAUMWELCH_H_
#define _CLOGBAUMWELCH_H_
#pragma once
#include "cInParam.h"
#include "cHmm.h"
#include "LogProb.h"

class cLogBaumWelch
{       private :
                        uint    mvNSample               ;
                        uint*   mvT                             ;
        public :                                                
                        cDMatrix*      mLogAlpha       ;
                        cDMatrix*      mLogBeta        ;
                        cDVector*      mLogRho         ;
                        cDMatrix*      mLogGamma       ;
                        cDMatrix**     mLogXsi         ;
                        cDMatrix*      mSumLogXsi      ;
                        cDVector       mLogVrais       ;
        public :
                cLogBaumWelch(uint theNSample, uint* theT, uint theNClass) ;
                cLogBaumWelch(const cInParam &theInParam) ;
                void LogForwardBackward(cDMatrix* theCondProba, cHmm& theHMM) ;
                uint GetSampleSize(uint theN){ return mvT[theN] ;}
                virtual ~cLogBaumWelch() ;
} ;
#endif //_CLOGBAUMWELCH_H_
