/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cViterbi.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CVITERBI_H_
#define _CVITERBI_H_
#pragma once
#include "cInParam.h"
#include "cHmm.h"
#include "cDistribution.h"

class cViterbi
{       public :
                uint            **mSeq          ;
                cDVector       mLogProb        ;
        public :
                cViterbi(cInParam &theInParam) ;
                ~cViterbi() ;
                void ViterbiPath(cInParam& theInParam, cHmm& theHMM) ;
} ;
#endif //_CVITERBI_H_
