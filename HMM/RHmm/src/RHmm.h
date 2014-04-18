/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: RHmm.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _RHMM_H_
#define _RHMM_H_
#pragma once

#ifdef _RDLL_

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


typedef enum ParamHMMEnum
{       eNClasses=0,
        eObsDim,
        eNMixt,
        eNProba,
        eDistrType
}ParamHMMEnum ;

typedef enum ParamAlgoBWEnum
{       eInitType=0,
        eNMaxIter,
        eTol,
        eVerbose,
        eNInitIter,
        eNMaxIterinit,
        eInitPoint
}ParamAlgoBWEnum ;

typedef enum HMMEnum
{       fInitProba=0,
        fTransMat,
        fDistr
}HMMEnum ;

typedef enum DistEnum
{       gType=0, 
        gNClasses
}DistEnum ;

#ifdef WIN32
        #define DECL_DLL_EXPORT __declspec(dllexport) 
#else
        #define DECL_DLL_EXPORT 
#endif // WIN32

#ifndef BEG_EXTERN_C
        #define BEG_EXTERN_C extern "C" {
        #define END_EXTERN_C }
#endif //BEG_EXTERN_C

#endif // _RDLL_
#endif //_RHMM_H_
