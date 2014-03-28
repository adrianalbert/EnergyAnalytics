/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cOTError.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _COTERROR_H_
#define _COTERROR_H_
#pragma once

#include <cstdlib>
#include <iostream>
#ifdef _RDLL_
        #include "R_ext/Error.h"
#endif // _RDLL_

#ifndef NULL
        #define NULL 0
#endif // NULL

#ifndef uint
        typedef unsigned int uint ;
#endif // uint

class cOTError
{
        public :
                cOTError(const char *theMess) ;
} ;

#endif // _COTERROR_H
