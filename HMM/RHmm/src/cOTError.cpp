/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cOTError.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "cOTError.h"

cOTError::cOTError(const char *theMess)
{
        if (theMess != (char *)NULL) 
#ifndef _RDLL_
                        std::cout << theMess << std::endl ;
                        exit(0) ;
#else
                        error(theMess) ;
#endif //_RDLL_

}
