/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: Debogage.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _DEBOGAGE_H_
#define _DEBOGAGE_H_
#pragma once
#ifdef _DEBOGAGE_
extern uint ourNProtect ;

#define SETVECTSEXP(theVect, theSEXP) {\
        PROTECT(theSEXP=allocVector(REALSXP, theVect.mSize)) ;\
        ourNProtect++ ; \
        for (register uint i = 0 ; i < theVect.mSize ; i++)\
                REAL(theSEXP)[i] = theVect[i] ; \
}
#define SETMATSEXP(theMat, theSEXP) {\
        PROTECT(theSEXP = allocMatrix(REALSXP, theMat.mNRow, theMat.mNCol)) ; \
        ourNProtect++ ; \
        for (register uint i = 0 ; i < theMat.mNRow ; i++) \
                for (register uint j = 0 ; j <  theMat.mNCol ; j++) \
                        REAL(theSEXP)[i+j*theMat.mNRow] = theMat.mMat[i][j] ; \
}

#define SETLISTVALSEXP(theVal, theSEXP) {\
        PROTECT(theSEXP = allocVector(VECSXP, theVal.mSize)) ; \
        ourNProtect++ ; \
        for (register uint i = 0 ; i < theVal.mSize ; i++) \
        {       SEXP myAux ; \
                SETVALSEXP(theVal[i], myAux) ; \
                SET_VECTOR_ELT(theSEXP, i, myAux) ; \
        } \
}

#define SETLISTVECTSEXP(theVect_, theNElt_, theSEXP_) {\
        SEXP myAux_ ; \
        PROTECT(theSEXP_= allocVector(VECSXP, theNElt_)) ;\
        ourNProtect++ ; \
        for (register uint i_= 0 ; i_< theNElt_; i_++) \
        {       SETVECTSEXP(theVect_[i_], myAux_) ; \
                SET_VECTOR_ELT(theSEXP_, i_, myAux_) ; \
        } \
}
#endif //_DEBOGAGE_

#endif //_DEBOGAGE_H_
