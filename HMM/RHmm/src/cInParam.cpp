/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cInParam.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

cInParam::cInParam(uint theNSample, uint theDimObs, cDVector* theY, distrDefinitionEnum theDistrType, uint theNClass, uint theNMixt, uint theNProba)
{       MESS_CREAT("cInParam")
        mDistrType = theDistrType ;
        mNClass = theNClass ;
        mNMixt = theNMixt ;
        mNProba = theNProba ;
        mNSample = theNSample ;
        mDimObs = theDimObs ;
        if (mNSample > 0)
        {       mY = new cDVector[mNSample] ;
                for (register uint i = 0 ; i < mNSample ; i++)
                        mY[i] = theY[i] ;
        }
        else
                mY = (cDVector *)NULL ;
}
cInParam::~cInParam()
{       MESS_DESTR("cInParam")
        if (mNSample != 0)
        {       for (register uint i = 0 ; i < mNSample ; i++)
                        mY[i].Delete() ;
                delete [] mY ;
                mNSample = 0 ;
        }
}

cInParam &cInParam::operator =(const cInParam &theSrc)
{       
        mDistrType = theSrc.mDistrType ;                
        mNClass = theSrc.mNClass ;
        if (mNSample > 0)
        {       for (register uint i = 0 ; i < mNSample ; i++)
                        mY[i].Delete() ;
                delete mY ;
        }
        mNSample = theSrc.mNSample ;
        mY = new cDVector[mNSample] ;
        
        mDimObs = theSrc.mDimObs ;
        mNProba = theSrc.mNProba ;
        mNMixt = theSrc.mNMixt ;
        
        for (register uint i = 0 ; i < mNSample ; i++)
                mY[i] = theSrc.mY[i] ;
        
        return(*this) ;
}

void cInParam::Print(void)
{
        Rprintf("NbSample = %d\n", mNSample) ;
        for (register uint n = 0 ; n < mNSample ; n++)
                Rprintf("mT[%d]=%d\n", n, mY[n].mSize) ;
}
