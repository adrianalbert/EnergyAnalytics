/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cHmm.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

cHmm::cHmm(distrDefinitionEnum theDistrType, uint theNClass, uint theDimObs, uint theNMixture, uint theNProba)
{       MESS_CREAT("cHmm")
        mDistrType = theDistrType ;
        mInitProba.ReAlloc(theNClass) ;

        cDMatrix *transMat = new cDMatrix(theNClass,theNClass,0.0);
        mTransMatVector.push_back(*transMat);
//      mTransMat.ReAlloc(theNClass, theNClass) ;

        switch(mDistrType)
        {       case eNormalDistr :
                        mDistrParam = new cUnivariateNormal(theNClass) ;
                break ;
                case eMultiNormalDistr :
                        mDistrParam = new cMultivariateNormal(theNClass, theDimObs) ;
                break ;
                case eMixtUniNormalDistr :
                        mDistrParam = new cMixtUnivariateNormal(theNClass, theNMixture) ;
                break ;
                case eMixtMultiNormalDistr :
                        mDistrParam = new cMixtMultivariateNormal(theNClass, theNMixture, theDimObs) ;
                break ;
                case eDiscreteDistr :
                        mDistrParam = new cDiscrete(theNClass, theNProba) ;
                break ;
                case eUnknownDistr :
                        mDistrParam = (cDistribution *)NULL ;
                break ;
        }
}

cHmm::cHmm(const cInParam &theInParam)
{       MESS_CREAT("cHmm")
        mInitProba.ReAlloc(theInParam.mNClass);
//      mTransMat.ReAlloc(theInParam.mNClass, theInParam.mNClass) ;

        cDMatrix *transMat = new cDMatrix(theInParam.mNClass,theInParam.mNClass,0.0);
        mTransMatVector.push_back(*transMat);
        delete transMat;

        mDistrType = theInParam.mDistrType ;
        switch(mDistrType)
        {       case eNormalDistr :
                        mDistrParam = new cUnivariateNormal(theInParam.mNClass) ;
                break ;
                case eMultiNormalDistr :
                        mDistrParam = new cMultivariateNormal(theInParam.mNClass, theInParam.mDimObs) ;
                break ;
                case eMixtUniNormalDistr :
                        mDistrParam = new cMixtUnivariateNormal(theInParam.mNClass, theInParam.mNMixt) ;
                break ;
                case eMixtMultiNormalDistr :
                        mDistrParam = new cMixtMultivariateNormal(theInParam.mNClass, theInParam.mNMixt, theInParam.mDimObs) ;
                break ;
                case eDiscreteDistr :
                        mDistrParam = new cDiscrete(theInParam.mNClass, theInParam.mNProba) ;
                break ;
                case eUnknownDistr :
                        mDistrParam = (cDistribution *)NULL ;
                break ;
        }
}

cHmm::~cHmm()
{       MESS_DESTR("cHmm")

        std::vector<cDMatrix>::iterator it;
        for (it = mTransMatVector.begin(); it < mTransMatVector.end(); it++ )
                it->Delete();

        // Probably more to free here
        mInitProba.Delete() ;

        delete mDistrParam ;
        mDistrParam = (cDistribution *)NULL ;
}

cHmm & cHmm::operator = (cHmm &theSrc)
{       
        mInitProba = theSrc.mInitProba ;
        this->mTransMatVector = theSrc.mTransMatVector ;

                /*      for (register uint i = 0 ; i < mvQ ; i++)
        {       mInitProba[i] = theSrc.mInitProba[i] ;
                for (register uint j = 0 ; j < mvQ ; j++)
                        mTransMat[i][j] = theSrc.mTransMat[i][j] ;
        }
*/
        mDistrParam->CopyDistr(theSrc.mDistrParam) ;
        return(*this) ;
}

void cHmm::Print(void)
{
register uint   i, j, k;

        Rprintf("ProbInit :\n") ;
        for (i = 0 ; i < mInitProba.mSize ; i++)
                Rprintf("\t%f", mInitProba[i]) ;


        for (k = 0 ; k < mTransMatVector.size(); k++)
        {
                Rprintf("\nMatrice de transition %u: \n", k) ;

                for (i = 0 ; i < mInitProba.mSize ; i++)
                {
                        for (j = 0 ; j < mInitProba.mSize ; j++)
                                Rprintf("\t%f", mTransMatVector[k][i][j]) ;
                        Rprintf("\n") ;
                }
        }
        Rprintf("\nParameters:\n") ;
        mDistrParam->Print() ;
}

uint cHmm::GetNParam(void)
{
uint myNClass = mInitProba.mSize ;
        return( myNClass * (1 + myNClass + mDistrParam->GetNParam() ) ) ;
}

uint cHmm::GetNFreeParam(void)
{
uint myNClass = mInitProba.mSize ;
	return((myNClass -1)*(myNClass + 1) + myNClass*mDistrParam->GetNFreeParam()) ;
}

void cHmm::SetParam(cDVector& theParam) 
{
uint myNClass = mInitProba.mSize ;
register uint k = 0 ;
        
        mInitProba[myNClass-1] = 1.0 ;
        for (register uint n = 0 ; n < myNClass - 1 ; n++)
        {       mInitProba[n] = theParam[k++] ;
                mInitProba[myNClass-1] -= mInitProba[n] ;
        }

        for (register uint n = 0 ; n < myNClass ; n++)
        {       mTransMatVector[0][n][myNClass-1] = 1.0 ;
                for (register uint p = 0 ; p < myNClass - 1 ; p++)
                {       mTransMatVector[0][n][p] = theParam[k++] ;
                        mTransMatVector[0][n][myNClass-1] -= mTransMatVector[0][n][p] ; // FIXME
                }
        }
        mDistrParam->SetParam(k, theParam) ;
}

void cHmm::GetParam(cDVector& theParam) 
{
uint myNClass = mInitProba.mSize ; 
register uint k = 0 ;
        for (register uint n = 0 ; n < myNClass - 1 ; n++)
                theParam[k++] = mInitProba[n] ;
        for (register uint n = 0 ; n < myNClass ; n++)
                for (register uint p = 0 ; p < myNClass - 1 ; p++)
                        theParam[k++] = mTransMatVector[0][n][p] ; // FIXME
        mDistrParam->GetParam(k, theParam) ;
}

