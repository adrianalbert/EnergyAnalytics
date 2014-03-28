/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cMixtUnivariateNormal.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

static void MixtUnivariateNormalDensity(cDVector& theY, uint theNMixt, cDVector& theMean, cDVector& theVar, cDVector& thep, double* theDens)
{
cDVector mySigma(theNMixt) ;
        for (register uint i = 0 ; i < theNMixt ; i++)
                mySigma[i] = sqrt(theVar[i]) ;

        for (register uint t = 0 ; t < theY.mSize ; t++)
        {       theDens[t] = 0.0 ;
                for (register uint i = 0 ; i < theNMixt ; i++)
                {
                double myCR = (theY[t] - theMean[i])/mySigma[i] ;
                        theDens[t] += thep[i]/(SQRT_TWO_PI*mySigma[i])*exp(-0.5*myCR*myCR) ;
                }
                if (theDens[t] < 1e-30)
                        theDens[t] = 1e-30 ;
        }
}

cMixtUnivariateNormal::cMixtUnivariateNormal(uint theNClass, uint theNMixt)
{       MESS_CREAT("cMixtUnivariateNormal")
        mvNClass = theNClass ;
        mvNMixt = theNMixt ;
        if ( (theNClass > 0) && (theNMixt > 0) )
        {       mMean = new cDVector[theNClass] ;
                mVar = new cDVector[theNClass] ;
                mp = new cDVector[theNClass] ;
                for (register uint i = 0 ; i < mvNClass ; i++)
                {       mMean[i].ReAlloc(theNMixt) ;
                        mVar[i].ReAlloc(theNMixt) ;
                        mp[i].ReAlloc(theNMixt) ;
                }
        }
        else
        {       mMean = mVar = mp = NULL ;
                mvNClass = mvNMixt = 0 ;
        }
}

cMixtUnivariateNormal::~cMixtUnivariateNormal()
{       MESS_DESTR("cMixtUnivariateNormal")
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       mMean[i].Delete() ;
                mVar[i].Delete() ;
                mp[i].Delete() ;
        }
        mMean = mVar = mp = NULL ;
        mvNClass = mvNMixt = 0 ;
}

void cMixtUnivariateNormal::ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)
{
        for (register uint n = 0 ; n < theNSample ; n++)
                for (register uint i = 0 ; i < mvNClass ; i++)
                        MixtUnivariateNormalDensity(theY[n], mvNMixt, mMean[i], mVar[i], mp[i], theCondProba[n][i]) ;
}

void cMixtUnivariateNormal::ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess)
{
uint myT = theY.GetSize() ;
	for (register uint t = 0 ; t < myT ; t++)
	{
	uint k = (mvNClass - 1)*(mvNClass + 1) ; // premier indice
		for (register uint j = 0 ; j < mvNClass ; j++)
		{	theGrad[j][t] = 0.0 ;
			theHess[j][t] = 0.0 ;
	
		double mySigma = sqrt(mVar[j][mvNMixt-1]) ;
		double myU = (theY[t] - mMean[j][mvNMixt-1])/mySigma ;
		double myFnMixt = exp(-myU*myU/2.0)/(SQRT_TWO_PI*mySigma) ;
			for (register uint n = 0 ; n < mvNMixt ; n++)
			{
			double mySigma = sqrt(mVar[j][n]) ;
			double myAux = (theY[t] - mMean[j][n])/mySigma ;
			double myAux2 = myAux*myAux ;
			double myDensity = exp(-myAux2/2.0)/(SQRT_TWO_PI*mySigma) ;
		
			// d/dm[j][n]
				theGrad[j][t][k] =  mp[j][n]*myAux/mySigma*myDensity ; 
			double myAux3 = (myAux2-1)/mVar[j][n]*myDensity ;
			// d/dVar[j][n]	
				theGrad[j][t][k+1] = mp[j][n]*myAux3/2.0 ;
			// d/dp[j][n]
				if (n < mvNMixt - 1)
					theGrad[j][t][k+2] = myDensity - myFnMixt ;
			// d^2/dm[j][n]^2	
				theHess[j][t][k][k] = mp[j][n]*myAux3 ;
			// d^2/dm[j][n]*dVar[j][n]
				theHess[j][t][k][k+1] = theHess[j][t][k+1][k] = mp[j][n]*(myAux*(myAux2-3.0))/(2.0*mySigma)*myDensity ;
			// d^2/dm[j][n]*dp[j][n]
				if ( n < mvNMixt - 1)
					theHess[j][t][k][k+2] = theHess[j][t][k+2][k] = myAux/mySigma*myDensity ;
			// d^2/dVar[j][n]^2
				theHess[j][t][k+1][k+1] = mp[j][n]*(myAux2*myAux2*-6*myAux2 + 3)/(4.0*mVar[j][n]*mVar[j][n])*myDensity ;
			// d^2/dVar[nj][n]*dp[j][n] 
				if (n < mvNMixt-1)
					theHess[j][t][k+1][k+2] = theHess[j][t][k+2][k+1] = myAux3/2.0 ;
				if (n < mvNMixt - 1)
					k += 3 ;
				else
					k += 2 ;
			}
		}
	}
}

void cMixtUnivariateNormal::ComputeCov(cDMatrix& theCov)
{
uint myBegIndex = (mvNClass - 1) * (mvNClass + 1) ;
uint myNFreeMixt = 3*mvNMixt - 1 ;
uint mySizeCour = theCov.GetNCols() ;
cDVector myU(mySizeCour, 0.0) ;

	for (register uint n = 0 ; n < mvNClass ; n++)
	{
		for (register uint i = myBegIndex + 2 ; i < myBegIndex + myNFreeMixt ; i+=3)
			myU[i] = -1.0 ;
		theCov = AddOneVariable(theCov, myU) ;
		mySizeCour++ ;
		myU.ReAlloc(mySizeCour, 0.0) ;
		myBegIndex += myNFreeMixt ;
	}

}

cDVector cMixtUnivariateNormal::GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd)
{
uint myNFreeParam = mvNMixt * 3 - 1 ;
cDVector myNumDistrParam ;
cDVector myNumMixt(myNFreeParam) ;
uint myIndCour = 0 ;
	for (register uint j = 0 ; j < mvNClass ; j++)
	{	GetSubVector(theNumDistrParam, myIndCour, myNFreeParam, myNumMixt) ;
		myNumDistrParam = cat(myNumDistrParam, myNumMixt) ;
		myNumDistrParam = cat(myNumDistrParam, (double)theNextInd) ;
		theNextInd++ ;
		myIndCour += myNFreeParam ;
	}
	return myNumDistrParam ;

}

void cMixtUnivariateNormal::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba)
{       
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       
        double mySumGammai = 0.0 ;
        register uint   n,
                                                t       ;
                for (n = 0 ; n < theInParam.mNSample ; n++)
                        for (t = 0 ; t < theInParam.mY[n].mSize  ; t++)
                                mySumGammai += theBaumWelch.mGamma[n][i][t] ;
                        
                for (register uint l = 0 ; l < mvNMixt ; l++)
                {                               
                        double myGammail ;
                        double mySumGammail = 0.0 ;
                        double myMoy = 0.0 ;
                        double myVar = 0.0 ;
                        for (n = 0 ; n < theInParam.mNSample ; n++)
                                for (t = 0 ; t < theInParam.mY[n].mSize ; t++)
                                {       
                                double  myStd = sqrt(mVar[i][l]),
                                                myCR = (theInParam.mY[n][t] - mMean[i][l])/myStd ;
                                        myGammail = theBaumWelch.mGamma[n][i][t] * mp[i][l] * 1/(SQRT_TWO_PI*myStd) * exp(-0.5*myCR*myCR) / theCondProba[n][i][t] ;
                                        mySumGammail += myGammail ; 
                                        myMoy += myGammail * theInParam.mY[n][t] ;
                                        myVar += myGammail * theInParam.mY[n][t] * theInParam.mY[n][t] ;
                                }
                        mp[i][l] = mySumGammail / mySumGammai ;
                        mMean[i][l] = myMoy/mySumGammail ;
                        mVar[i][l] = myVar/mySumGammail ;
                        mVar[i][l] -= mMean[i][l] * mMean[i][l] ;
                }
        }
}               

void cMixtUnivariateNormal::InitParameters(cBaumWelchInParam &theInParam)
{
#ifdef _RDLL_
        GetRNGstate();
#endif //_RDLL_
double  myMoy = 0, 
                myVar = 0,
                mystdev         ;
register uint s = 0 ;
        for (register uint n = 0 ; n < theInParam.mNSample ; n++)
                for (register uint t = 0 ; t < theInParam.mY[n].mSize  ; t++)
                        {       myMoy = ((double)s*myMoy + theInParam.mY[n][t])/(double)(s+1) ;
                                myVar = ((double)s*myVar + theInParam.mY[n][t]*theInParam.mY[n][t])/(double)(++s) ;
                        }
        myVar -= myMoy*myMoy ;
        mystdev = sqrt(myVar) ;

        for (register uint i = 0 ; i < mvNClass ; i++)
                {       double mySomme = 0.0 ;
                        register uint l ;
                        for (l = 0 ; l < mvNMixt ; l++)
                        {       mMean[i][l] =  -2*mystdev + myMoy + 2*mystdev * unif_rand() ;
                                mVar[i][l] = 0.5*myVar + 3*myVar * unif_rand() ;         
                                mp[i][l] = unif_rand() ;        
                                mySomme += mp[i][l] ;
                        }
                        for (l = 0 ; l < mvNMixt ; l++)
                                mp[i][l] /= mySomme ;
                }
                        
#ifdef _RDLL_
        PutRNGstate() ;
#endif //_RDLL_
}

void cMixtUnivariateNormal::CopyDistr(cDistribution* theSrc)
{
cMixtUnivariateNormal* mySrc = dynamic_cast<cMixtUnivariateNormal *>(theSrc) ;
	if (mySrc)
	{	mvNClass = mySrc->mvNClass ;
        mvNMixt = mySrc->mvNMixt ;
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       mMean[i] = mySrc->mMean[i] ;
                mVar[i] = mySrc->mVar[i] ;
                mp[i] = mySrc->mp[i] ;
        }
	}
	else
		cOTError("Wrong distribution in cMixtUnivariateNormal") ;
}

cMixtUnivariateNormal::cMixtUnivariateNormal(cDistribution& theSrc)
{
	CopyDistr(&theSrc) ;
}

void cMixtUnivariateNormal::Print()
{
        Rprintf("Parameters\n") ;
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       Rprintf("State %d\n", i) ;
                for (register uint j = 0 ; j < mvNMixt ; j++)
                        Rprintf("\tm[%d]=%lf - s[%d]=%lf - p[%d]=%lf\n", j, mMean[i][j], 
                                j, sqrt(mVar[i][j]), j, mp[i][j]) ;
                Rprintf("\n") ;
        }
}

void cMixtUnivariateNormal::GetParam(uint theDeb, cDVector& theParam)
{
register uint k = theDeb ;
        for (register uint n = 0 ; n < mvNClass ; n++)
                for (register uint p = 0 ; p < mvNMixt ; p++)
                {       theParam[k++] = mMean[n][p] ;
                        theParam[k++] = mVar[n][p] ;
                        if (p < mvNMixt-1)
                                theParam[k++] = mp[n][p] ;
                }       
}

void cMixtUnivariateNormal::SetParam(uint theDeb, cDVector& theParam)
{
register uint k = theDeb ;
        for (register uint n = 0 ; n < mvNClass ; n++)
        {       mp[n][mvNMixt-1] = 1.0 ;
                for (register uint p = 0 ; p < mvNMixt ; p++)
                {       mMean[n][p] = theParam[k++] ;
                        mVar[n][p] = theParam[k++] ;
                        if (p < mvNMixt-1)
                        {       mp[n][p] = theParam[k++] ;
                                mp[n][mvNMixt-1] -= mp[n][p] ;
                        }
                }
        }
}

