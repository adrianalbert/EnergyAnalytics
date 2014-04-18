/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cBaumWelch.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"


cBaumWelch::cBaumWelch(uint theNSample, uint* theT, uint theNClass)
{	MESS_CREAT("cBaumWelch") 
	mtNSample = theNSample ;
	if (mtNSample == 0)
	{       mtT = NULL ;
		mLogVrais.Delete() ;
		mAlpha = NULL ;
		mBeta = NULL ;
		mGamma = NULL ;
		mXsi = NULL ;
		mSumXsi = NULL ;
		mRho = NULL ;
		return ;
	}
	mtT = new uint[mtNSample] ;
	mLogVrais.ReAlloc(mtNSample) ;
	
	mAlpha = new cDMatrix[mtNSample] ;
	mBeta = new cDMatrix[mtNSample] ;
	mDelta = new cDMatrix[mtNSample] ;
	mGamma = new cDMatrix[mtNSample] ;
	mXsi = new cDMatrix*[mtNSample] ;
	mSumXsi = new cDMatrix[mtNSample] ;
	mRho = new cDVector[mtNSample] ;
	for (register uint n = 0 ; n < mtNSample ; n++)
	{	mtT[n] = theT[n] ;
		mAlpha[n].ReAlloc(theNClass, mtT[n]) ;
		mDelta[n].ReAlloc(theNClass, mtT[n]) ;
		mBeta[n].ReAlloc(theNClass, mtT[n]) ;
		mGamma[n].ReAlloc(theNClass, mtT[n]) ;
		mXsi[n] = new cDMatrix[mtT[n]] ;
		for (register uint t = 0 ; t < mtT[n] ; t++)
			mXsi[n][t].ReAlloc(theNClass, theNClass) ;
		mSumXsi[n].ReAlloc(theNClass, theNClass) ;
		mRho[n].ReAlloc(mtT[n]) ;
	}       
}

cBaumWelch::cBaumWelch(const cInParam &theInParam)
{	MESS_CREAT("cBaumWelch") 
	mtNSample = theInParam.mNSample ;
	if (mtNSample == 0)
	{	mtT = NULL ;
		mLogVrais.Delete() ;
		mAlpha = NULL ;
		mDelta = NULL ;
		mBeta = NULL ;
		mGamma = NULL ;
		mXsi = NULL ;
		mRho = NULL ;
		return ;
	}       
	mtT = new uint[mtNSample] ;
	mLogVrais.ReAlloc(mtNSample) ;
	
	mAlpha = new cDMatrix[mtNSample] ;
	mBeta = new cDMatrix[mtNSample] ;
	mDelta = new cDMatrix[mtNSample] ;
	mGamma = new cDMatrix[mtNSample] ;
	mXsi = new cDMatrix*[mtNSample] ;
	mSumXsi = new cDMatrix[mtNSample] ;
	mRho = new cDVector[mtNSample] ;
	for (register uint n = 0 ; n < mtNSample ; n++)
	{       mtT[n] = (theInParam.mY[n].mSize)/theInParam.mDimObs ;
		mAlpha[n].ReAlloc(theInParam.mNClass, mtT[n]) ;
		mDelta[n].ReAlloc(theInParam.mNClass, mtT[n]) ;
		mBeta[n].ReAlloc(theInParam.mNClass, mtT[n]) ;
		mGamma[n].ReAlloc(theInParam.mNClass, mtT[n]) ;
		mXsi[n] = new cDMatrix[mtT[n]] ;
		for (register uint t=0 ; t < mtT[n] ; t++)
			mXsi[n][t].ReAlloc(theInParam.mNClass, theInParam.mNClass) ;
		mSumXsi[n].ReAlloc(theInParam.mNClass, theInParam.mNClass) ;
		mRho[n].ReAlloc(mtT[n]) ;
	}       
}

cBaumWelch::~cBaumWelch()
{	MESS_DESTR("cBaumWelch") 
	if (mtNSample > 0)
	{	for (register uint n = 0 ; n < mtNSample ; n++)
		{	mAlpha[n].Delete() ;
			mDelta[n].Delete() ;
			mBeta[n].Delete() ;
			mGamma[n].Delete() ;
			for (register uint t = 0 ; t < mtT[n] ; t++)
				mXsi[n][t].Delete() ;
			delete [] mXsi[n] ;
			mSumXsi[n].Delete() ;
			mRho[n].Delete() ;
		}
		delete [] mtT ;
		delete [] mRho ;
		delete [] mXsi ;
		delete [] mSumXsi ;
		delete [] mGamma ;
		delete [] mBeta ;
		delete [] mDelta ;
		delete [] mAlpha ;
	}
}

void cBaumWelch::ForwardBackward(cDMatrix* theCondProba, cHmm& theHMM)
{
register uint i,j ;
register int t ;
double myAux, mySum   ;
uint myNClass = theHMM.mInitProba.mSize ;
double myLLH  ;	
	for (register uint n = 0 ; n < mtNSample ; n++)
	{
	int myT = (int)mtT[n] ;
		mRho[n][0] = 0.0 ;

		for (i = 0 ; i < myNClass ; i++)
		{	mAlpha[n][i][0] = theHMM.mInitProba[i] * theCondProba[n][i][0] ;
			mRho[n][0] += mAlpha[n][i][0] ; 
		}
		for ( i = 0 ; i < myNClass ; i++)
		{	mAlpha[n][i][0] /= mRho[n][0] ; // Normalisation
			mDelta[n][i][0] = mAlpha[n][i][0] ;
		}

		myLLH = log(mRho[n][0]) ;

	//forward
		for (t = 0 ; t < myT-1 ; t++)
		{	mRho[n][t+1] = 0.0 ;
			for (j = 0 ; j < myNClass ; j++)
			{	myAux = 0.0 ;
			for (i = 0 ; i < myNClass ; i++)
				myAux += mAlpha[n][i][t] * theHMM.mTransMatVector[t][i][j] ; /* FIXME: Is t correct or should we shift? */
			mAlpha[n][j][t+1] = myAux * theCondProba[n][j][t+1] ;
			mRho[n][t+1] += mAlpha[n][j][t+1] ;
			}
			for (j = 0 ; j < myNClass ; j++)
			{	mAlpha[n][j][t+1] /= mRho[n][t+1] ;
				mDelta[n][j][t+1] = mAlpha[n][j][t+1] ;
			}

			myLLH += log(mRho[n][t+1]) ;
		}

	// backward
		for (i = 0 ; i < myNClass ; i++)
			mBeta[n][i][myT-1] = 1.0 ;
		for (t = myT-2 ; t >= 0 ; t--)
		{	for (i = 0 ; i < myNClass ; i++)
			{   myAux = 0.0 ;
				for (j = 0 ; j < myNClass ; j++)
					myAux +=  theHMM.mTransMatVector[t][i][j] * theCondProba[n][j][t+1] * mBeta[n][j][t+1] ; /* FIXME: Is t correct or should we shift? */
				mBeta[n][i][t] = myAux/mRho[n][t] ;
			}
		}
	
		// Calcul des Gamma et des Xsi et LogVrais
		mLogVrais[n] = myLLH ;
		for (t = 0 ; t < myT ; t++)
		{	mySum = 0.0 ;
			for (i = 0 ; i < myNClass ; i++)
			{	mGamma[n][i][t] = mAlpha[n][i][t] * mBeta[n][i][t] ;
			mySum += mGamma[n][i][t] ;
			}
			for (i = 0 ; i < myNClass ; i++)
			mGamma[n][i][t] /= mySum ;
		}

		// Calcul des Xsi
		for (i = 0 ; i < myNClass ; i++)
			for (j = 0 ; j < myNClass ; j++)
			{	mSumXsi[n][i][j] = 0.0 ;
				for (t = 0 ; t < myT - 1 ; t++)
				{	mXsi[n][t][i][j] = mGamma[n][i][t] * theHMM.mTransMatVector[t][i][j] * theCondProba[n][j][t+1] * mBeta[n][j][t+1]  / (mRho[n][t] * mBeta[n][i][t]) ; 
					mSumXsi[n][i][j] += mXsi[n][t][i][j] ;
				}
			}
	}
}

void cBaumWelch::OutForwardBackward(cDMatrix* theCondProba, cHmm& theHMM, bool theLogData)
{
	ForwardBackward(theCondProba, theHMM) ;
uint myNClass = theHMM.mInitProba.mSize ;
	
// Compute "true" log(Rho), log(Alpha), log(Beta) or "true" Rho, Alpha, Beta
	for (register uint n = 0 ; n < mtNSample ; n++)
	{
	int myT = (int)mtT[n] ;
		
		if (theLogData)
		{	
		double myLogRho = 0.0 ;
			for (register uint i = 0 ; i < myNClass ; i++)
				mBeta[n][i][myT-1] = 0.0 ;
		
			for (register int t = myT - 2 ; t >= 0 ; t--)
			{	myLogRho += log(mRho[n][t]) ;
				for (register uint i = 0 ; i < myNClass ; i++)
					mBeta[n][i][t] = log(mBeta[n][i][t]) + myLogRho ;
			}
			
			myLogRho = 0.0 ;
			for (register int t = 0 ; t < myT ; t++)
			{	myLogRho += log(mRho[n][t]) ;
				mRho[n][t] = myLogRho ;
				for (register uint i = 0 ; i < myNClass ; i++)
					mAlpha[n][i][t] = log(mAlpha[n][i][t]) + myLogRho ;
			}
		}
		else
		{	
		double myNewRho = 1.0 ;
			
			for (register int t = myT - 2 ; t >= 0 ; t--)
			{	myNewRho *= mRho[n][t] ;
				for (register uint i = 0 ; i < myNClass ; i++)
					mBeta[n][i][t] *= myNewRho ;
			}

			myNewRho = 1.0 ;
			for (register int t = 0 ; t < myT ; t++)
			{	myNewRho *= mRho[n][t] ;
				mRho[n][t] = myNewRho ;
				for (register uint i = 0 ; i < myNClass ; i++)
					mAlpha[n][i][t] *= myNewRho ;
			}
		}
	}
}

