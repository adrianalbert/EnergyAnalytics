/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: MultivariateNormalUtil.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

void SymetricInverseAndDet(cDMatrix& theMat, double& theDet, cDMatrix& theInvMat)
{
        LapackInvAndDet(theMat, theInvMat, theDet) ;
}

void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, double*  theDens)
{
register uint i, j, t ;
double  myAux, myRapport ;

uint myDimObs = theMu.mSize ;   
        myRapport = pow(SQRT_TWO_PI, (int)myDimObs)*sqrt(theDet) ;

uint myT = thex.mSize / myDimObs ;

        for ( t = 0 ; t < myT ; t++)
        {       myAux = 0.0 ;
                for (i = 0 ; i < myDimObs ; i++)
                        for (j = 0 ; j < myDimObs ; j++)
                                myAux += (thex[t+i*myT]-theMu[i]) * theInvCov[i][j] * (thex[t+j*myT]-theMu[j]) ;
                theDens[t] = exp(-0.5*myAux)/myRapport ;
        }
}

void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, cDVector&  theDens)
{
register uint i, j, t ;
double  myAux, myRapport ;

uint myDimObs = theMu.mSize ;   
        myRapport = pow(SQRT_TWO_PI, (int)myDimObs)*sqrt(theDet) ;

uint myT = thex.mSize / myDimObs ;

        for ( t = 0 ; t < myT ; t++)
        {       myAux = 0.0 ;
                for (i = 0 ; i < myDimObs ; i++)
                        for (j = 0 ; j < myDimObs ; j++)
                                myAux += (thex[t+i*myT]-theMu[i]) * theInvCov[i][j] * (thex[t+j*myT]-theMu[j]) ;
                theDens[t] = exp(-0.5*myAux)/myRapport ;
        }
}

void InvCovMatDeriv(cDMatrix& theInvCov, cDMatrix* theGrad, cDMatrix** theHess)
{
uint mySize = theInvCov.GetNCols() ;
cDMatrix myCovPrime = Zeros(mySize, mySize) ;
cDMatrix myCovSeconde = Zeros(mySize, mySize) ;
uint myNParam = mySize*(mySize+1)/2 ;
uint k = 0 ;
	for (register uint i = 0 ; i < mySize ; i++)
		for (register uint j = i ; j < mySize ; j++)
		{	myCovPrime[i][j] = myCovPrime[j][i] = 1.0 ;
		cDMatrix myAuxMat = myCovPrime * theInvCov ;
			theGrad[k] = -1.0 * theInvCov * myAuxMat ;
		uint l = 0 ;
			for (register uint p = 0 ; p < mySize ; p++)
				for (register uint q = p ; q < mySize ; q++)
				{	myCovSeconde[p][q] = myCovSeconde[q][p] = 1.0 ;
					theHess[k][l] = theHess[l][k] = -1 * theInvCov * myCovSeconde * theGrad[k] - theGrad[k] * myCovSeconde * theInvCov ;
					l++ ;
					myCovSeconde[p][q] = myCovSeconde[q][p] = 0.0 ;
				}
			k++ ;
			myCovPrime[i][j] = myCovPrime[j][i] = 0.0 ;
		}
}

void SymDetDeriv(cDMatrix& theMat, cDVector& theGrad, cDMatrix& theHess)
{
// Gradient
cDMatrix myAuxMat = theMat ;
uint myNCol = theMat.GetNCols() ;
uint k = 0 ;
	for (register uint i = 0 ; i < myNCol ; i++)
	{	for (register uint j = i ; j < myNCol ; j++)
		{	myAuxMat[i][j] = myAuxMat[j][i] = 0 ;
		double myG0 = LapackDet(myAuxMat) ;
			myAuxMat[i][j] = myAuxMat[j][i] = 1.0 ;
		double myG1 = LapackDet(myAuxMat) ;
			myAuxMat[i][j] = myAuxMat[j][i] = -1.0 ;
		double myGm1 = LapackDet(myAuxMat) ;
		double myA = (myG1 + myGm1)/2.0 - myG0 ;
		double myB = (myG1 - myGm1)/2.0 ;
			theGrad[k++] = 2.0*myA*theMat[i][j] + myB ;
			myAuxMat[i][j] = myAuxMat[j][i] = theMat[i][j] ;
		}
	}
// Hessian
cDMatrix myInvZ = Zeros(9,9) ;	
      myInvZ[0][0] = 0.25 ;
      myInvZ[0][1] = -0.5 ;
      myInvZ[0][2] = 0.25 ;
      myInvZ[0][3] = -0.5 ;
      myInvZ[0][4] = 1.0 ;
      myInvZ[0][5] = -0.5 ;
      myInvZ[0][6] = 0.25 ;
      myInvZ[0][7] = -0.5 ;
      myInvZ[0][8] = 0.25 ;
      myInvZ[1][0] = -0.25 ;
      myInvZ[1][2] = 0.25 ;
      myInvZ[1][3] = 0.5 ;
      myInvZ[1][5] = -0.5 ;
      myInvZ[1][6] = -0.25 ;
      myInvZ[1][8] = 0.25 ;
      myInvZ[2][0] = -0.25 ;
      myInvZ[2][1] = 0.5 ;
      myInvZ[2][2] = -0.25 ;
      myInvZ[2][6] = 0.25 ;
      myInvZ[2][7] = -0.5 ;
      myInvZ[2][8] = 0.25 ;
      myInvZ[3][1] = 0.5 ;
      myInvZ[3][4] = -1.0 ;
      myInvZ[3][7] = 0.5 ;
      myInvZ[4][0] = 0.25 ;
      myInvZ[4][2] = -0.25 ;
      myInvZ[4][6] = -0.25 ;
      myInvZ[4][8] = 0.25 ;
      myInvZ[5][3] = 0.5 ;
      myInvZ[5][4] = -1.0 ;
      myInvZ[5][5] = 0.5 ;
      myInvZ[6][1] = -0.5 ;
      myInvZ[6][7] = 0.5 ;
      myInvZ[7][3] = -0.5 ;
      myInvZ[7][5] = 0.5 ;
      myInvZ[8][4] = 1.0 ;

cDMatrix myInd(9,2) ;
	k = 0 ;
		for (register int i =-1 ; i < 2 ; i++)
			for (register int j = -1 ; j < 2 ; j++)
			{	myInd[k][0] = (double)i ;
				myInd[k][1] = (double)j ;
				k++ ;
			}

	k = 0 ;
	myAuxMat = theMat ;
cDVector myG(9) ;
	for (register uint i = 0 ; i < myNCol ; i++)
	{	for (register uint j = i ; j < myNCol ; j++)
		{	
		uint l = 0 ;
		double myx = theMat[i][j] ;
			for (register uint p = 0 ; p < myNCol ; p++)
			{	for (register uint q = p ; q < myNCol ; q++)
				{	
				double myy = theMat[p][q] ;	
					for (register uint r = 0 ; r < 9 ; r++)
					{	myAuxMat[i][j]= myAuxMat[j][i] = myInd[r][0] ;
						myAuxMat[p][q]=myAuxMat[q][p] = myInd[r][1] ;
						myG[r] = LapackDet(myAuxMat) ;
						myAuxMat[i][j] = myAuxMat[j][i] = myx ;
						myAuxMat[p][q] = myAuxMat[q][p] = myy ;
					}
				cDVector myCoef = myInvZ * myG ;
					if ( (i == p) && (j == q) )
						theHess[k][l] = theHess[l][k] = 2.0*(6.0*myCoef[0]*myx*myx+3.0*(myCoef[1]+myCoef[2])*myx+myCoef[3]+myCoef[4]+myCoef[5]) ;
					else
						theHess[k][l] = theHess[l][k] = 4.0*myCoef[0]*myx*myy+2.0*(myCoef[1]*myx+myCoef[2]*myy)+myCoef[4] ;
					l++ ;
				} // for q
			} // for p
			k++ ;
		} // for j
	} // for i
}

void MultivariateNormalDensityDeriv(cDVector& thex, cDVector& theMu, cDMatrix& theCov, cDMatrix& theInvCov, double theDet, cDVector* theGrad, cDMatrix* theHess)
{
uint myDimObs = theMu.mSize ;   
uint myT = thex.mSize / myDimObs ;	

cDVector myDens(myT) ;
	MultivariateNormalDensity(thex, theMu, theInvCov, theDet, myDens) ;
	MultivariateNormalDensityDeriv(thex, theMu, theCov, theInvCov, theDet, myDens, theGrad, theHess) ;
}

void MultivariateNormalDensityDeriv(cDVector& thex, cDVector& theMu, cDMatrix& theCov, cDMatrix& theInvCov, double theDet, cDVector& theDensity, cDVector* theGrad, cDMatrix* theHess)
{
uint myDimObs = theMu.mSize ;   
uint myNCovParam = myDimObs*(myDimObs + 1)/2 ;
uint myT = theDensity.GetSize() ;
	
cDVector myGradDet(myNCovParam) ;
cDMatrix myHessDet(myNCovParam, myNCovParam) ;
	SymDetDeriv(theCov, myGradDet, myHessDet) ;
cDMatrix* myGradInvCov = new cDMatrix[myNCovParam] ;
cDMatrix** myHessInvCov = new cDMatrix*[myNCovParam] ;
	for (register uint k = 0 ; k < myNCovParam ; k++)
	{	myGradInvCov[k].ReAlloc(myNCovParam, myNCovParam) ;
		myHessInvCov[k] = new cDMatrix[myNCovParam] ;
		for (register uint l = 0 ; l < myNCovParam ; l++)
			myHessInvCov[k][l].ReAlloc(myNCovParam, myNCovParam) ;
	}
	InvCovMatDeriv(theInvCov, myGradInvCov, myHessInvCov) ;

//
	for (register uint t = 0 ; t < myT ; t++)
	{
	// Derivée par rapport à Mu
	cDVector myx(myDimObs) ;
		for (register uint i = 0 ; i < myDimObs ; i++)
			myx[i] = (thex[t+i*myT]-theMu[i]) ;
	
	cDVector myGradMu = theDensity[t]*theInvCov*myx ;
		SetSubVector(myGradMu, 0, theGrad[t]) ;
	
	// Dérivée par rapport à Cov
	cDVector myGradLambda(myNCovParam) ;
		for (register uint k = 0 ; k < myNCovParam ; k++)
			myGradLambda[k]  = AsDouble(Transpose(myx) * myGradInvCov[k] * myx) /-2.0 ;
	
		myGradLambda -= 0.5*myGradDet/theDet ;
		myGradLambda *= theDensity[t] ;
		SetSubVector(myGradLambda, myDimObs, theGrad[t]) ;

	// dérivée seconde par rapport à Mu
	cDMatrix myHessMu2 = theInvCov*myx*Transpose(myx)*theInvCov ;
		myHessMu2 -= theInvCov ;
		myHessMu2 *= theDensity[t] ;
		SetSubMatrix(myHessMu2, 0, 0, theHess[t]) ;

   // dérivée seconde par rapport à Mu et à Cov
	cDMatrix myHessMuLambda(myNCovParam, myDimObs) ;
	cDVector myAuxVect1 = theInvCov * myx ;
		for (register uint k = 0 ; k < myNCovParam ; k++)
		{
		cDVector myAuxVect2 = myGradInvCov[k] * myx * theDensity[t] ;
			for (register uint i = 0 ; i < myDimObs ; i++)
				myHessMuLambda[k][i] = myAuxVect2[i] + myAuxVect1[i]*myGradLambda[k] ;
		}
		SetSubMatrix(myHessMuLambda, myDimObs, 0, theHess[t]) ;
	cDMatrix myHessMuLambdaPrime = Transpose(myHessMuLambda) ;
		SetSubMatrix(myHessMuLambdaPrime, 0, myDimObs, theHess[t]) ;

	// dérivée seconde par rapport à Cov
	cDMatrix myHessLambda2 = -0.5*theDensity[t]*myHessDet/theDet ;
		myHessLambda2 += 0.5*theDensity[t]*myGradDet*Transpose(myGradDet)/(theDet*theDet) ;
		if (theDensity[t] != 0.0)
			myHessLambda2 += myGradLambda * Transpose(myGradLambda)/theDensity[t] ;
		for (register uint k = 0 ; k < myNCovParam ; k++)
			for (register uint l = k ; l < myNCovParam ; l++)
			{
			double myDoubleAux =  AsDouble(Transpose(myx)*myHessInvCov[k][l] * myx) * theDensity[t] ;
				myHessLambda2[k][l] -= 0.5*myDoubleAux ;
				if (l != k)
					myHessLambda2[l][k] -= 0.5*myDoubleAux ;
			}
		SetSubMatrix(myHessLambda2, myDimObs, myDimObs, theHess[t]) ;
	}

	for (register uint i = 0 ; i < myNCovParam ; i++)
	{	myGradInvCov[i].Delete() ;
		for (register uint j = 0 ; j < myNCovParam ; j++)
			myHessInvCov[i][j].Delete() ;
		delete [] myHessInvCov[i] ;
	}

	delete [] myGradInvCov ;
	delete [] myHessInvCov ;
}
