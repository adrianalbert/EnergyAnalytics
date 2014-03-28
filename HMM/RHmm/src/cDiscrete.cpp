/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDiscrete.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

cDiscrete::cDiscrete(uint theNClass, uint theNProba) 
{
        MESS_CREAT("cDiscrete")
        if ( (theNClass > 0) && (theNProba > 0) )
    {
                mvNClass = theNClass ;
                cDMatrix *emissionMat = new cDMatrix(theNClass,theNProba,0.0);
        mProbaMatVector.push_back(*emissionMat);
        delete emissionMat;
    } else
    {
        mvNClass = 0 ;
    }
}

cDiscrete::~cDiscrete()
{
        MESS_DESTR("cDiscrete")
}

uint cDiscrete::GetNProba(void)
{
        if (mvNClass > 0)
                        return mProbaMatVector[0].mNCol ;
        else
                return 0 ;
}

void cDiscrete::Print()
{
	for (uint h = 0; h<mProbaMatVector.size();h++ )
	{	Rprintf("Position %d\n",h);
		for (register uint i = 0 ; i < mvNClass ; i++)
		{	Rprintf("State %d :\t", i) ;
			for (register uint j = 0 ; j < GetNProba() ; j++)
				Rprintf("P[%d]=%lf\t", j, mProbaMatVector[h][i][j]) ;
			Rprintf("\n") ;
		}
	}
}

void cDiscrete::ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)
{
register uint   i,
                                n,
                                t       ;

        for (n = 0 ; n < theNSample ; n++)
        {
                for (i = 0 ; i < mvNClass ; i++)
                {
                        for (t = 0 ; t < theY[n].mSize ; t++)
                                theCondProba[n][i][t] = mProbaMatVector[t][i][(uint)theY[n][t]];
                }
        }
}

void cDiscrete::ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess)
{
uint myT = theY.GetSize() ;

uint myNFreeProba = GetNProba() - 1 ;
	for (register uint t = 0 ; t < myT ; t++)
	{	
	uint myBegIndex = (mvNClass - 1) * (mvNClass + 1) ;
		for (register uint i = 0 ; i < mvNClass ; i++)
		{	theHess[i][t] = 0.0 ;
			theGrad[i][t] = 0.0 ;
			for (register uint j = 0 ; j < myNFreeProba ; j++)
			{	if ( (uint)theY[t] == j)
					theGrad[i][t][myBegIndex + j] = 1.0L ;
				else
				{	if ( (uint)theY[t] == myNFreeProba) 
						theGrad[i][t][myBegIndex + j] = -1.0L ;
				}
			}
			myBegIndex += myNFreeProba ;
		}
	}
}

void cDiscrete::ComputeCov(cDMatrix& theCov)
{
uint myBegIndex = (mvNClass - 1) * (mvNClass + 1) ;
uint myNFreeProba = GetNProba() - 1 ;
uint mySizeCour = theCov.GetNCols() ;
cDVector myU(mySizeCour, 0.0) ;

	for (register uint n = 0 ; n < mvNClass ; n++)
	{
		for (register uint i = myBegIndex ; i < myBegIndex + myNFreeProba ; i++)
			myU[i] = -1.0 ;
		theCov = AddOneVariable(theCov, myU) ;
		mySizeCour++ ;
		myU.ReAlloc(mySizeCour, 0.0) ;
		myBegIndex += myNFreeProba ;
	}
}

cDVector cDiscrete::GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd)
{
uint myNFreeProba = GetNProba() - 1 ;
cDVector myNumDistrParam ;
cDVector myNumProba(myNFreeProba) ;
uint myIndCour = 0 ;
	for (register uint j = 0 ; j < mvNClass ; j++)
	{	GetSubVector(theNumDistrParam, myIndCour, myNFreeProba, myNumProba) ;
		myNumDistrParam = cat(myNumDistrParam, myNumProba) ;
		myNumDistrParam = cat(myNumDistrParam, (double)theNextInd) ;
		theNextInd++ ;
		myIndCour += myNFreeProba ;
	}
	return myNumDistrParam ;
}

/* FIXME: This doesn't work for variable emissions */
void cDiscrete::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba)
{
        uint i;
        uint myNProba = GetNProba() ;

        for (i = 0 ; i < mvNClass ; i++)
        {
                double myDenominateur = 0.0 ;
        uint n,t;

        for (n = 0 ; n < theInParam.mNSample ; n++)
                for (t = 0 ; t < theInParam.mY[n].mSize; t++)
                        myDenominateur += theBaumWelch.mGamma[n][i][t];

        for (uint k = 0 ; k < myNProba ; k++)
        {
                for (t=0;t<mProbaMatVector.size();t++)
                        mProbaMatVector[t][i][k] = 0.0;

                for (n = 0 ; n < theInParam.mNSample ; n++)
                for (t = 0 ; t < theInParam.mY[n].mSize ; t++)
                        mProbaMatVector[t][i][k] += theBaumWelch.mGamma[n][i][t]*(theInParam.mY[n][t]==k);

                /* FIXME: for variable emissions */
                if (myDenominateur > MIN_DBLE)
                mProbaMatVector[0][i][k] /= myDenominateur;
            else
                mProbaMatVector[0][i][k] = 0.0 ;
        }
        }
}

void cDiscrete::InitParameters(cBaumWelchInParam& theInParam)
{
        register uint   i, t ;
        uint myNProba = GetNProba() ;

#ifdef _RDLL_
        GetRNGstate();
#endif //_RDLL_

        for (t = 0 ; t < mProbaMatVector.size();t++)
        {
                for (i = 0 ; i < mvNClass ; i++)
                {
                        register uint j;
                        double mySum = 0.0 ;

                        for(j = 0 ; j < myNProba ; j++)
                        {
                                mProbaMatVector[t][i][j] =  unif_rand() ;
                                mySum += mProbaMatVector[t][i][j];
                        }

                        /* Make it a probability measure */
                        for (j=0;j<myNProba;j++)
                                mProbaMatVector[t][i][j] /= mySum;
                }
        }
#ifdef _RDLL_
        PutRNGstate() ;
#endif //_RDLL_
}

void cDiscrete::CopyDistr(cDistribution* theSrc)
{
cDiscrete* mySrc = dynamic_cast<cDiscrete *>(theSrc) ;
	if (mySrc)
	{	mvNClass = mySrc->mvNClass;
        mProbaMatVector = mySrc->mProbaMatVector;
	}
	else
		cOTError("Wrong distribution parameter") ;
}

cDiscrete::cDiscrete(cDistribution& theSrc) 
{
	CopyDistr(&theSrc) ;
}

void cDiscrete::GetParam(uint theDeb, cDVector& theParam)
{
        uint myNProba = this->GetNProba();
        uint k = theDeb ;

        for (uint t = 0 ; t < mProbaMatVector.size(); t++)
                for (uint n = 0 ; n < mvNClass ; n++)
                        for (uint p = 0 ; p < myNProba - 1 ; p++)
                                theParam[k++] = mProbaMatVector[t][n][p];
}

void cDiscrete::SetParam(uint theDeb, cDVector& theParam)
{
        uint myNProba = GetNProba();
        uint k = theDeb;

        for (uint t = 0 ; t < mProbaMatVector.size(); t++)
        {
                for (uint n = 0 ; n < mvNClass ; n++)
                {
                        mProbaMatVector[t][n][myNProba-1] = 1.0 ;
                        for (register uint p = 0 ; p < myNProba - 1 ; p++)
                        {
                                mProbaMatVector[t][n][p] = theParam[k++] ;
                                mProbaMatVector[t][n][myNProba-1] -= mProbaMatVector[t][n][p] ;
                        }
                }
        }
}
