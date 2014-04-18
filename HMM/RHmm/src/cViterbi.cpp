/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cViterbi.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

cViterbi::cViterbi(cInParam& theInParam)
{       MESS_CREAT("cViterbi")
        if (theInParam.mNSample > 0)
        {
        register uint n ;
                mSeq = new uint*[theInParam.mNSample] ;
                for (n = 0 ; n < theInParam.mNSample ; n++)
                        mSeq[n] = new uint[theInParam.mY[n].mSize] ;
                
                mLogProb.ReAlloc(theInParam.mNSample) ;
                for (n = 0 ; n < theInParam.mNSample ; n++)
                        mLogProb[n] = -1e100 ;
        } 
        else
        {       mSeq = NULL ;
                mLogProb.Delete() ;
        }
}

cViterbi::~cViterbi()
{       MESS_DESTR("cViterbi")
        if (mLogProb.mSize > 0)
        {       for (register uint n = 0 ; n < mLogProb.mSize ; n++)
                        delete [] mSeq[n] ;
                delete [] mSeq ;
                mLogProb.Delete() ;
        }
}

void cViterbi::ViterbiPath(cInParam &theInParam, cHmm &theHMM)
{
int             myIndAux                                                ;
double  myAux,
                myAux1                                                  ;
uint    myNSample = theInParam.mNSample ;

cDMatrix* myProbaCond = new cDMatrix[myNSample] ;
        for (register uint n = 0 ; n < myNSample ; n++)
        {       
        uint mySize = theInParam.mY[n].mSize/theInParam.mDimObs ;
                myProbaCond[n].ReAlloc(theInParam.mNClass, mySize) ;
        }

cDVector* myDelta = new cDVector[theInParam.mNClass] ; 
int** myPsi = new int*[theInParam.mNClass] ;
        theHMM.mDistrParam->ComputeCondProba(theInParam.mY, myNSample, myProbaCond) ;
        for (register uint n = 0 ; n < myNSample ; n++)
        {       
        uint mySize = theInParam.mY[n].mSize/theInParam.mDimObs ;
                for (register uint j = 0 ; j < theInParam.mNClass ; j++)
                {       myPsi[j] = new int[mySize] ;
                        myDelta[j].ReAlloc(mySize) ;
                }
        // Initialization
                for (register uint i = 0 ; i < theInParam.mNClass  ; i++)
                {       myDelta[i][0] = log(theHMM.mInitProba[i]) + log(myProbaCond[n][i][0]) ;
                        myPsi[i][0] = 0 ;
                }

        // Recursion
                for (register int t = 0 ; t < (int)mySize -1 ; t++)
                {
                        for (register uint j = 0 ; j < theInParam.mNClass  ; j++)
                        {       myAux = myDelta[0][t] + log(theHMM.mTransMatVector[t][0][j]) ;
                                myIndAux = 0 ;
                                for (register uint i = 1 ; i < theInParam.mNClass ; i++)
                                {
                                        if ((myAux1 = myDelta[i][t] + log(theHMM.mTransMatVector[t][i][j])) > myAux)
                                        {       myAux = myAux1 ;
                                                myIndAux = i ;
                                        }
                                }
                                myDelta[j][t+1] = myAux + log(myProbaCond[n][j][t+1]) ;
                                myPsi[j][t+1] = myIndAux ;
                        }
                }
        // Terminaison
                mLogProb[n] = myDelta[0][mySize-1] ;
                mSeq[n][mySize-1] = 0 ;
                for (register uint i = 1 ; i < theInParam.mNClass ; i++)
                {       if (myDelta[i][mySize-1] > mLogProb[n])
                        {       mLogProb[n] = myDelta[i][mySize-1] ;
                                mSeq[n][mySize-1] = i ;
                        }
                }

                for (register int t =  (int)(mySize-2) ; t >= 0 ; t--)
                        mSeq[n][t] = myPsi[mSeq[n][t+1]][t+1] ;

/*              for (register uint j = 0 ; j < theInParam.mNClass ; j++)
                {       myPsi[j] = new int[theInParam.mY[n].mSize] ;
                        myDelta[j].ReAlloc(theInParam.mY[n].mSize) ;
                }
*/              
                for (register uint j = 0 ; j < theInParam.mNClass ; j++)
                {
                        delete [] myPsi[j] ;
                        myDelta[j].Delete() ;
                }
        }
        for (register uint n = 0 ; n < myNSample ; n++)
                myProbaCond[n].Delete() ;
//      delete myPsi ;
//      delete myDelta ;
//      delete myProbaCond ;
}

