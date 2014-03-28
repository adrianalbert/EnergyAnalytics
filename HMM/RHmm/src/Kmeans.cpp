/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: Kmeans.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

#ifndef _RDLL_
void mkmeans(cDVector& theYt, uint theNClass, int* theSeq)
{ 
cDVector myKmeans(theNClass) ;
uint myT = theYt.mSize ;
double myMoy = 0.0,
         myVar = 0.0 ;
register uint t ;
        for (t = 0 ; t < myT ; t++)
        {       myMoy = ((double)t * myMoy + theYt[t])/(double)(t+1) ;
                myVar = ((double)t * myVar + theYt[t] * theYt[t])/((double)(t+1)) ;
        }
        myVar -= myMoy*myMoy ;

        register uint k ;
double mystdev = sqrt(myVar) ;
        for (k = 0 ; k < theNClass ; k++)
                myKmeans[k] =  -2*mystdev + myMoy + 2*mystdev * unif_rand() ;

int myNbChangement  ;
int *myNbObs = new int[theNClass] ;

int myIter = 0 ;
        for (t = 0 ; t < myT ; t++)
                theSeq[t] = 0 ;
        do
        {       myIter++ ;
                myNbChangement = 0 ;
                for (t = 0 ; t < myT ; t++)
                {       double myMin ;
                        int kmin = 0 ;
                        myMin = fabs(theYt[t] - myKmeans[0]) ;
                        for (k = 1 ; k < theNClass ; k++)
                        { double myAux ;
                                myAux = fabs(theYt[t] - myKmeans[k]) ;
                                if ( myAux < myMin  ) 
                                {       myMin = myAux ;
                                        kmin = k ;
                                }
                        }
                        if (theSeq[t] != kmin)
                        {       theSeq[t] = kmin ;
                                myNbChangement++ ;
                        }       
                }
                if (myNbChangement > 0)
                {       register uint i ;
                        for (i = 0 ; i < theNClass ; i++)
                        {       myNbObs[i] = 0 ;
                                myKmeans[i] = 0.0 ;
                        }
                        for (t = 0 ; t < (int)myT ; t++)
                        {       i = theSeq[t] ;
                                myKmeans[i] = ((double)myNbObs[i] * myKmeans[i] + theYt[t]) / (double)(myNbObs[i] + 1) ;
                                myNbObs[i]++ ;
                        }
                }
        }
        while ( (myNbChangement > 0) && (myIter < 1000) ) ;
        myKmeans.Delete() ;
        delete myNbObs ;
}


void mkmeans(cDVector& theYt, uint theNClass, uint theDimObs, int* theSeq)
{ 
uint myT = (theYt.mSize/theDimObs) ;
cDMatrix myKmeans(theNClass, theDimObs) ;

cDVector       myMoy(theDimObs),
                        myVar(theDimObs) ;

register uint   i,
                                t ;
        for (t = 0 ; t < myT ; t++)
                for (i = 0 ; i < theDimObs ; i++)
                {       myMoy[i] = ((double)t * myMoy[i] + theYt[t+myT*i])/(double)(t+1) ;
                        myVar[i] = ((double)t * myVar[i] + theYt[t+myT*i] * theYt[t+myT*i])/((double)(t+1)) ;
                }
        for(register uint i = 0 ; i < theDimObs ; i++)
                myVar[i] -= myMoy[i]*myMoy[i] ;


register uint k ;
cDVector myStDev(theDimObs)  ;
        for (register uint i = 0 ; i < theDimObs ; i++)
                myStDev[i] = sqrt(myVar[i]) ;

        for (k = 0 ; k < theNClass ; k++)
                for(i = 0 ; i < theDimObs ; i++)
                        myKmeans[k][i] =  -2*myStDev[i] + myMoy[i] + 2*myStDev[i] * unif_rand() ;


int myNbChangement  ;
int* myNbObs = new int[theNClass] ;

        for (t = 0 ; t < myT ; t++)
                theSeq[t] = 0 ;
int myIter = 0 ;
        do
        {       myNbChangement = 0 ;
                myIter++ ;
                for (t = 0 ; t < myT ; t++)
                {       double myMin ;
                        int kmin = 0 ;
                        myMin = 0.0 ;
                        for (i = 0 ; i < theDimObs ; i++)
                        {       double myDist = theYt[t+myT*i] - myKmeans[0][i] ;
                                myMin += myDist * myDist  ;
                        }
                        for (k = 1 ; k < theNClass ; k++)
                        {       double myAux = 0.0 ;
                                for (i = 0 ; i < theDimObs ; i++)
                                {       double myDist = theYt[t+myT*i] - myKmeans[k][i] ;
                                        myAux +=  myDist * myDist ;
                                }
                                if ( myAux < myMin  ) 
                                {       myMin = myAux ;
                                        kmin = k ;
                                }
                        }
                        if (theSeq[t] != kmin)
                        {       theSeq[t] = kmin ;
                                myNbChangement++ ;
                        }       
                }
                if (myNbChangement > 0)
                {       for (k = 0 ; k < theNClass ; k++)
                        {       myNbObs[k] = 0 ;
                                for (i = 0 ; i < theDimObs ; i++)
                                        myKmeans[k][i] = 0.0 ;
                        }
                        for (t = 0 ; t < (int)myT ; t++)
                        {       k = theSeq[t] ;
                                for (i = 0  ; i < theDimObs ; i++)
                                        myKmeans[k][i] = ((double)myNbObs[k] * myKmeans[k][i] + theYt[t+myT*i]) / (double)(myNbObs[k] + 1) ;
                                myNbObs[k]++ ;
                        }
                }
        }
        while ( (myNbChangement > 0) && (myIter < 1000) );
        delete myNbObs ;
//      myKmeans.Delete() ;

}
#endif // _RDLL_
