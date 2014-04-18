/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cRUtils.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"
#ifdef _RDLL_
/*
 *      Récupérer une seule valeur ï¿½ partir d'une liste SEXP à la place n° theNum
 */
void cRUtil::GetValSexp(SEXP theSEXP, uint theNum, uint &theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        theVal = INTEGER(myAux)[0] ;
}

void cRUtil::GetValSexp(SEXP theSEXP, uint theNum, int &theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        theVal = INTEGER(myAux)[0] ;
}

void cRUtil::GetValSexp(SEXP theSEXP, uint theNum, double &theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        theVal = REAL(myAux)[0] ;
}

void cRUtil::GetValSexp(SEXP theSEXP, uint theNum, char* theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        sprintf(theVal, CHAR(STRING_ELT(myAux, 0))) ;
}
void cRUtil::GetValSexp(SEXP theSEXP, uint theNum, SEXP &theVal)
{
        theVal = VECTOR_ELT(theSEXP, theNum) ;
}

/*
 *      Récupérer une vecteur ï¿½ partir d'une liste SEXP à la place n° theNum
 */
void cRUtil::GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, int* theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theDim ; i++)
                theVal[i] = INTEGER(myAux)[i] ;
}
void cRUtil::GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, uint* theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theDim ; i++)
                theVal[i] = INTEGER(myAux)[i] ;
}

void cRUtil::GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, double* theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theDim ; i++)
                theVal[i] = REAL(myAux)[i] ;
}
void cRUtil::GetVectSexp(SEXP theSEXP, uint theNum, cDVector& theVal)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theVal.mSize ; i++)
                theVal[i] = REAL(myAux)[i] ;
}
/*
 *      Récupérer une matrice á partir d'une liste SEXP à la place n° theNum
 */
void cRUtil::GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, int** theMat)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theLigne ; i++)
                for (register uint j = 0 ; j < theCol ; j++)
                        theMat[i][j] = INTEGER(myAux)[i+j*theLigne] ;
}
void cRUtil::GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, uint** theMat)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theLigne ; i++)
                for (register uint j = 0 ; j < theCol ; j++)
                        theMat[i][j] = INTEGER(myAux)[i+j*theLigne] ;
}
void cRUtil::GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, double** theMat)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theLigne ; i++)
                for (register uint j = 0 ; j < theCol ; j++)
                        theMat[i][j] = REAL(myAux)[i+j*theLigne] ;
}

void cRUtil::GetMatSexp(SEXP theSEXP, uint theNum, cDMatrix& theMat)
{
SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;
        for (register uint i = 0 ; i < theMat.mNRow ; i++)
                for (register uint j = 0 ; j < theMat.mNCol ; j++)
                {       theMat[i][j] = REAL(myAux)[i+j*theMat.mNRow] ;
                }
}

/**
 * Retrieves the list of matrix and stores it in theList.
 * We assume that at least one element has been stored before.
 */
void cRUtil::GetMatListSexp(SEXP theSEXP, uint theNum, std::vector<cDMatrix> &theList)
{
        SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;

        if (isMatrix(myAux))
        {
                /* Fill as first matrix */
                GetMatSexp(theSEXP,theNum,theList[0]);
        } else
        {
                uint i;
                uint nrow = theList.at(0).mNRow;
                uint ncol = theList.at(0).mNCol;

                for (i=0;i<(uint)length(myAux);i++)
                {
                        if (theList.size() <= i)
                        {
                                cDMatrix *mat = new cDMatrix(nrow,ncol,0.0);
                                theList.push_back(*mat);
                        }
                        GetMatSexp(myAux, i, theList.at(i) );
                }
        }
}

/**
 * Retrieves the emission probabilities which either can be stored as list of vectors
 * or as a list of matrices (time-dependent then).
 *
 * We assume that at least one element has been stored before.
 */
void cRUtil::GetEmissionSexp(SEXP theSEXP, uint theNum, std::vector<cDMatrix> &theList)
{
        SEXP myAux = VECTOR_ELT(theSEXP, theNum) ;

        uint nrow = theList.at(0).mNRow;
        uint ncol = theList.at(0).mNCol;
        uint i,j;

        if (!isVector(myAux))
                return;

        if (!isMatrix(VECTOR_ELT(myAux,0)))
        {
                /* Parameter is a list of vectors, as the first elements is no matrix */
                cDVector vec;
                vec.ReAlloc(ncol);

                for (i=0;i<nrow;i++)
                {
                        GetVectSexp(myAux, i, vec);

                        for (j=0;j<ncol;j++)
                                theList.at(0)[i][j] = vec[j];
                }
        } else
        {
                /* Parameter is a list of matrices */

                for (i=0;i<(uint)length(myAux);i++)
                {
                        if (theList.size() <= i)
                        {
                                cDMatrix *mat = new cDMatrix(nrow,ncol,0.0);
                                theList.push_back(*mat);
                        }
                        GetMatSexp(myAux, i, theList.at(i) );
                }
        }
}

/*
 *      Récupérer l'ensemble des nombres dans une liste de nombres
 */
void cRUtil::GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, int* theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetValSexp(myAux, i, theVal[i]) ;
}
void cRUtil::GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, uint* theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetValSexp(myAux, i, theVal[i]) ;
}

void cRUtil::GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, double* theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetValSexp(myAux, i, theVal[i]) ;
}

/*
 * Récupérer l'ensemble des vecteurs dans une liste de vecteur
 */
void cRUtil::GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, int** theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetVectSexp(myAux, i, theDim, theVal[i]) ;
}
void cRUtil::GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, uint** theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetVectSexp(myAux, i, theDim, theVal[i]) ;
}
void cRUtil::GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, double** theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetVectSexp(myAux, i, theDim, theVal[i]) ;
}
void cRUtil::GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, cDVector* theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       GetVectSexp(myAux, i, theVal[i]) ;
        }
}
/*
 *      Récupérer l'ensemble des matrices d'une liste de matrices
 */
void cRUtil::GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, int*** theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetMatSexp(myAux, i, theLigne, theCol, theVal[i]) ;
}

void cRUtil::GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, uint*** theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetMatSexp(myAux, i, theLigne, theCol, theVal[i]) ;
}

void cRUtil::GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, double*** theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetMatSexp(myAux, i,theLigne, theCol, theVal[i]) ;
}

void cRUtil::GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, cDMatrix* theVal)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNElt ; i++)
                GetMatSexp(myAux, i, theVal[i]) ;
}

/*
 * Store the given SEXP vectors in an array of vectors
 */
void cRUtil::GetListListVectSexp(SEXP theSEXP, uint theNum, uint theNList1, uint theNList2, cDVector** theVect)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNList1 ; i++)
        {       GetListVectSexp(myAux, i, theNList2, theVect[i]) ;
        }
}

/*
 * Store the given SEXP matrices in an array of matrices.
 */
void cRUtil::GetListListMatSexp(SEXP theSEXP, uint theNum, uint theNList1, uint theNList2, cDMatrix** theMat)
{
SEXP myAux ;
        GetValSexp(theSEXP, theNum, myAux) ;
        for (register uint i = 0 ; i < theNList1 ; i++)
                GetListMatSexp(myAux, i, theNList2, theMat[i]) ;

}

/*
 *      Remplit une seule valeur dans un SEXP à la place n° theNum
 */
void cRUtil::set_val_sexp(int theVal, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(REALSXP, 1)) ;
        INTEGER(theSEXP)[0] = theVal ;
}
void cRUtil::set_val_sexp(uint theVal, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(REALSXP, 1)) ;
        INTEGER(theSEXP)[0] = theVal ;
}
void cRUtil::set_val_sexp(double theVal, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(REALSXP, 1)) ;
        REAL(theSEXP)[0] = theVal ;
}
/*
 *      Remplit un vecteur de taille theDim dans un SEXP 
 */
void cRUtil::SetVectSexp(int *theVect, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP=allocVector(INTSXP, theDim)) ;
        for (register uint i = 0 ; i < theDim ; i++)
                INTEGER(theSEXP)[i] = theVect[i] ;
}

void cRUtil::SetVectSexp(uint *theVect, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP=allocVector(INTSXP, theDim)) ;
        for (register uint i = 0 ; i < theDim ; i++)
                INTEGER(theSEXP)[i] = theVect[i] ;
}

void cRUtil::SetVectSexp(double *theVect, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP=allocVector(REALSXP, theDim)) ;
        for (register uint i = 0 ; i < theDim ; i++)
                REAL(theSEXP)[i] = theVect[i] ;
}
void cRUtil::SetVectSexp(cDVector& theVect, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP=allocVector(REALSXP, theVect.mSize)) ;
        for (register uint i = 0 ; i < theVect.mSize ; i++)
                REAL(theSEXP)[i] = theVect[i] ;
}


/*
 *      Remplit une matrice de taille theLigne x theCol dans un SEXP
*/
void cRUtil::SetMatSexp(int **theMat, uint theLigne, uint theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocMatrix(INTSXP, theLigne, theCol)) ;
        for (register uint i = 0 ; i < theLigne ; i++)
                for (register uint j = 0 ; j < theCol ; j++)
                        INTEGER(theSEXP)[i+j*theLigne] = theMat[i][j] ;
}

void cRUtil::SetMatSexp(uint **theMat, uint theLigne, uint theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocMatrix(INTSXP, theLigne, theCol)) ;
        for (register uint i = 0 ; i < theLigne ; i++)
                for (register uint j = 0 ; j < theCol ; j++)
                        INTEGER(theSEXP)[i+j*theLigne] = theMat[i][j] ;
}

void cRUtil::SetMatSexp(double** theMat, uint theLigne, uint theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocMatrix(REALSXP, theLigne, theCol)) ;
        for (register uint i = 0 ; i < theLigne ; i++)
                for (register uint j = 0 ; j < theCol ; j++)
                        REAL(theSEXP)[i+j*theLigne] = theMat[i][j] ;
}

void cRUtil::SetMatSexp(cDMatrix& theMat, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocMatrix(REALSXP, theMat.mNRow, theMat.mNCol)) ;
        for (register uint i = 0 ; i < theMat.mNRow ; i++)
                for (register uint j = 0 ; j <  theMat.mNCol ; j++)
                        REAL(theSEXP)[i+j*theMat.mNRow] = theMat[i][j] ;
}

/*
 * Remplit une liste de theDim Nombres dans un SEXP
 */
void cRUtil::SetListValSexp(int* theVal, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theDim)) ;
        for (register uint i = 0 ; i < theDim ; i++)
        {       SEXP myAux ;
                set_val_sexp(theVal[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListValSexp(uint* theVal, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theDim)) ;
        for (register uint i = 0 ; i < theDim ; i++)
        {       SEXP myAux ;
                set_val_sexp(theVal[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListValSexp(double* theVal, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theDim)) ;
        for (register uint i = 0 ; i < theDim ; i++)
        {       SEXP myAux ;
                set_val_sexp(theVal[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}
void cRUtil::SetListValSexp(cDVector& theVal, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theVal.mSize)) ;
        for (register uint i = 0 ; i < theVal.mSize ; i++)
        {       SEXP myAux ;
                set_val_sexp(theVal[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

/*
 * Remplit une liste de theNElt vecteur de taille theDim dans un SEXP
 */
void cRUtil::SetListVectSexp(int** theVal, uint theNElt, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetVectSexp(theVal[i], theDim, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListVectSexp(uint** theVal, uint theNElt, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetVectSexp(theVal[i], theDim, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListVectSexp(double** theVal, uint theNElt, uint theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux  ;
                SetVectSexp(theVal[i], theDim, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}
/*
 *  Remplit une liste de theNElt vecteur de tailles deifférentes theDim[i] dans un SEXP
 */
void cRUtil::SetListVectSexp(int** theVal, uint theNElt, uint *theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetVectSexp(theVal[i], theDim[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListVectSexp(uint** theVal, uint theNElt, uint *theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetVectSexp(theVal[i], theDim[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListVectSexp(double** theVal, uint theNElt, uint *theDim, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux  ;
                SetVectSexp(theVal[i], theDim[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}


void cRUtil::SetListVectSexp(cDVector* theVal, uint theNElt, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux  ;
                SetVectSexp(theVal[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListVectSexp(cDMatrix& theVal, SEXP &theSEXP)
{
        uint theNElt;

        theNElt = theVal.mNRow;

        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        mvNbProtect++ ;
        for (uint i = 0 ; i < theNElt; i++)
        {
               SEXP myAux  ;
               SetVectSexp(theVal[i], theVal.mNCol, myAux) ;
               SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }

}

/*
 * Remplit une liste de theNElt matrice de taille theLigne x theCol dans un SEXP
 */
void cRUtil::SetListMatSexp(int*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetMatSexp(theVal[i], theLigne, theCol, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}
void cRUtil::SetListMatSexp(uint*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetMatSexp(theVal[i], theLigne, theCol, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}
void cRUtil::SetListMatSexp(double*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetMatSexp(theVal[i], theLigne, theCol, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}
/*
 * Remplit une liste de theNElt matrices de tailles différentes theLigne[i] x theCol[i] dans un SEXP
 */
void cRUtil::SetListMatSexp(int*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetMatSexp(theVal[i], theLigne[i], theCol[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}
void cRUtil::SetListMatSexp(uint*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetMatSexp(theVal[i], theLigne[i], theCol[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}
void cRUtil::SetListMatSexp(double*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetMatSexp(theVal[i], theLigne[i], theCol[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}


void cRUtil::SetListMatSexp(cDMatrix* theVal, uint theNElt, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNElt)) ;
        for (register uint i = 0 ; i < theNElt ; i++)
        {       SEXP myAux ;
                SetMatSexp(theVal[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

/*
* Remplit une liste de theNList1 elements de listes de theNList2 elements de vecteurs dans un SEXP
*/
void cRUtil::SetListListVectSexp(cDVector** theVect, uint theNList1, uint theNList2, SEXP &theSEXP)
{
        mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNList1)) ;
        for (register uint i = 0 ; i < theNList1 ; i++)
        {       SEXP myAux ;
                SetListVectSexp(theVect[i], theNList2, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

/*
* Remplit une liste de theNList1 elements de listes de theNList2 elements de matrices dans un SEXP
*/
                
void cRUtil::SetListListMatSexp(cDMatrix** theMat, uint theNList1, uint theNList2, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNList1)) ;
        for (register uint i = 0 ; i < theNList1 ; i++)
        {       SEXP myAux ;
                SetListMatSexp(theMat[i], theNList2, myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

void cRUtil::SetListListMatSexp(cDMatrix** theMat, uint theNList1, uint* theNList2, SEXP &theSEXP)
{       mvNbProtect++ ;
        PROTECT(theSEXP = allocVector(VECSXP, theNList1)) ;
        for (register uint i = 0 ; i < theNList1 ; i++)
        {       SEXP myAux ;
                SetListMatSexp(theMat[i], theNList2[i], myAux) ;
                SET_VECTOR_ELT(theSEXP, i, myAux) ;
        }
}

#endif /* _RDLL_ */
