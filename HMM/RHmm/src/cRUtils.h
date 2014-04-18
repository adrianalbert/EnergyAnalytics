/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cRUtils.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CRUTILS_H_
#define _CRUTILS_H_
#pragma once
#ifdef _RDLL_
#include "OTMathUtil.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <vector>

#ifndef uint
        typedef unsigned int uint ;
#endif // uint

class cRUtil
{       private :
                int     mvNbProtect     ;
        public :
                cRUtil(){mvNbProtect = 0 ;} ;
                void EndProtect(void){if (mvNbProtect > 0) {UNPROTECT(mvNbProtect); mvNbProtect = 0 ; }} ;
                 ~cRUtil(){mvNbProtect = 0 ;};
#ifdef _DEBUG
				int GetNProtect(void){ return mvNbProtect ;} ;
#endif // _DEBUG
				/*
                 *      R�cup�rer une seule valeur � partir d'une liste SEXP � la place n� theNum
                 */
                void GetValSexp(SEXP theSEXP, uint theNum, uint &theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, int &theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, double &theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, char* theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, SEXP &theVal) ;
                /*
                *       R�cup�rer une vecteur � partir d'une liste SEXP � la place n� theNum
                */
                void GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, double* theVal) ;
                void GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, int* theVal) ;
                void GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, uint* theVal) ;
                void GetVectSexp(SEXP theSEXP, uint theNum, cDVector& theVal) ;
                /*
                *       R�cup�rer une matrice � partir d'une liste SEXP � la place n� theNum
                */
                void GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, int** theMat) ;
                void GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, uint** theMat) ;
                void GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, double** theMat) ;
                void GetMatSexp(SEXP theSEXP, uint theNum, cDMatrix& theMat) ;
                void GetMatListSexp(SEXP theSEXP, uint theNum, std::vector<cDMatrix> &theList);
                void GetEmissionSexp(SEXP theSEXP, uint theNum, std::vector<cDMatrix> &theList);
                /*
                *       R�cup�rer l'ensemble des nombres dans une liste de nombres
                */
                void GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, int* theVal) ;
                void GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, uint* theVal) ;
                void GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, double* theVal) ;
                /*
                * R�cuperer l'ensemble des vecteurs dans une liste de vecteur
                */
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, int** theVal) ;
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, uint** theVal) ;
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, double** theVal) ;
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, cDVector* theVal) ;
                /*
                *       R�cup�rer l'ensemble des matrices d'une liste de matrices
                */
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, int*** theVal) ;
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, uint*** theVal) ;
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, double*** theVal) ;
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, cDMatrix* theVal) ;
                
                /*
                 * R�cup�rer l'ensemble des vecteurs dans une liste de liste de vecteurs
                 */
                void GetListListVectSexp(SEXP theSEXP, uint theNum, uint theNList1, uint theNList2, cDVector** theVect) ;

                /*
                 * R�cup�rer l'ensemble des vecteurs dans une liste de liste de matrices
                 */
                void GetListListMatSexp(SEXP theSEXP, uint theNum, uint theNList1, uint theNList2, cDMatrix** theVect) ;

                /*
                *       Remplit une seule valeur dans un SEXP � la place n� theNum
                */
                void set_val_sexp(int theVal, SEXP &theSEXP) ;
                void set_val_sexp(uint theVal, SEXP &theSEXP) ;
                void set_val_sexp(double theVal, SEXP &theSEXP) ;
                /*
                *       Remplit un vecteur de taille theDim dans un SEXP 
                */
                void SetVectSexp(int *theVect, uint theDim, SEXP &theSEXP) ;
                void SetVectSexp(uint *theVect, uint theDim, SEXP &theSEXP) ;
                void SetVectSexp(double *theVect, uint theDim, SEXP &theSEXP) ;
                void SetVectSexp(cDVector& theVect, SEXP &theSEXP) ;
          /*
                *       Remplit une matrice de taille theLigne x theCol dans un SEXP
                */
                void SetMatSexp(int **theMat, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetMatSexp(uint **theMat, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetMatSexp(double **theMat, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetMatSexp(cDMatrix& theMat, SEXP &theSEXP) ;
                /*
                * Remplit une liste de theDim Nombres dans un SEXP
                */
                void SetListValSexp(int* theVal, uint theDim, SEXP &theSEXP) ;
                void SetListValSexp(uint* theVal, uint theDim, SEXP &theSEXP) ;
                void SetListValSexp(double* theVal, uint theDim, SEXP &theSEXP) ;
                void SetListValSexp(cDVector& theVal,  SEXP &theSEXP) ;
                /*
                * Remplit une liste de theNElt vecteur de taille theDim dans un SEXP
                */
                void SetListVectSexp(int** theVal, uint theNElt, uint theDim, SEXP &theSEXP) ;
                void SetListVectSexp(uint** theVal, uint theNElt, uint theDim, SEXP &theSEXP) ;
                void SetListVectSexp(double** theVal, uint theNElt, uint theDim, SEXP &theSEXP) ;
                /*
                * Remplit une liste de theNElt vecteurs de tailles diff�rentes theDim[i] dans un SEXP
                */
                void SetListVectSexp(int** theVal, uint theNElt, uint* theDim, SEXP &theSEXP) ;
                void SetListVectSexp(uint** theVal, uint theNElt, uint *theDim, SEXP &theSEXP) ;
                void SetListVectSexp(double** theVal, uint theNElt, uint *theDim, SEXP &theSEXP) ;
                void SetListVectSexp(cDVector* theVal, uint theNElt, SEXP &theSEXP) ;
                void SetListVectSexp(cDMatrix& theVal, SEXP &theSEXP) ;

                /*
                * Remplit une liste de theNElt matrice de taille theLigne x theCol dans un SEXP
                */
                void SetListMatSexp(int*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetListMatSexp(uint*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetListMatSexp(double*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP) ;
                /*
                * Remplit une liste de theNElt matrice de tailles diffr�ntes theLigne[i] x theCol[i] dans un SEXP
                */
                void SetListMatSexp(int*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP) ;
                void SetListMatSexp(uint*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP) ;
                void SetListMatSexp(double*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP) ;
                void SetListMatSexp(cDMatrix* theVal, uint theNElt, SEXP &theSEXP) ;

                /*
                * Remplit une liste de theNList1 elements de listes de theNList2 elements de vecteurs dans un SEXP
                */
                void SetListListVectSexp(cDVector** theVect, uint theNList1, uint theNList2, SEXP &theSEXP) ;

                /*
                * Remplit une liste de theNList1 elements de listes de theNList2 elements de matrices dans un SEXP
                */
                void SetListListMatSexp(cDMatrix** theMat, uint theNList1, uint theNList2, SEXP &theSEXP) ;
                void SetListListMatSexp(cDMatrix** theMat, uint theNList1, uint* theNList2, SEXP &theSEXP) ;

        } ;

#endif // _RDLL_ 
#endif //_CRUTILS_H_
