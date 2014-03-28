/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDMatrix.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/


#ifndef CDMATRIX_H
#define CDMATRIX_H
#pragma once

#ifndef MIN_DBLE
        #define MIN_DBLE 1e-16L
#endif //MIN_DBLE


#include "cOTError.h"
#include "cDVector.h"
#include "R_ext/Lapack.h"
#ifndef D_PRECISION
        #define D_PRECISION 16
#endif //D_PRECISION


class cDMatrix 
{
private :
                uint mvNRow ;
                uint mvNCol ;
                uint mvSize ;      // total size
                double* mvV ;                  
                double** mvRow ;           
                double* mvVm1 ;       // these point to the same data, but are 1-based 
                double** mvRowm1 ;

    // internal helper function to create the array
    // of GetRow pointers

    void Initialize(uint theNRow, uint theNCol) ;  
    void Copy(const double* theFlatMat) ;
    void Set(const double& theVal) ;
public:

    operator double**() ;
        operator double**() const ;
    int Size() const ;
    void Delete() ;


    // constructors
    cDMatrix() ;
    cDMatrix(const cDMatrix& theMatrix) ;
    cDMatrix(uint theNRow, uint theNCol, const double& theValue = (double)0.0) ;
    cDMatrix(uint theNRow, uint theNCol, const double* theFlatMat) ;

    // destructor
    //
    virtual ~cDMatrix() ;

    // reallocating
    //
    cDMatrix& ReAlloc(uint theNRow, uint theNCol) ;

    // assignments
    //
    cDMatrix& operator =(const cDMatrix& theSrcMat) ;
        cDMatrix& operator =(const cDVector& theVect) ;
    cDMatrix& operator =(const double& theScalar) ;
    uint GetNRows(void) const ;
        uint GetNCols(void) const ;
        double* GetCol (uint theIndex) ;
        double& operator()(int i) ;
        const double& operator()(int i) const ;
        double& operator()(int i, int j) ;
        const double& operator() (int i, int j) const ;
} ;


extern std::ostream& operator << (std::ostream& theStream, const cDMatrix& theMatrix) ;
extern cDMatrix operator +(const cDMatrix& theLeftMat, const cDMatrix& theRightMat) ;
extern cDMatrix operator +=(cDMatrix& theLeftMat, const cDMatrix& theRightMat) ;
extern cDMatrix operator-(const cDMatrix& theLeftMat, const cDMatrix& theRightMat) ;
extern cDMatrix operator-=(cDMatrix& theLeftMat, const cDMatrix& theRightMat) ;
extern cDMatrix Transpose(const cDMatrix& theLeftMat) ;
extern cDMatrix Transpose(const cDVector& theVect) ;
extern cDMatrix operator*(const cDMatrix& theLeftMat, const cDMatrix& theRightMat) ;
extern cDMatrix operator*(const cDMatrix& theLeftMat, const double& theVal) ;
extern cDMatrix operator*=(cDMatrix& theLeftMat, const double& theVal) ;
extern cDMatrix operator/(const cDMatrix& theLeftMat, const double& theVal) ;
extern cDMatrix operator/=(cDMatrix& theLeftMat, const double& theVal) ;
extern cDMatrix operator*(const double& theVal, const cDMatrix& theLeftMat) ;
extern cDMatrix operator*(const cDVector &theVect, const cDMatrix& theLeftMat) ;
extern cDVector operator*(const cDMatrix& theMat, const cDVector& theVect) ;
extern cDVector AsVector(const cDMatrix& theMat) ;
extern double AsDouble(const cDMatrix& theMat) ;
extern cDMatrix Zeros(uint theN, uint theP) ;
extern cDMatrix Identity(uint theN) ;
extern cDMatrix Diag(cDVector& theVect) ;
extern cDMatrix Inv(cDMatrix& theMatrix) ;
extern void GetSubMatrix(cDMatrix& theSrcMatrix, uint theSize, cDMatrix& theDestMatrix) ;
extern void GetSubMatrix(cDMatrix& theSrcMatrix,  uint theNRow, uint theNCol, cDMatrix& theDestMatrix) ;
extern void SetSubMatrix(cDMatrix& theSrcMatrix, uint theFirtRow, uint theFirstCol, cDMatrix& theDestMatrix) ;
extern void LapackInvAndDet(cDMatrix& theMatrix, cDMatrix& theInvMatrix, double& theDet) ;
extern double LapackDet(cDMatrix& theMatrix) ;
extern void AddColRow(const cDVector& theColRow, cDMatrix& theMat) ;

#endif // CDMATRIX_H

