/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDVector.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/


#ifndef CDVECTOR_H_
#define CDVECTOR_H_

#include <cmath> 
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "cOTError.h"
#include "ExitError.h"

#define D_PRECISION 16

class cDMatrix ;

class cDVector 
{
private:
    double* mvV ;                  
    double* mvVm1 ;        // pointer adjustment for optimzied 1-offset indexing
    uint mvN ;

    // internal helper function to create the array
    // of row pointers  

    void Initialize(uint theN) ;
    void Copy(const double*  theV) ;
    void Set(const double& theVal) ;
    // access

  public:
        void Delete() ;

    // destructor
    virtual ~cDVector() ;

    // constructors
    cDVector() ;
    cDVector(const cDVector& theVect) ;
    cDVector(uint theN, const double& theValue = 0.0) ;
    cDVector(uint theN, const double* theVect) ;

    // methods
    // 
    void ReAlloc(uint theN) ;
    void ReAlloc(uint theSize, double theVal) ;
    void ReAlloc(uint theSize, double* theVect) ;
   
    // assignments
    //
    cDVector& operator=(const cDVector& theVect) ;
    cDVector& operator=(const double* theVect) ;
        cDVector& operator=(const double& theScalar) ;
    uint GetSize(void) const ;
        double* GetVect(void) const ;
        uint Dim(void) const ;
    uint Size(void) const ;

    inline double& operator()(int theIndex) ;
    inline const double& operator() (int theIndex) const ;
    double& operator[](int theIndex) ;
    const double& operator[](int theIndex) const ;

        friend bool operator ==(const cDVector& theSrcVect, const cDVector& theCompVect) ;
        friend bool operator <(const cDVector& theSrcVect, const cDVector& theCompVect) ;
        friend bool operator <=(const cDVector& theSrcVect, const cDVector& theCompVect) ;
        friend bool operator >(const cDVector& theSrcVect, const cDVector& theCompVect) ;
        friend bool operator >=(const cDVector& theSrcVect, const cDVector& theCompVect) ;
        friend bool operator !=(const cDVector& theSrcVect, const cDVector& theCompVect) ;
} ;


extern bool operator ==(const cDVector& theSrcVect, const cDVector& theCompVect) ;
extern bool operator !=(const cDVector& theSrcVect, const cDVector& theCompVect) ;
extern std::ostream& operator<<(std::ostream& theStream, const cDVector& theVect) ;
extern cDVector operator*(const cDVector& theVect, const double& theValue) ;
extern cDVector operator *=(cDVector& theVect, const double& theValue) ;
extern cDVector operator*(const double& theVal, const cDVector& theVect) ;
extern cDVector operator/(const cDVector& theVect, const double& theVal) ;
extern cDVector operator /=(cDVector& theVect, const double& theVal) ;
extern cDVector operator+(const cDVector& theLeftVect, const cDVector& theRightVect) ;
extern cDVector operator += (cDVector& theVect, const cDVector& theRightVect) ;
extern cDVector operator-(const cDVector& theVect, const cDVector& theRightVect) ;
extern cDVector operator-=(cDVector& theVect, const cDVector& theRightVect) ;
extern double AsDouble(const cDVector& theVect) ;
extern cDVector Zeros(uint theN) ;
extern void GetSubVector(const cDVector& theSrcVector, uint theFirstIndex, uint theSize, cDVector& theDestVector) ;
extern void SetSubVector(const cDVector& theSrcVector, uint theFirtIndex, cDVector& theDestVector) ;
extern cDVector CopyDouble(double* theVect, uint theSize) ;
extern double ScalarProduct(cDVector& theVect1, cDVector& theVect2) ;
extern cDVector cat(const cDVector& theVect1, const cDVector& theVect2) ;
extern cDVector cat(const cDVector& theVect, double theDouble) ;
extern cDVector cat(double theDouble, const cDVector& theVect) ;
#endif // CDVECTOR_H_



