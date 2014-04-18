/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDMatrix.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "cDMatrix.h"

void cDMatrix::Initialize(uint theNRow, uint theNCol)
{       mvSize = theNRow*theNCol ;
    mvNRow = theNRow ;
        mvNCol = theNCol ;
        if ( (mvV = new double[mvSize]) == NULL) 
                throw cOTError("Memory allocation failed") ;
        if ( (mvRow = new double*[theNRow]) == NULL)
                throw cOTError("Memory allocation failed") ;
        if ( (mvRowm1 = new double*[theNRow]) == NULL)
                throw cOTError("Memory allocation failed") ;
        double* myPointer = mvV ;             
        mvVm1 = mvV - 1 ;
        for (uint i = 0 ; i < theNRow ; i++)
        {       mvRow[i] = myPointer ;
            mvRowm1[i] = myPointer-1 ;
            myPointer += theNCol ;
        }
        mvRowm1 -- ;
}
   
void cDMatrix::Copy(const double* theFlatMat)
{
uint mySize = mvNRow * mvNCol ;
        for (uint i = 0 ; i < mySize ; i++)
            mvV[i] = theFlatMat[i] ;      
}

void cDMatrix::Set(const double& theVal)
{
uint mySize = mvNRow * mvNCol ;
        for (uint i = 0 ; i < mySize ; i++)
        mvV[i] = theVal ;
}
    
void cDMatrix::Delete()
{     
        if (mvV == NULL) 
                return ;
        else
                delete [] (mvV) ;     
        if (mvRow != NULL) 
                delete [] (mvRow) ;
        mvRowm1 ++ ;
        if (mvRowm1 != NULL ) 
                delete [] (mvRowm1) ;
        mvV = NULL ;
        mvRow = NULL ;
        mvRowm1 = NULL ;
}

cDMatrix::operator double**()
{       return  mvRow ; 
}

cDMatrix::operator double**() const 
{       return mvRow ; 
}

int cDMatrix::Size() const
{       return mvSize ; 
}

cDMatrix::cDMatrix() : mvNRow(0), mvNCol(0), mvSize(0), mvV(0), mvRow(0), mvVm1(0), mvRowm1(0) 
{
}

cDMatrix::cDMatrix(const cDMatrix& theMatrix)
{
        Initialize(theMatrix.mvNRow, theMatrix.mvNCol) ;
        Copy(theMatrix.mvV) ;
}

cDMatrix::cDMatrix(uint theNRow, uint theNCol, const double& theValue)
{
        Initialize(theNRow, theNCol) ;
        Set(theValue) ;
}

cDMatrix::cDMatrix(uint theNRow, uint theNCol, const double* theFlatMat)
{
        Initialize(theNRow, theNCol) ;
        Copy(theFlatMat) ;
}

cDMatrix::~cDMatrix()
{
        Delete() ;
}

cDMatrix& cDMatrix::ReAlloc(uint theNRow, uint theNCol)
{       if (GetNRows() == theNRow && GetNCols() == theNCol)
            return *this ;
        Delete() ;
        Initialize(theNRow, theNCol) ;       
        return *this ;
}

cDMatrix& cDMatrix::operator =(const cDMatrix& theSrcMat)
{       if (mvV == theSrcMat.mvV)
            return *this ;
        if (mvNRow == theSrcMat.mvNRow  && mvNCol == theSrcMat.mvNCol)
            Copy(theSrcMat.mvV) ;
        else
        {
            Delete() ;
            Initialize(theSrcMat.mvNRow, theSrcMat.mvNCol) ;
            Copy(theSrcMat.mvV) ;
        }
        return *this ;
}

cDMatrix& cDMatrix::operator =(const cDVector& theVect)
{       
uint myNRow = theVect.GetSize() ;
        if (mvNRow == myNRow  && mvNCol == 1)
            Copy(theVect.GetVect()) ;
        else
        {
            Delete() ;
            Initialize(myNRow, 1) ;
            Copy(theVect.GetVect()) ;
        }
        return *this ;
}

cDMatrix& cDMatrix::operator =(const double& theScalar)
{ 
        Set(theScalar) ; 
        return *this ;
}

uint cDMatrix::GetNRows(void) const 
{       
        return mvNRow ; 
}

uint cDMatrix::GetNCols(void) const 
{
        return mvNCol ;
}

double* cDMatrix::GetCol (uint theIndex)
{
double* myTemp = new double(mvNRow) ;
        for(uint i = 0 ; i < mvNRow ; i++) 
                myTemp[i] = mvRow[i][theIndex] ;
        return myTemp ;
}

double& cDMatrix::operator()(int i)
{ 
        return mvV[i] ; 
}

const double& cDMatrix::operator()(int i) const
{ 
        return mvV[i] ; 
}

double& cDMatrix::operator()(int i, int j)
{ 
        return  mvRow[i][j] ; 
}
    
const double& cDMatrix::operator() (int i, int j) const
{
        return  mvRow[i][j] ; 
}

std::ostream& operator << (std::ostream& theStream, const cDMatrix& theMatrix)
{
uint theNRow = theMatrix.GetNRows() ;
uint theNCol = theMatrix.GetNCols() ;
    
    for (uint i = 0 ; i < theNRow ; i++)
    {
        for (uint j = 0 ; j < theNCol - 1 ; j++)
        {
          theStream << theMatrix[i][j] << " " ;
        }
        theStream << theMatrix[i][theNCol-1] << "\n" ;
    }

    return theStream ;
}

cDMatrix operator +(const cDMatrix& theLeftMat, const cDMatrix& theRightMat)
{
uint myNRow = theLeftMat.GetNRows() ;
uint myNCol = theLeftMat.GetNCols() ;
    myAssert(myNRow == theRightMat.GetNRows(), "operator +: Matrices must have the same dimensions") ;
    myAssert(myNCol == theRightMat.GetNCols(), "operator +: Matrices must have the same dimensions") ;
cDMatrix myTmpMat(myNRow, myNCol) ;
    for (uint i = 0 ; i < myNRow ; i++)
        for (uint j = 0 ; j < myNCol ; j++)
            myTmpMat[i][j] = theLeftMat[i][j] + theRightMat[i][j] ;
    return myTmpMat ;
}

cDMatrix operator +=(cDMatrix& theLeftMat, const cDMatrix& theRightMat)
{
uint myNRow = theLeftMat.GetNRows() ;
uint myNCol = theLeftMat.GetNCols() ;
    myAssert(myNRow == theRightMat.GetNRows(), "operator +=: Matrices must have the same dimensions") ;
    myAssert(myNCol == theRightMat.GetNCols(), "operator +=: Matrices must have the same dimensions") ;

        for (uint i = 0 ; i < myNRow ; i++)
        for (uint j = 0 ; j < myNCol ; j++)
            theLeftMat[i][j] += theRightMat[i][j] ;
    return theLeftMat ;
}

cDMatrix operator-(const cDMatrix& theLeftMat, const cDMatrix& theRightMat)
{
uint myNRow = theLeftMat.GetNRows() ;
uint myNCol = theLeftMat.GetNCols() ;
    myAssert(myNRow == theRightMat.GetNRows(), "operator -: Matrices must have the same dimensions") ;
    myAssert(myNCol == theRightMat.GetNCols(), "operator -: Matrices must have the same dimensions") ;
cDMatrix myTmpMat(myNRow, myNCol) ;
    for (uint i = 0 ; i < myNRow ; i++)
        for (uint j = 0 ; j < myNCol ; j++)
                        myTmpMat[i][j] = theLeftMat[i][j] - theRightMat[i][j] ;
    return myTmpMat ;
}

cDMatrix operator-=(cDMatrix& theLeftMat, const cDMatrix& theRightMat)
{
uint myNRow = theLeftMat.GetNRows() ;
uint myNCol = theLeftMat.GetNCols() ;
    myAssert(myNRow == theRightMat.GetNRows(), "operator -=: Matrices must have the same dimensions") ;
	myAssert(myNCol == theRightMat.GetNCols(), "operator -=: Matrices must have the same dimensions") ;
    for (uint i = 0 ; i < myNRow ; i++)
        for (uint j = 0 ; j < myNCol ; j++)
                        theLeftMat[i][j] -= theRightMat[i][j] ;
    return theLeftMat ;
}

cDMatrix Transpose(const cDMatrix& theLeftMat)
{
uint myNRow = theLeftMat.GetNRows() ;
uint myNCol = theLeftMat.GetNCols() ;
cDMatrix myTmpMat(myNCol, myNRow) ;
   for (uint i = 0 ; i < myNRow ; i++)
        for (uint j = 0 ; j < myNCol ; j++)
            myTmpMat[j][i] = theLeftMat[i][j] ;

    return myTmpMat ;
}

cDMatrix Transpose(const cDVector& theVect)
{
uint myNCol = theVect.GetSize() ;
cDMatrix myTmpMat(1, myNCol) ;
    for (uint i = 0 ; i < myNCol ; i++)
                myTmpMat[0][i] = theVect[i] ;
    return myTmpMat ;
}

cDVector AsVector(const cDMatrix& theMat)
{
uint n ;
	myAssert((theMat.GetNCols() == 1) || (theMat.GetNRows() == 1), "AsVector: Matrix must have one row or one column") ;
	if ( (n = theMat.GetNCols()) == 1 )
	{	
	cDVector myVect(n) ;
		for ( register uint i = 0 ; i < n ; i++)
			myVect[i] = theMat[i][0] ;
		return myVect ;
	}
	else
	{	
	cDVector myVect(n) ;
		for ( register uint i = 0 ; i < n ; i++)
			myVect[i] = theMat[0][1] ;
		return myVect ;
	}
}

static cDMatrix MatMult(const cDMatrix  &theLeftMat, const cDMatrix& theRightMat)
{
uint myNRow = theLeftMat.GetNRows() ;
uint myNCol = theLeftMat.GetNCols() ;
uint myRightNCol = theRightMat.GetNCols() ;
cDMatrix myTmpMat(myNRow, myRightNCol) ;
double mySum ;

    for (uint i = 0 ; i < myNRow ; i++)
                for (uint k = 0 ; k < myRightNCol ; k++)
                {       mySum = 0.0 ;
                        for (uint j = 0 ; j < myNCol ; j++)
                                mySum = mySum +  theLeftMat[i][j] * theRightMat[j][k] ;
                        myTmpMat[i][k] = mySum ; 
                }
    return myTmpMat ;
}

cDMatrix operator*(const cDMatrix& theLeftMat, const cDMatrix& theRightMat)
{
  return MatMult(theLeftMat, theRightMat) ;
}

static cDMatrix ScalMult(const cDMatrix  &theLeftMat, const double& theVal)
{
  cDMatrix myTmpMat = theLeftMat ;
  uint myNRow = theLeftMat.GetNRows() ;
  uint myNCol = theLeftMat.GetNCols() ;

  for(uint i = 0 ;i < myNRow ;i++)
    for(uint j = 0 ;j < myNCol ;j++)
      myTmpMat[i][j] *= theVal ;
  return myTmpMat ;
}

cDMatrix operator*(const cDMatrix& theMat, const double& theVal)
{
  return ScalMult(theMat, theVal) ;
}

cDMatrix operator*=(cDMatrix& theMat, const double& theVal)
{
  theMat = ScalMult(theMat, theVal) ;
  return theMat ;
}

cDMatrix operator/(const cDMatrix& theMat, const double& theVal)
{
  if (theVal == 0.0)
          throw cOTError("division by 0") ;
  return ScalMult(theMat, 1.0/theVal) ;
}

cDMatrix operator/=(cDMatrix& theMat, const double& theVal)
{
  if (theVal == 0.0)
          throw cOTError("division by 0") ;
  theMat = ScalMult(theMat, 1.0/theVal) ;
  return theMat ;
}

cDMatrix operator*(const double& theVal, const cDMatrix& theMat)
{
  return ScalMult(theMat, theVal) ;
}

static cDMatrix MatMult(const cDVector& theVect, const cDMatrix& theMat)
{
uint myNRow = theVect.GetSize() ;
uint myNCol = theMat.GetNCols() ;
cDMatrix myTmpMat(myNRow, myNCol) ;
    
    for (uint i = 0 ; i < myNCol ; i++)
    {   for (uint j = 0 ; j < myNRow ; j++)
            myTmpMat[i][j] = theVect[i] * theMat[0][j] ;
    }
    return myTmpMat ;
}

cDMatrix operator*(const cDVector &theVect, const cDMatrix& theMat)
{
  return MatMult(theVect, theMat) ;
}

static cDVector MatMult(const cDMatrix& theMat, const cDVector& theVect)
{
uint myNRow = theMat.GetNRows() ;
uint myNCol = theMat.GetNCols() ;
cDVector myTmpVect(myNRow) ;
double mySum ;

    for (uint i = 0 ; i < myNRow ; i++)
    {   mySum = 0.0 ;
    const double* myRow = theMat[i] ;
        for (uint j = 0 ; j < myNCol ; j++)
            mySum = mySum +  myRow[j] * theVect[j] ;
        myTmpVect[i] = mySum ; 
    }
    return myTmpVect ;
}

cDVector operator*(const cDMatrix& theMat, const cDVector& theVect)
{
    return MatMult(theMat, theVect) ;
}

cDMatrix Zeros(uint theN, uint theP)
{
cDMatrix myTempMatrix(theN, theP) ;
        return myTempMatrix ;
}

double AsDouble(const cDMatrix& theMat)
{
	return theMat[0][0] ;
}

cDMatrix Identity(uint theN)
{
cDMatrix myTempMatrix(theN, theN) ;
        for (register uint i=0 ; i < theN ; i++)
                myTempMatrix[i][i] = 1.0L ;
        return myTempMatrix ;
}

cDMatrix Diag(cDVector& theVect)
{
uint mySize = theVect.GetSize() ;
cDMatrix myTempMatrix(mySize, mySize) ;
        for (register uint i = 0 ; i <mySize ; i++)
                myTempMatrix[i][i] = theVect[i] ;

        return myTempMatrix ;
}

cDMatrix Inv(cDMatrix& theMatrix)
{
uint myNCol = theMatrix.GetNCols() ;
cDMatrix myTempMatrix(myNCol, myNCol) ;
double myDet ;

        LapackInvAndDet(theMatrix, myTempMatrix, myDet) ;
        if (std::fabs(myDet) < MIN_DBLE)
                throw cOTError("Non inversible matrix") ;
        return myTempMatrix ;
}

void GetSubMatrix(cDMatrix& theSrcMatrix, uint theSize, cDMatrix& theDestMatrix)
{
	GetSubMatrix(theSrcMatrix, theSize, theSize, theDestMatrix) ;
}

void GetSubMatrix(cDMatrix& theSrcMatrix, uint theNRow, uint theNCol, cDMatrix& theDestMatrix)
{
	if ( (theNRow > theSrcMatrix.GetNRows()) || (theNCol > theSrcMatrix.GetNCols()) )
		throw cOTError("Wrong matrix size in GetSubMatrix") ;
	theDestMatrix.ReAlloc(theNRow, theNCol) ;
	
	for (register uint i = 0 ; i < theNRow ; i++)
		for (register uint j = 0 ; j < theNCol ; j++)
			theDestMatrix[i][j] = theSrcMatrix[i][j] ;
}

void SetSubMatrix(cDMatrix& theSrcMatrix, uint theFirstRow, uint theFirstCol, cDMatrix& theDestMatrix)
{
uint myNRows = theSrcMatrix.GetNRows() ;
uint myNCols = theSrcMatrix.GetNCols() ;

	if ( (theDestMatrix.GetNRows()  < myNRows + theFirstRow) || (theDestMatrix.GetNCols() < myNCols + theFirstCol) )
		throw cOTError("Wrong matrix size in SetSubMatrix") ;

	for (register uint i = 0 ; i < myNRows ; i++)
		for (register uint j = 0 ; j < myNCols ; j++)
			theDestMatrix[i+theFirstRow][j+theFirstCol] = theSrcMatrix[i][j] ;
}

void AddColRow(const cDVector& theColRow, cDMatrix& theMat) 
{
uint myNRow = theMat.GetNRows() ;
uint myNCol = theMat.GetNCols() ;
uint mySize = theColRow.GetSize() ;
	if ( (myNRow != myNCol) || (myNRow + 1 != mySize) )
		throw cOTError("Wrong sizes in AddColRow") ;
cDMatrix mySrcMatrix = theMat ;
	theMat.ReAlloc(mySize, mySize) ;
	SetSubMatrix(mySrcMatrix, 0, 0, theMat) ;
	for (register uint i = 0 ; i < mySize  ; i++)
		theMat[i][mySize-1] = theMat[mySize-1][i] = theColRow[i] ;
}

void LapackInvAndDet(cDMatrix& theMatrix, cDMatrix& theInvMatrix, double& theDet)
{
uint myNCol = theMatrix.GetNCols() ;

double  *myAP = new double[myNCol*(myNCol + 1)/2],
                *myW = new double[myNCol],
                *myZ = new double[myNCol*myNCol],
                *myWork = new double[myNCol * 3] ;
int myInfo,
        myN = (int)(myNCol),
        myldz = (int)(myNCol) ;

        for (register int i = 0 ; i < myN ; i++)
                for (register int j = i ; j < myldz ; j++)
                        myAP[i+(j+1)*j/2]  = theMatrix[i][j] ;

        F77_NAME(dspev)("V", "U", &myN, myAP, myW, myZ, &myldz, myWork, &myInfo) ;

        if (myInfo != 0)
                throw cOTError("Non inversible matrix") ;
        theDet = 1.0L ;
cDVector myInvEigenValue = cDVector(myNCol) ;

cDMatrix myEigenVector(myNCol, myNCol) ;
        for (register uint i = 0 ; i < myNCol ; i++)
        {       theDet *= myW[i] ;
                myInvEigenValue[i] = 1.0 /myW[i] ;
                for (register int j = 0 ; j < myN ; j++)
                        myEigenVector[i][j] = myZ[i + j*myN] ;
        }
        theInvMatrix =  myEigenVector ;
cDMatrix myAuxMat1 = Diag(myInvEigenValue), myAuxMat2 = Transpose(myEigenVector) ;
cDMatrix myAuxMat = myAuxMat1 * myAuxMat2 ;
        theInvMatrix = theInvMatrix * myAuxMat ;
        
        delete myAP ;
        delete myW ;
        delete myZ ;
        delete myWork ;
}

double LapackDet(cDMatrix& theMatrix)
{
uint myNCol = theMatrix.GetNCols() ;

double *myAP = new double[myNCol*(myNCol + 1)/2],
       *myW = new double[myNCol],
       *myZ = new double[myNCol*myNCol],
       *myWork = new double[myNCol * 3] ;
int myInfo,
    myN = (int)(myNCol),
    myldz = (int)(myNCol) ;

	for (register int i = 0 ; i < myN ; i++)
		for (register int j = i ; j < myldz ; j++)
			myAP[i+(j+1)*j/2]  = theMatrix[i][j] ;

    F77_NAME(dspev)("V", "U", &myN, myAP, myW, myZ, &myldz, myWork, &myInfo) ;

double myDet ;
	if (myInfo != 0)
		myDet = 0.0 ;
	else
	{	myDet = 1.0L ;
        for (register uint i = 0 ; i < myNCol ; i++)
        {	myDet *= myW[i] ;
        }
	}        
    delete myAP ;
    delete myW ;
    delete myZ ;
    delete myWork ;

	return myDet ;
}
