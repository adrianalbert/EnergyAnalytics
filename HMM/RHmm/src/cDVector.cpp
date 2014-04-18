/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDVector.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "cDVector.h"

// cDVector methods

void cDVector::Initialize(uint theN)
{	myAssert(mvV == NULL, "Problem") ;
	if ( (mvV = new double[theN]) == NULL)
		throw cOTError("Memory allocation failed") ;
	mvVm1 = mvV-1 ;
	mvN = theN ;
}
    
void cDVector::Copy(const double*  theV)
{	for (uint i = 0 ; i < mvN ; i++)
		mvV[i] = theV[i] ;    
}
    
void cDVector::Set(const double& theVal)
{
	for (uint i = 0 ; i< mvN ; i++)
	    mvV[i] = theVal ;
}
    
void cDVector::Delete()
{     
/* do nothing, if no memory has been previously allocated */
	if (mvV == NULL) 
		return ;

/* if we are here, then matrix was previously allocated */
	delete [] (mvV) ;     

	mvV = NULL ;
	mvVm1 = NULL ;
	mvN = 0 ;
}

cDVector::~cDVector() 
{
	Delete() ;
}

cDVector::cDVector() : mvV(0), mvVm1(0), mvN(0)  
{
} 

cDVector::cDVector(const cDVector& theVect) : mvV(0), mvVm1(0), mvN(0)
{
	Initialize(theVect.mvN) ;
	Copy(theVect.mvV) ;
}

cDVector::cDVector(uint theN, const double& theValue) :  mvV(0), mvVm1(0), mvN(0)
{
	Initialize(theN) ;
	Set(theValue) ;
}

cDVector::cDVector(uint theN, const double* theVect) :  mvV(0), mvVm1(0), mvN(0)
{
	Initialize(theN) ;
	Copy(theVect) ;
}

void cDVector::ReAlloc(uint theN)
{	if (mvN != theN) 
	{
		Delete() ;
		Initialize(theN) ;
	}
}

void cDVector::ReAlloc(uint theSize, double theVal)
{	if (mvN != theSize) 
	{	Delete() ;
		Initialize(theSize) ;
	}
	for (register uint i = 0 ; i < theSize ; i++)
		mvV[i] = theVal ; 
}
 
void cDVector::ReAlloc(uint theSize, double* theVect)
{
	if (mvN != theSize) 
	{	Delete() ;
		Initialize(theSize) ;
	}
	for (register uint i = 0 ; i < theSize ; i++)
		mvV[i] = theVect[i] ; 
}

    cDVector& cDVector::operator=(const cDVector& theVect)
    {
	if (mvV == theVect.mvV)
	    return *this ;

	if (mvN == theVect.mvN)	 // no need to re-alloc
	    Copy(theVect.mvV) ;
	else
	{
	    Delete() ;
	    Initialize(theVect.mvN) ;
	    Copy(theVect.mvV) ;
	}

	return *this ;
    }       
 
	cDVector& cDVector::operator=(const double* theVect)
    {
	    Copy(theVect) ;
	return *this ;
    }       
 
	cDVector& cDVector::operator=(const double& theScalar)
    { 
	Set(theScalar) ;  
	return *this ;
    }

	uint cDVector::GetSize(void) const 
    {
	return  mvN ; 
    }

	double* cDVector::GetVect(void) const
	{
		return mvV ;
	}

	uint cDVector::Dim(void) const 
    {
	return  mvN ; 
    }
    
	uint cDVector::Size(void) const 
    {
	return  mvN ; 
    }
	
	bool operator ==(const cDVector& theSrcVect, const cDVector& theCompVect) 
	{
    bool myBool=true ;
		if(theSrcVect.mvN != theCompVect.mvN) 
			return false ;

		for(uint i=0 ; myBool && i < theSrcVect.mvN ;i++)
		{
		  if (theSrcVect.mvV[i] != theCompVect.mvV[i]) 
			  myBool=false ;
		}
		return myBool ;    
	}

	bool operator <(const cDVector& theSrcVect, const cDVector& theCompVect) 
	{
    bool myBool=true ;
		if(theSrcVect.mvN != theCompVect.mvN) 
			return false ;

		for(uint i=0 ; myBool && i < theSrcVect.mvN ;i++)
		{
		  if (theSrcVect.mvV[i] >= theCompVect.mvV[i]) 
			  myBool = false ;
		}
		return myBool ;    
	}

	bool operator <=(const cDVector& theSrcVect, const cDVector& theCompVect) 
	{
    bool myBool=true ;
		if(theSrcVect.mvN != theCompVect.mvN) 
			return false ;

		for(uint i=0 ; myBool && i < theSrcVect.mvN ;i++)
		{
		  if (theSrcVect.mvV[i] > theCompVect.mvV[i]) 
			  myBool = false ;
		}
		return myBool ;    
	}

	bool operator >(const cDVector& theSrcVect, const cDVector& theCompVect) 
	{
    bool myBool=true ;
		if(theSrcVect.mvN != theCompVect.mvN) 
			return false ;

		for(uint i=0 ; myBool && i < theSrcVect.mvN ;i++)
		{
		  if (theSrcVect.mvV[i] <= theCompVect.mvV[i]) 
			  myBool = false ;
		}
		return myBool ;    
	}

	bool operator >=(const cDVector& theSrcVect, const cDVector& theCompVect) 
	{
    bool myBool=true ;
		if(theSrcVect.mvN != theCompVect.mvN) 
			return false ;

		for(uint i=0 ; myBool && i < theSrcVect.mvN ;i++)
		{
		  if (theSrcVect.mvV[i] < theCompVect.mvV[i]) 
			  myBool = false ;
		}
		return myBool ;    
	}

	bool operator !=(const cDVector& theSrcVect, const cDVector& theCompVect) 
	{
		return !(theSrcVect == theCompVect) ;
	}
    
	double& cDVector::operator()(int theIndex)
    { 
	return mvVm1[theIndex] ; 
    }
 
	const double& cDVector::operator() (int theIndex) const
    {
	return mvVm1[theIndex] ; 
    }

	double& cDVector::operator[](int theIndex)
    { 
	return mvV[theIndex] ; 
    }

    const double& cDVector::operator[](int theIndex) const
    {
	return mvV[theIndex] ; 
    }

	std::ostream& operator<<(std::ostream& theStream, const cDVector& theVect)
	{
	uint myN=theVect.GetSize() ;
		for (uint i=0 ; i<myN ; i++)
		  theStream  << theVect[i] << " " << std::endl ;
		return theStream ;
	}

	static cDVector ScalMult(const cDVector& theVect, const double& theValue)
	{
	uint myN = theVect.GetSize() ; 
	cDVector myTmpVect(myN) ;
		for (uint i = 0 ; i < myN ; i++)
				myTmpVect[i] = theVect[i] * theValue ;
		return myTmpVect ; 
	} 

	cDVector operator*(const cDVector& theVect, const double& theValue)
	{
	  return ScalMult(theVect, theValue) ;
	}

	cDVector operator *=(cDVector& theVect, const double& theValue)
	{
		theVect = ScalMult(theVect, theValue) ;
		return theVect ;
	}

	cDVector operator*(const double& theVal, const cDVector& theVect)
	{
	  return ScalMult(theVect, theVal) ;
	}

	cDVector operator/(const cDVector& theVect, const double& theVal)
	{
		if (theVal == 0.0)
			throw cOTError("Division by 0 (cDVector operator '/')") ;
		return ScalMult(theVect, 1.0/theVal) ;
	}

	cDVector operator /=(cDVector& theVect, const double& theVal)
	{
		if (theVal == 0.0)
			throw cOTError("Division by 0 (cDVector operator '/=')") ;
		theVect = ScalMult(theVect, 1/theVal) ;
		return theVect ;
	}

	cDVector operator+(const cDVector& theLeftVect, const cDVector& theRightVect)
	{
	uint myN = theLeftVect.GetSize() ;
		myAssert(myN==theRightVect.GetSize(), "operator +: vectors must have the same dimensions") ;
	cDVector myTmpVect(myN) ;
		for (uint i = 0 ; i < myN ; i++)
			myTmpVect[i] = theLeftVect[i] + theRightVect[i] ;
		return myTmpVect ;
	}

	cDVector operator += (cDVector& theVect, const cDVector& theRightVect)
	{
	uint myN = theVect.GetSize() ;
		myAssert(myN == theRightVect.GetSize(), "operator +=: vectors must have the same dimensions") ;
		for (uint i = 0 ; i < myN ; i++)
				theVect[i] += theRightVect[i] ;
		return theVect ;
	}

	cDVector operator-(const cDVector& theVect, const cDVector& theRightVect)
	{
	uint myN = theVect.GetSize() ;
		myAssert( myN == theRightVect.GetSize(), "operator -: vectors must have the same dimensions") ;
	cDVector myTmpVect(myN) ;
		for (uint i = 0 ; i < myN ; i++)
			myTmpVect[i] = theVect[i] - theRightVect[i] ;
		return myTmpVect ;
	}

	cDVector operator-=(cDVector& theVect, const cDVector& theRightVect)
	{
	uint myN = theVect.GetSize() ;
		myAssert(myN == theRightVect.GetSize(), "operator -=: vectors must have the same dimensions") ;
		for (uint i = 0 ; i < myN ; i++)
			theVect[i] -= theRightVect[i] ;
		return theVect ;
	}

cDVector Zeros(uint theN)
	{
	cDVector myTmpVect(theN, (double)0.0) ;
		return myTmpVect ; 
	}

double AsDouble(const cDVector& theVect)
{ return theVect[0] ;
}
    
cDVector CopyDouble(double* theVect, uint theSize)
{
cDVector myTmpVect(theSize, theVect) ;
	return myTmpVect ; 
}

double ScalarProduct(cDVector& theVect1, cDVector& theVect2)
{
uint mySize = theVect1.GetSize() ;
	if (mySize != theVect2.GetSize())
		throw cOTError("Wrong vector sizes in ScalarProduct") ;

double myRes = 0.0 ;
	for (register uint i = 0 ; i < mySize ; i++)
		myRes += theVect1[i]*theVect2[i] ;
	return myRes ;
}

void GetSubVector(const cDVector& theSrcVector, uint theFirstIndex, uint theSize, cDVector& theDestVector)
{
	if (theSrcVector.GetSize() < theSize + theFirstIndex)
		throw cOTError("Wrong vector size in GetSubVector") ;
	theDestVector.ReAlloc(theSize) ;
	for (register uint i = 0 ; i < theSize ; i++)
		theDestVector[i] = theSrcVector[i+theFirstIndex] ;
}

void SetSubVector(const cDVector& theSrcVector, uint theFirstIndex, cDVector& theDestVector)
{	
uint mySize = theSrcVector.GetSize() ;
	if (theFirstIndex + mySize > theDestVector.GetSize())
		throw cOTError("Wrong vector size in SetSubVector") ;

	for (register uint i = 0 ; i < mySize ; i++)
		theDestVector[i+theFirstIndex] = theSrcVector[i] ;
}

cDVector cat(const cDVector& theVect1, const cDVector& theVect2)
{
uint mySize1 = theVect1.GetSize() ;
uint mySize2 = theVect2.GetSize() ;
cDVector myVect(mySize1+mySize2) ;
	SetSubVector(theVect1, 0, myVect) ;
	SetSubVector(theVect2, mySize1, myVect) ;
	return myVect ;
}

cDVector cat(const cDVector& theVect, double theDouble)
{
uint mySize = theVect.GetSize() ;

cDVector myVect(mySize+1) ;
	SetSubVector(theVect, 0, myVect) ;
	myVect[mySize] = theDouble ;
	return myVect ;

}

cDVector cat(double theDouble, const cDVector& theVect)
{
uint mySize = theVect.GetSize() ;

cDVector myVect(mySize + 1) ;
	SetSubVector(theVect, 1, myVect) ;
	myVect[0] = theDouble ;
	return myVect ;
}
