/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: StatUtil.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

cDMatrix AddOneVariable(cDMatrix& theCovCour, cDVector& theVect)
{
uint	mySize = theVect.GetSize() ;
cDVector myU(mySize, 0.0) ;
cDVector myRes(mySize+1, 0.0) ;
cDMatrix myAuxMat = Transpose(theVect) ;
	for (register uint i= 0; i < mySize ; i++)
	{	myU[i] = 1.0 ;
		myRes[i] = AsDouble(myAuxMat * theCovCour * myU) ;
		myU[i] = 0.0 ;
	}
	myU = 1.0 ;
	myRes[mySize] = AsDouble(myAuxMat * theCovCour * theVect) ;
	AddColRow(myRes, theCovCour) ;
	return theCovCour ;
}

