/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: RHmm.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

#ifdef _RDLL_



BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RBaumWelch (SEXP theParamHMM, SEXP theYt, SEXP theParamBW)
{
initEnum myTypeInit ;
distrDefinitionEnum myDistrType ;
uint myNbClasses, myDimObs, myNbMixt, myNbProba  ;
cRUtil myRUtil ;
// Retrieve parameters for given HMM
	myRUtil.GetValSexp(theParamHMM, eNClasses, myNbClasses) ;
	myRUtil.GetValSexp(theParamHMM, eObsDim, myDimObs) ;
	myRUtil.GetValSexp(theParamHMM, eNMixt, myNbMixt) ;
	myRUtil.GetValSexp(theParamHMM, eNProba, myNbProba) ;

char myStr[255] ;
char *myString = (char *)myStr ;
	myRUtil.GetValSexp(theParamHMM, eDistrType, myString) ;
	myDistrType = eUnknownDistr ;
	if (strcmp(myString, "NORMAL") == 0)
	{	if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myString, "DISCRETE") == 0)
			myDistrType = eDiscreteDistr ;
		else
		{	if (strcmp(myString, "MIXTURE") == 0)
				if (myDimObs==1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
		}
	}

	myRUtil.GetValSexp(theParamBW, eInitType, (char *)myString) ;
	if (strcmp(myString, "RANDOM") == 0)
		myTypeInit = eRandom ;
	else
	{	if (strcmp(myString, "KMEANS") == 0)
			myTypeInit = eKMeans ;
		else
			if (strcmp(myString, "USER") == 0)
				myTypeInit = eUser ;
	}

uint myNbIterMax ;
	myRUtil.GetValSexp(theParamBW, eNMaxIter, myNbIterMax) ;
double myTol ;
	myRUtil.GetValSexp(theParamBW, eTol, myTol) ;
uint myVerbose ;
	myRUtil.GetValSexp(theParamBW, eVerbose, myVerbose) ;
uint myNbIterInit ;
	myRUtil.GetValSexp(theParamBW, eNInitIter, myNbIterInit) ;
uint myNbIterMaxInit ;
	myRUtil.GetValSexp(theParamBW, eNMaxIterinit, myNbIterMaxInit) ;

uint myNbSample = length(theYt) ;

cDVector* myY = new cDVector[myNbSample] ;

	for (register uint n = 0 ; n < myNbSample ; n++)
	{       
	SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myY[n].ReAlloc(length(myAux), REAL(myAux)) ;
	}

cBaumWelchInParam myParamEntree = cBaumWelchInParam(myNbSample, myDimObs, myY, myDistrType, myNbClasses, myNbMixt, myNbProba) ;

	myParamEntree.mNMaxIter = myNbIterMax ;
	myParamEntree.mTol = myTol ;
	myParamEntree.mVerbose = myVerbose ;
	myParamEntree.mNInitIter = myNbIterInit ;
	myParamEntree.mNMaxIterInit = myNbIterMaxInit ;
	myParamEntree.mInitType = myTypeInit ;  

cHmmFit myParamSortie = cHmmFit(myParamEntree) ;

	if (myTypeInit == eUser)
	{	cHmm myHMM = cHmm(myParamEntree) ;
	SEXP myAux1 ;
		myRUtil.GetValSexp(theParamBW, eInitPoint, myAux1) ;
		myRUtil.GetVectSexp(myAux1, 0, myHMM.mInitProba) ;
		myRUtil.GetMatSexp(myAux1, 1, myHMM.mTransMatVector[0]) ; /* FIXME */
		if (myHMM.mTransMatVector.size() > 1)
			warning("Time-inhomogeneous Markov chain not supported yet for BaumWelch algorithm.");

	SEXP myAux ;
		myRUtil.GetValSexp(myAux1, 2, myAux) ; // $distribution
		switch (myDistrType)
		{	case eNormalDistr :
			{	cUnivariateNormal* myParam = dynamic_cast<cUnivariateNormal *>(myHMM.mDistrParam) ;                     
				myRUtil.GetVectSexp(myAux, 3, myParam->mMean) ;
				myRUtil.GetVectSexp(myAux, 4, myParam->mVar) ;
			}
			break ;
			case eMultiNormalDistr :
			{	cMultivariateNormal* myParam = dynamic_cast<cMultivariateNormal *>(myHMM.mDistrParam) ;
				myRUtil.GetListVectSexp(myAux, 3, myNbClasses, myParam->mMean) ; 
				myRUtil.GetListMatSexp(myAux, 4, myNbClasses, myParam->mCov) ;
			}
			break ;
			case eMixtUniNormalDistr :
			{	cMixtUnivariateNormal* myParam = dynamic_cast<cMixtUnivariateNormal *>(myHMM.mDistrParam) ;
				myRUtil.GetListVectSexp(myAux, 4, myNbClasses, myParam->mMean) ;
				myRUtil.GetListVectSexp(myAux, 5, myNbClasses, myParam->mVar) ;
				myRUtil.GetListVectSexp(myAux, 6, myNbClasses, myParam->mp) ;
			}
			break ;
			case eMixtMultiNormalDistr :
			{	cMixtMultivariateNormal* myParam = dynamic_cast<cMixtMultivariateNormal *>(myHMM.mDistrParam) ;
				myRUtil.GetListListVectSexp(myAux, 4, myNbClasses, myNbMixt, myParam->mMean) ;
				myRUtil.GetListListMatSexp(myAux, 5, myNbClasses, myNbMixt, myParam->mCov) ;
				myRUtil.GetListVectSexp(myAux, 6, myNbClasses, myParam->mp) ;
			}
			break ;
			case eDiscreteDistr :
			{	cDiscrete* myParam = dynamic_cast<cDiscrete *>(myHMM.mDistrParam) ;
				myRUtil.GetEmissionSexp(myAux, 3, myParam->mProbaMatVector);

				if (myParam->mProbaMatVector.size() > 1)
					warning("Variable discrete emission probabilities not supported yet for BaumWelch algorithm.");
			}
			break ;
		}
		myParamSortie.CopyHmm(myHMM) ;
	}
	else
		myParamSortie.BaumWelchAlgoInit(myParamEntree) ;        
	myParamSortie.BaumWelchAlgo(myParamEntree) ;
	
	for (register uint n = 0 ; n < myNbSample ; n++)
		myY[n].Delete() ;
	delete[] myY ;


SEXP myRes,
	 myAux[6]       ;

myRUtil.SetVectSexp(myParamSortie.mInitProba, myAux[0]) ;
	myRUtil.SetMatSexp(myParamSortie.mTransMatVector[0], myAux[1]) ; // TODO: Do it for the real list

	switch (myDistrType)
	{	case eNormalDistr :
		{	cUnivariateNormal* myParam =  dynamic_cast<cUnivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetVectSexp(myParam->mMean, myAux[2]) ;
			myRUtil.SetVectSexp(myParam->mVar, myAux[3]) ;
		}
		break ;
		case eMultiNormalDistr :
		{	cMultivariateNormal* myParam = dynamic_cast<cMultivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListVectSexp(myParam->mMean, myNbClasses, myAux[2]) ;
			myRUtil.SetListMatSexp(myParam->mCov, myNbClasses, myAux[3]) ;
		}
		break ;
		case eDiscreteDistr :
		{	cDiscrete* myParam = dynamic_cast<cDiscrete *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListVectSexp(myParam->mProbaMatVector[0], myAux[2]); // TODO: Do it for the real list
		}
		break ;
		case eMixtUniNormalDistr :
		{	cMixtUnivariateNormal* myParam = dynamic_cast<cMixtUnivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListVectSexp(myParam->mMean, myNbClasses, myAux[2]) ;
			myRUtil.SetListVectSexp(myParam->mVar, myNbClasses, myAux[3]) ;
			myRUtil.SetListVectSexp(myParam->mp, myNbClasses,myAux[4]) ;
		}
		break ;
		case eMixtMultiNormalDistr :
		{	cMixtMultivariateNormal* myParam = dynamic_cast<cMixtMultivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListListVectSexp(myParam->mMean, myNbClasses, myNbMixt, myAux[2]) ;
			myRUtil.SetListListMatSexp(myParam->mCov, myNbClasses, myNbMixt, myAux[3]) ;
			myRUtil.SetListVectSexp(myParam->mp, myNbClasses, myAux[4]) ;
		}
		break ;
		default :
		break ;
	}

	PROTECT(myAux[5] = allocVector(REALSXP, 5)) ;
	REAL(myAux[5])[0] = myParamSortie.mLLH ;
	REAL(myAux[5])[1] = myParamSortie.mBic ;
	REAL(myAux[5])[2] = (double)( myParamSortie.mNIter) ;
	REAL(myAux[5])[3] =  myParamSortie.mTol ;
	REAL(myAux[5])[4] = myParamSortie.mAic ;

	uint myNAlloc ;
	switch (myDistrType)
	{	case eNormalDistr :
		case eMultiNormalDistr :
			myNAlloc = 5 ;    
		   break ;
		case eMixtUniNormalDistr :
		case eMixtMultiNormalDistr :
			myNAlloc = 6 ;    
		  break ;
		case eDiscreteDistr :
		   myNAlloc = 4 ;    
		break ;
		case eUnknownDistr :
		default :
			myNAlloc = 0 ;
		break ;
	}
	PROTECT(myRes = allocVector(VECSXP, myNAlloc)) ;
 
	for (register uint i = 0 ; i < 3 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	for ( register uint i = 3 ; i < myNAlloc-1 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	SET_VECTOR_ELT(myRes, myNAlloc-1, myAux[5]) ;

	UNPROTECT(2) ;
	myRUtil.EndProtect() ;
	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RViterbi(SEXP theHMM, SEXP theYt)
{
distrDefinitionEnum myDistrType ;
uint myDimObs=1, myNbClasses, myNbProba=0, myNbMixt=0 ;
cRUtil myRUtil ;

SEXP myDistSEXP ;
	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba        
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNbClasses) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{	myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
		if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myStr, "DISCRETE") == 0)
		{	myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNbProba) ;
		}
		else
		{	if (strcmp(myStr, "MIXTURE") == 0)
			{	myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
				myRUtil.GetValSexp(myDistSEXP, 2, myNbMixt) ;
			}
		}
	}
uint    myNbSample = length(theYt) ;    
uint*   myT = new uint[myNbSample]      ;
//double        **myY   ;
cDVector* myY = new cDVector[myNbSample] ;
	for (register uint n = 0 ; n < myNbSample ; n++)
	{
	SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm myHMM = cHmm(myDistrType, myNbClasses, myDimObs, myNbMixt, myNbProba) ;

	myRUtil.GetVectSexp(theHMM, fInitProba, myHMM.mInitProba) ;
	myRUtil.GetMatListSexp(theHMM, fTransMat, myHMM.mTransMatVector) ;
	
	switch (myDistrType)
	{	case eNormalDistr :
		{	
		cUnivariateNormal* myLoi = (cUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{
		cMultivariateNormal* myLoi = (cMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNbClasses, myLoi->mCov) ;
		}
		break ;
		case eMixtUniNormalDistr :
		{
		cMixtUnivariateNormal* myParam = (cMixtUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNbClasses, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNbClasses, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;
		case eMixtMultiNormalDistr :
		{
		cMixtMultivariateNormal* myParam = (cMixtMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNbClasses, myNbMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNbClasses, myNbMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{
		cDiscrete* myParam = (cDiscrete *)(myHMM.mDistrParam) ;
			myRUtil.GetEmissionSexp(myDistSEXP, 3, myParam->mProbaMatVector);
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}
				
cInParam myParamEntree(myNbSample, myDimObs, myY) ;
	myParamEntree.mDimObs = myDimObs ;
	myParamEntree.mNMixt = myNbMixt ;
	myParamEntree.mNProba = myNbProba ;
	myParamEntree.mNClass = myNbClasses ;
	myParamEntree.mDistrType = myDistrType ;
cViterbi myViterbi = cViterbi(myParamEntree) ;
	myViterbi.ViterbiPath(myParamEntree, myHMM) ;

SEXP myAux[2] ;
	myRUtil.SetListVectSexp(myViterbi.mSeq, myNbSample, myT, myAux[0]) ;
	myRUtil.SetListValSexp(myViterbi.mLogProb, myAux[1]) ;

	SEXP myRes ;
	PROTECT(myRes = allocVector(VECSXP, 2)) ;
	for (register uint i = 0 ; i < 2 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	myRUtil.EndProtect() ;
	UNPROTECT(1) ;
	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP Rforwardbackward(SEXP theHMM, SEXP theYt, SEXP theLogData)
{
distrDefinitionEnum myDistrType ;
uint  myDimObs=1,
      myNbClasses,
      myNbProba=0,
      myNbMixt=0 ;
cRUtil myRUtil ;


SEXP myDistSEXP ;

int myValBool = *INTEGER(theLogData) ;
	
bool myLogData = (myValBool != 0) ;

	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba        
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNbClasses) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
		if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{if (strcmp(myStr, "DISCRETE") == 0)
		{	myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNbProba) ;
		}
		else
		{	if (strcmp(myStr, "MIXTURE") == 0)
			{	myRUtil.GetValSexp(myDistSEXP, 2, myNbMixt) ;
				myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
			}
		}
	}
uint myNbSample = length(theYt) ;    
uint* myT = new uint[myNbSample] ;

cDVector* myY = new cDVector[myNbSample] ;

	for (register uint n = 0 ; n < myNbSample ; n++)
	{
	SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm myHMM = cHmm(myDistrType, myNbClasses, myDimObs, myNbMixt, myNbProba) ;
	myRUtil.GetVectSexp(theHMM, fInitProba, myHMM.mInitProba) ;
	myRUtil.GetMatListSexp(theHMM, fTransMat, myHMM.mTransMatVector) ;

	switch (myDistrType)
	{	case eNormalDistr :
		{
		cUnivariateNormal* myLoi = (cUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{
		cMultivariateNormal* myLoi = (cMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNbClasses, myLoi->mCov) ;
		}
		break ;
		case eMixtUniNormalDistr :
		{
		cMixtUnivariateNormal* myParam = (cMixtUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNbClasses, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNbClasses, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eMixtMultiNormalDistr :
		{
		cMixtMultivariateNormal* myParam = (cMixtMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNbClasses, myNbMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNbClasses, myNbMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{
		cDiscrete* myParam = (cDiscrete *)(myHMM.mDistrParam) ;
			myRUtil.GetEmissionSexp(myDistSEXP, 3, myParam->mProbaMatVector);
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

cDMatrix* myProbaCond = new cDMatrix[myNbSample] ;

	for (register uint n = 0 ; n < myNbSample ; n++)
		myProbaCond[n].ReAlloc(myNbClasses, myT[n]) ;

	myHMM.mDistrParam->ComputeCondProba(myY, myNbSample, myProbaCond) ;

cBaumWelch myBaumWelch=cBaumWelch(myNbSample, myT, myNbClasses) ;

	myBaumWelch.OutForwardBackward(myProbaCond, myHMM, myLogData) ;

	for (register uint n = 0 ; n < myNbSample ; n++)
	{       myProbaCond[n].Delete() ;
		myY[n].Delete() ;
	}

	delete [] myY ;

	delete [] myProbaCond ;

SEXP    myAux[7] ;



	myRUtil.SetListMatSexp(myBaumWelch.mAlpha, myNbSample,myAux[0]) ;
	myRUtil.SetListMatSexp(myBaumWelch.mBeta, myNbSample, myAux[1]) ;
	myRUtil.SetListMatSexp(myBaumWelch.mDelta, myNbSample, myAux[2]) ;
	myRUtil.SetListMatSexp(myBaumWelch.mGamma, myNbSample, myAux[3]) ;
	myRUtil.SetListListMatSexp(myBaumWelch.mXsi, myNbSample, myT, myAux[4]) ;
	myRUtil.SetListVectSexp(myBaumWelch.mRho, myNbSample, myAux[5]) ;
	myRUtil.SetListValSexp(myBaumWelch.mLogVrais, myAux[6]) ;

	delete [] myT ;
SEXP myRes ;
	PROTECT(myRes = allocVector(VECSXP, 7)) ;
	for (register int i = 0 ; i < 7 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	myRUtil.EndProtect() ;

UNPROTECT(1) ;
	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RLogforwardbackward(SEXP theHMM, SEXP theYt)
{
distrDefinitionEnum myDistrType ;
uint	myDimObs=1, 
		myNbClasses,
	myNbProba=0,
	myNbMixt=0 ;
cRUtil  myRUtil ;


SEXP myDistSEXP ;

	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba	
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNbClasses) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{	myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
			if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{       if (strcmp(myStr, "DISCRETE") == 0)
		{       myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNbProba) ;
		}
		else
		{       if (strcmp(myStr, "MIXTURE") == 0)
			{       myRUtil.GetValSexp(myDistSEXP, 2, myNbMixt) ;
				myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
			}
		}
	}
uint    myNbSample = length(theYt) ;    
uint*   myT = new uint[myNbSample] ;

cDVector* myY = new cDVector[myNbSample] ;


	for (register uint n = 0 ; n < myNbSample ; n++)
	{       SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm myHMM = cHmm(myDistrType, myNbClasses, myDimObs, myNbMixt, myNbProba) ;
	myRUtil.GetVectSexp(theHMM, fInitProba, myHMM.mInitProba) ;
	myRUtil.GetMatListSexp(theHMM, fTransMat, myHMM.mTransMatVector) ;

	switch (myDistrType)
	{       case eNormalDistr :
		{       cUnivariateNormal* myLoi = (cUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{       cMultivariateNormal* myLoi = (cMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNbClasses, myLoi->mCov) ;
		}
		break ;
		case  eMixtUniNormalDistr :
		{       cMixtUnivariateNormal* myParam = (cMixtUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNbClasses, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNbClasses, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eMixtMultiNormalDistr :
		{       cMixtMultivariateNormal* myParam = (cMixtMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNbClasses, myNbMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNbClasses, myNbMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{       cDiscrete* myParam = (cDiscrete *)(myHMM.mDistrParam) ;
				myRUtil.GetEmissionSexp(myDistSEXP, 3, myParam->mProbaMatVector);
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

cDMatrix* myProbaCond = new cDMatrix[myNbSample] ;

	for (register uint n = 0 ; n < myNbSample ; n++)
		myProbaCond[n].ReAlloc(myNbClasses, myT[n]) ;

		myHMM.mDistrParam->ComputeCondProba(myY, myNbSample, myProbaCond) ;

cLogBaumWelch myLogBaumWelch=cLogBaumWelch(myNbSample, myT, myNbClasses) ;
	myLogBaumWelch.LogForwardBackward(myProbaCond, myHMM) ;

	for (register uint n = 0 ; n < myNbSample ; n++)
	{       myProbaCond[n].Delete() ;
		myY[n].Delete() ;
	}

	delete [] myY ;

	delete [] myProbaCond ;

SEXP    myAux[6] ;
uint*   myLigne = new uint[myNbSample] ;


	for (register uint n = 0 ; n < myNbSample ; n++)
		myLigne[n] = myNbClasses ;

	myRUtil.SetListMatSexp(myLogBaumWelch.mLogAlpha, myNbSample,myAux[0]) ;
	myRUtil.SetListMatSexp(myLogBaumWelch.mLogBeta, myNbSample, myAux[1]) ;
	myRUtil.SetListMatSexp(myLogBaumWelch.mLogGamma, myNbSample, myAux[2]) ;
	myRUtil.SetListListMatSexp(myLogBaumWelch.mLogXsi, myNbSample, myT, myAux[3]) ;
	myRUtil.SetListVectSexp(myLogBaumWelch.mLogRho, myNbSample, myAux[4]) ;
	myRUtil.SetListValSexp(myLogBaumWelch.mLogVrais, myAux[5]) ;

	delete [] myLigne ;
	delete [] myT ;
SEXP myRes ;
	PROTECT(myRes = allocVector(VECSXP, 6)) ;
	for (register int i = 0 ; i < 6 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	myRUtil.EndProtect() ;

UNPROTECT(1) ;
	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RComputeCov(SEXP theHMM, SEXP theYt)
{

distrDefinitionEnum myDistrType ;
uint	myDimObs=1, 
		myNClass,
		myNProba=0,
		myNMixt=0 ;
cRUtil  myRUtil ;



SEXP myDistSEXP ;
	
	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba	
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNClass) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{	myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
			if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myStr, "DISCRETE") == 0)
		{       myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNProba) ;
		}
		else
		{	if (strcmp(myStr, "MIXTURE") == 0)
			{       myRUtil.GetValSexp(myDistSEXP, 2, myNMixt) ;
				myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
			}
		}
	}


uint    myNSample = length(theYt) ;    
uint*   myT = new uint[myNSample] ;

cDVector* myY = new cDVector[myNSample] ;


	for (register uint n = 0 ; n < myNSample ; n++)
	{       SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm myHmm = cHmm(myDistrType, myNClass, myDimObs, myNMixt, myNProba) ;
	myRUtil.GetVectSexp(theHMM, fInitProba, myHmm.mInitProba) ;
	myRUtil.GetMatListSexp(theHMM, fTransMat, myHmm.mTransMatVector) ;

		switch (myDistrType)
	{       case eNormalDistr :
		{       cUnivariateNormal* myLoi = (cUnivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{       cMultivariateNormal* myLoi = (cMultivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNClass, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNClass, myLoi->mCov) ;
		}
		break ;
		case  eMixtUniNormalDistr :
		{       cMixtUnivariateNormal* myParam = (cMixtUnivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNClass, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNClass, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNClass, myParam->mp) ;
		}
		break ;

		case eMixtMultiNormalDistr :
		{       cMixtMultivariateNormal* myParam = (cMixtMultivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNClass, myNMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNClass, myNMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNClass, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{       cDiscrete* myParam = (cDiscrete *)(myHmm.mDistrParam) ;
				myRUtil.GetEmissionSexp(myDistSEXP, 3, myParam->mProbaMatVector);
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

cInParam myInParam(myNSample, myDimObs, myY, myDistrType, myNClass, myNMixt, myNProba) ;
uint myNFreeParam = myHmm.GetNFreeParam() ;
cDerivative myDerivative(myInParam, myNFreeParam) ;

	myDerivative.ComputeDerivative(myHmm, myInParam) ;

cDMatrix myCov ;
	myDerivative.ComputeCov(myHmm, myCov) ;


	for (register uint n = 0 ; n < myNSample ; n++)
	{	myY[n].Delete() ;
		
	}

	delete [] myY ;

	delete [] myT;

SEXP myRes ;
	myRUtil.SetMatSexp(myCov, myRes) ;
	myRUtil.EndProtect() ;

	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RScoreAndInformation(SEXP theHMM, SEXP theYt)
{

distrDefinitionEnum myDistrType ;
uint	myDimObs=1, 
		myNClass,
		myNProba=0,
		myNMixt=0 ;
cRUtil  myRUtil ;



SEXP myDistSEXP ;
	
	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba	
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNClass) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{	myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
			if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myStr, "DISCRETE") == 0)
		{       myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNProba) ;
		}
		else
		{	if (strcmp(myStr, "MIXTURE") == 0)
			{       myRUtil.GetValSexp(myDistSEXP, 2, myNMixt) ;
				myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
			}
		}
	}


uint    myNSample = length(theYt) ;    
uint*   myT = new uint[myNSample] ;

cDVector* myY = new cDVector[myNSample] ;


	for (register uint n = 0 ; n < myNSample ; n++)
	{       SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm myHmm = cHmm(myDistrType, myNClass, myDimObs, myNMixt, myNProba) ;
	myRUtil.GetVectSexp(theHMM, fInitProba, myHmm.mInitProba) ;
	myRUtil.GetMatListSexp(theHMM, fTransMat, myHmm.mTransMatVector) ;

	switch (myDistrType)
	{       case eNormalDistr :
		{       cUnivariateNormal* myLoi = (cUnivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{       cMultivariateNormal* myLoi = (cMultivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNClass, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNClass, myLoi->mCov) ;
		}
		break ;
		case  eMixtUniNormalDistr :
		{       cMixtUnivariateNormal* myParam = (cMixtUnivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNClass, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNClass, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNClass, myParam->mp) ;
		}
		break ;

		case eMixtMultiNormalDistr :
		{       cMixtMultivariateNormal* myParam = (cMixtMultivariateNormal *)(myHmm.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNClass, myNMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNClass, myNMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNClass, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{       cDiscrete* myParam = (cDiscrete *)(myHmm.mDistrParam) ;
				myRUtil.GetEmissionSexp(myDistSEXP, 3, myParam->mProbaMatVector);
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

cInParam myInParam(myNSample, myDimObs, myY, myDistrType, myNClass, myNMixt, myNProba) ;
uint myNFreeParam = myHmm.GetNFreeParam() ;
cDerivative myDerivative(myInParam, myNFreeParam) ;

	myDerivative.ComputeDerivative(myHmm, myInParam) ;
cDVector myScore(myNFreeParam) ;
cDMatrix myInformation(myNFreeParam, myNFreeParam) ;
	myDerivative.ComputeScoreAndInformation(myScore, myInformation) ;

	for (register uint n = 0 ; n < myNSample ; n++)
	{	myY[n].Delete() ;
		
	}

	delete [] myY ;

	delete [] myT;


SEXP myAux[2] ;
	myRUtil.SetVectSexp(myScore, myAux[0]) ;
	myRUtil.SetMatSexp(myInformation, myAux[1]) ;
SEXP myRes ;
	PROTECT(myRes = allocVector(VECSXP, 2)) ;
	for (register int i = 0 ; i < 2 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	myRUtil.EndProtect() ;

UNPROTECT(1) ;
	return(myRes) ;}
END_EXTERN_C

#endif //_RDLL_
