/* ************************************************************************ *
 * ************************************************************************ *

   File: gaussianObsProb.h

   The class CGaussianObsProb defines operations 
   for the gaussian observation distribution B

  * ************************************************************************ *

   Authors: Daniel DeMenthon & Marc Vuilleumier
   Date:  2-18-99 

 * ************************************************************************ *

   Modification Log:
	

 * ************************************************************************ *
   Log for new ideas:
 * ************************************************************************ *
               Language and Media Processing
               Center for Automation Research
               University of Maryland
               College Park, MD  20742
 * ************************************************************************ *
 * ************************************************************************ */
 
//===============================================================================

class CGaussianObsProb: public CObsProb{
  public:
	CGaussianObsProb(int nbSymbols, int nbStates);
        CGaussianObsProb(int nbSymbols, int nbStates, double *means, double* stds);
        CGaussianObsProb(ifstream &hmmFile, int nbStates);
        ~CGaussianObsProb(void);
	void Start(void);
	void StartIter(void);
	void BWSum(double *gamma, CObs *obs);
        void SKMSum(int state, CObs *obs);
	double EndIter();
	void End();
	CObs* PickObservation(int state);
	void Print(ostream &outFile);
        CObs** MapStateToObs(void);
	double at(int state, CObs *obs);
  //	int GetM(void){return M};
#if 1
	CObs* ReadObsFrom(ifstream &inFile);
	void ReadFileHeader(ifstream &inFile);
	void PrintFileHeader(ostream &outFile);
#endif

  private:
	void InitParameters(void);
        void InitParameters(double *means, double* stds);
  private:
       double *mMean;
       double *mStd;
       double *mSumMean, *mSumVar, *mDiv;
};

//===============================================================================
//===============================================================================

//===============================================================================
