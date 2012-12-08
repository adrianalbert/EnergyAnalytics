/* ************************************************************************ *
 * ************************************************************************ *

   File: hmm.h
   The class CHMM defines operations for HMM

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

//===============================================================================
//===============================================================================

//===============================================================================

const double DELTA = 5*1E-3;// threshold used for exiting Baum-Welch or K-Means loop
//===============================================================================

class CHMM{

public:

	CHMM(CStateTrans *a, CObsProb *b, CInitStateProb *pi);
//      CHMM(char* hmmFileName);
	~CHMM(void);
	double Forward(double **alpha, double *scale, CObs **obs, long T, boolean doLog);
	void Backward(double **beta, double *scale, CObs **obs, long T);
        double Viterbi(CObs **obs, long T, int *q);
        double ViterbiLog(CObs **obs, long T, int *q);
	double BaumWelchCore(CObs **obs, long T, double *gamma, double **xi, boolean doLog);
        double IterBaumWelch(CObsSeq *obsSeq, double *gamma, double **xi);
	void LearnBaumWelch(CObsSeq *obsSeq);
        double SegmentalKMeansCore(CObs **obs, long T);
        double IterSegmentalKMeans(CObsSeq *obsSeq);
        void LearnSegmentalKMeans(CObsSeq *obsSeq);
        void LearnHybridSKM_BW(CObsSeq *obsSeq);
        CObsSeq* GenerateSequences(long nbSequences, long nbObs, int seed);
        CObs** GenerateObservations(long T);
 #if 1
        CObsSeq* ReadSequences(ifstream &inputSeqFile);
 #endif
	void PrintStatesAndExpectedObs(CObsSeq *obsSeq, 
						 ostream& stateFile, ostream& bestObsFile);
        double FindDistance(CObsSeq *obsSeq, ostream &outFile);
        double FindViterbiDistance(CObsSeq *obsSeq, ostream &outFile);
        double FindCrossEntropyDistance(CObsSeq *obsSeq, ostream &outFile);
        double FindQToPiProb(long T, int *q);
        void Print(ostream &outFile);
        int GetN(void){return mN;};

protected:
        int mN;// nb of states
	CStateTrans *mA;
	CObsProb *mB;
	CInitStateProb *mPi;
};

//===============================================================================
//===============================================================================
