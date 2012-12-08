/* ************************************************************************ *
 * ************************************************************************ *

   File: obsSeq.h
   The class CObsSeq defines operations on observation sequences for HMM

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

class CObsSeq{

public:
  CObsSeq(void){};
//  CObsSeq(long nbSequences, long nbObs);
  CObsSeq(CObs *obsType, long nbSequences, long nbObs);
  CObsSeq(CObs *obsType, ifstream &inFile);
  ~CObsSeq(void);

public:
  void PrintHeader(ostream &sequenceFile);
  void Print(ostream &sequenceFile);
  void Make1DArray(void);
  inline long GetNbSequences(void){return mNbSequences;};

public:
  long mNbSequences;
  long *mNbObs;// array of nb of observations for each sequence
  CObs ***mObs;
  CObs **m1DArray;
  
public:
	CObs* mObsType;
	long mObsCount;
};

//===============================================================================
//===============================================================================
