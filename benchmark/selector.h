// Copyright Andy Singleton, 1993,1994
// This code is released for non-commercial use only
// For questions or upgrades contact:
// Andy Singleton, Creation Mechanics Inc.
// PO Box 248, Peterborough, NH 03458
// Internet: p00396@psilink.com
// Compuserve 73313,757
// Phone: (603) 563-7757

// W.Langdon@cs.ucl.ac.uk $Revision: 1.2 $
 
//Modifications (reverse order):
// WBL  8 Apr 1995  Add find
 
/////////////////////////////////////////////////////////////////////
// SELECTOR.H          Generic selector for GA use
//
//  SELECTOR OBJECT - selects individuals by scaled scores
//   also provides statistics
//	 also sorts
//  To use:
//			CSelector CSelectorOb(size,scaling_method, ratio_of_max_to_min)
//			CSelectorOb.reset()
//			CSelector.sort()
//			For each individual in selection window
//				CSelectorOb.add(fitness_score, individual_number)
//			Selected_individual = CSelectorOb.roulette()
/////////////////////////////////////////////////////////////////////

#ifndef _SELECTOR_H
#define _SELECTOR_H


#define FHANDLE int			// individual handle.  Could be int or pointer

class FitPair {
public:
	float rawscore;		// Fitness Score
	FHANDLE handle;				// Individual handle
};

#define IMAX 200			// maximum number of individuals to score

enum {EStraight,ERank,EMinToMax,EAveToMax,ETournament};
// EStraight - scaled fitness = raw score if >0, exp(raw score) if <0.  Weak on negatives
// ERank - Scaled fitness according to rank = 1+(maxratio-1)*rank/(icount-1)
// EMinToMax - Scaled fitness line from min to max = 1+((maxratio-1)/(max-min))*(raw_score-min)
// EAveToMax - Scaled fitness line from ave to max.  Cutoff negative values at .01
// ETournament - tournament fitness.  Return the best.

// return positive fitness from +- value
//#define POSVAL(v) (v>0? 1+v : exp(v))
#define POSVAL(v) (v>0? v : 0)

class CSelector {				// Class to score fitness and select individuals
	FitPair* pop;		// population so far
	float* scalescore;	// scaled scores
    int size;					// max number of pairs to insert
	int didscale;				// flag whether we scaled already
	int icount;					// count of individuals in array
	int method;					// scaling method, enumerated above
	int pmin;					// population minimum (number of individual)
	int pmax;					// population maximum (number of individual)
	float rawsum;				// population sum (raw)
	float scalesum;			// population sum (scaled)
	float maxratio;			// ratio of max to ave or min

public:
	CSelector(int m = EStraight,int s=IMAX);		// scaling method, size
	~CSelector() {delete[] pop;delete[] scalescore;};
	void reset() {icount = 0;pmin=0;pmax=0;rawsum=0;didscale=0;};
	void add(float score,FHANDLE h);				// add an individual
	FHANDLE roulette();							// roulette wheel selection on scaled scores
	float average() {return icount>0? rawsum/icount:0;};
	float vmin() {return pop[pmin].rawscore;}; // minimum value
	FHANDLE imin()  {return pop[pmin].handle;};		// minimum individual
	float vmax() {return pop[pmax].rawscore;};
	FHANDLE imax() {return pop[pmax].handle;};
	void sort();								// sort in score order
	FHANDLE GetHandle(int i) {return pop[i].handle;};
	int count() {return icount;};
	int find(FHANDLE h);
};

#endif		
/////////////////////////////////////////////////////////////////////
// EOF SELECTOR.H
/////////////////////////////////////////////////////////////////////
