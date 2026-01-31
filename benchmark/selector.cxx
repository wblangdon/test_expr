// Copyright Andy Singleton, 1993,1994
// This code is released for non-commercial use only
// For questions or upgrades contact:
// Andy Singleton, Creation Mechanics Inc.
// PO Box 248, Peterborough, NH 03458
// Internet: p00396@psilink.com
// Compuserve 73313,757
// Phone: (603) 563-7757

// W.Langdon@cs.ucl.ac.uk $Revision: 1.3 $
 
//Modifications (reverse order):
// WBL  8 Apr 1995  Add find
 
// Generic SELECTOR object for GA use

#include "pch.h"
#pragma hdrstop
#include "selector.h"

// *********************** CSelector Object methods **********************


CSelector::CSelector(int m,int s)
{
	method = m;
	maxratio=100;
	size=s;
	pop= new FitPair[s];
	scalescore= new float[s];
	reset();
}

void CSelector::add(float score,FHANDLE h) {					// add an individual
	pop[icount].rawscore=score;
	pop[icount].handle=h;
	rawsum+=score;
	if (icount!=0) {
		pmin = (pop[pmin].rawscore < score? pmin : icount);
		pmax = (pop[pmax].rawscore > score? pmax : icount);
	}
	icount++;
}

int CSelector::find(FHANDLE h) {
	for(int i=0; i<icount; i++) {
		if(pop[i].handle==h) return i;
		}
	return -1;
}

void CSelector::sort()	 		// sort in score order, ascending
{
	int i,j,mn;
	float tempscore;
        FHANDLE temph;

	for (i=0;i<icount;i++)
	{	mn=i;
		for (j=i+1;j<icount;j++)
		{
			if (pop[j].rawscore < pop[mn].rawscore)
				mn=j;
		}
		tempscore=pop[i].rawscore;
		temph=pop[i].handle;
		pop[i].rawscore=pop[mn].rawscore;
		pop[i].handle=pop[mn].handle;
		pop[mn].rawscore=tempscore;
        pop[mn].handle=temph;
	}    	
}

FHANDLE CSelector::roulette() {	// scale and do selection on scaled scores
	int i,omin,winner;
	float slot,wedge,slope;
	if (!didscale) {
					// scale, according to method
			switch(method) {
			case EStraight:
				for (i=0;i<icount;i++) scalescore[i]=POSVAL(pop[i].rawscore);
				break;
			case ERank:
				sort();
				slope=(icount - 1 != 0? (maxratio-1)/(icount-1):0);
                // put in a score by rank
				for (i=0;i<icount;i++) scalescore[i]=1+slope*i;
				break;
			case EMinToMax:
				slope=(vmax()-vmin()!=0?(maxratio-1)/(vmax()-vmin()):0);
				for (i=0;i<icount;i++) scalescore[i]=1+slope*(pop[i].rawscore-vmin());
				break;
			case EAveToMax:
				slope=(vmax()-average()!=0? (maxratio-1)/(vmax()-average()): 0);
				for (i=0;i<icount;i++) {
					slot=1+slope*(pop[i].rawscore-average());
					scalescore[i]= slot>0? slot:.01;
				}
				break;
			}

		scalesum=0;
		for (i=0;i<icount;i++) scalesum+=scalescore[i];
		didscale=TRUE;
	}
	if (method!=ETournament)
        {
		slot=scalesum*rand_0to1();
		wedge=scalescore[0];
		winner=0;
		for (i=1;slot>=wedge && i<icount;i++) {
			winner++;
			wedge+=scalescore[i];
		}
	}
	else
	{					// ETOURNAMENT - return best of out of maxratio
		winner=rnd(icount);
		for (i=1;i<maxratio;i++)
		{
        	omin=rnd(icount);
			if (pop[omin].rawscore > pop[winner].rawscore)
				winner=omin;
        }
	}    	
	return pop[winner].handle;
}
