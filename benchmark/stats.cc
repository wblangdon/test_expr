// stats.ccc  version to support queue2.cc, from GPQUICK-2.1 chrome.cc
// W.Langdon cs.ucl.ac.uk $Revision: 1.91x $

//Modifications (reverse order):
// WBL 31 Jan 2026  stub for bench_rand.make

#include <math.h>	//for DisplayStats
#include <assert.h>
#include "pch.h"
#include "chrome.h"
#include "prob.h"
#ifndef STACK_LIB
//#include "Primitives.h"
#ifdef QUEUE_LIB
#include "queue2.h"
#else
#include "gp.h"
#endif /*QUEUE_LIB*/
#endif

int Chrome::same(Chrome* other) const
{
return 0;

}//end same
