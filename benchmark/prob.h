#ifndef POLY_LIB
#define POLY_LIB true
// avx.h   Include file for Intel AVX 512 float symbolic regression problems based on bool.h r1.14
// W.B.Langdon cs.bham.ac.uk 3 March 1998

// version "$Revision: 1.11 $"

//Modifications (in reverse order)
// WBL  4 Jan 2019  Add CACHE_LINE_SIZE for PTHREADS
// WBL 15 Dec 2018  Add PTHREADS
// WBL 22 Nov 2018  based on bool.h r1.14
// WBL  5 Aug 2016  Replace fit64 with raw value
// WBL  6 Jun 2016  try adding CROSSFILE, INTRON64
// WBL 24 May 2016  for gcc 4.8.5
// WBL 25 Jan 1999  add ALLFUNC_LIB (order 4 max)
// WBL  9 Jan 1999  use SMC_GP
// WBL 11 Feb 1998  use non-macro definition of setmaxfit
// WBL 11 Feb 1998  Support requiring food to be eaten in order


//#define INTRON64 1
//#define CROSSFILE 1

//#define ALLFUNC_LIB true
typedef float retval;
typedef void reteval; //Values held on EvalSP evalstack retval[MAXTESTCASES] 
#define GENERATIONAL true
//PTHREADS depends on GENERATIONAL
#define PTHREADS true
#define npthreads 96
#define CACHE_LINE_SIZE 64

#define NUM_TREES 1

//#define TIMINGS true
//RAND_TREE assumes MULTREE inuse
#define RAND_TREE true
#define RAND_TREE_ARITIES true

class ChromeParams;
extern float setmaxfit(const ChromeParams* p); //create evalstack, keep going for ever


extern int Reset();

//#define MULTREE
#define MAXTREEDEPTH 0 /*not in use*/

//MAXTESTCASES must be aligned
#define MAXTESTCASES 48
//requires FASTEVAL
#define FASTCONST
//Override POLY_LIB in chrome.h
//NOT WANTED #define FASTTRAVERSE

#define NUM_OPERATIONS 1
#define NUM_OTHERS 0
#define num_pareto_components  1
#define num_ellite_pareto      1
typedef float scoretype;

//the following not used.. but need to keep compier happy

//#define MAX_BREAKS         1
#define NUM_TEST_SEQUENCES 1
#define NUM_TEST_PHASES    1
/*

#define last_op  0

retval random_value(int& seed); //seed > 0
*/

extern      float max_fitness;

extern      int end_gens(); //see if need to move to the next test phase
                             //return <> 0 if wish to save a bench point

extern      int num_sol_found; //for save and restore to dumpfile

extern      int problem_order;
#ifdef INTRON64
extern      unsigned char intron64 [100000000];
//extern    unsigned char fit64 [1000000];
extern      retval raw[100000000];
extern      INT32 chrome64_age;
//#define CLEAR_INTRON assert(EXPRLEN==1000000); memset(intron64,0,EXPRLEN*sizeof(unsigned char)); memset(fit64,0,EXPRLEN*sizeof(unsigned char)); chrome64_age = gencount
  #define CLEAR_INTRON assert(EXPRLEN==100000000); memset(intron64,0,EXPRLEN*sizeof(unsigned char)); memset(raw,  0,EXPRLEN*sizeof(retval));        chrome64_age = gencount
#endif /*INTRON64*/

#endif /*POLY_LIB*/
