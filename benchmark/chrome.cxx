// chrome.cxx  version to support stack.cc, from GPQUICK-2.1
// W.Langdon cs.ucl.ac.uk $Revision: 1.243 $

#define debug

//Modifications (reverse order):

// WBL 12 Jan 2019  Since gp.cc does not use DisplaySolStats remove oldpop
// WBL  5 Jan 2019  avx.cc nolonger using SetupEval, remove expr_table
// WBL  1 Jan 2019  Init treeinfo, matesizeavail only if used, 
//                  add NUM_TREES+1 index to tree_starts and use to replace many traverse
// WBL 30 Dec 2018  Tidyup PTHREADS debug code
// WBL 15 Dec 2018  Add PTHREADS
// WBL 23 Nov 2018  for avx split retval into reteval and retval
// WBL 26 Aug 2016  experiment with reporting size and maxdepth
// WBL 24 Jun 2016  To count children with perfect score add oldpop

// WBL 28 Oct 1999  Add displaydepthlimit (cf r1.155) as option on Write
//                  Nb r1.155 may be buggy

// WBL 29 Sep 1999  Add JoinTree

// WBL 22 Sep 1999  Reimplement MutateC, MutateShrink
//                  Turn off MutateNode on terminals if MutateC is active.

// WBL 21 Sep 1999  Significant speed in MutateNode by using RndOp

// WBL 20 Aug 1999  Improve debug in chrome::Load, remove displaydepthlimit

// WBL 17 Aug 1999  For Quintic poly re-runs
//                  Add displaydepthlimit and remove debug
//                  Add negative values for pCrossFairLimit
//                  Add (Boolean) pInitSparse

// WBL 29 Jun 1999  Add Homologous diference monitoring code, assumes simple
//                  fitness function, like BOOL_LIB
//                  NOT SUITABEL FOR GENERAL USE

// WBL 30 Apr 1999  Conversion for SiliconGraphics MIPS4, avoid EXPRLEN arrays

// WBL 27 Apr 1999  Allow chromes to exceed 32767

// WBL 13 Jan 1999  Add pUniInitWt, pMinInitSizeExpr and pInitSizeExpr

// WBL  8 Dec 1998  Change algorthim for CrossFair to ensure expected size 
//                  change on all crossovers is zero and increase range

// WBL 10 Aug 1998  Make iscomment dump comment if #debug
//                  Add display of depth to write_trace

// WBL 31 Jul 1998  Implement constants 

// WBL  7 Jul 1998  Implement FASTTRAVERSE, add arityop

// WBL 26 Mar 1997  BUGFIX 1ptXO (include roots of non-matching subtrees as
//                  potential crossover points, cf discussion with rmp)

// WBL 24 Mar 1997  Make adjustable point mutation rate ensure prog is changed

// WBL 24 Mar 1997  Add pMuteCrossWt and adjustable point mutation rate
//                  replace newchrome by calling same

// WBL 22 Mar 1997  Add 1-point crossover plus its variants
//                  MULTREE required by these changes

// WBL 22 Mar 1997  Ensure calls clear_rank at start of each GENERATION

// WBL 30 Dec 1997  Make Pop::Fitness display rval if pTraceTest is set

// WBL  1 Sep 1997  Add CrossSubTree_t_chnglen and CrossSubTree_calls

// WBL 31 Aug 1997  Add pCrossFairWt and pCrossFairLimit

// WBL 18 Aug 1997  Fix MutateSubTree crossfile output for pMuteFairLimit<>0

// WBL 18 Aug 1997  BUGFIX pMuteFairLimit<0 code to avoid parsimony bias

// WBL 16 Aug 1997  Simplify code of MutateSubTree by completly
//                  separating RAND_TREE conditional code from the rest and
//                  removing other conditionals within it 
//                  Make negative pMuteFairLimit mean 1st rand_tree pMuteFair
//                  Remove rand_tree size exceeded displays

// WBL 12 Aug 1997  To eliminate parsimony bias, tweak new version (RAND_TREE)
//                  of pMuteFair so reselects both mutation points if either
//                  indicates a tree that its too big for rand_tree to
//                  generates 
//                  add pMuteFairLimit, if==0 then old version of pMuteFair

// WBL  9 Aug 1997  Add support for RAND_TREE

// WBL  1 Aug 1997  Add pMaxHaltExpr and Temperature

// WBL 30 Jul 1997  make pStartTemp and pEndTemp doubles

// WBL 28 Jul 1997  Allow pMuteFair to grow single terminal trees

// WBL 28 Jul 1997  BUGFIX ensure Pop::Pop sets initeval ok if popsize=1
//                  Dont try and delete newpop unless it has been created

// WBL 25 Jul 1997  Ensure doanneal and docrossanneal are working, 
//                  rewrite GetAnnealTarget to use selectParant, add accept

// WBL 24 Jun 1997  BUGFIX ensure go_until does not stop early, this ensures
//                  newpop is copied to pop

// WBL 23 Jun 1997  Disable crossfile output

// WBL  8 Jun 1997  Restrict pMuteFair to doubling tree size etc.

// WBL 10 May 1997  Add crossfile output to MutateSubTree

// WBL 10 May 1997  Add pMuteFair

// WBL 30 Apr 1997  Add first level of PDGP support

// WBL  2 Apr 1997  Add ORIGINAL_DEME to allow queue2 to operate identically

// WBL 14 Mar 1997  BUGFIX Ensure only one individual per crossfile line 
//                  Ie copy clones to appear on own line.

// WBL 15 Feb 1997  BUGFIX (GENERATIONAL) ensure answer returned by generate 
//                  refers to chrome just created. Was going wrong if last
//                  chrome in the population was a copy. Improve layout of
//                  generate while about it.

// WBL 13 Feb 1997  Add pTournOff

// WBL 31 Dec 1996  Reduce number of warnings with newer C++ compilers,
//                  Make deleting newpop cleaner, 
//                  ensure lastchrome is initialised

// WBL 14 Dec 1996  Add pMaxDepth

// WBL 23 May 1996  Add pCPU

// WBL 16 May 1996  BUGFIX ensure MutateSubTree dont crash with depth==0

// WBL 26 Apr 1996  Add display of location of parents and offspring
//                  to crossfile. Trace deme bug
//                  BUGFIX correct calculation of offset in select_candidates
//                  further testing required

// WBL 12 Apr 1996  BUGFIX move update_rank to InsertMember so mutation and
//                  seeding work with pElitist. Nb remove other calls to it

// WBL 11 Apr 1996  Allow multiple seeds in gp.load
//                  Allow comments in gp.load (based on isscoment from test.cc)

// WBL  8 Apr 1996  Get mutation working

// WBL 21 Mar 1996  Add chrome differs from mum or dad to crossover log

// WBL 28 Feb 1996  Add logfile of crossover details

// WBL 27 Feb 1996  Disable static_check when used in stack problem
//                  Restore oldversion of GetTarget for STACK_LIB

// WBL 12 Dec 1995  Make Chrome::Chrome more robust to errors in the dumpfile

// WBL 10 Dec 1995  Ensure Pop::Pop calls clear_rank

// WBL 13 Jun 1995  Add binary save mode to Pop::Write 
//                  and input stream to Pop::Pop and Chrome::Chrome

// WBL 27 Jun 1995  Add reinitialise and reinit_trees

// WBL 23 May 1995  Add reporting which tree is used for crossover

// WBL  6 May 1995  Remove MaxExpr, ExprBytes, MaxTreeExpr, funccount,
//                  from Chrome. Replace with params pMaxExpr, pMaxTreeExpr
//                  and probl->funccount. Replace depth with new SubWrite arg

// WBL  4 May 1995  Restore copying fitnessvalue on crossover (see r1.29,1.28)
//                  in order to support STORE_FIT. Add call to update_rank
//                  Add pStoreFit

// WBL  2 May 1995  Add pDirCrossWt

// WBL 26 Apr 1995  Reduce size of chrome by removing idx from node
//                  since it is never used. Mainly done by redefining
//                  macros in chrome.h but minor change to load

// WBL 18 Apr 1995  Improve layout of debug info

// WBL 17 Apr 1995  Add optional last_tree argument to chrome::write
//                  Add and use pMaxTreeExpr

// WBL  8 Apr 1995  Add check function exists in tree. This also
//                  ensures we get right function where they are ambigious

// WBL  7 Apr 1995  Make gp.load case sensitive by chaning FindFunc

// WBL 23 Mar 1995  Move pElitist from Worstof to select_candidates

// WBL 23 Mar 1995  Add Chrome* optional parameter to GetFitnessObj

// WBL 22 Mar 1995  Dont call Copy when creating a new Chrome by crossover

// WBL 15 Mar 1995  Add pElitist

// WBL 12 Mar 1995  Add pTournComp and pNicheSize

// WBL 19 Feb 1995  Add passing target to Bestof and Worstof

// WBL 19 Feb 1995  Add call to static_check in chrome::Load

// WBL 28 Jan 1995  Add some debug traces to chrome::Load() and SubLoad()

// WBL 29 Nov 1994  Implement pnames, change Load and Saveto use it.
//                  add Load(char*)

// WBL 15 Oct 1994  Make chrome::Load() work with MULTREE 
//                  NB I have not implemented constants

// WBL 14 Oct 1994  Call static_check when creating entirely new chrome

// WBL 13 Sep 1994  Make compilable by solaris gcc

// WBL  9 Sep 1994  UNDO part of 17 Aug 1994: ``InsertMember nolonger
//                  places worse child into pop'' 

// WBL  8 Sep 1994  Make crosstree call static_check

// WBL 25 Aug 1994  Make funcbag and array, add funcbagcount.
//                  Replace getfuncbag by addtofuncbag

// WBL 24 Aug 1994  Add Problem::Bestof and Worstof and use it in place
//                  of IsBetter tournaments. Add Pop::select_candidates

// WBL 23 Aug 1994  Remove copying hits and write_hits() (to problem
//                  specific file). Replace call of write_hits by
//                  write(nfitness*). Use this inplace of fout<<fvalue.
//                  Make Pop::Pop call problem specific NewChrome.
//                  Add debug output to stub Problem::fitness

// WBL 17 Aug 1994  TRY AND COMMENT OUT GetTarget() avoid selecting BestMember
//                  InsertMember nolonger places worse child into pop (nb
//                  appears to be very CPU inensive, illusion of prog length?

// WBL 15 Aug 1994  Add write_hits(). 
//                  Make  Chrome::Dup() set mum and write_trace() trap xtree<0

// WBL  6 Jul 1994  Add two dimensional torridal demes. Add pPopWidth,
//		    pDemeWidth

// WBL 1  Jul 1994  Remove class stats, Pop::DisplayStats and
//                  Pop::DisplayDStats to new file stats.cc

// WBL 28 Jun 1994  Make nfit take average of cube rather than linear average
//                  Tidy displays in DisplayDStats

// WBL 27 Jun 1994  Correct display of gencount in DisplayStats
//                  Add class stats and make DisplayStats use it.
//                  Complete impementation of DisplayDStats

// WBL 24 Jun 1994  Direct crossover to a particular tree
//                  Protect sqrt from rounding errors

// WBL 14 Jun 1994  Add support for pPopSeed, pTestSeed and PTrace

// WBL 14 Jun 1994  Add ChromeParams::Save() and ChromeParams::Load()
//                  based upon gpcplus3.0 cpv.cc.
//                  Add Pop::DisplayDStats()
//		    Add Pop::Write()

// WBL 16 May 1994  Add Pop::DisplayStats(). Fix MULTREE CrossOver bug

// WBL 12 May 1994  fix % by zero error in rnd() on SETNODE if varnum==0,
//		    Probably only arises if you have the wrong compiler switch

// WBL  5 May 1994  Make pMaxAge = 0 disable age limit,
//		    Add support for P-M random numbers

// GPQUICK
// C++ prefix stack implementation
// Divides GP system into classes:
//    Chrome-Function  - representation of the genetic program
//    Pop - Runs the GA
//    Problem - fitness and primitive functions

//#include <math.h>	//for DisplayStats
#include <assert.h>
#include "pch.h"
#include "chrome.h"
#ifdef PTHREADS
#include <unistd.h>
#include <pthread.h>
#endif

//#define CROSSFILE true
#ifdef CROSSFILE
#define crossfile cout
#endif

// Global variables
char scratch[100];

int gpquick_seed = 1;

#ifdef FASTEVAL
        // Global evaluation variables.  Eliminates parameter passing
        // If you can get IpGlobal into a register, you will run near the speed
        // limit for S-expression evaluation
        Chrome * ChromeGlobal;                  // currently evaluating Chrome
        evalnode * ExprGlobal = NULL; //[EXPRLEN];   // expression expanded to evalnodes
        evalnode * IpGlobal;            // Current instruction pointer
#ifdef FASTEVAL
        evalnode* fast_tree_starts [NUM_TREES];     // start of each tree 
#endif
#endif


//**************************************************************************
//////////////////////////////////// Utility functions
int min(int value1, int value2)
   {
          return ( (value1 < value2) ? value1 : value2);
   }
int max(int value1, int value2)
   {
          return ( (value1 > value2) ? value1 : value2);
   }

// Comment out if included in your implementation (e.g. Borland)
        // upper case string
char* strupr(char* s)
{
        int i=0;
        while (s[i])
                if (s[i]>='a' && s[i]<='z') s[i]=s[i]-32;
        i++;
        return s;
}

        // compare strings without case sensitivity
int strcmpi(char* s1, char* s2)         // only checks equality
{
        int mismatch=0;
        while (*s1 && *s2 && !mismatch)
        {
                if (*s1>='a' && *s1<='z' && *s2<='Z')
                        mismatch = *s1-32 != *s2;               
                else if (*s1>='A' && *s1<='Z' && *s2>='a')
                mismatch = *s1+32 !=*s2;

                else
                        mismatch = *s1!=*s2;
                s1++;s2++;
        }
        return mismatch || *s1 || *s2;
}


//**************************************************************************
//////////////Function member function (separated for use in a DLL)
char* Function::getprint(Chrome* chrome)
{
        // NUMBER has its own getprint to write the correct number
        // This is fancy footwork for other functions with operands
        // these are written as <name>_<operand number>
        if (varnum) sprintf(scratch,"%s_%-d",name,VARNUM(chrome->expr[chrome->ip]));
        return (varnum? scratch : name);
}
#ifndef INTRON64
//WBL 7 Sep 2016 ignore new parameters for bool.cc intron64
char* Function::getprint(Chrome* chrome,const int dummy_ip,char* dummy_scratch)
  {return getprint(chrome);}
#endif /*INTRON64*/



//**************************************************************************
/// Chrome functions ///////////////////////////////////////////

#ifdef PTHREADS
int DoCross = 0; //ugly hack to see if in initial generation or not
#endif
Chrome::Chrome(ChromeParams* p,Problem* pr,Function** f,retval* c,BOOL doinit,
			     istream* fin )
#ifdef TRACE
                  :num_children (0)
#endif
#ifdef PTHREADS
		  ,expr(NULL)
#endif
{
// initialize a expr using the functions in array f

//cout << "Chrome::Chrome\n";
        int depth = p->params[pInitExpr];
        const int md = p->params[pMinInitExpr]-1;
//        funcbag=cf;
        params = p;
        probl = pr;
        funclist = f;
        constlist = c;
#ifdef PTHREADS
        if(DoCross == 0 || params->params[pThreads]==0) //delay till thread
#endif
        expr=new node[params->params[pMaxExpr]];
        ip=0;
        birth = 0;
        nfitness  = probl->GetFitnessObj(this);
        // allocates a FitnessValue object;  GetFitnessObj  calls new

//cout<<"Size of chrome "<<sizeof(*this)<<" + "
//    <<sizeof(node)*params->params[pMaxExpr]<<endl;

        if (doinit)
        {
        if (fin!=NULL && *fin)
        {//load each line into a separate tree
	 //psuedo binary. 1 byte per opcode but in ascii
#ifdef PDGP
	  if(params->params[pPDGP]) assert(0==1);//not implemented yet
#endif /*PDGP*/
		int ip = 0;
#ifdef MULTREE
		tree_starts[NUM_TREES] = -1;
		for (int t = 0; t < NUM_TREES; t++ ) {
			tree_starts[t] = ip;
#endif
			char c;
			for(c=fin->get(); (c != '\n') && (c != EOF); ip++) {
//cout<<"`"<<c<<"'";debug
				assert(ip<params->params[pMaxExpr]);
				expr[ip].op = c-'A';
				c = fin->get();
			};//endfor each line
		if(c == EOF) {                                //sanity checks
			cout<<"Unexpected EOF on loading dumpfile\n";
			assert(c != EOF); }
		if( (ip-tree_starts[t]) <= 0) {
			cout<<"Empty tree in dumpfile\n";
			assert( (ip-tree_starts[t]) > 0 );}
		if( (ip-tree_starts[t]) != SubLen(tree_starts[t]) ) {
			cout<<"Corrupt tree in dumpfile\n";
			assert( (ip-tree_starts[t])==SubLen(tree_starts[t]));}
#ifdef MULTREE
//cout<<"end tree "<<t<<endl; debug
		int valid = probl->static_check(this,t);
		if (valid>0) cout<<"static_check reports "<<valid
				 <<" on tree "<<t<<endl;
		};//endfor all trees
		tree_starts[NUM_TREES] = ip;
#endif
	}
	else {
#ifdef PDGP
	  if(params->params[pPDGP]) initpdgp();
	  else
#endif /*PDGP*/
#ifdef MULTREE
	tree_starts[NUM_TREES] = -1;
	for (int t = 0; t < NUM_TREES; t++ ){
	        tree_starts[t] = ip;
#ifndef STACK_LIB
	const BOOL UniInit = (params->params[pUniInitWt] > 0 && 
			      rnd(100) < params->params[pUniInitWt]);
	do {
#endif
	        ip = tree_starts[t]; //discard earlier attempt at tree
#endif
		if(UniInit) {
                // DO uniform random of size between pMinInitSizeExpr and pInitSizeExpr
		  int limit=0;
		  do {//loop because there may be no tree of the chosen length
		    assert(limit++<1000); //beter than infinite loop?
		  } while(!rand_tree(params->params[pMinInitSizeExpr] +
				     rnd(1 +
					 params->params[pInitSizeExpr] -
					 params->params[pMinInitSizeExpr]),
					 t));
		}
		else {
                // DO "ramped half and half from pMinInitExpr to depth"
                if (depth > 0)
                {
#ifdef MULTREE
                  SubInit(probl->funcbag[t],1,rnd(depth-md)+md,rnd(2), md, t);
#else
                  SubInit(probl->funcbag[0],1,rnd(depth-md)+md,rnd(2), md);       // go to some depth less than max, full or half full
#endif
                }
                else
                        SETNODE(expr[ip],0,0);    // stick in a stub constant
		}//endif UniInit
#ifdef MULTREE
#ifndef STACK_LIB
	}while (probl->static_check(this,t)>rnd(100));
#endif
	}//end for each tree
	tree_starts[NUM_TREES] = ip; //ip set by SubInit
#endif
        }//end else fil
	}//end doinit
}

#ifndef PDGP
#ifdef MULTREE
void Chrome::reinit_trees(BOOL flag[NUM_TREES])
{
node* save = new node[params->params[pMaxExpr]];
memcpy(save,expr,params->params[pMaxExpr]*sizeof(node));

int save_starts[NUM_TREES+1];
save_starts[NUM_TREES] = tree_starts[NUM_TREES];
{for (int t = 0; t < NUM_TREES; t++) {
	save_starts[t] = tree_starts[t]; }}

        const int depth = params->params[pInitExpr];
        const int md = params->params[pMinInitExpr]-1;

        ip=0;

	tree_starts[NUM_TREES] = -1;
	{for (int t = 0; t < NUM_TREES; t++ ){
		tree_starts[t] = ip;
//debug	if(tree_starts[t] != save_starts[t])
//		cout<<"New tree_starts["<<t<<"]="<<tree_starts[t]
//end debug	    <<" original="<<save_starts[t]<<endl;
	if(flag[t]) {
	const BOOL UniInit = (params->params[pUniInitWt] > 0 && 
			      rnd(100) < params->params[pUniInitWt]);
	do {
	        ip = tree_starts[t]; //discard earlier attempt at tree

		if(UniInit) {
		// DO uniform random of size between pMinInitSizeExpr and pInitSizeExpr
		int limit=0;
		do {//loop because there may be no tree of the chosen length
		  assert(limit++<1000); //beter than infinite loop?
		} while(!rand_tree(params->params[pMinInitSizeExpr] +
				   rnd(1 +
				       params->params[pInitSizeExpr] -
				       params->params[pMinInitSizeExpr]),
				   t));
		}
		else {
                // DO "ramped half and half from pMinInitExpr to depth"
                if (depth > 0)
                {
                  SubInit(probl->funcbag[t],1,rnd(depth-md)+md,rnd(2), md, t);
                }
                else
                        SETNODE(expr[ip],0,0);    // stick in a stub constant
		}
	}while (probl->static_check(this,t)>rnd(100));
	if(t == NUM_TREES-1) tree_starts[NUM_TREES] = (depth<=0)? 1 : ip; //ip set by SubInit
	}//endif flag
	else {
		ip = tree_starts[t] + save_starts[t+1]-save_starts[t];
		if(ip>=params->params[pMaxExpr]) {
			cout<<"reinit_trees OUTOFSPACE on tree "<<t<<endl;
			assert(0==1);
		}
		memcpy(&(expr[tree_starts[t]]),&(save[save_starts[t]]),
		       save_starts[t+1]-save_starts[t]                 );
		if(t == NUM_TREES-1) tree_starts[NUM_TREES] = ip;
	}//end else flag
	}}//end for each tree

delete[] save;

}//end reinit_trees
#endif
#endif /*PDGP*/

Chrome* Chrome::Copy()
{
        Chrome* nw=new Chrome(params,probl,funclist,constlist,FALSE);
        nw->Dup(this);
        return nw;
}

void Chrome::Dup(Chrome* source)                // make a copy nondestructively
{
#ifdef TRACE
	source->num_children++;
	mum.birth = source->birth;
#endif
#ifdef MULTREE
	memcpy (tree_starts, source->tree_starts, sizeof(tree_starts));
#endif
//        funcbag=source->funcbag;
        params=source->params;
        funclist=source->funclist;
        constlist=source->constlist;
        
#ifdef PTHREADS
        if(DoCross == 0 || params->params[pThreads]==0) //delay till thread
#endif
        memcpy(expr,source->expr,params->params[pMaxExpr]*sizeof(node));
        nfitness->Copy(source->nfitness);
#ifdef ORIGINAL_RANK
        nfitness->update_rank();
#endif
        ip=0;
}

#ifdef PTHREADS
reteval Chrome::evalAll(int tree) {assert(strlen("PTHREADS Chrome::evalAll not implemented ")==0);}
#else
#ifdef MULTREE
reteval Chrome::evalAll(int tree) // eval the whole expression anew
#else
reteval Chrome::evalAll() // eval the whole expression anew
#endif
{
#ifndef FASTEVAL
#ifdef MULTREE
        ip=tree_starts[tree]-1;
#else
        ip=-1;
#endif
        return eval();
#else
        
#ifdef MULTREE
	IP=fast_tree_starts[tree];
#else
        IP=ExprGlobal;
#endif
        //cout<<"FUNC="<< IP->op;
        return (IP->ef)();
#endif
}
#endif /*PTHREADS*/

#if 0
//debug code. inefficient but fewer side effects
void Chrome::SetupEval(const int tree /*= -1*/,evalnode* xip /*NULL*/) {
  assert(tree_starts[1]>0);
  const int length = tree_starts[1]; //SubLen(0);
  if(xip==NULL && ExprGlobal==NULL) {
    ExprGlobal = new evalnode[params->params[pMaxExpr]];
  }
  evalnode* ip = (xip!=NULL)? xip : ExprGlobal;
  int args=1; //sanity check only
  for(int i=0;i<length;i++,ip++){
    const node* spp = &expr[i];
    SETEVALNODE(ip,spp->op,constlist[spp->op]);
    args += PARGNUM(spp)-1;
  }
  assert(args == 0);
}
#else /*end debug code*/
void Chrome::SetupEval(const int tree /*= -1*/,evalnode* xip /*NULL*/)
{
                // expand stack into evalnodes in global eval stack
#ifdef FASTEVAL
        if(xip==NULL && ExprGlobal==NULL) {
	  ExprGlobal = new evalnode[params->params[pMaxExpr]];
	}
        int args;
        evalnode* ip = (xip!=NULL)? xip : ExprGlobal;
#ifdef MULTREE
        node* spp;
	int start; int end;
	if(tree==-1) {	spp=expr;
			start = 0;    end = NUM_TREES;}
	else         {	spp=expr+tree_starts[tree];
			start = tree; end = tree+1;}
	for (int t = start; t < end; t++ ){
		fast_tree_starts[t] = ip;
#else
        node* spp=expr;?
#endif
#ifdef PDGP
	if(params->params[pPDGP]) {
	  int d = params->net_limits[tree];
	  ip = SetupEvalpdgp(ip,d,pdgp(d,0));
	}
	else {
#endif
#ifdef FASTTRAVERSE
	typedef struct {int args; evalnode* ip;} s;
	s* stack = new s[params->params[pMaxExpr]]; //[EXPRLEN];
	s* sp = stack;
//	memset(ExprGlobal,0,sizeof(ExprGlobal));
//	memset(stack,0,sizeof(stack));
#endif /*FASTTRAVERSE*/
        args=1;
        while (args > 0)
        {
                SETEVALNODE(ip,spp->op,constlist[spp->op]);
#ifdef FASTTRAVERSE
		if(PARGNUM(spp)>0) { //save info so can recognise end of args
		  //push(ip,args);
		  ++sp;
//		  assert(sp < &stack[EXPRLEN]);
		  sp->args = args;
		  sp->ip   = ip;
		}
		else {
		  SETEVALJUMP(ip,ip+1);    //terminal so jump 1
		  while(sp>stack && sp->args==args) {  //fill in closed forward jumps
		    SETEVALJUMP(sp->ip,ip+1);  //jump past closing terminal
		    sp--;
		  }
		}
#endif /*FASTTRAVERSE*/
                args=args+PARGNUM(spp)-1;
                ip++;
                spp++;
        }
#ifdef FASTTRAVERSE
	assert(sp==stack);
	delete[] stack;
#endif /*FASTTRAVERSE*/
	if(tree!=-1) {
#ifdef FASTTRAVERSE
		ip->ef   = 0; //clear next opcode Beware array bounds
		ip->jump = 0; //and set other pointers ??????TRAP for now
#else
		ip->op = 0; //clear next opcode Beware array bounds
                            //and set other pointers
#endif /*FASTTRAVERSE*/
		for(t = 0; t<NUM_TREES; t++) {
			if(t!=tree) fast_tree_starts[t] = ip;} 
	}
                // set global eval pointers
#ifdef PDGP
       }//endif pPDGP
#endif /*PDGP*/
#ifdef MULTREE
	}//end for each tree
#endif
        ChromeGlobal=this;
#endif
}
#endif /*debug code*/

#ifdef MULTREE
void Chrome::SubInit(CSelector* funcbag, int argstogo, int depth,
	BOOL isfull, int mindepth, int tree)
#else
void Chrome::SubInit(CSelector* funcbag, int argstogo,int depth,
	BOOL isfull, int mindepth)     // Initialize a subtree half full
#endif
{
//cout <<"SubInit( " << argstogo << ", " << depth << ", " << isfull <<
//", " << isstart << " )\n";

        int i;
        int maxargs,thisargs;
//calculate max number of args there is space to add. 
        maxargs = params->params[pMaxExpr] -(ip+argstogo);
	assert(maxargs>=0); //can sometimes trip this WBL    
#ifdef MULTREE
	{int m = params->params[pMaxTreeExpr]-(ip-tree_starts[tree]+argstogo);
	if((m>=0) && (m<maxargs)) maxargs = m;
	}
#endif
//cout <<"SubInit Calling roulette() 1\n";
        i=funcbag->roulette();
//cout <<"SubInit roulette() 1, returned "<<i<<endl;
        // terminal required
        if (maxargs == 0 || depth <= 0 )
                while (funclist[i]->argnum>0) i=funcbag->roulette();

        // node required
        else if (isfull || mindepth>0)		//support for MULTREE
                if (maxargs > 5)
                        while (funclist[i]->argnum == 0) i=funcbag->roulette();
                else                    // not enough room.
					// Take what you can get. 
					// NB may cause prog to loop WBL 
					// Hence added assert 25/8/94WBL
                        while (funclist[i]->argnum>maxargs) i=funcbag->roulette();

        // terminal allowed 50% of the time
        else
                if (rnd(2))             // terminal required
                        while (funclist[i]->argnum>0) i=funcbag->roulette();
                else            // node required
                        if (maxargs > 5)
                                while (funclist[i]->argnum == 0) i=funcbag->roulette();
                        else                    // not enough room.  Take what you can get.
                                while (funclist[i]->argnum>maxargs)
i=funcbag->roulette();
	if (funclist[i]->varnum == 0)
        	{SETNODE(expr[ip],i,0);} 	//WBL bug fix 001 12-May-94
	else
        	{SETNODE(expr[ip],i,rnd(funclist[i]->varnum));}
        ip++;
        thisargs=funclist[i]->argnum;
        argstogo += funclist[i]->argnum-1;


        const int trunk = (thisargs>0 && params->params[pInitSparse])? 
	                   rnd(thisargs) : -1;
        for (i=0;i<thisargs;i++)
        {
	  if(trunk == -1 || trunk == i )
#ifdef MULTREE
                SubInit(funcbag,argstogo,depth-1,isfull, mindepth-1, tree);
	  else
                SubInit(funcbag,argstogo,0,      isfull, 0,          tree);
#else
                SubInit(funcbag,argstogo,depth-1,isfull, mindepth-1);
	  else
                SubInit(funcbag,argstogo,0,      isfull, 0);
#endif
                argstogo--;
        }
//cout <<"ENDSubI( " << argstogo << ", " << depth << ", " << isfull << 
//", " << isstart << " )\n";
}

void Chrome::Traverse() {             // skip the next expression
        int args = 1;
        while (args > 0)
        {
             args=args+ARGNUM()-1;
             ip++;
        }
}


void Chrome::TraverseGlobal() {             // skip the next expression
#ifndef FASTTRAVERSE
        int args = 1;
        while (args > 0)
        {
             args=args+PARGNUM(IP)-1;
             IP++;
        }
#endif /*FASTTRAVERSE*/
}

int Chrome::Depth(int ip, int end)    //Relative depth of end from ip
{                                     //assume end>=start
//        int tree_stack[end-ip+2];
        int* tree_stack = new int[end-ip+2];
//      int tree_stack [EXPRLEN];
        tree_stack[0] = -1;
	int sp = 0;
        while (ip<=end)
        {
		while (tree_stack[sp] == 0) { //remove finished branches
		        --sp; //-O3 bug in mips4 cc?
			tree_stack[sp] = tree_stack[sp] - 1;
		} 
		tree_stack[++sp] = ARGNUM();
		ip++;
        }
	delete[] tree_stack;
	return sp;
}//end Depth

int Chrome::DepthMax(int ip, int end) //Max Relative depth 
{                                     //assume end>=start
//        int tree_stack[end-ip+2];
        int* tree_stack = new int[end-ip+2];
//      int tree_stack [EXPRLEN];
        tree_stack[0] = -1;
	int sp = 0;
	int max = 0;
        while (ip<=end)
        {
		while (tree_stack[sp] == 0) { //remove finished branches
		        --sp; //-O3 bug in mips4 cc?
			tree_stack[sp] = tree_stack[sp] - 1;
		} 
		tree_stack[++sp] = ARGNUM();
		if(sp>max) max = sp;
		ip++;
        }
	delete[] tree_stack;
	return max;
}//end DepthMax

int Chrome::DepthMax(const int in) {ip=in; return DepthMax();}
int Chrome::DepthMax() //Max Relative depth 
{
        int max = 1;
        ip++;
        const int args = ARGNUM();
	for(int i=0;i<args;i++) {
	  const int d = DepthMax()+1;
	  if(d>max) max = d;
	}
	return max;
}//end DepthMax


                                // Write the next subexpression to a stream
                                // pretty flag controls parentheses, indenting
                                // PRETTY_NONE = no indenting, just a stream, with parens
                                // PRETTY_NOPARENS = indented, one function per line, no parens
                                // PRETTY_PARENS = indented, one function per line, parens
void Chrome::SubWrite(ostream& ostr,const int pretty, const int depth, 
		      const int displaydepthlimit,
		      int* outsize, int* outdepth) {
        int args,i;
	const int sd = params->params[pTrace] & pTraceDynamic;
        int size     = (outsize !=NULL)? *outsize  : 0;
	int maxdepth = (outdepth!=NULL)? *outdepth : 0;
        ip++;
        args = ARGNUM();
        if(depth<=displaydepthlimit) {
        if (pretty)                                     // indent?
        {
                ostr<<endl;
                for (i=0;i<depth;i++){
                                ostr << " : ";
                }
        }
        else
                ostr << ' ';
        if (args>0)
        {                               // write (FUNC args)
	  if(!sd) {
                if (pretty != PRETTY_NOPARENS)
                        ostr << '(';
                ostr << funclist[FUNCNUM(expr[ip])]->getprint(this,ip,scratch) <<flush;
	  }
                while (args-- > 0) {
		  int s = 1;
		  int d = 0;
		  SubWrite(ostr,pretty,depth+1,displaydepthlimit,&s,&d);
		  size += s;
		  if(d > maxdepth) maxdepth = d;
		}
	  if(!sd) {
                if (pretty != PRETTY_NOPARENS)
                        ostr << ")";
	  } else 
                ostr << " " << size << "," << maxdepth;
        }
        else {                          // write an atom
	  if(!sd) ostr << funclist[FUNCNUM(expr[ip])]->getprint(this,ip,scratch) <<flush;
        }
	} else {
	  if((depth-1)==displaydepthlimit) 
                ostr << "..";
                while (args-- > 0) {
		  int s = 1;
		  int d = 0;
		  SubWrite(ostr,pretty,depth+1,displaydepthlimit,&s,&d);
		  size += s;
		  if(d > maxdepth) maxdepth = d;
		}
		//ostr << "[" << *size << "," << *maxdepth << "]";
	}//end if displaydepthlimit
        if(outsize !=NULL) *outsize  = size;
	if(outdepth!=NULL) *outdepth = 1+maxdepth;
}

/*<>*/
#ifndef MULTREE
MULTREE required by following routines
#endif
#define XFUNCNUM(c) c->expr[c->ip].op
#define XARGNUM(c)  c->funclist[XFUNCNUM(c)]->argnum
BOOL ipmatch(const BOOL strict, const Chrome* p1, const Chrome* mate) {
  return (strict)? 
            (XFUNCNUM(p1) == XFUNCNUM(mate)) :
             (XARGNUM(p1) == XARGNUM(mate));
}//end ipmatch

void findonepoint(const int tree, const BOOL strict, Chrome* p1, Chrome* mate, 
		  int& thissubptr,int& matesubptr) {	
  int* matchp1   = new int[p1->params->params[pMaxExpr]];
  int* matchmate = new int[p1->params->params[pMaxExpr]];
  assert(  p1->tree_starts[NUM_TREES]>0);
  assert(mate->tree_starts[NUM_TREES]>0);
  int size = 0;
  //p1->ip   =   p1->tree_starts[tree];
  //p1->Traverse();
  const int end = p1->tree_starts[tree+1];
  p1->ip       =   p1->tree_starts[tree];
  mate->ip     = mate->tree_starts[tree];

  do {
    assert(size<p1->params->params[pMaxExpr]);
    matchp1[size]   = p1->ip;   //include roots of non-matching subtrees
    matchmate[size] = mate->ip;
    size++;
    if(ipmatch(strict,p1,mate)) {
      p1->ip++; 
      if(p1!=mate) mate->ip++;
    }
    else {
      p1->Traverse();
      if(p1!=mate) mate->Traverse();
    }
  } while(p1->ip < end);

  const int t = rnd(size);
  thissubptr = matchp1[t];
  matesubptr = matchmate[t];
  delete matchmate;
  delete matchp1;
}//end findonepoint

void findtwopoints(const BOOL internal, const BOOL fair, const int tree,
		   Chrome* p1, Chrome* mate, 
		   int& thissubptr, int& matesubptr) {	
  if(internal) {
    thissubptr=p1->GetIntNode(tree);
    if(!fair)
      matesubptr=mate->GetIntNode(tree);
  }
  else {
    thissubptr=p1->GetAnyNode(tree);
    if(!fair)
      matesubptr=mate->GetAnyNode(tree);
  }
}//end findtwopoints
/*<>*/
/*<>*/

void InsertNodes(const int newlen, 
		 node* newexpr,        const int thissublen, 
		 const node* mateexpr, const int matesublen ) {
  const int remlen = newlen-thissublen;
  node* temp = new node[remlen];

  memcpy(temp,&newexpr[thissublen],remlen*sizeof(node));
  memcpy(newexpr,mateexpr,matesublen*sizeof(node));
  memcpy(&newexpr[matesublen],temp,remlen*sizeof(node));

  delete[] temp;
}//end InsertNodes

void Chrome::InsertTree(const int tree, const int thislen,
                        const int thissubptr, const int thissublen,
  const node* mateexpr, const int matesubptr, const int matesublen) {
#ifdef PTHREADS
  if(params->params[pThreads]==0) //delay till thread
#endif
  InsertNodes(thislen-thissubptr,
	      &expr[thissubptr], thissublen,
	      &mateexpr[matesubptr],matesublen);
#ifdef MULTREE
  for (int i=tree+1; i<NUM_TREES+1; i++) {
    tree_starts[i] = tree_starts[i] + matesublen-thissublen;}
#endif
}//end InsertTree

#ifdef TRACE
//also used for PTHREADS
void birthstats(Chrome::parent* p, 
		const int tree, const int birth, const int thissubptr, const int thissublen) {
  p->birth  = birth;
  p->xpoint = thissubptr;
#ifdef MULTREE
  p->xtree  = tree;
#endif
#ifdef PTHREADS
  p->xsize  = thissublen; //subtree size
#endif
}//end birthstats
#endif //TRACE
/*<>*/




extern int CrossSubTree_t_chnglen = 0;
extern int CrossSubTree_calls = 0;
extern int CrossFH_calls = 0;
extern int CrossFH_choice = 0;
extern int CrossFH_differ = 0;
extern int CrossFH_fitdiffer = 0;

int Debug = FALSE;
int Chrome::treeinfo(const int ip, const int root, const int depth, 
		     const int argnum,
		           nodeinfo subtreeinfo[]) const
{
  subtreeinfo[ip].backip = root;
  subtreeinfo[ip].depth  = depth;
  subtreeinfo[ip].argnum = argnum;
  int size = 1;
  if(Debug) cout<<"("<<flush;
  for(int i = 0; i<ARGNUM(); i++) 
    size += treeinfo(ip+size,ip,depth+1,i,subtreeinfo);
  subtreeinfo[ip].size = size;
  if(Debug) cout<<")="<<size<<" "<<flush;
  return size;
}

inline void Chrome::route(const nodeinfo n[], const int nip, int route[])
{
  int ip = nip;
  for(int i=n[nip].depth;i>0;i--) {
    route[i] = n[ip].argnum;
    ip       = n[ip].backip;
  }
}
inline int Chrome::compareroutes(const nodeinfo a[], const int aip, 
				 const nodeinfo b[], const int bip )
{
  int* ra = new int[a[aip].depth+1];
  int* rb = new int[b[bip].depth+1];

  route(a,aip,ra);
  route(b,bip,rb);
  int i;
  for(i=1; ra[i]==rb[i] && i<=a[aip].depth && i<=b[bip].depth; i++) {;};
    //    cout<<"ra["<<i<<"] "<<ra[i]<<" "
    //	<<"rb["<<i<<"] "<<rb[i]<<" "<<flush;

  delete[] rb;
  delete[] ra;
  //cout<<" mismatch "<<i-1<<endl;
  return i-1;
}

/*----
int Chrome::treesizes(const int ip, 
		            int subtreesizes[]) const
{
  int size = 1;
  if(Debug) cout<<"("<<flush;
  for(int i = 0; i<ARGNUM(); i++) 
    size += treesizes(ip+size,subtreesizes);
  subtreesizes[ip] = size;
  if(Debug) cout<<")="<<size<<" "<<flush;
  return size;
}----*/

extern BOOL smallertaken = FALSE; //for reporting effect of depth or size limits
#ifdef MULTREE
Chrome* Chrome::CrossTree(Chrome* mate, int tree, const BOOL fair, 
	const BOOL onepoint, const BOOL strict )  //work on one tree only
#else
Chrome* Chrome::CrossTree(Chrome* mate, const BOOL fair)
#endif
// Return a new Chrome that is a cross with mate

{
#ifdef PDGP
  if(params->params[pPDGP]) assert(0==1);// should not get here
#endif
//debug//cout<<"\nSTART CrossTree, this=";
//this->write(PRETTY_NONE,cout);
//cout<<"\nCrossTree, mate=";
//mate->write(PRETTY_NONE,cout);
        Chrome* newexpr = new Chrome(params,probl,funclist,constlist,FALSE);
#ifdef STACK_LIB
	tree = rnd(NUM_TREES);  //work on one tree only
#endif
#ifdef TRACE
	num_children++;
	newexpr->mum.birth = birth;
#endif
        newexpr->nfitness->Copy(this->nfitness); //to suuport STORE_FIT
	//load of stuff done by Dup but dont need as also done by new Chrome()
        newexpr->ip=0;
#ifdef TRACE
	mate->num_children++; //this incremented by Copy()
#endif
        int thislen;        // total length of expression
        int thissubptr;     // pointer to subexpression
        int thissublen;
        int matelen;
        int matetreelen;
        int matesubptr;
        int matesublen;

#ifdef MULTREE
	assert(      tree_starts[NUM_TREES]>0);
	assert(mate->tree_starts[NUM_TREES]>0);
        thislen=      tree_starts[NUM_TREES];
        matelen=mate->tree_starts[NUM_TREES];
#else
        thislen=SubLen(0);
	matelen=mate->SubLen(0);
#endif
       	nodeinfo* thistreeinfo  = NULL;
       	nodeinfo* matetreeinfo  = NULL;
	BOOL*     matesizeavail = NULL;
	int chnglen;
#ifdef MULTREE
#ifndef STACK_LIB
	do { //while static_check says not ok
#endif
#endif
	do { //while fair && size checks fail
	if(onepoint) 
	  findonepoint(tree,strict,this,mate,thissubptr,matesubptr);
	else 
	  findtwopoints(rnd(101)>params->params[pUnRestrictWt], //should this be 100<??
			fair,tree,this,mate,thissubptr,matesubptr);
        thissublen = SubLen(thissubptr);
	if(fair){
	      if(Debug) cout<<"\nthissublen "<<thissublen<<" "<<endl;
	      if((params->params[pCrossFairLimit] >0 ) &&
		 (2*thissublen+1 > params->params[pCrossFairLimit]))
		matesublen = 0;
	      else { //choose matesublen so mean size change=0 if possible
		/* init treeinfo only once 1 Jan 2019 */
		if(thistreeinfo==NULL) {
		  thistreeinfo  = new nodeinfo[thislen];
		  matetreeinfo  = new nodeinfo[matelen];
		  const int maxsize = (2*thislen+1>matelen)? 2*thislen+1 : matelen;
		  matesizeavail = new BOOL[maxsize+1];
		  memset(matesizeavail,0,(maxsize+1)*sizeof(BOOL));
#ifdef MULTREE
		  treeinfo(tree_starts[tree],-1,1,0,thistreeinfo);
		  matetreelen=mate->treeinfo(mate->tree_starts[tree],-1,1,0,matetreeinfo);
		  assert(matelen==matetreelen); //POLY_LIB debug code
		  {for(int i=mate->tree_starts[tree];i<mate->tree_starts[tree]+matetreelen;i++)
		    matesizeavail[matetreeinfo[i].size] = TRUE;}
#else
		  treeinfo(0,-1,1,0,thistreeinfo);
		  mate->treeinfo(0,-1,1,0,matetreeinfo);
		  matetreelen=matelen;
		  {for(int i=0;i<matetreelen;i++) 
		    matesizeavail[matetreeinfo[i].size] = TRUE;}
#endif
		}//end init treeinfo
		
		const int maxmatesublen = (params->params[pCrossFairLimit]<0)? 
		                          matelen : 2*thissublen+1;
		int l, numneg=0, sumneg=0, numpos=0, sumpos=0;
		for(l=1;l<thissublen;l++)
		  if(matesizeavail[l]) {numneg++; sumneg += thissublen-l;}
		for(l=thissublen+1;l<=maxmatesublen;l++)
		  if(matesizeavail[l]) {numpos++; sumpos += l-thissublen;}
		//cout<<" sumneg "<<sumneg<<" sumpos "<<sumpos<<" "<<flush;
		if(sumpos==0 || sumneg==0)
		  matesublen = (matesizeavail[thissublen]==0)? 0 : thissublen;
		else {
		  const float pzero= matesizeavail[thissublen]? 1.0/thissublen:
		                                                0.0;
		  const float ppos = (1.0 - pzero)/
		                     (numpos + float(numneg*sumpos)/sumneg);
		  const float pneg = (1.0 - pzero)/
		                     (numneg + float(numpos*sumneg)/sumpos);
		  const float total = numneg*pneg + pzero + numpos*ppos;
		  assert(total<1.00001&&total>.99999);

		  CSelector* wheel = new CSelector(EStraight,maxmatesublen);
		  for(l=1;l<thissublen;l++)
		    if(matesizeavail[l])        wheel->add(pneg,l);
		  if(matesizeavail[thissublen]) wheel->add(pzero,thissublen);
		  for(l=thissublen+1;l<=maxmatesublen;l++)
		    if(matesizeavail[l])        wheel->add(ppos,l);
		  matesublen = wheel->roulette();
		  delete wheel;
		}
	      }
	      if(Debug) cout<<" matesublen "<<matesublen<<" "<<flush;

	      //Chose one of the subtrees of size=matesublen
	      //either at random or choose closes to thissubptr
	      const BOOL homologous = (/*params->params[pCrossFHWt]>0 && */
				       rnd(100) < params->params[pCrossFHWt]) ?
		                          TRUE : FALSE;
	      BOOL Hactive = FALSE;
	      BOOL Hchoice = FALSE;
	      int maxmisdepth = 0;
	      matesubptr = -1;
	      const int r = rnd(matetreelen);
	      for(int i=0;i<matetreelen;i++) {
		int ip = mate->tree_starts[tree] + (i+r)%matetreelen;
		if(Debug) cout<<ip<<" "<<flush;
		if(matetreeinfo[ip].size==matesublen) {
		  if(homologous) {
		    if(Hactive) Hchoice = TRUE;
		    Hactive = TRUE;
		    const int misdepth = compareroutes(thistreeinfo,thissubptr,
						       matetreeinfo,ip       );
		    if(misdepth>maxmisdepth) {
		      maxmisdepth = misdepth;
		      matesubptr=ip;
		    }
		  }
		  else {
		    matesubptr=ip;
		    break;
		  }}
	      }
	      if(Hactive) CrossFH_calls++;
	      if(Hchoice) CrossFH_choice++;
	      if(Debug) cout<<" matesubptr "<<matesubptr
			    <<" matesublen "<<matesublen<<endl;
	      assert(matesublen<=3 || matesubptr>=0);//POLY_LIB sanity check
	    }//endif fair
	else
        matesublen=mate->SubLen(matesubptr);

	if(matesubptr>= 0) {//!fair or !matesizeavail[thissublen]
//cout<<"\ntree="<<tree;
//cout<<" thislen="<<thislen<<" matelen="<<matelen;
//cout<<" thissubptr="<<thissubptr<<" thissublen="<<thissublen;
//cout<<" matesubptr="<<matesubptr<<" matesublen="<<matesublen<<endl;
#ifdef CROSSFILE
	if(crossfile) {
		crossfile<<"cross "<<tree<<" ";
		crossfile<<"mum ";
		crossfile<<this->birth<<" ";
		crossfile<<this->nfitness->fvalue<<" ";
		crossfile<<thislen<<" ";	//size of program
		crossfile<<SubLen(this->tree_starts[tree])<<" "; //size of tree
		crossfile<<thissubptr<<" ";	//root of crossover point
		crossfile<<thissublen<<" ";	//size of crossover tree
		crossfile<<Depth(tree_starts[tree],thissubptr)<<" "<<flush;
						//depth of crossover point
		crossfile<<DepthMax(thissubptr,thissubptr+thissublen-1)<<" ";
						//depth of crossover fragment
		crossfile<<"dad "<<flush;
		crossfile<<mate->birth<<" ";
		crossfile<<mate->nfitness->fvalue<<" ";
		crossfile<<matelen<<" ";	//size of program
		crossfile<<mate->SubLen(mate->tree_starts[tree])<<" ";
		crossfile<<matesubptr<<" ";	//root of crossover point
		crossfile<<matesublen<<" ";	//size of crossover tree
		crossfile<<mate->Depth(mate->tree_starts[tree],matesubptr)<<" ";
						//depth of crossover point
		crossfile<<mate->DepthMax(matesubptr,matesubptr+matesublen-1)
			 <<" "<<flush;		//depth of crossover fragment
	}
#endif
	if(((params->params[pMaxDepth] != 0) &&
	   ((Depth(tree_starts[tree],thissubptr) +
	     mate->DepthMax(matesubptr,matesubptr+matesublen-1)) -1 >
	    params->params[pMaxDepth])) ||
	   (thislen+matesublen-thissublen >params->params[pMaxExpr])){
	// take smaller side of swap
	       assert(matelen+thissublen-matesublen<=params->params[pMaxExpr]);
		chnglen = thissublen-matesublen;
		if(!smallertaken) {
		  smallertaken=TRUE;
		  cout<<"CrossTree(*,"<<tree<<","<<fair<<","<<onepoint<<","
		      <<strict<<") taking smaller "<<flush;
		  cout<<"("<<params->params[pMaxDepth]<<"!=0 && "
		      <<Depth(tree_starts[tree],thissubptr)<<"+"
		      <<"mate->DepthMax("<<matesubptr<<","
		      <<matesubptr+matesublen-1<<")"
		      <<mate->DepthMax(matesubptr,matesubptr+matesublen-1)
		      <<" > "
		      <<params->params[pMaxDepth]<<") || ("
		      <<thislen<<"+"<<matesublen<<"-"<<thissublen
		      <<" > "<<params->params[pMaxExpr]<<")\n";
		}//endif firsttime smaller taken
#ifdef CROSSFILE
		if(crossfile) crossfile<< ", 1 ";
#endif
#ifdef PTHREADS
                if(params->params[pThreads]==0) //delay till thread
#endif
                memcpy(newexpr->expr,mate->expr,matelen*sizeof(node));
#ifdef MULTREE
	        memcpy(newexpr->tree_starts,mate->tree_starts,sizeof(tree_starts));
		assert(newexpr->tree_starts[NUM_TREES]>0);
#endif
		newexpr->InsertTree(tree,matelen, matesubptr,matesublen,
				    expr,         thissubptr,thissublen);
#ifdef TRACE
		birthstats(&newexpr->mum,tree,mate->birth,matesubptr,matesublen);
		birthstats(&newexpr->dad,tree,this->birth,thissubptr,thissublen);
#endif
        }
        else {
		chnglen = matesublen-thissublen;
#ifdef CROSSFILE
		if(crossfile) crossfile<< ", 0 ";
#endif
#ifdef PTHREADS
                if(params->params[pThreads]==0) //delay till thread
#endif
                memcpy(newexpr->expr,expr,thislen*sizeof(node));
#ifdef MULTREE
	        memcpy(newexpr->tree_starts,tree_starts,sizeof(tree_starts));
		assert(newexpr->tree_starts[NUM_TREES]>0);
#endif
		newexpr->InsertTree(tree,thislen, thissubptr,thissublen,
				    mate->expr,   matesubptr,matesublen);
#ifdef TRACE
		birthstats(&newexpr->mum,tree,this->birth,thissubptr,thissublen);
		birthstats(&newexpr->dad,tree,mate->birth,matesubptr,matesublen);
#endif
        }
        }//endif matesubptr ok
#ifdef MULTREE
#ifndef STACK_LIB
	}while (matesubptr<0);
	}while (probl->static_check(newexpr,tree)>rnd(100));
#endif
#endif
#ifdef CROSSFILE
if(crossfile) {
crossfile<<same(newexpr);       //same as mum?
crossfile<<" ";
crossfile<<mate->same(newexpr); //same as dad?
crossfile<<" ";
crossfile<<same(mate);          //mum same as dad?
crossfile<<" ";
}
#endif
//cout<<"CrossTree, newexpr=";
//newexpr->write(PRETTY_NONE,cout);
//cout<<endl;

	CrossSubTree_t_chnglen += chnglen;
	CrossSubTree_calls++;

	if(matesizeavail) delete[] matesizeavail;
	//delete[] matetreesizes;
	if(matetreeinfo) delete[] matetreeinfo;
	if(thistreeinfo) delete[] thistreeinfo;

        return newexpr;
}

Chrome* Chrome::JoinTree(Chrome* mate, int tree)
// Return a new Chrome that is a cross at root with mate
{
#ifdef PTHREADS
        assert(params->params[pThreads]==0);   // not implemented
#endif
        Chrome* newexpr = new Chrome(params,probl,funclist,constlist,FALSE);
#ifdef TRACE
	num_children++;
	newexpr->mum.birth = birth;
#endif
        newexpr->nfitness->Copy(this->nfitness); //to suuport STORE_FIT
        newexpr->ip=0;
#ifdef TRACE
	mate->num_children++; //this incremented by Copy()
#endif	
	assert(      tree_starts[NUM_TREES]>0);
	assert(mate->tree_starts[NUM_TREES]>0);
	const int thislen    = tree_starts[NUM_TREES];// total length
	const int thissubptr = tree_starts[tree]; // pointer to subexpression
        const int thissublen = SubLen(thissubptr);
        const int matesubptr = mate->tree_starts[NUM_TREES-1];
        const int matesublen = mate->SubLen(matesubptr);

#ifdef PTHREADS
        if(params->params[pThreads]==0) //delay till thread
#endif
	memcpy(newexpr->expr,expr,thislen*sizeof(node));
	memcpy(newexpr->tree_starts,tree_starts,sizeof(tree_starts));
	assert(newexpr->tree_starts[NUM_TREES]>0);

	if(((params->params[pMaxDepth] != 0) &&
	   ((DepthMax(tree_starts[tree]) + 1 > params->params[pMaxDepth])||
	     mate->DepthMax(mate->tree_starts[tree]) + 1 > 
	                                       params->params[pMaxDepth])) ||
	   (thislen+matesublen+1 >params->params[pMaxExpr])){
	// no space for joint
		if(!smallertaken) {
		  smallertaken=TRUE;
		  cout<<"JointTree(*,"<<tree<<") too big or too deep "<<flush;
		  cout<<"("<<params->params[pMaxDepth]<<"!=0 && "
		      <<"DepthMax{"<<DepthMax(tree_starts[tree])<<"}+1 or "
		      <<"mate->DepthMax{"
		      <<mate->DepthMax(mate->tree_starts[tree])<<"}+1"<<" > "
		      <<params->params[pMaxDepth]<<") || ("
		      <<thislen<<"+"<<matesublen<<"+1"
		      <<" > "<<params->params[pMaxExpr]<<"))\n";
		}//endif firsttime smaller taken
		return newexpr;
        }
        else {
	  node root;
	  SETNODE(root,probl->RndOp(tree,2),0); //any binary function

	  newexpr->InsertTree(tree,  thislen, thissubptr,thissublen,
			      &root, 0,                   1);

	  newexpr->InsertTree(tree, 1+thislen-thissublen, 1+thissubptr, 0,
			      expr, thissubptr,           thissublen);

	  newexpr->InsertTree(tree, 1+thislen, 1+thissubptr+thissublen, 0,
			      mate->expr, matesubptr,     matesublen);

#ifdef TRACE
		birthstats(&newexpr->mum,tree,this->birth,thissubptr,thissublen);
		birthstats(&newexpr->dad,tree,mate->birth,matesubptr,matesublen);
#endif
        }

	//sanity check
	const int newlen = newexpr->tree_starts[NUM_TREES-1] + 
	                   newexpr->SubLen(newexpr->tree_starts[NUM_TREES-1]);
	assert(newlen == 1 + thislen + matesublen);
	assert(newexpr->tree_starts[NUM_TREES]>0);
	assert(newlen == newexpr->tree_starts[NUM_TREES]);
        return newexpr;
}//end JoinCross

#ifdef MULTREE
int Chrome::GetAnyNode(int tree)        // get a node for crossover
{
	assert(tree_starts[NUM_TREES]>0);
	assert(tree<NUM_TREES);
        return tree_starts[tree] + rnd(tree_starts[tree+1]-tree_starts[tree]);
}
#else
int Chrome::GetAnyNode()                // get a node for crossover
{
        return rnd(SubLen(0));
}
#endif

#ifdef MULTREE
int Chrome::GetIntNode(int tree)         // get an internal node for crossover
{
	assert(tree_starts[NUM_TREES]>0);
	assert(tree<NUM_TREES);
        int rval;
        int l=tree_starts[tree+1]-tree_starts[tree];
        if (l>2)
        {
                rval=tree_starts[tree] + rnd(l);
                while (funclist[FUNCNUM(expr[rval])]->argnum <1)
                        rval=tree_starts[tree] + rnd(l);
        }
        else
                rval=tree_starts[tree];
#else
int Chrome::GetIntNode()                // get an internal node for crossover
{
        int rval;
        int l=SubLen(0);
        if (l>2)
        {
                rval=rnd(l);
                while (funclist[FUNCNUM(expr[rval])]->argnum <1)
                        rval=rnd(l);
        }
        else
                rval=0;
#endif
        return rval;
}

inline BOOL MutateNode(const int ip, Function** funclist, 
		       Problem* probl, const int tree, 
	   node* expr)
{
  const int old_f = expr[ip].op;
  const int f = probl->RndOp(tree,ARGNUM());
/*RndOp is faster
  const int args=ARGNUM();
  do { //find a function with same number of  arguments
    f=funcbag->roulette();
  } while (funclist[f]->argnum!=args);
*/

#ifdef CROSSFILE
  if(crossfile) {
    crossfile<<ip<<" ";
    crossfile<<funclist[expr[ip].op]->name<<" ";
}
#endif /*CROSSFILE*/
//  if (funclist[f]->varnum == 0)
    {SETNODE(expr[ip],f,0);} //WBL bug fix 001
//  else {
//    assert(1==0); //return doesnt support varnum change
//    SETNODE(expr[ip],f,rnd(funclist[f]->varnum));}
#ifdef CROSSFILE
  if(crossfile) {
    crossfile<<funclist[expr[ip].op]->name<<" ";
}
#endif /*CROSSFILE*/
  return (f != old_f);
}//end MutateNode
                                // Mutate the current Chrome
void Chrome::Mutate()                   // mutate nodes with rate r
{
#ifdef PTHREADS
        assert(params->params[pThreads]==0);   // not implemented
#endif
//replace about pMuteRate/1024 of nodes in program with a function
//that takes the same number of arguments
        int end;
        int rate=params->params[pMuteRate];
#ifdef CROSSFILE
#ifdef MULTREE 
	assert(tree_starts[NUM_TREES]>0);
        int len=tree_starts[NUM_TREES];
#else
        int len=SubLen(0);
#endif
	node* temp = new node[len];
	memcpy(temp,expr,len*sizeof(node));
#endif /*CROSSFILE*/
#ifdef MULTREE
for (int t = 0; t < NUM_TREES; t++ ) {
#ifdef PDGP
  if(params->params[pPDGP]) {
        ip  = params->pdgpindex[params->net_limits[t]];
	end = params->pdgpindex[params->net_limits[t+1]];
      }
  else
#endif /*PDGP*/
      {
	/*the following takes about 5% of runtime in big programs and 100% Mutate
	ip = tree_starts[t];
	Traverse();
        end=ip;
	ip = tree_starts[t];
	speed up Jan 2019 */
	assert(t<NUM_TREES);
	end=tree_starts[t+1];
	assert(end>0);
      }
#else
      {
        ip=0;
        Traverse();
        end=ip;
        ip=0;
      }
#endif
#ifdef CROSSFILE
	if(crossfile) {
		crossfile<<"point "<<t<<" ";
		crossfile<<"mum ";
		crossfile<<this->birth<<" "; //requires TRACE
		crossfile<<this->nfitness->fvalue<<" ";
		crossfile<<len<<" ";	//size of program (not used!)
		crossfile<<end-ip<<" "; //size of tree
		crossfile<<ip<<" ";	//root of crossover point
		crossfile<<"M ";	//size of crossover tree
}
#endif /*CROSSFILE*/
        if(rate>0) {
        while (ip<end)
        {
//rnd takes about 12% of runtime in big programs and 100% Mutate
	        //for VAN leave variables alone and constants to MutateC
                if((ARGNUM()>0 || params->params[pMuteConstWt]==0) && 
		   rnd(1024) < rate) 
                {
		  MutateNode(ip,funclist,probl,t,expr);
                }
                ip++;
        }
      }
  else { //do -rate mutations (no check multiple changes to same node)
    for(int n=0; n<(-rate);n++) { 
      //choose a node to mutate but avoid opcodes if only one of its arity
      {int args;
      do { ip   = tree_starts[t] + rnd(end-tree_starts[t]); 
	   args = ARGNUM();
      } while(probl->Arityc(t,args)<=1); }

      do{} while(!MutateNode(ip, funclist, probl, t, expr));
    }
  }//end rate positive
#ifdef MULTREE
        }//endfor NUM_TREES
#endif
#ifdef TRACE
	mum.xpoint = -2;
#endif
#ifdef CROSSFILE
	if(crossfile) {
		crossfile<< ", 0 ";

crossfile<<(memcmp(expr,temp,sizeof(node)*len)==0);//same as mum?
crossfile<<" ";
crossfile<<"0";//mate->same(newexpr); //same as dad?
crossfile<<" ";
crossfile<<"0";//same(mate);          //mum same as dad?
crossfile<<" ";
}
	delete[] temp;
#endif /*CROSSFILE*/
}


void Chrome::MutateC()                  // jiggle constants with a random walk 
{
  if(probl->funcconstmax==probl->funcconstmin) return; //only one constant
//replace about pMuteRate/1024 of constant in program with another
  const int rate=params->params[pMuteRate];
#ifdef MULTREE
  assert(tree_starts[NUM_TREES]>0);
for (int t = 0; t < NUM_TREES; t++ )
{ ip = tree_starts[t];
#else
{ ip=0;
#endif
#ifdef PDGP
  if(params->params[pPDGP]) assert(0==1);// not implemented
#endif
  int args=1;
  while (args > 0) {
    if (rnd(1024) < rate) {
    const int opcode = expr[ip].op;
    if(opcode >= probl->funcconstmin && opcode <= probl->funcconstmax) {
      retval newconst = constlist[opcode];
      //near gaussian noise, sigma = 2 for ints otherwise 5%, mean=0;
      double noise = 0;
      for (int i=0;i<12;i++) noise += drand_0to1()-0.5;
      newconst += noise*((int(newconst)==newconst)? 2.0 : newconst/20.0);

      //find closet match in constlist
      int bot = probl->funcconstmin; //binary chop search, constlist will be inorder
      int top = probl->funcconstmax;
      do {
	int i = (top+bot)/2;
	if(newconst < constlist[i]) top = i - 1;
	else                        bot = i + 1;
      } while(top>=bot);
      /*const TBD*/ int newop =(bot > probl->funcconstmax)? top :
                       (top < probl->funcconstmin)? bot :
	               ((newconst-constlist[top])<= (constlist[bot]-newconst))?
	                                            top : bot;
      /*************************************************************
#define dif(x,y) (x<y)? (y-x):(x-y)
      retval delta = FLT_MAX;
      int newop2;
      for(int i=probl->funcconstmin;i<=probl->funcconstmax;i++) {
	const retval d = dif(constlist[i],newconst);
	if(d<delta) { newop2 = i; delta = d; }
      }
#undef dif

      if(newop!=newop2) 
	cout<<"ERROR newconst "<<newconst
	    <<" funcconstmin "<<probl->funcconstmin
	    <<" funcconstmax "<<probl->funcconstmax
	    <<" top "<<top<<" bot "<<bot
	    <<" newop "<<newop<<" newop2 "<<newop2<<endl;
      *************************************************************/

      if(newop==opcode) { //force some change
	if( newop == probl->funcconstmin || 
           (newop <  probl->funcconstmax && rnd(2)==0)) newop++;
	else                                            newop--;}
      //cout<<constlist[opcode]<<"("<<newconst<<")"<<constlist[newop]<<" ";

      SETNODE(expr[ip],newop,);

    }}//end mutation required
    args=args+ARGNUM()-1;
    ip++;
  }
  //cout<<endl;
}//endfor each tree
#ifdef TRACE
 mum.xpoint = -3;
#endif
}//end MutateC

void Chrome::MutateShrink()             // shrink by promoting a subtree
{
#ifdef PDGP
        if(params->params[pPDGP]) assert(0==1);// not implemented
#endif
  assert(tree_starts[NUM_TREES]>0);

  const int tree       = 0;//sort out later
  const int thislen    = tree_starts[NUM_TREES];
  const int thissubptr = GetIntNode(tree);
  const int thissublen = SubLen(thissubptr);
  if(thissublen<=1) return;//no subtree to promote
  const int matesubptr = thissubptr+1+rnd(thissublen-1);//not own root
  const int matesublen = SubLen(matesubptr);
  InsertTree(tree,thislen, thissubptr,thissublen,
	     expr,         matesubptr,matesublen);
#ifdef TRACE
  birthstats(&mum,tree,birth,thissubptr,thissublen);
#endif
}

int MutateSubTree_t_chnglen = 0;
int MutateSubTree_g_chnglen = 0;
int MutateSubTree_calls = 0;
//int MutateSubTree_smaller = 0;
//float MutateSubTree_bias = 0;

void DisplayMStats(ostream & fout ) //WBL
{
fout<<" MutateSubTree "<<MutateSubTree_calls;
fout<<" length change "<<MutateSubTree_g_chnglen<<flush;
fout<<" total "<<MutateSubTree_t_chnglen<<flush;
//fout<<" forced reductions "<<MutateSubTree_smaller
fout<<endl;
}//end DisplayMStats

void Chrome::MutateSubTree(BOOL fair)             // Replace a sub tree
{
#ifdef PDGP
        if(params->params[pPDGP]) assert(0==1);// not implemented
#endif
#ifdef PTHREADS
        assert(params->params[pThreads]==0);   // not implemented
#endif
        int tree,treelen,subtreelen;
	int chnglen;
	int t;
        int depth = params->params[pInitExpr];
#ifdef MULTREE
	assert(tree_starts[NUM_TREES]>0);
        int len=tree_starts[NUM_TREES];
#else
        int len=SubLen(0);
#endif
	node* temp = new node[len];
	memcpy(temp,expr,len*sizeof(node));
#ifdef MULTREE
	t = rnd(NUM_TREES); //choose tree
#endif
#ifndef RAND_TREE /*********************************************************/
#ifdef MULTREE
	assert(t<NUM_TREES);
        tree=tree_starts[t]+rnd(tree_starts[t+1]-tree_starts[t];
#else
        tree=rnd(len);
#endif
        treelen=SubLen(tree);

#ifdef CROSSFILE
	if(crossfile) {
		crossfile<<"cross "<<t<<" ";
		crossfile<<"mum ";
		crossfile<<this->birth<<" "; //requires TRACE
		crossfile<<this->nfitness->fvalue<<" ";
		crossfile<<len<<" ";	//size of program
		crossfile<<SubLen(this->tree_starts[t])<<" "; //size of tree
		crossfile<<tree<<" ";	//root of crossover point
		crossfile<<treelen<<" ";	//size of crossover tree
		crossfile<<Depth(tree_starts[t],tree)<<" "<<flush;
						//depth of crossover point
		crossfile<<DepthMax(tree,tree+treelen-1)<<" ";
						//depth of crossover fragment
}
#endif /*CROSSFILE*/
	if(fair)
	  depth=DepthMax(tree,tree+treelen-1);

	if(params->params[pMaxDepth] != 0) {
#ifdef MULTREE
	  const int tree_depth = Depth(tree_starts[tree],tree);
#else
	  const int tree_depth = Depth(0,tree);
#endif
	  if(depth+tree_depth > params->params[pMaxDepth])
	    depth = params->params[pMaxDepth] - tree_depth;
	}

//write(PRETTY_NONE, cout);
//cout<<"tree "<<tree<<" depth "<<depth<<" "<<"remove "<<treelen<<flush;

	int chnglen;
	do {
        // DO "ramped half and half from 2 to depth"
	// go to some depth less than max, full or half full
	do {
        // if fair keep retrying if length is increased too much
	ip = tree;
	if(fair) {
	  // DO "ramped half and half from 0 to existing depth +/- 1"
#ifdef MULTREE
	  SubInit(probl->funcbag[t],1,depth+rnd(3)-2,rnd(2), 0, t);
//	  SubInit(probl->funcbag[t],1,rnd(depth+1),rnd(2), FALSE, t);
//	  SubInit(probl->funcbag[t],1,depth,rnd(2), FALSE, t);
#else
	  SubInit(probl->funcbag[0],1,rnd(depth+1),rnd(2), 0);
#endif /*MULTREE*/
	}
	else {
#ifdef MULTREE
	SubInit(probl->funcbag[t],1,rnd(depth-1)+1,rnd(2), 1, t);
#else
	SubInit(probl->funcbag[0],1,rnd(depth-1)+1,rnd(2), 1);
#endif
        }
	subtreelen=SubLen(tree); //length of new tree

	chnglen = subtreelen-treelen;
        } while(fair && chnglen > len && len > 15 );
	--depth; //if tree too big try creating it again but smaller
	} while (((len+chnglen) > params->params[pMaxExpr]) && (depth > 1));
#else /*RAND_TREE*********************************************************/
	BOOL ok;
	int tree_depth;
        ostrstream* optr = NULL;
	  do { //if pMuteFairLimit<>0 keep going until accept mutation

	if(tree_starts[0]!=0) cout<<"tree_starts[0] "<<tree_starts[0]<<endl;
	assert(tree_starts[0]==0);
        tree=tree_starts[t]+rnd(tree_starts[t+1]-tree_starts[t]);
	if(tree>=len) cout<<"tree "<<tree<<" len "<<len<<endl;
	assert(tree<len);
        treelen=SubLen(tree);

#ifdef CROSSFILE
	if(crossfile) {
	  if(optr!=NULL) delete optr;
	  optr = new ostrstream();
#define crossfile *optr
		crossfile<<"cross "<<t<<" ";
		crossfile<<"mum ";
		crossfile<<this->birth<<" "; //requires TRACE
		crossfile<<this->nfitness->fvalue<<" ";
		crossfile<<len<<" ";	//size of program
		crossfile<<SubLen(this->tree_starts[t])<<" "; //size of tree
		crossfile<<tree<<" ";	//root of crossover point
		crossfile<<treelen<<" ";	//size of crossover tree
		crossfile<<Depth(tree_starts[t],tree)<<" "<<flush;
						//depth of crossover point
		crossfile<<DepthMax(tree,tree+treelen-1)<<" ";
						//depth of crossover fragment
#define crossfile cout
}
#endif /*CROSSFILE*/
	if(fair)
	  depth=DepthMax(tree,tree+treelen-1);

	if(params->params[pMaxDepth] != 0) {
	  assert(1==0);//bug?
	  int tree_depth = Depth(tree_starts[t],tree);
	  //const int tree_depth = Depth(tree_starts[tree],tree);
	  if(depth+tree_depth > params->params[pMaxDepth])
	    depth = params->params[pMaxDepth] - tree_depth;
	}
	do {//if pMaxExpr exceeded keep retrying until depth<=1
	do {//if fair & pMuteFairLimit keep retry if length increased too much
	//restore original in case changed by previous iteration of this loop
	memcpy(expr,temp,len*sizeof(node)); 
	ip = tree;
	if(fair) {
	  if(params->params[pMuteFairLimit]!=0) {
	  if(ip==tree_starts[t] && (ARGNUM()==0 //create tree as initial pop
//#ifdef ANT_LIB /* cope with special case that tree size 2 is illegal */
//use same code for POLY_LIB, sizes 2 and 4 illegal
	     || (treelen==3 && params->params[pMuteFairLimit]<0)
//#endif /*ANT_LIB*/
	     )) { 
	    cout<<"RAND_TREE Replacing tree "<<treelen<<endl;
	    const int depth = params->params[pInitExpr];
	    const int md    = params->params[pMinInitExpr]-1;
	    //ignore pUniInitWt for timebeing too complicated
                // DO "ramped half and half from pMinInitExpr to depth"
                if (depth > 0) {
                  SubInit(probl->funcbag[t],1,rnd(depth-md)+md,rnd(2), md, t);
		  // go to some depth less than max, full or half full
                } else
                        SETNODE(expr[ip],0,0);    // stick in a stub constant
	    ok = TRUE;
	  }
	  else {
	    if(params->params[pMuteFairLimit]<0) {
//#ifdef ANT_LIB /* cope with special case that tree size 2 is illegal */
//use same code for POLY_LIB, sizes 2 and 4 illegal
//NB do not bury code with side effects inside assert
	      if(treelen==3)      {ok = rand_tree(3,t);        assert(ok);}
	      else if(treelen==4) {ok = rand_tree(5-rnd(3),t); assert(ok);}
	      else
//#endif /*ANT_LIB*/
	      if(((3*treelen)/2) > -params->params[pMuteFairLimit])
		ok = FALSE;
	      else {
		int limit=0;
		do {
		  ok = rand_tree((3*treelen)/2-rnd(1+(treelen/2)*2),t);
		  assert(limit++<1000); //beter than infinite loop?
		} while(!ok);
	      }
	    }
	    else {
	    if(treelen>params->params[pMuteFairLimit]) {
	      ok = FALSE;
	      //cout<<"rand_tree "<<treelen<<","<<t<<" too big"<<endl;
	    }
	    else {
	    //choose another subtree at random and make mutation tree same size
	      const int size=SubLen(tree_starts[t]+rnd(tree_starts[t+1]-tree_starts[t]));
	    ip = tree;//restore ip, used by SubLen
	    ok = rand_tree(size,t);
	    //if(!ok) cout<<"rand_tree("<<size<<","<<t<<") failed"<<endl;
	    }}
	  }
	  } else { //pMuteFairLimit
	  // DO "ramped half and half from 0 to existing depth +/- 1"
	  SubInit(probl->funcbag[t],1,depth+rnd(3)-2,rnd(2), 0, t);
	  }//endif pMuteFairLimit
	}
	else {
        // DO "ramped half and half from 2 to depth" ignore pMinInitExpr too complicated
	// go to some depth less than max, full or half full
	SubInit(probl->funcbag[t],1,rnd(depth-1)+1,rnd(2), 1, t);
        }//endif fair
	subtreelen=SubLen(tree); //length of new tree

	chnglen = subtreelen-treelen;
        } while(fair && params->params[pMuteFairLimit] ==0 && chnglen > len && len > 15 );
	--depth; //if tree too big try creating it again but smaller
	} while (((len+chnglen) > params->params[pMaxExpr]) && (depth > 1));
      } while (fair && params->params[pMuteFairLimit] !=0 && ((!ok) || 
	       (params->params[pMaxDepth] != 0 &&
	        (tree_depth+DepthMax(tree,tree+subtreelen-1)-1 > 
		 params->params[pMaxDepth]                       ))));
#ifdef CROSSFILE
	if(crossfile) {
	  cout<<optr->str()<<flush;
	  if(optr!=NULL) delete optr;
	}
#endif
#endif /*RAND_TREE*/
//if(fair&&chnglen!=0)cout<<" MutateSubTree(fair) change "<<chnglen<<endl;
//else cout<<endl;

	if(len+chnglen>params->params[pMaxExpr]) { //may trip this?
#ifdef CROSSFILE
	  if(!crossfile)
#endif	    
	      {
		cout<<"MutateSubTree too big: len+chnglen=";
		cout<<len<<"+"<<chnglen<<">"<<params->params[pMaxExpr];
		cout<<endl;
	      }
		memcpy(expr,temp,len*sizeof(node)); //restore original
		chnglen = 0;
	}
	else {
	memcpy(expr+tree+subtreelen,temp+tree+treelen,
	       (len-(tree+treelen))*sizeof(node));         // add the rest
	}
#ifdef CROSSFILE
	if(crossfile) {
		crossfile<<"dad "<<flush;
		crossfile<<"0 "; //mate->birth<<" ";
		crossfile<<"0 "; //mate->nfitness->fvalue<<" ";
		crossfile<<"0 ";//matelen<<" ";	//size of program
		crossfile<<"0 ";//mate->SubLen(mate->tree_starts[tree])<<" ";
		crossfile<<"0 ";//matesubptr<<" ";	//root of crossover point
		crossfile<<chnglen+treelen<<" ";//size of crossover tree
		crossfile<<"0 ";//depth of crossover point
		crossfile<<DepthMax(tree,tree+chnglen+treelen-1)
			 <<" "<<flush;		//depth of crossover fragment

		crossfile<< ", 0 ";

crossfile<<(chnglen==0 && memcmp(expr,temp,sizeof(node)*len)==0);//same as mum?
crossfile<<" ";
crossfile<<"0";//mate->same(newexpr); //same as dad?
crossfile<<" ";
crossfile<<"0";//same(mate);          //mum same as dad?
crossfile<<" ";
}
#endif /*CROSSFILE*/
	delete[] temp;
#ifdef MULTREE
	assert(tree_starts[NUM_TREES]>0);
        for(int i=t+1;i<NUM_TREES+1;i++) tree_starts[i] += chnglen;
#endif
#ifdef TRACE
	mum.xpoint = tree;
#endif
//}while (probl->static_check(this,t)>rnd(100)); implemnt sometime?
	MutateSubTree_g_chnglen += chnglen;
	MutateSubTree_t_chnglen += chnglen;
	MutateSubTree_calls++;
	if((MutateSubTree_calls % params->params[pPopSize])==0){
	  MutateSubTree_g_chnglen = 0;
	}
}//end MutateSubTree

#ifdef MULTREE
void Chrome::write(const int pretty,    ostream& ofile, 
		   const int last_tree, const int depthlimit) //WBL
{
assert(last_tree<NUM_TREES);
for (int t = 0; t <= last_tree; t++ )
	   {
	      ip = tree_starts[t]-1;
	      int depth = 0;
	      ofile<<endl; probl->WriteTreeName(t,ofile); ofile<<"\t=";
#ifdef PDGP
	      if(params->params[pPDGP]) writepdgp(t,ofile);
	      else
#endif /*PDGP*/
	     {int size=1; int maxdepth=0;
	      SubWrite(ofile,pretty, depth, depthlimit, &size,&maxdepth); // write the expression
	      assert(tree_starts[t+1]==ip+1);
	     }
	   }
assert(tree_starts[NUM_TREES]>0);
}
#endif

#ifdef TRACE
void Chrome::write_trace(ostream& fout) //WBL
{
	fout << " Created " << birth;
	fout << " fitness ";
	nfitness->write (fout);
#ifdef PDGP
        if(!params->params[pPDGP]) //all programs are same length!
#endif
	fout << " len " << tree_starts[NUM_TREES-1] + 
			   SubLen(tree_starts[NUM_TREES-1]);
	fout << " children " << num_children;
	if ( mum.birth >= 0 )
	      { 
		fout << "\tmum=" << mum.birth;
		fout << "x" << mum.xpoint;
              }
	if ( dad.birth >= 0 )
	      { fout << ",dad=" << dad.birth;
		fout << "x" << dad.xpoint << ",";
              }
#ifdef MULTREE
	if ( mum.birth >= 0 && mum.xtree >= 0 )
		probl->WriteTreeName(mum.xtree);

#endif
	fout << " depth ";
#ifdef MULTREE
	assert(tree_starts[NUM_TREES]>0);
	assert(tree_starts[NUM_TREES] == tree_starts[NUM_TREES-1] + SubLen(tree_starts[NUM_TREES-1]));
	for(int t=0;t<NUM_TREES;t++) 
	  fout << DepthMax(tree_starts[t],
			   tree_starts[t]+SubLen(tree_starts[NUM_TREES-1])-1);
#else
	fout << DepthMax(0,SubLen(0)-1);
#endif
}//end Chrome::write_trace()
#endif //TRACE



BOOL Chrome::iscomment(const int c,istream& istr) const
{       
	if(c=='%') {
#ifdef debug
		do{char t=istr.get();cout<<t;} 
		while((istr.peek()!=EOF)&&(istr.peek()!='\n'));
		cout<<endl;
#else
		do{istr.get();} while((istr.peek()!=EOF)&&(istr.peek()!='\n'));
#endif
		return TRUE;
	}
	else
	       return FALSE;
}

#define EATSPACE c=istr.peek();while(isspace(c)||c==':'||iscomment(c,istr)) {istr.get();c=istr.peek();}
#ifndef debug
#define GETTOK tp=0;c=istr.peek();while(c!=EOF && !isspace(c) && c!=':' && c!='(' && c!=')' && c!='%' &&tp<79) {scratch[tp++]=istr.get(); c=istr.peek();} scratch[tp]='\0'
#else
#define GETTOK tp=0;c=istr.peek(); \
  while(c!=EOF && (!isspace(c)||(tp>0&&c=='\n')) \
	&&c!=':' && c!='(' && c!=')' && c!='%' &&tp<79) \
  { c=istr.get(); \
    if(c!='\n') scratch[tp++]=c; \
    c=istr.peek(); } \
  scratch[tp]='\0';cout<<" `"<<scratch<<"'"<<flush;//debug
#endif //debug

int Chrome::Load(istream& istr,int issource)
// Load an expression from a stream
// issource indicates readable or binary input

{
        int x,tempop; //, tempidx;
        int rval=LOAD_OK;
        node* buf;

        //Get a new clean FitnessValue object
        delete nfitness;
        nfitness  = probl->GetFitnessObj(this);
        birth = 0;

        if (!issource)
        {
#ifdef PDGP
                if(params->params[pPDGP]) assert(0==1);// not implemented
#endif
                for(x=0;x<params->params[pMaxExpr];x++)
                {
                        istr >> tempop;
                        expr[x].op =(unsigned) tempop;
//                        istr >> tempidx;
//                        expr[x].idx=(unsigned) tempidx;
                }
        }
        else            // load a Source file
        {
                ip=0;
                buf=new node[params->params[pMaxExpr]];
		memset(buf,0,sizeof(node)*params->params[pMaxExpr]);
#ifdef PDGP
		node* bufptr = buf;
#endif /*PDPP*/
#ifdef MULTREE
	        char c;
	        int tp;
		rval = LOAD_OK;
	        tree_starts[NUM_TREES] = -1;
              	for (int t = 0; t < NUM_TREES && rval == LOAD_OK; t++ ){
	        tree_starts[t] = ip;
		EATSPACE;
                GETTOK;
                int temp = probl->TreeNameMatch(t,scratch);
#ifdef dbug
cout <<"tree "<<scratch<<((temp==0)? " name ok" : "BAD_TREE_NAME")<<endl;
#endif //debug
		if (temp != 0)
                        rval=BAD_TREE_NAME;
		else{
		EATSPACE;
		GETTOK;
		if (scratch[0] != '=')
                        rval=BAD_TREE_FORMAT;
		else{
#ifdef PDGP
                if(params->params[pPDGP]) rval=SubLoadpdgp(istr,&bufptr,t);
		else
#endif
		rval=SubLoad(istr,buf,t);
		}}
#ifdef debug
		cout<<endl;
#endif
		}//end for each tree
	        tree_starts[NUM_TREES] = ip;
#else
                rval=SubLoad(istr,buf);
#ifdef debug
		cout<<endl;
#endif
#endif //MULTREE
                if (rval==LOAD_OK)
                {
                        delete[] expr;
                        expr=buf;
#ifdef GEN_MEM
			rval = nfitness->Load(istr);
#endif
                } else
                        delete[] buf;
        }

#ifdef MULTREE
        for (int t = 0; t < NUM_TREES && rval == LOAD_OK; t++ ){
	        ip = tree_starts[t];
		int valid = probl->static_check(this,t);
		if (valid>0) cout<<"static_check reports "<<valid
				<<" on tree "<<t<<endl;
	}
#endif
        return rval;
}

#ifdef MULTREE
int Chrome::SubLoad(istream& istr,node* buf, int tree)
#else
int Chrome::SubLoad(istream& istr,node* buf)
#endif
// load a sub-expression from source format
// Return LOAD_OK or an error value
// text expression is in istr
// nodes go into buf

{
        int rval = LOAD_OK;
        char c;
//      char token[80];
        int tp;
        int func;
        int args;
        int vnum;
        char* s;
        scratch[0]='\0';
        EATSPACE;
        if (istr.peek()==')' || istr.peek()==EOF)
        {
                rval=LOAD_TOOFEW;                       // missing expression body
        }
        else if (ip>=params->params[pMaxExpr])
        {
                rval=LOAD_TOOLONG;
        }
        else if (istr.peek()=='(')              // Get an expression and the close paren
        {
                istr >> c;
#ifdef debug
cout<<c<<flush;
#endif //debug
#ifdef MULTREE
                rval = SubLoad(istr,buf,tree);
#else
                rval = SubLoad(istr,buf);
#endif
                if (rval==LOAD_OK)
        {
                        EATSPACE;
                        istr >> c;
                        if (c!=')') {
			        cout<<"LOAD_TOOMANY `"<<c<<"'"<<endl;
                                rval=LOAD_TOOMANY;}
#ifdef debug
			else
				cout<<c<<flush;
#endif //debug
        }
        }
        else                            // Get an expression.  Return if you hit a close paren.
        {
                GETTOK;
                if (strlen(scratch)==0)
                        rval=LOAD_TOOFEW;
                else
                {
//DIRTY FRIG - I DONT USE CONSTANTS (EXCEPT 0 and 1)
//                        if (isdigit(scratch[0]) || scratch[0]=='-')     // it's a number.  Function 0
//                        {
//                                func=atoi(scratch);
//                                if (func>=0-CONSTARRAYSIZE/2 && func<CONSTARRAYSIZE/2)
//                                {
//                                        SETNODE(buf[ip],0,func+CONSTARRAYSIZE/2);
//                                        ip++;
//                                        rval=LOAD_OK;
//                                }
//                                else rval=LOAD_BADCONST;
//                        }
//                        else       // look up a function name
                        {
                                vnum=0;
                                if (strchr(scratch,'_')!=NULL)
                                {
                                        // parse it to take off the variable number?
                                        // This is fancy footwork for functions with operands
                                        // Except for NUMBER (function 0) these functions are
                                        // written as <name>_operand
                                        s=strrchr(scratch,'_');
                                        if (strchr("0123456789",s[1]))  // it is an underscore followed by numbers
                                        {
                                                vnum=atoi(s+1);
                                                if (vnum>=0 && vnum <= CONSTARRAYSIZE)
                                                {
                                                        s[0]='\0';              // found a valid function with variable number appended.  Keep function name.
                                                }
                                                else
                                vnum=0;
                                        }
                                }
#ifdef MULTREE
                                func=FindFunc(scratch,tree);
#else
                                func=FindFunc(scratch);
#endif
                                if (func<0) {
				  cout<<"LOAD_BADFUNC `"<<scratch<<"'"<<endl;
                                        rval=LOAD_BADFUNC;
				}
                                else
                                {
                                        SETNODE(buf[ip],func,vnum);
                                        ip++;
                    rval=LOAD_OK;
                                        // get the arguments
                    args=funclist[func]->argnum;
                                        while(args>0 && rval==LOAD_OK)
                    {
#ifdef MULTREE
                                                rval=SubLoad(istr,buf,tree);
#else
                                                rval=SubLoad(istr,buf);
#endif
                                                args--;
                                        }
//                                      if (rval == LOAD_TOOFEW)
//                      ;               // restore the token for the error message
                                }
                        }
                }
        }
        return rval;
}


//**************************************************************************

#ifdef MULTREE
int Chrome::FindFunc(const char* funcname, const int tree) const// find a function index by name
#else
int Chrome::FindFunc(char* funcname)            // find a function index by name, or return -1
#endif
{
        int rval=-1;
        int i;
        for (i=0;i<probl->funccount && rval<0;i++)
                if (strcmp(funcname,funclist[i]->name)==0)
#ifdef MULTREE
			if (probl->funcbag[tree]->find(i) >= 0)
#endif
                        rval=i;

        return rval;
}

Chrome::~Chrome()
{
        delete[] expr;
        delete nfitness;
}

//**************************************************************************
//// Generic Problem Functions /////////////////////////////////////

Problem::Problem()
// Set up the tables
// You add the primitive functions in your subclass constructor

{
        int i;

        funclist = new Function*[FUNCARRAYSIZE];
	for (i=0;i<NUM_TREES;i++)
	{
	        funcbag[i] = new CSelector(EStraight,256);
		funcbag[i]->reset();
	}
        varlist = new retval[CONSTARRAYSIZE];
        constlist = new retval[CONSTARRAYSIZE]; //expected to be order
        funcconstmin = 0;
        funcconstmax = -1;
        funccount = 0;
        // set up constant table in problem subclass constructor
	memset(constlist,0,sizeof(constlist));
	memset(arityc,0,sizeof(arityc));
	memset(arityop,0,sizeof(arityop));
}

Problem::~Problem()
{
        int i;
        delete[] constlist;
        delete[] varlist;

        for (i=0;i<funccount;i++) delete funclist[i];
	for (i=0;i<NUM_TREES;i++) delete funcbag[i];
#ifdef RAND_TREE
//fails to compile under gcc
//	for (i=0;i<NUM_TREES;i++) delete rand_tree[i];
#endif
        delete[] funclist;
}

FitnessValue* Problem::GetFitnessObj(Chrome* chrome)
{
    FitnessValue* fv = new FitnessValue;
    fv->fvalue=0;
    // fv has to be deleted by the Chrome Object
    return fv;
}

//CSelector* Problem::getfuncbag()         // update the function selector and return its pointer
//{
//        int i;
//        funcbag->reset();
//        for (i=0;i<funccount;i++)
//                funcbag->add(funclist[i]->weight/(2+funclist[i]->argnum),i);
//        return funcbag;
//}

void Problem::addtofuncbag(int tree, Function* f )
//NB must be called in correct sequence, ie after AddF has funccount++
{
	funcbag[tree]->add(f->weight/(2+f->argnum),funccount-1);
	if(funccount>1) {
	  const int a = f->argnum;
	  assert(a>=0 && a<=max_arity);
	  arityc[tree][a]++;
	  if(arityop[tree][a]==0) arityop[tree][a]=funccount-1;
	  assert(arityc[tree][a]+arityop[tree][a]==funccount);//must group functions
      }//endif ignore ConstFunc(0)
}

float Problem::fitness(Chrome* chrome)//, evalnode* xip)
// This is just a stub.  You must add your own virtual fitness function

{
	cout <<"Problem::fitness (stub)\n";//debug
        float f=0;
        // clear variables if required
        // call the installed fitness function
        chrome->nfitness->fvalue = f;
        return f;
}

extern int BestofCount = 0;
extern int BestofRandomCount = 0;
extern int WorstofCount = 0;
extern int WorstofRandomCount =0 ;

Chrome* Problem::Bestof(const PtrChrome list[], const int listsize, 
			int* bestinlist, const int target	)
{ //WBL
#ifdef PARETO
assert (1==19); //should be using non-virtual function in problem class!
#endif /*PARETO*/
if ( bestinlist == NULL ) {int dummy; bestinlist = &dummy;}//discard
*bestinlist   = 0;
Chrome * best = list[0];
for (int i=1; i<listsize; i++)
{	if (list[i]->nfitness->IsBetter(best->nfitness))
	{	*bestinlist = i;
		best        = list[i];
	}
}
//for reporting
if(*bestinlist == 0) {
  BOOL random = true;
  for (int i=1; random && (i<listsize); i++) {	
    if (best->nfitness->IsBetter(list[i]->nfitness)) random = FALSE;
  }
  if(random) BestofRandomCount++;
}
BestofCount++;
return best;

}//end Problem::Bestof


Chrome* Problem::Worstof(const PtrChrome list[], const int listsize,
			 const int timenow, int* worstinlist,
			 const int target ) //WBL
{
#ifdef PARETO
assert (1==19); //should be using non-virtual function in problem class!
#endif /*PARETO*/
// Will select a target that is too old, even if more fit.

*worstinlist = 0;
for (int i=1; i<listsize; i++)
{
	if (list[i]->params->params[pMaxAge] != 0 &&
	     (timenow - list[i]->birth) > list[i]->params->params[pMaxAge])
	     {	*worstinlist = i;
		break;
	     }
	    
	if (list[*worstinlist]->nfitness->IsBetter(list[i]->nfitness))
	     {	*worstinlist = i;
	     }
}
//for reporting
if(*worstinlist == 0) {
  BOOL random = true;
  for (int i=1; random && (i<listsize); i++) {	
    if (list[i]->nfitness->IsBetter(list[0]->nfitness)) random = FALSE;
  }
  if(random) WorstofRandomCount++;
}
WorstofCount++;
return list[*worstinlist];

//////////////////////////////////////////////////////////////////////
//
// Nice Idea, but leads to even more dominance of population by best
// WBL 17-Aug-94. Might be marginaly quicker?
//  		if (target == BestMember) 
//			target = winner; //avoid discarding best if can
//
//		else 
//////////////////////////////////////////////////////////////////////
}//end Problem::Worstof



//**************************************************************************
//// Pop Functions /////////////////////////////////////////////////

Pop::Pop(Problem* prob,ChromeParams* par,const int size, istream* fin)
#ifdef MULTREE
	        :num_crosses_cleared_at(-1)
#endif
// set up a population for a particular problem
// Creates <size> new Chromes
// evaluates the first Chrome.  The rest will be evaluated in the first
// Size-1 calls to generate

{
        int i;
        isAborted=FALSE;
        problem = prob;
        gencount = 0;
        start=time(NULL);
#ifdef PARETO
	clear_rank();
#else
        BestMember=0;
        BestFitness=NULL;
#endif

        params=par;
        par->funccount = prob->funccount;
#ifndef PTHREADS
        if (par->params[pMaxExpr] > EXPRLEN)
                par->params[pMaxExpr]=EXPRLEN;
#endif /*disable limit, EXPRLEN only used in stats.cc*/
	if((par->params[pMaxDepth] !=0) && 
	   (par->params[pInitExpr] > par->params[pMaxDepth]))
                par->params[pInitExpr] = par->params[pMaxDepth];

        popsize=size;
        pop = new Chrome*[size];
#ifdef GENERATIONAL
	if(params->params[pGenerational])
	  newpop = new Chrome*[size];
	//oldpop = new Chrome*[size]; memset(oldpop,0,size*sizeof(Chrome*));
#endif
	int seeds_loaded = 0;

	if(params->params[pSeeds]==0)
	for (i=0;i<size;i++) pop[i]=prob->NewChrome(par,fin);
	else {
	        if(fin->peek()!=EOF) 
		do{int c=fin->get();if(c!='$'){cout<<char(c);};} 
	        	while((fin->peek()!=EOF)&&(fin->peek()!='\n'));
		cout<<endl;
		for (i=0;i<size;i++) { //load all seeds in fin
			pop[i]=prob->NewChrome(par);
			int rval = pop[i]->Load(*fin, TRUE);
			if (rval != 0) {
				if(i==0) cout<<"Load error - "<<rval<<endl;
				else     cout<<"Loaded "<<i<<" seeds"<<endl;
				break;
			}
			seeds_loaded = i+1;
		}
		for (;i<(size-params->params[pSeeds]);i++) pop[i]=prob->NewChrome(par);
	}
//cout << "Created Pop, Now testing pop[0]\n";
        // now eval 1 guy
        initeval=TRUE;
        nexteval=1;
        Fitness(pop[0]);
        InsertMember(0,pop[0],TRUE);
	UpdateBest(  0,pop[0]);
	if (nexteval>=popsize)
	  initeval=FALSE;

	for (;i<size;i++) {
		pop[i] = pop[i%seeds_loaded]->Copy();
		if(rnd(100) < params->params[pMuteSeedWt]) {
		        int prob,wheel;
	                prob=rnd(params->params[pMuteNodeWt]+params->params[pMuteConstWt]+params->params[pMuteShrinkWt]+params->params[pMuteSubTreeWt]+params->params[pMuteFairWt]);
			wheel=params->params[pMuteNodeWt];
			if (prob<wheel)
				pop[i]->Mutate();
			else if (prob<(wheel += params->params[pMuteConstWt]))
				pop[i]->MutateC();
			else if (prob<(wheel += params->params[pMuteShrinkWt]))
				pop[i]->MutateShrink();
			else if (prob<(wheel +=params->params[pMuteSubTreeWt]))
				pop[i]->MutateSubTree(FALSE);
			else 
				pop[i]->MutateSubTree(TRUE);
		}
	}
}

int spin=0;

pthread_mutex_t mutex;
float Pop::Fitness(Chrome* chrome /*, evalnode* xip /*= NULL*/)
// uses the problem to evalue the fitness of member X
// performs setup on the Chrome
// returns the fitness value

{
        //chrome->SetupEval(-1,xip);
#ifdef CLEAR_INTRON
        CLEAR_INTRON;
#endif /*CLEAR_INTRON*/
#ifdef TRACE
        retval rval = problem->fitness(chrome);//,xip);
	if ( params->params[pTrace] & pTraceTest ) {  
#ifdef PTHREADS
//NB do not bury code with side effects inside assert
	  {const int e = pthread_mutex_lock(&mutex);   assert(e==0);}
	  cout << "\nTested ";
	  const int thread = chrome->thread;
	  if(thread>=0) cout << "by thread "<<thread<<" ";
	  chrome->write(PRETTY_NONE, cout);
	  cout<<" fitness "<<rval<<"\n";
	  {const int e = pthread_mutex_unlock(&mutex); assert(e==0);}
#endif /*PTHREADS*/
          }
        return rval;
#else
        return problem->fitness(chrome)//,xip);
#endif //TRACE
}


        // 

#ifndef PARETO
//to support pthreads
void Pop::UpdateBest(const int slot, const Chrome* NewChrome) {
        if (BestFitness==NULL)
        {
               BestMember=slot;
               BestFitness=NewChrome->nfitness;
        }
        else if(NewChrome->nfitness->IsBetter(BestFitness))
        {
cout<<"Improved solution found "<<NewChrome->nfitness->fvalue
    <<" ("<<BestFitness->fvalue<<")"
    <<" in slot "<<slot<<" at "<<gencount<<endl;
                BestMember=slot;
                BestFitness=NewChrome->nfitness;
        }
        else if(slot==BestMember)
        {
	  //if generational BestMember will always be less than slot
		const UINT endpop = (initeval? nexteval : popsize);
cout<<"Scanning whole pop"<<flush;
		problem->Bestof(pop,endpop,&BestMember);
cout<<" Updating BestFitness"<<endl;
                BestFitness=pop[BestMember]->nfitness;
        }
}
#endif /*not PARETO*/

void Pop::InsertMember(int slot,Chrome* NewChrome,int nodelete)
// uses the problem to evalue the fitness of member X
// performs setup on the Chrome
// updates the BestMember and returns the fitness value

{
#ifndef ORIGINAL_RANK
        NewChrome->nfitness->update_rank();
#endif
        NewChrome->birth=gencount;

//        if (!nodelete && pop[slot]->nfitness->IsBetter(NewChrome->nfitness) )
//        {
//	// Check NewChrome is not worse than existing
//#ifdef TRACE
//	        if ( params->params[pTrace] & pTraceDynamic )
//		    {	cout << endl;
//			cout << "Discarding new ";
//                	NewChrome->write_trace(); 
//			cout << "\nRather than " << slot; 
//                	pop[slot]->write_trace(); 
//			cout << endl;
//		    }
//#endif
//		delete NewChrome;
//	}//end discard NewChrome (Nb parents credited with child anyway
//	else
	{
	Chrome* oldchrome = pop[slot];  
        if (!nodelete)
        {
#ifdef GENERATIONAL
	  if(params->params[pGenerational])
	    newpop[slot]=NewChrome;
	  else
#else
	// replace the chrome in this slot
	// remember to delete oldchrome later
#endif
           pop[slot]=NewChrome;
        }// end discard old chrome in this slot

#ifdef PARETO
	if(NewChrome->nfitness->ParetoBest()) {
	  cout<<"Improved solution found "
	      <<" in slot "<<slot<<" at "<<gencount<<endl;
	  NewChrome->write_trace(cout);
	  NewChrome->write(PRETTY_NONE, cout);
	  cout << endl;
	}
#else /*not PARETO*/
        // Update Best Member
	//UpdateBest(slot,NewChrome);
#endif
#ifdef TRACE
	if ( params->params[pTrace] & pTraceDynamic )
	  {  cout << "\nInserting " << slot << " ";
	     NewChrome->write_trace(cout);
	     NewChrome->write(PRETTY_NONE, cout);
	     cout << endl;
          }
#endif //TRACE

#ifdef GENERATIONAL
	if((!nodelete) && (!params->params[pGenerational]))
#else
        if (!nodelete)
#endif
        {
#ifdef TRACE
	if ( params->params[pTrace] & pTraceDynamic )
	     {  cout << endl;
		cout << "Deleting " << slot; oldchrome->write_trace(); }
#endif
            delete oldchrome;
        }// end discard old chrome in this slot

	}//end if add NewChrome to pop
}

Pop::Pop(Problem* prob,ChromeParams* par)
#ifdef MULTREE
	        :num_crosses_cleared_at(-1)
#endif
// set up just the core of a population;; let a subclass allocate the members
{
        isAborted=FALSE;
        problem = prob;
        gencount = 0;
        params=par;
        popsize=0;
        pop=NULL;
#ifndef PARETO
        BestMember=-1;
        BestFitness=prob->GetFitnessObj();    // allocates FitnessValue object on heap
           
        BestFitness->fvalue=1-MAXFLOAT;
#endif
}


Pop::~Pop()
{
      int i;
      for (i=0;i<popsize;i++) delete pop[i];
      delete[] pop;
#ifdef GENERATIONAL
      if(params->params[pGenerational]) {
      delete[] newpop; //when called, population should be in pop only 
      //for (i=0;i<popsize;i++) if(oldpop[i]) delete oldpop[i];
      //delete[] oldpop;
      }
#endif
}

//#ifndef GENERATIONAL
#ifndef PARETO
Chrome* Pop::best()
{
        return pop[BestMember];
}
#endif
//#endif

void Pop::select_candidates(	int target, int tourn_size,
				Chrome* candidates[], int slots[],
				BOOL best) const
{
int  dw, pw, offset;
UINT region;

region = params->params[pMateRadius];
dw = params->params[pDemeWidth];
pw = params->params[pPopWidth];

if (region == 0 || region > popsize) region = popsize;

if (dw==0 || pw==0)
{	dw = 1;
	pw = 1;
}

#ifdef ORIGINAL_DEME
offset=popsize+target-dw/2-region*pw/dw/2;
#else
offset=popsize+target-dw/2-pw*(region/dw/2);
#endif /*ORIGINAL_DEME*/

const int Elitist = params->params[pElitist];

int endstop = 0;
for (int i=0;  i<tourn_size; endstop++)
{	assert(endstop<1000);		//could happend with small demes
	int r = rnd(region);		//select from within deme
	int x = r/dw*pw + r%dw;		//convert from deme to pop index
	slots [i]     = (x+offset)%popsize;
	if ( best || (Elitist==0) || 
	     (!pop[slots[i]]->nfitness->unique_best())) {
		candidates[i] = pop[slots[i]];
		i++;
	}
}
}//end Pop::select_candidates

enum {docross,do1pt,dost1pt,domutate,doanneal,docrossanneal,docopy,dojoin};

Chrome* Pop::selectParent(const int target, int* slotno) const
{                                                                               // target identifies selection region in pop
// select one parent. Return pointer to selected parent
        int ts;
        Chrome* best;

        ts = params->params[pTournSize];
        // Only tournament selection is implemented in GPQUICK.
        // Add your own methods to taste
        if (params->params[pSelectMethod] == ETournament)
        {
//		Chrome* candidates [ts];
//		int     slot      [ts];
		Chrome* candidates [10];
		int     slot      [10];
		int	index;
		select_candidates(target,ts,candidates,slot,TRUE);
		best = problem->Bestof(candidates,ts,&index,target);
#ifdef CROSSFILE
		if(crossfile) {	crossfile<<"P "<<slot[index]<<" "<<flush;}
#endif
		if(slotno!=NULL) *slotno = slot[index];
        }
        return best;
}

#ifdef PTHREADS
void Dummy(Chrome* chrome) {
  const int length = chrome->SubLen(0);
  //do something time consuming
  float F=0;
  for(int i=0;i<length;i++){
    for(int t=0;t<48;t++){
      for(int j=0;j<100;j++) F += sqrt(t+1.0/(1+i));//F++;//
    }
  } 
  chrome->nfitness->fvalue = -(1.0f + 1.0f/F); //bloat
}

//https://computing.llnl.gov/tutorials/pthreads/samples/condvar.c
int slot;
//GCC 7.3.1 objects to int end
Pop* thispop = NULL;
int ids[npthreads];
//evalnode* expr_table[npthreads]; //= calloc(npthreads,sizeof(evalnode));
//pthread_mutex_t mutex;
//unsigned long int active_threads = 0l;

void* thread_fitness(void *t) {
  const int my_id = *(int*)t;
  assert(my_id >= 0);
  assert(my_id < npthreads);
  const int popsize = thispop->popsize;
  int Mum = -1; //avoid second used of lock/unlock
  int Dad;
  Chrome* mum;
  Chrome* dad;
  for (;;) {
    {const int e = pthread_mutex_lock(&mutex);   assert(e==0);}
    /*if(Mum != -1) { //STUFF todo from previous loop
      mum->num_children--;
      if(mum->num_children==0) {delete mum; thispop->pop[Mum] = NULL;}
      dad->num_children--;
      if(dad->num_children==0) {delete dad; thispop->pop[Dad] = NULL;}
    }*/

    //if(my_id<64) active_threads |= (1l<<my_id);
    //const unsigned long int initial_threads = active_threads;
    const int i = slot++;
    {const int e = pthread_mutex_unlock(&mutex); assert(e==0);}
    if(i >= popsize) {
      /**
      if(my_id<64) active_threads &= (~0l)^(1l<<my_id);
      char buff[80];
      sprintf(buff,"%012lx %012lx",initial_threads,active_threads);
      cout << "Last "<<i<<" thread "<<my_id
	   <<" active "<<buff<<"\n";
      pthread_mutex_unlock(&mutex);      /**/
      break;
    }

    Chrome* newguy = (DoCross)? thispop->newpop[i] : thispop->pop[i];
    if(DoCross) { //100% crossover on later generations
    assert(newguy->expr != NULL);
    //assert(newguy->expr == NULL);
    //newguy->expr = new node[thispop->params->params[pMaxExpr]];

    assert(newguy->mum.birth >= 0);
    assert(newguy->dad.birth >= 0);
    Mum = newguy->mum.birth % popsize;
    Dad = newguy->dad.birth % popsize;
    mum = thispop->pop[Mum];
    dad = thispop->pop[Dad];
    assert(mum->num_children>0);
    assert(dad->num_children>0);

    const int mumlen    = mum->tree_starts[NUM_TREES]; //Sublen(0)
    const int mumsubptr = newguy->mum.xpoint;
    const int mumsublen = newguy->mum.xsize;
    const int dadsubptr = newguy->dad.xpoint;
    const int dadsublen = newguy->dad.xsize;

    /*emulate
    memcpy(newexpr->expr,mate->expr,matelen*sizeof(node));
    InsertNodes(thislen-thissubptr,
	      &expr[thissubptr], thissublen,
	      &mateexpr[matesubptr],matesublen);
    */
    const int remlen = mumlen-mumsubptr-mumsublen;

    assert(mumsubptr                  >= 0 && mumsubptr           <  newguy->tree_starts[NUM_TREES]);
    assert(mumsubptr+dadsublen        >  0 && mumsubptr+dadsublen <= newguy->tree_starts[NUM_TREES]);
    assert(mumsubptr+dadsublen+remlen ==                             newguy->tree_starts[NUM_TREES]);

    memcpy( newguy->expr,                      mum->expr,                     mumsubptr*sizeof(node));
    memcpy(&newguy->expr[mumsubptr],          &dad->expr[dadsubptr],          dadsublen*sizeof(node));
    memcpy(&newguy->expr[mumsubptr+dadsublen],&mum->expr[mumsubptr+mumsublen],remlen   *sizeof(node));

    /*do at start of second loop to save lock/unlock
    {const int e = pthread_mutex_lock(&mutex);   assert(e==0);}
      mum->num_children--;
      if(mum->num_children==0) {delete mum; thispop->pop[Mum] = NULL;}
      dad->num_children--;
      if(dad->num_children==0) {delete dad; thispop->pop[Dad] = NULL;}
    {const int e = pthread_mutex_unlock(&mutex); assert(e==0);}
    */

    }//end do rest of crossover

    //assert(expr_table[my_id]);
    //printf("thread_fitness(): thread %d, calling Dummy(thispop->(pop[%d])\n",my_id, i);
    newguy->thread = my_id;
    //time_t start; time (&start);
    thispop->Fitness(newguy);//,expr_table[my_id]);
    /** //Dummy(thispop->pop[i]);
    time_t now; time (&now);
    char begin[80]; memcpy(begin,ctime(&start)+11,8);begin[11]='\0';
    char done[80];  memcpy(done, ctime(&now)  +11,8); done[11]='\0';
    **
    pthread_mutex_lock(&mutex);
    /** //if(my_id==0)
    if(my_id<64) active_threads &= (~0l)^(1l<<my_id);
    char buff[80];
    sprintf(buff,"%012lx %012lx",initial_threads,active_threads);
    cout << "Tested "<<i<<" by thread "<<my_id
	 <<" from "<<begin<<" to "<<done
	 <<" fitness "<<thispop->pop[i]->nfitness->fvalue
	 <<" active "<<buff<<"\n";
    **
    thispop->UpdateBest(i,thispop->pop[i]);
    pthread_mutex_unlock(&mutex);
    **/
  }
  pthread_exit(NULL);
}//end thread_fitness

//gdb print *(Chrome* (*)[4000])pop

void Pop::generation_fitness(const int start) {
  time_t now; time (&now);
  const int nthreads = params->params[pThreads];
  slot = start;
  cout<<"generation_fitness("<<start<<") "<<gencount;
  cout<<" slot "<<slot<<" pThreads "<<nthreads<<" "<<ctime(&now);

  if(nthreads > npthreads){
    cout<<"pThreads "<<nthreads<<" exceeds "<<npthreads<<endl;
    exit(1);
  }
  if(thispop==NULL) {
    thispop = this;
    for(int i=0;i<nthreads;i++) {
      ids[i] = i;
      //expr_table[i] = new evalnode[thispop->params->params[pMaxExpr]];
    }
    /* Initialize mutex (only once) and condition variable objects */
    {const int e = pthread_mutex_init(&mutex, NULL); assert(e==0);}
  }
  if(DoCross){for(int i=0;i<popsize;i++) assert(newpop[i]->birth==(gencount+1)-popsize+i);} //sanity check

  if(DoCross) { //100% crossover on later generations
    for(int s=0;s<popsize;s++) {
      if(pop[s]->num_children==0) {delete pop[s]; pop[s] = NULL;}
      
      assert(thispop->newpop[s]->expr == NULL);
      thispop->newpop[s]->expr = new node[thispop->params->params[pMaxExpr]];
  }}
  pthread_t threads[nthreads];
  /* For portability, explicitly create threads in a joinable state **
  not needed https://stackoverflow.com/questions/11806793/what-is-the-usage-of-pthread-create-joinable-in-pthread
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  */
  for(int i=0;i<nthreads;i++) {
    {const int e = pthread_create(&threads[i], NULL, thread_fitness, (void *)&ids[i]); assert(e==0);}
  }
  time (&now);
  cout<<"created "<<nthreads<<" threads. "<<ctime(&now);

  /* Wait for all threads to complete */
  for(int i=0;i<nthreads;i++) {
    {const int e = pthread_join(threads[i], NULL); assert(e==0);}
  }

  if(DoCross) { //100% crossover on later generations
  for(int s=0;s<popsize;s++) {
    //if(pop[s]) assert(pop[s]->num_children==0);
    if(pop[s]) delete pop[s];
    pop[s] = newpop[s];
  }}

  for(int i=start;i<popsize;i++) UpdateBest(i,pop[i]); //avoid 2nd use of mutex

  time (&now);
  cout<<"generation_fitness Waited and joined with "<<nthreads<<" threads. "<<ctime(&now);
  DoCross = 1;

}//end Pop::generation_fitness
#endif /*PTHREADS*/


#ifdef MULTREE
Chrome* Pop::generate(int tree) //tree to crossover on, < 0 do at random
#else
Chrome* Pop::generate()
#endif
// generate a single new chrome.  Return fitness.
// This implements a one processor, steady stage GA
// Its reproductive operators are Copy, Subtree Crossover, and Node Mutation
// It uses tournament selection to select parents
// It selects in a one dimensional local region
// It uses global tournament anti-selection to select a member for replacement

// Virtual - add your own GA here if desired

{
#ifndef STACK_LIB
#ifdef MULTREE
	if ( tree < 0 ) tree = rnd ( NUM_TREES ); //defualt
#endif //MULTREE
#endif

        gencount++;

        if (initeval)                 // still on initial evaluations?
                                                        // don't generate. Finish evaluating initial pop.
        {
          InsertMember(nexteval,pop[nexteval],TRUE);
#ifdef PTHREADS
          if(params->params[pThreads]==0)
#endif
	  {
            Fitness(pop[nexteval]);
            UpdateBest(nexteval,pop[nexteval]);
	  } //end if pthreads wait until fitness has been calculated

	  int target=nexteval;
          nexteval++;
          if (nexteval>=popsize)
                initeval=FALSE;
	  return pop[target];
          }

        int target;
//        int i,region,ts,offset;
        int prob,wheel;
        int dowhat;

        Chrome* best;
        Chrome* secondbest;
        Chrome* trial;
        Chrome* trialNoH = NULL;
                // decide what to do - Cross, Mutate, Copy
        
              prob=rnd(params->params[pCrossWt]+params->params[pCross1ptWt]+params->params[pCrossSt1ptWt]+params->params[pMuteWt]+params->params[pCopyWt]+params->params[pAnnealMuteWt]+params->params[pAnnealCrossWt]+params->params[pJoinWt]);
                wheel=params->params[pCrossWt];
                if (prob<wheel)
                        dowhat=docross;
                else if (prob<(wheel += params->params[pCross1ptWt]))
                        dowhat=do1pt;
                else if (prob<(wheel += params->params[pCrossSt1ptWt]))
                        dowhat=dost1pt;
                else if (prob<(wheel += params->params[pMuteWt]))
                        dowhat=domutate;
                else if (prob<(wheel += params->params[pCopyWt]))
                        dowhat=docopy;
                else if (prob<(wheel += params->params[pAnnealMuteWt]))
                        dowhat=doanneal;
                else if (prob<(wheel += params->params[pAnnealCrossWt]))
                        dowhat=docrossanneal;
                else
                        dowhat=dojoin;
           
                // Perform reproduction
                switch (dowhat)
                {
                case docross:
                case docrossanneal:
                case do1pt:
                case dost1pt:
                case dojoin:
		  if (dowhat!=docrossanneal) {
                        target=GetTarget();                     // Find a member to replace
		  }
		  else {
                        target=GetAnnealTarget();
		  }
                        best=selectParent(target);
                        secondbest=selectParent(target);
			{BOOL fair = FALSE;
		        if(params->params[pCrossFairWt]!=0)
			  fair = (rnd(100) < params->params[pCrossFairWt]);
#ifdef MULTREE
#ifdef DIRECTEDXOVER
        		if ((params->params[pDirCrossWt] != 0)      &&
			    (rnd(100) < params->params[pDirCrossWt])   )
				tree = best->ChooseCrossTree(secondbest);
#endif
			num_crosses[tree]++;
#ifdef PDGP
                        if(params->params[pPDGP])
			  trial=best->CrossNet(secondbest, tree);
			else
#endif /*PDGP*/
switch(dowhat) {
case docrossanneal:
case docross: 
  {
  const int S2_gpquick_seed = gpquick_seed;
  trial=best->CrossTree(secondbest, tree, fair, FALSE, FALSE);
/***********************************************************************
  if(fair && params->params[pCrossFHWt]>0) {
    const int save_gpquick_seed           = gpquick_seed;
    const int save_pCrossFHWt             = params->params[pCrossFHWt];
    const int save_CrossSubTree_t_chnglen = CrossSubTree_t_chnglen;
    const int save_CrossSubTree_calls     = CrossSubTree_calls;
    const int save_CrossFH_calls          = CrossFH_calls;
    const int save_CrossFH_choice         = CrossFH_choice;
    gpquick_seed               = S2_gpquick_seed;
    params->params[pCrossFHWt] = 0;
    trialNoH=best->CrossTree(secondbest, tree, fair, FALSE, FALSE);
    if(gpquick_seed!=save_gpquick_seed) cout<<"trialNoH changed seed\n";
    gpquick_seed               = save_gpquick_seed;
    params->params[pCrossFHWt] = save_pCrossFHWt;
    CrossSubTree_t_chnglen     = save_CrossSubTree_t_chnglen;
    CrossSubTree_calls         = save_CrossSubTree_calls;
    CrossFH_calls              = save_CrossFH_calls;
    CrossFH_choice             = save_CrossFH_choice;
  }
***********************************************************************/
  }
break;
case do1pt  :trial=best->CrossTree(secondbest, tree, FALSE, TRUE, FALSE);break;
case dost1pt:trial=best->CrossTree(secondbest, tree, FALSE, TRUE, TRUE );break;
case dojoin :trial=best->JoinTree (secondbest, tree);break;
default     :assert(1==0);
}//endswitch
#else
                        trial=best->CrossTree(secondbest);
#endif
			if((params->params[pMuteCrossWt] != 0) &&
			   (rnd(100) < params->params[pMuteCrossWt]))
			  goto domutation;
		       }
                break;
                case domutate:
		case doanneal:
		  if(dowhat==domutate) {
                        target=GetTarget();                     // Find a member to replace
		  }
		  else {
                        target=GetAnnealTarget();
		  }
                        best=selectParent(target);
                        trial = best->Copy();
domutation:
                        prob=rnd(params->params[pMuteNodeWt]+params->params[pMuteConstWt]+params->params[pMuteShrinkWt]+params->params[pMuteSubTreeWt]+params->params[pMuteFairWt]);
                        wheel=params->params[pMuteNodeWt];
                        if (prob<wheel)
//Mutate with RndOp takes about 37% of runtime in big programs and 100% Mutate
                             trial->Mutate();
                        else if (prob<(wheel += params->params[pMuteConstWt]))
                             trial->MutateC();
                        else if (prob<(wheel += params->params[pMuteShrinkWt]))
                             trial->MutateShrink();
                        else if (prob<(wheel +=params->params[pMuteSubTreeWt]))
                             trial->MutateSubTree(FALSE);
			else 
                             trial->MutateSubTree(TRUE);
                 break;

                  case docopy:
                        target=GetTarget();                     // Find a member to replace
                        best=selectParent(target);
                        trial=best->Copy();
                  }//endswitch

                      // Update the pop array
                       // fitness?
	               if(params->params[pRepeatEval] || (!trial->same(best)))
                       {
#ifdef PTHREADS
                          if(params->params[pThreads]==0)
#endif
                            Fitness(trial);
                        }
		       if(trialNoH!=NULL) {
#ifdef PTHREADS
			 assert(params->params[pThreads]==0); //not supported by PTHREADS
#endif
			 if(params->params[pRepeatEval] || 
			         (!trialNoH->same(best))) {
			   if(!trialNoH->same(trial))
			     Fitness(trialNoH);
			   else
			    trialNoH->nfitness->fvalue=trial->nfitness->fvalue;
			 }
			 if(!trialNoH->same(trial))
			   CrossFH_differ++;
			 if(trialNoH->nfitness->fvalue != 
			    trial->nfitness->fvalue)
			   CrossFH_fitdiffer++;
			 delete trialNoH;
		       }

                  if ((dowhat==doanneal) || (dowhat==docrossanneal)) {
#ifdef PTHREADS
		      assert(params->params[pThreads]==0); //not supported by PTHREADS
#endif
                                // test whether mutated chome is fit enough
		      if(!accept(trial->nfitness,pop[target]->nfitness)) {
#ifdef CROSSFILE
		       if(crossfile) {
			 crossfile<<"!";
			 crossfile<<trial->nfitness->fvalue<<" "<<flush;
		       }
#endif
		       delete trial;
                       trial = pop[target]->Copy(); //restore original chrome
		      }//endif not fit enough
		  }//endif anneal

#ifdef CROSSFILE
#ifdef PTHREADS
		       assert(params->params[pThreads]==0); //not supported by PTHREADS
#endif
		       if(crossfile) {
			 crossfile<<"off ";
			 crossfile<<gencount<<" ";
			 crossfile<<trial->nfitness->fvalue<<endl;
		       }
#endif
                        InsertMember(target,trial);
#ifdef PTHREADS
                      if(params->params[pThreads]==0) //wait until fitness is calculated
#endif
                        UpdateBest(  target,trial);
        return trial;
}

#ifdef MULTREE
Chrome* Pop::go_until(time_t max_time, int maxevals, float maxfitness,
	              int tree) //default random tree
#else
Chrome* Pop::go_until(time_t max_time, int maxevals, float maxfitness)
#endif
// generate until time max_time, evals maxevals, or fitness maxfitness
// Do this to run the GA in a bigger event loop

{
        //int done = FALSE;
        int didevals = 0;
#ifdef PTHREADS
	if(params->params[pThreads]) {
	  const int allwheel = params->params[pCrossWt]+params->params[pCross1ptWt]+params->params[pCrossSt1ptWt]+params->params[pMuteWt]+params->params[pCopyWt]+params->params[pAnnealMuteWt]+params->params[pAnnealCrossWt]+params->params[pJoinWt];
	  if(params->params[pCrossWt] != allwheel) {
	    cout<<"ERROR pThreads: "<<params->params[pThreads]<<" only supports crossover "
		<<"pCrossWt: "<<params->params[pCrossWt]<<" of "<<allwheel<<endl;
	    exit(1);
	  }
	  if(params->params[pGenerational]==0 ) {
	    cout<<"ERROR pThreads: "<<params->params[pThreads]<<" needs pGenerational"<<endl;
	    exit(1);
	  }
	  if(params->params[pRepeatEval] == 0) {
	    cout<<"ERROR pThreads: "<<params->params[pThreads]<<" needs pRepeatEval"<<endl;
	    exit(1);
	  }
	}
	int start = -1;
#endif
#ifndef GENERATIONAL
#ifdef MULTREE
        Chrome* lastchrome=generate(tree);
#else
        Chrome* lastchrome=generate();
#endif //MULTREE
        didevals++;
#else
        Chrome* lastchrome=NULL;
#endif //GENERATIONAL

        //while(maxevals == 0 || didevals<maxevals)
        while ((max_time == 0 || time(NULL) < max_time) &&
        (maxevals == 0 || didevals<maxevals) &&
        ( lastchrome==NULL || (lastchrome->nfitness->fvalue < maxfitness)
#ifdef GENERATIONAL 
	 || (params->params[pGenerational] && ((gencount+1)%popsize!=0))
#endif
	)
        &&!Aborted())
        {
#ifdef GENERATIONAL 
	       if((gencount+1)%popsize==0 && 
		  (!initeval)             && params->params[pGenerational]) {
#ifdef PARETO
		 clear_rank();
#else
		 BestFitness = NULL;
		 if(params->params[pElitist] != 0) {
		   cout<<"Copying BestMember of pop ("<<BestMember;
		   cout<<") to newpop 0 "<<gencount<<endl;
		   gencount++;
		   didevals++;
		   Chrome* newchrome = pop[BestMember]->Copy();
		   InsertMember(0,newchrome);
		   UpdateBest(  0,newchrome);
		 }
#endif
     {time_t now; time (&now); cout <<"go_until end if new gen "<< ctime(&now);}
	       }
#endif

#ifdef MULTREE
               lastchrome = generate(tree);
#else
               lastchrome = generate();
#endif
#ifdef PTHREADS
               if(start<0) start = gencount;
#endif
               didevals++;
#ifdef GENERATIONAL
	       if((gencount+1)%popsize==0) {
	       if(gencount > popsize  && params->params[pGenerational] ) {
		 time_t now; time (&now);
		 cout<<"Copying newpop to pop, at "<<gencount<<" "<< ctime(&now);
#ifdef PTHREADS
		 if(params->params[pThreads]==0) //wait until crossovers done
#endif
		 for(int s=0; s<popsize; s++) {
		   //if(oldpop[s]) delete oldpop[s];
		   //oldpop[s] = pop[s]; //save for children stats only
		   delete pop[s];
		   pop[s] = newpop[s];
		 }
	}
#ifdef PTHREADS
	       if(params->params[pThreads] && start>=0)
		 generation_fitness(start%popsize);
#endif
	       }
#endif //GENERATIONAL
        }
        return lastchrome;
}

#ifndef PDGP
#ifdef MULTREE
void Pop::reinitialise(BOOL change_flags[NUM_TREES])
{
#ifdef PARETO
clear_rank();
#else
BestFitness=NULL;
#endif

for ( UINT i = 0; i < popsize; i++ ) {
	pop[i]->reinit_trees(change_flags);
	Fitness(pop[i]);
	InsertMember(i,pop[i],TRUE);
	UpdateBest(  i,pop[i]);
	}
}//end reinitialise
#endif
#endif /*PDGP*/

#ifdef STACK_LIB
int Pop::GetTarget()
// Tournament anti-selection for replacement.  Usually size 1 (random selection) or 2
// Will select a target that is too old, even if more fit.
{
        int target = rnd(popsize);
        int i,winner;
        if (params->params[pKillTourn]>1 && 
	    (params->params[pMaxAge] == 0 ||
	     (gencount - pop[target]->birth) <= params->params[pMaxAge]))
        // pick a target to replace
        for (i=1;i<params->params[pKillTourn];i++)
        {
                winner=rnd(popsize);
                if (params->params[pMaxAge] != 0 &&
		    (gencount - pop[winner]->birth) > params->params[pMaxAge])
                	return winner;

		if (pop[target]->nfitness->IsBetter(pop[winner]->nfitness))
                        target = winner;
        }
#ifdef CROSSFILE
	if(crossfile) {crossfile<<"T "<<target<<" "<<flush;}
#endif
        return target;
}
#else
int Pop::GetTarget()
// Tournament anti-selection for replacement.  Usually size 1 (random selection) or 2
{
#ifdef GENERATIONAL
  if(params->params[pGenerational]) return gencount % popsize;
#endif
//	Chrome* candidates [params->params[pKillTourn]];
//	int     slot       [params->params[pKillTourn]];
	Chrome* candidates [10];
	int     slot       [10];

        // pick a target to replace

        int s=rnd(popsize); //chose any deme in pop
#ifdef QUEUE_LIB
	slot[0] = s;
	candidates[0] = pop[s];

	select_candidates(s, params->params[pKillTourn] - 1,
			  &candidates[1], &slot[1], FALSE ); //chose from same deme
#else /*QUEUE_LIB*/
	select_candidates(s, params->params[pKillTourn],  //chose from same
			  candidates, slot, FALSE );      //deme
#endif /*QUEUE_LIB*/
	int index;
	problem->Worstof(candidates, params->params[pKillTourn], gencount, 
                         &index, s);
#ifdef CROSSFILE
	if(crossfile) {crossfile<<"T "<<slot[index]<<flush;}
#endif
        return slot[index];
}
#endif

double Pop::Temperature() const
{
const int generate = params->params[pGenerational]?
                         (gencount/popsize)*popsize :
                         gencount;

const double fract_done = double(generate)/double(params->params[pGenerateLimit]);

if (params->params_dbl[pStartTemp]==params->params_dbl[pEndTemp])
  return params->params_dbl[pStartTemp];
else if (params->params_dbl[pEndTemp]==0)
  return params->params_dbl[pStartTemp]*(1-fract_done);
else {
  const double ratio = log(params->params_dbl[pEndTemp]/
			   params->params_dbl[pStartTemp]);
  return params->params_dbl[pStartTemp]*exp(fract_done*ratio);
}

}//end Temperature

BOOL Pop::accept(const FitnessValue* newfv, const FitnessValue* oldfv) const
{
double temperature = Temperature();

// cout<<"\ntemperature "<<temperature<<" at "<<gencount<<", "
//     <<newfv->fvalue<<" <- "<<oldfv->fvalue<<endl;

#ifdef PARETO
assert (0==1); //cant do Simulated Annealing with multi-objective ffn
#endif /*PARETO*/

return ( newfv->fvalue > oldfv->fvalue ||
	(temperature==0 && newfv->fvalue == oldfv->fvalue) ||
	(temperature >0 && 
	 exp((newfv->fvalue - oldfv->fvalue)/temperature) > rand_0to1() ));

}//end accept

int Pop::GetAnnealTarget()
{
#ifdef GENERATIONAL
  if(params->params[pGenerational]) return gencount % popsize;
#endif
  int best;
  selectParent(rnd(popsize), &best);
  return best;
}

void Pop::Write(ostream & fout, BOOL binary) //WBL
{
for ( UINT i = 0; i < popsize; i++ )
    if(!binary)
      { fout << endl << i ;
#ifdef TRACE
	pop[i]->write_trace(); //fix some time. Problem with const
#endif //TRACE
	pop[i]->write(PRETTY_NONE, fout);
	fout << endl;
      }
    else { //psuedo binary. 1 byte per opcode but in ascii
#ifdef PDGP
        if(params->params[pPDGP]) assert(0==1);//not implemented yet
#endif /*PDGP*/
#ifdef MULTREE
	for (int t = 0; t < NUM_TREES; t++ ) {
		int start = pop[i]->tree_starts[t];
		int len   = pop[i]->SubLen(pop[i]->tree_starts[t]);
		assert(len == pop[i]->tree_starts[t+1]);
#else
		int start = 0;
		int len   = pop[i]->SubLen(0);
#endif
		for(int x = 0; x<len; x++)
			fout<<char(pop[i]->expr[start+x].op+'A');
		fout<<endl;
#ifdef MULTREE
		assert(pop[i]->tree_starts[t+1]-pop[i]->tree_starts[t]==len);
	   }
#endif
      }
}//end Pop::Write

//************************************************************************
//////////////////////ChromeParams methods

ChromeParams::ChromeParams()
//Set default parameters here in the constructor
{
        varcount = 0;
        funccount = 0;

	memset(params,0,sizeof(params));
	memset(params_dbl,0,sizeof(params_dbl));
	memset(params_str,0,sizeof(params_str));
	memset(pnames,0,sizeof(pnames));
#ifdef PDGP
	memset(pdgpwidth,0,sizeof(pdgpwidth));
	memset(pdgpindex,0,sizeof(pdgpindex));
	memset(net_limits,0,sizeof(net_limits));
	pdgpheight = 0;
#endif /*PDGP*/

      params[pMaxExpr] = 50;         // Maximum expression size in nodes
      params[pMaxTreeExpr] = 50;     // Soft initial Maximum tree size in nodes
      params[pMaxHaltExpr] = 0;      // Soft limit on tree size in nodes
      params[pUniInitWt] = 0;         // Uniform (or ramped half) (out of 100)
      params[pInitSparse] = 0;        // ramped half or sparse (cant mix them)
      params[pInitExpr] = 6;          // Maximum initial expression depth
      params[pMinInitExpr] = 2;       // Minimum initial expression depth
      params[pMinInitSizeExpr] = 3;   // Minimum initial expression size
      params[pInitSizeExpr] = 64;     // Maximum initial expression size
      params[pMaxDepth] = 0;          // Maximum expression depth
      params[pMuteRate] = 100;        // Node Mutation rate per 1000
      params[pCrossSelf] = 1;         // enable cross with self
      params[pUnRestrictWt] = 70;        // any point crossover out of 100
      params[pCrossWt] = 100;       // Crossover weight
      params[pCross1ptWt] = 0;      // One point Crossover weight
      params[pCrossSt1ptWt] = 0;    // Strict one point Crossover weight
      params[pCrossFairWt] = 0;     // Fair Crossover weight (out of 100)
      params[pCrossFairLimit] = 0;  // Fair Crossover insert size limit
                                    // if <0 no limit on Fair XO insert subtree
      params[pCrossFHWt] = 0;       // Fair Homologous weight (out of 100)
      params[pJoinWt] = 0;          // Join two trees with new binary node
      params[pMuteCrossWt] = 0;     // Crossover folowed by mutate (out of 100)
      params[pMuteWt] = 30;        // overall Mutation weight
      params[pAnnealCrossWt] = 0;   // Crossover Annealing Weight
      params[pAnnealMuteWt] = 0;   // Mutate Annealing weight
      params[pCopyWt] = 10;        // Copy weight
      params[pDirCrossWt] = 0;        // Directed choice of tree (given cross)
      params[pMuteNodeWt] = 100;   // Node Mutate weight
      params[pMuteConstWt] = 100;    // MutateC weight
      params[pMuteShrinkWt] = 100;   // MutateShrink weight
      params[pMuteSubTreeWt] = 100;  // MutateSubTree weight
//    params[pMuteFairWt];           // MutateSubTree weight (equal heights)
      params[pMuteFairLimit] = 0;    // Size of rand_tree tables, 0 => none
	                             // -ve +-50%, +ve fairmute size limit
      params[pSeeds] = 0;	     // Number seeds to load into init pop
      params[pMuteSeedWt] = 100;     // Chance of mutating seed out of 100
      params[pSelectMethod] = ETournament;  // Selection method = Tournament
      params[pGenerational] = 0;      // Steady state or Generational
      params[pThreads] = 0;           // If generational, how many pthreads to use
      params[pTournSize] = 7;         // Tournament size
      params[pTournComp] = 0;         // Used to separate ties if have pareto
      params[pTournOff] = 0;         // Turn off selection
      params[pNicheSize] = 0;         // Used to spread population if pareto
      params[pElitist] = 0;           // <>0 keep best. Only for pareto so far
      params[pMateRadius] = 500;       // Mating radius
      params[pGaussRegion] = 0;         // Guassian mate selection 0 disable 1 enable
      params[pRepeatEval] = 0;         // Repeat eval on copy? 0 no 1 yes
      params[pKillTourn] = 2;         // Size of "Kill Tournament" for replacement
      params[pMaxAge] = 2000;      // Maximum age before replaced even if better
      params[pParsimony] = 0;         // Parsimony factor
      params[pFitnessCases] = 20;        // # fitness cases
      params[pStoreFit] = 0;          //Dont test fitness of unchanged code
      params[pCPU] = 0;
      params[pNeighbours] = 20000;    // number of neighbours analysed by ALL
      params[pSAANWt] = 100;

      params_dbl[pStartTemp] = 0;     //Accept same or better fitness
      params_dbl[pEndTemp] = 0;       //Linear temp change
      params_dbl[pSamePen] = -1;      // -1 >= subtract 0.5 else fitness*=x
      params_dbl[pAllMinFit] = 24;    // min fitness to be analysed by ALL

// set pnames. Could be done as preset but this avoid order change problems
      strcpy(pnames[pPopSize],		"pPopSize" );
      strcpy(pnames[pGenerateLimit],	"pGenerateLimit" );
      strcpy(pnames[pPopSeed],		"pPopSeed" );
      strcpy(pnames[pTestSeed],		"pTestSeed" );
      strcpy(pnames[pTrace],	 	"pTrace" );
      strcpy(pnames[pMaxExpr],	 	"pMaxExpr" );
      strcpy(pnames[pMaxTreeExpr],	"pMaxTreeExpr" );
      strcpy(pnames[pMaxHaltExpr],	"pMaxHaltExpr" );
      strcpy(pnames[pUniInitWt],	"pUniInitWt" );
      strcpy(pnames[pInitSparse],	"pInitSparse" );
      strcpy(pnames[pInitExpr],	 	"pInitExpr" );
      strcpy(pnames[pMinInitExpr],	"pMinInitExpr" );
      strcpy(pnames[pMinInitSizeExpr],	"pMinInitSizeExpr" );
      strcpy(pnames[pInitSizeExpr],	"pInitSizeExpr" );
      strcpy(pnames[pMaxDepth],	 	"pMaxDepth" );
      strcpy(pnames[pMuteRate],	 	"pMuteRate" );
      strcpy(pnames[pCrossSelf],	"pCrossSelf" );
      strcpy(pnames[pUnRestrictWt],	"pUnRestrictWt" );
      strcpy(pnames[pCrossWt],	 	"pCrossWt" );
      strcpy(pnames[pCross1ptWt],	"pCross1ptWt" );
      strcpy(pnames[pCrossSt1ptWt],	"pCrossSt1ptWt" );
      strcpy(pnames[pCrossFairWt],	"pCrossFairWt" );
      strcpy(pnames[pCrossFairLimit],	"pCrossFairLimit" );
      strcpy(pnames[pCrossFHWt],	"pCrossFHWt" );
      strcpy(pnames[pJoinWt],	        "pJoinWt" );
      strcpy(pnames[pMuteCrossWt],	"pMuteCrossWt" );
      strcpy(pnames[pMuteWt],	 	"pMuteWt" );
      strcpy(pnames[pMuteNodeWt],	"pMuteNodeWt" );
      strcpy(pnames[pMuteConstWt],	"pMuteConstWt" );
      strcpy(pnames[pMuteShrinkWt],	"pMuteShrinkWt" );
      strcpy(pnames[pMuteSubTreeWt],	"pMuteSubTreeWt" );
      strcpy(pnames[pMuteFairWt],	"pMuteFairWt" );
      strcpy(pnames[pMuteFairLimit],	"pMuteFairLimit" );
      strcpy(pnames[pAnnealMuteWt],	"pAnnealMuteWt" );
      strcpy(pnames[pAnnealCrossWt],	"pAnnealCrossWt" );
      strcpy(pnames[pSeeds],		"pSeeds" );
      strcpy(pnames[pMuteSeedWt],	"pMuteSeedWt" );
      strcpy(pnames[pCopyWt],	 	"pCopyWt" );
      strcpy(pnames[pDirCrossWt],	"pDirCrossWt" );
      strcpy(pnames[pSelectMethod],	"pSelectMethod" );
      strcpy(pnames[pGenerational],	"pGenerational" );
      strcpy(pnames[pThreads],		"pThreads" );
      strcpy(pnames[pTournSize],	"pTournSize" );
      strcpy(pnames[pTournComp],	"pTournComp" );
      strcpy(pnames[pTournOff],         "pTournOff" );
      strcpy(pnames[pNicheSize],	"pNicheSize" );
      strcpy(pnames[pElitist],		"pElitist" );
      strcpy(pnames[pPopWidth],	 	"pPopWidth" );
      strcpy(pnames[pDemeWidth],	"pDemeWidth" );
      strcpy(pnames[pMateRadius],	"pMateRadius" );
      strcpy(pnames[pGaussRegion],	"pGaussRegion" );
      strcpy(pnames[pRepeatEval],	"pRepeatEval" );
      strcpy(pnames[pKillTourn],	"pKillTourn" );
      strcpy(pnames[pMaxAge],	 	"pMaxAge" );
      strcpy(pnames[pParsimony],	"pParsimony" );
      strcpy(pnames[pFitnessCases],	"pFitnessCases" );
      strcpy(pnames[pStoreFit],		"pStoreFit" );
      strcpy(pnames[pCPU],		"pCPU" );
      strcpy(pnames[pNeighbours],	"pNeighbours" );
      strcpy(pnames[pPDGPinit],      	"pPDGPinit" );
      strcpy(pnames[pSAANWt],		"pSAANWt" );
      strcpy(pnames[pSIANWt],		"pSIANWt" );
      strcpy(pnames[pSAINWt],		"pSAINWt" );
      strcpy(pnames[pSIINWt],		"pSIINWt" );
      strcpy(pnames[pSSAANWt],		"pSSAANWt" );
      strcpy(pnames[pSSIANWt],		"pSSIANWt" );
      strcpy(pnames[pSSAINWt],		"pSSAINWt" );
      strcpy(pnames[pSSIINWt],		"pSSIINWt" );
      strcpy(pnames[pPDGP],      	"pPDGPwidth" );

      strcpy(pnames[pStartTemp+PARAM_DBL_OFF], "pStartTemp" );
      strcpy(pnames[pEndTemp  +PARAM_DBL_OFF], "pEndTemp" );
      strcpy(pnames[pSamePen  +PARAM_DBL_OFF], "pSamePen" );
      strcpy(pnames[pAllMinFit+PARAM_DBL_OFF], "pAllMinFit" );

      strcpy(pnames[pDumpfile1+PARAM_STR_OFF], "pDumpfile1" );
      strcpy(pnames[pDumpfile2+PARAM_STR_OFF], "pDumpfile2" );
      strcpy(pnames[pZipCmd   +PARAM_STR_OFF], "pZipCmd" );
      strcpy(pnames[pPopfile  +PARAM_STR_OFF], "pPopfile" );
}


ChromeParams::~ChromeParams()
{
}

void ChromeParams::Load(istream& in )
{
const	int len = 22+param_str_size;//maximum line length+1 of gp.ini
	char buff[len];

	for(; in.getline(buff,len), strlen(buff)>0; ) {
//		cout<<"getline returned "<<len<<" ,"<<flush; //debug
//		cout<<buff<<endl;
		Load(buff);
	}
}//end ChromeParams::Load()


void ChromeParams::Load(char* cs)
{
        int x;
        for (x=0;x<PARAM_TOTAL;x++) {
		char* s = strstr(cs,&pnames[x][0]);
              	if ( s != NULL ) {
//			cout<<"'"<<s<<"' params["<<x<<"]= "<<params[x]<<flush;
			char* s1 = strstr(s,":");
              		if ( s1 != NULL ) {
			        if(x<PARAM_COUNT)
#ifdef PDGP
				  if(x==pPDGP) {
				    loadpdgp(++s1);
				    params[x] = TRUE;
				  }
				  else
#endif /*PDGP*/
					sscanf(++s1,"%d",&params[x]);
			        else if(x<PARAM_STR_OFF) 
					sscanf(++s1,"%le",
					        &params_dbl[x-PARAM_DBL_OFF]);
				else
					sscanf(++s1,"%s",
						&params_str[x-PARAM_STR_OFF][0]);
//				cout<< "'"<<s1<<"' params["<<x<<"]= "
//				    <<params[x]<<endl;
			}
		}
	}
}//end ChromeParams::Load()


void ChromeParams::Save(ostream& out)
{
	char buf [21+param_str_size];
	int x;

        for (x=0;x<PARAM_COUNT;x++) {
		sprintf(buf, "%-15s:", &pnames[x][0]);
		out<<buf;
#ifdef PDGP
		if(x==pPDGP)
		  savepdgp(out);
		else 
#endif /*PDGP*/
		  out<<" "<<params[x]<<endl;
	}

        for (x=0;x<PARAM_DBL_COUNT;x++) {
		sprintf(buf, "%-15s: %g", &pnames[x+PARAM_DBL_OFF][0],
			params_dbl[x]);
		out<<buf<<endl;
	}

        for (x=0;x<PARAM_STR_COUNT;x++) {
		sprintf(buf, "%-15s:%s", &pnames[x+PARAM_STR_OFF][0],
			params_str[x]);
		out<<buf<<endl;
	}

}//end ChromeParams::Save()
