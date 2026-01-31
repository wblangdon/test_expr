// Chrome.h   version to support stack.cc, from GPQUICK-2.1
// W.Langdon cs.ucl.ac.uk $Revision: 1.172 $

//Modifications (reverse order):
// WBL  1 Jan 2019  Speed crossover by saving length in tree_starts
// WBL 15 Dec 2018  For PTHREADS.
// WBL 23 Nov 2018  for avx split retval into reteval and retval
// WBL 27 Aug 2016  Increase EXPRLEN for edge of bloat search
// WBL 26 Aug 2016  experiment with reporting size and maxdepth
// WBL 12 Jul 2016  Add exon
// WBL  5 Jul 2016  Add DisplayZStats()
// WBL 24 Jun 2016  Add DisplaySolStats()
// WBL 23 Jun 2016  Add stats::nmax
// WBL  6 Jun 2016  r1.125 for ~/gp/eden_mux add INTRON64 (cf r1.147)
// WBL 28 Oct 1999  Add displaydepthlimit to Write and SubWrite
// WBL  4 Oct 1999  Add pTraceLPStats

// WBL 29 Sep 1999  Add JoinTree and DepthMax (0 and 1 arguments)

// WBL 30 Apr 1999  Conversion for SiliconGraphics MIPS4
//                  Make FASTCONST independent of FASTTRAVERSE

// WBL 27 Apr 1999  Allow chromes to exceed 32767

// WBL  8 Aug 1998  Move costats from stats.cc

// WBL  7 Jul 1998  Implement FASTTRAVERSE, add arityop

// WBL 12 Jun 1998  Add gtree* optional argument to rand_tree

// WBL 24 Mar 1998  Add pMuteCrossWt, arityc

// WBL 22 Mar 1998  Support 1point and strict onepoint crossover

// WBL  7 Jan 1998  move class stats from stats.cc into this file

// WBL 31 Dec 1997  add classes ghashtable and ghashchain, cf stats.cc

// WBL  7 Sep 1997  Support pTraceAllTests

// WBL  1 Sep 1997  Support display of CrossSubTree_t_chnglen

// WBL 31 Aug 1997  Add pCrossFairWt and pCrossFairLimit

// WBL 12 Aug 1997  Add pMuteFairLimit

// WBL  9 Aug 1997  Add support for RAND_TREE

// WBL  1 Aug 1997  Add pMaxHaltExpr and length() and Temperature()

// WBL 30 Jul 1997  Add support for floating point parameters

// WBL 28 Jul 1997  Add support for Simulated Annealing and Hill Climbing

// WBL 30 Apr 1997  Add first level of support for PDGP (chrome only)

// WBL 13 Feb 1997  Add pTournOff

// WBL 17 Dec 1996  Add pGenerational, GENERATIONAL and newpop

// WBL 16 Dec 1996  Make GatherPstats support tree depth info

// WBL 14 Dec 1996  Add pMaxDepth

// WBL 14 Jun 1996  Boost EXPRLEN from 1000 to 4500 (for 42 tree S.Wales)

// WBL  8 Apr 1996  Add MutateSubTree

// WBL 15 Mar 1996  Add hashcode, Chrome::same and DisplayDupStats

// WBL  9 Mar 1996  Add FitnessValue::Load and Save for GEN_MEM

// WBL 28 Feb 1996  Add Depth and crossfile

// WBL 23 Oct 1995  Remove NUM_TREES, PARETO, TRACE, MULTREE, retval to prob.h
//                  And BOOL to pch.h

// WBL 20 Oct 1995  Increase NUM_TREES to 5 (for 2 more adfs)

// WBL 15 Oct 1995  Restore PARETO

// WBL 22 Sep 1995  Remove PARETO and set NUM_TREES to 3

// WBL 30 Jun 1995  Add ability to dump and restore whole population

// WBL 20 Jun 1995  Add ProbLTStat

// WBL  1 May 1995  Add DisplayLPstats and GatherPstats

// WBL  6 May 1995  Pack Chrome tighter, remove fields which are the same inall

// WBL  4 May 1995  Add pStoreFit.
//                  Add FitnessValue::update_rank, remove unused arg of Copy
//                  cf 23 March 1995

// WBL  2 May 1995  Add pDirCrossWt

// WBL  1 May 1995  Add DisplayPstats

// WBL 26 Apr 1995  Reduce size of population by removing idx from node
//                  and evalnode (in the latter boost op from unsigned char to
//                  int, as it takes no more space and might be faster)

// WBL 18 Apr 1995  Decrement NUM_TREES (for adf2 tree removal)

// WBL 17 Apr 1995  Add optional last_tree argument to chrome::write
//                  Add pMaxTreeExpr

// WBL  8 Apr 1995  Increase NUM_TREES to 16 (for 4 local adfs) etc

// WBL  7 Apr 1995  Increase NUM_TREES to 12 (for list and two adfs)

// WBL 23 Mar 1995  Move pElitist tests into chrome.cxx

// WBL 23 Mar 1995  Add Optional Chrome* parameters to GetFitnessObj and Copy

// WBL 15 Mar 1995  Add pElitist

// WBL 12 Mar 1995  add pTournComp and pNicheSize

// WBL 26 Feb 1995  Update definition of Bestof and Worstof

// WBL 19 Feb 1995  Add def PARETO and use it to remove BestMember etc

// WBL 13 Feb 1995  Add pTracePareto

// WBL 31 Jan 1995  make strcmpi accessible by test.cc

// WBL 25 Jan 1995  Increase NUM_TREES to 10 (for list)

// WBL 29 Nov 1994  Implment pnames and re-order chromeParams, add Load(char*)

// WBL 15 Oct 1994  Add support for loading populations

// WBL  8 Sep 1994  Add static_check

// WBL 25 Aug 1994  Make CSelector different for each tree, remove funcbag
//                  from chrome (use problem pointer instead)

// WBL 24 Aug 1994  Add Bestof and Worstof to Problem. Add select_candidates

// WBL 22 Aug 1994  remove (to problem files) hits and write_hits()
//                  add Problem::NewChrome

// WBL 15 Aug 1994  Add hits and write_hits()
//                  Initialise xtree to -1 rather than 0

// WBL  6 Jul 1994  Add pPopWidth and pDemeWidth

// WBL  6 Jul 1994  Put NUM_TREES to 6 (for adf1)

// WBL 30 Jun 1994  Add pTraceStatic+pTraceDynamic+pTraceDStats

// WBL 28 Jun 1994  Supply default for go_until tree param

// WBL 21 Jun 1994  Put back NUM_TREES to 5. Add parameter to generate()
//                  go_until() and CrossTree()

// WBL 21 Jun 1994  Increase NUM_TREES from 5 to 7 (for two adfs)

// WBL 16 Jun 1994  Add PopSeed and pTestSeed, pTrace. Add TRACE

// WBL 14 Jun 1994  Add Save and Load to ChromeParams (based on GPCPlus3.0
//                  gpv.cc). Add pGenerateLimit
//                  Add DisplayDStats() and Write()

// WBL 16 May 1994  Add DisplayStats() to Pop, Fix MULTREE bug

// GPQUICK
// C++ prefix stack implementation
// Divides GP system into classes:
//    Chrome-Function  - representation of the genetic program
//    Pop - Runs the GA
//    Problem - fitness and primitive functions


#ifndef _CHROME_H
#define _CHROME_H

#include "selector.h"
#include "prob.h"
#include <assert.h>

// Support for multiple trees within one individual in pop WBL
// Only tested with FASTEVAL so far
#define MULTREE

// keep and print genelogical trace info
//TRACE depends on PTHREADS
#define TRACE

// Define this if you do multiple evals, one Chrome at a time
// It speeds up the EVAL calls significantly
// Penalties are:  Time to expand the chrome once before evaluations
//                 Uses globals, so only one Chrome at a time
#define FASTEVAL
#ifndef SMC_GP
#ifndef POLY_LIB
#define FASTTRAVERSE
#endif /*POLY_LIB no ifs in poly*/
#endif

	class Problem;
	class Chrome;
        typedef Chrome* PtrChrome;
        class FitnessValue;
#ifdef RAND_TREE
        class gtree;
        class rand_tree_tables;
extern  gtree* fulltree(const int depth, const int arity, int& size);
#endif

#ifdef PTHREADS
typedef struct node* nodeptr;
//typedef struct evalnode* evalnodeptr;
#ifdef FASTEVAL
		// evaluation code with global pointers
//Jan 2019 pass IP and stack via arguments rather than globals
	typedef reteval (*EVALFUNC)(node** ip, retval** sp);
#else
		// evaluation code with pointer passing
	typedef reteval (*EVALFUNC)(Chrome*,const int my_id);
#endif
#else
#ifdef FASTEVAL
		// evaluation code with global pointers
	typedef reteval (*EVALFUNC)();
#else
		// evaluation code with pointer passing
	typedef reteval (*EVALFUNC)(Chrome*);
#endif
#endif /*PTHREADS*/


#ifdef PDGP
	// One byte instruction, 8 bits links
typedef struct {unsigned char op; unsigned char lk[MAX_ARITY];} node;
#else
	// One byte instruction node, 8 bits function index
typedef struct node {unsigned char op;} node;  // node type
#endif /*PDGP*/
	// eval node with pointers for direct function calls
typedef struct evalnode {
          EVALFUNC ef; 
#ifdef FASTTRAVERSE
          evalnode* jump; 
#else
          int op;
#endif
#ifdef FASTCONST
          retval value;
#endif
        } evalnode; // node type

	// Grab the function index
#define FUNCNUM(c) (c.op)
	// argument count for the current function
#define ARGNUM() (funclist[FUNCNUM(expr[ip])]->argnum)
#define PARGNUM(ip) (funclist[ip->op]->argnum)
	// Grab the operand
#define VARNUM(c) (0)
#define PVARNUM(ip) (0)
	// Build a node "c"
#define SETNODE(c,o,v) c.op=o
#ifdef FASTTRAVERSE //--------------------------------------------------
#ifdef FASTCONST
#define SETEVALNODE(ip,o,v) ip->ef=funclist[o]->evalfunc;\
                            ip->value=v
#else
#define SETEVALNODE(ip,o,v) ip->ef=funclist[o]->evalfunc
#endif /*FASTCONST*/
#define SETEVALJUMP(ip,j)   ip->jump=j
#else /*FASTTRAVERSE*/
#ifdef FASTCONST
#define SETEVALNODE(ip,o,v) ip->op=o;ip->ef=funclist[o]->evalfunc;\
                            ip->value=v
#else
#define SETEVALNODE(ip,o,v) ip->op=o;ip->ef=funclist[o]->evalfunc
#endif /*FASTCONST*/
#endif /*FASTTRAVERSE --------------------------------------------------*/

	// Function evaluation stuff
	// These macros may change internal form for faster evaluation.
	// Arguments will be removed.  Use them.
#ifdef FASTEVAL

		//Define an EVALFUNC with no arguments
#define OPDEF(Op) reteval Op()
		//Get a pointer to the chrome being evaluated
#define CHROMEP ChromeGlobal
		//current instruction pointer
#define IP IpGlobal
#ifdef FASTCONST
		// get its argument
#define GETVAL IP->value
#else
#define GETVAL (0)
#endif /*FASTCONST*/
		// traverse an unused argument in an eval
#ifdef FASTTRAVERSE
#define TRAVERSE() IP=IP->jump
#else
#define TRAVERSE() CHROMEP->TraverseGlobal()
extern evalnode* ExprGlobal; //for INTRON64
#endif /*FASTTRAVERSE*/
		//Evaluate the next expression
#define EVAL ((++IP)->ef)()
#else
		//Define an EVALFUNC
#define OPDEF(Op) reteval Op(Chrome* curchrome)
		//Get a pointer to the chrome being evaluated
#define CHROMEP curchrome
		//current instruction pointer
#define IP CHROMEP->ip
		// get its argument
#define GETIDX (CHROMEP->expr[IP].idx)
		// traverse an unused argument in an eval
#define TRAVERSE() CHROMEP->Traverse()
		//Evaluate the next expression
#define EVAL CHROMEP->eval()
#endif

		//function and memory arrays FUNCARRAYSIZE <= 256
#define FUNCARRAYSIZE 256
#define CONSTARRAYSIZE 256

// cut numeric overflows.  Bound returns to 10^15
#define BIGFLOAT ( (retval) 1.0e15 )
#define SMALLFLOAT ( (retval) 1.0e-15 )

#define BOUNDF(f) (f==0? f : (f>0 ?((f)>BIGFLOAT ? BIGFLOAT : ((f)<SMALLFLOAT? SMALLFLOAT : (f))) : ((f)<-BIGFLOAT ? -BIGFLOAT : ((f)>-SMALLFLOAT? -SMALLFLOAT : (f))) ))


// Compatibility stuff
#define UINT unsigned int


extern ofstream& crossfile;

// All primitives in a problem are subclasses of this FUNCTION object

class Function {
public:
	int serial;             // serial number in universal function list (not implemented)
	char name[30];                  // Function name
	int argnum;             // number of arguments
	int varnum;                             // number of variables in variable table of opcode
	int weight;             // selection frequency relative to other functions
	Function() {};
	Function(const int a,const int v,const int w,const EVALFUNC e,const char* n)
		{argnum=a;varnum=v;weight=w;evalfunc=e;strcpy(name,n);}; 
	virtual char* getprint(Chrome* st);       // Printed name may differ
	virtual char* getprint(Chrome* st, const int ip, char* scratch);
	char* getname() {return name;};
	EVALFUNC evalfunc;                      // pointer to evaluation code.  Not a virtual for speed reasons
#ifndef FASTEVAL
	reteval eval(Chrome* st) {return (evalfunc)(st);};       // active ingredient
#else
	//reteval eval() {return (evalfunc)();};
#endif
};

//******************* Parameters
enum {
pPopSize,         // population size
pGenerateLimit,         // Terminate GP if no solution
pPopSeed,         // Seed for Pop etc.
pTestSeed,         // Seed for application to generate tests
pTrace,         // display trace info
pMaxExpr,          // Maximum expression size in nodes
pMaxTreeExpr,      // Soft initial maximum tree size in nodes
pMaxHaltExpr,      // Stop run if continously exceeded
pUniInitWt,        // Uniform initial expression (v 1/2 and 1/2) (out of 100)
pInitSparse,       // half and half or sparse (0 or 1 only!)
pInitExpr,           // Maximum initial expression depth
pMinInitExpr,        // Minimum initial expression depth
pInitSizeExpr,       // Maximum initial expression size (if Uniform initial)
pMinInitSizeExpr,    // Minimum initial expression size (if Uniform initial)
pMaxDepth,           // Maximum expression depth
pMuteRate,          // Node mutation rate per 1000
pCrossSelf,         // Allow self crossover?
pUnRestrictWt,         // any point crossover per 100
pCrossWt,         // Crossover weight on generate
pCross1ptWt,      // One point Crossover weight on generate
pCrossSt1ptWt,    // Strict one point Crossover weight on generate
pCrossFairWt,     // Fair Crossover weight on generate (out of 100)
pCrossFairLimit,  // Fair Crossover limit on inserted tree, 0 => no limit
pCrossFHWt,       // Fair Homologous Crossover weight on generate (out of 100)
pJoinWt,          // Joint two trees with new binary root node
pMuteCrossWt,     // Chance of mutating child produced by crossover
pMuteWt,         // overall Mutation weight on generate
pMuteNodeWt,         // Normal-Mutation weight on generate
pMuteConstWt,         // C-Mutation weight on generate
pMuteShrinkWt,         // Shrink-Mutation weight on generate
pMuteSubTreeWt,        // Koza subtree replacement Mutation weight on generate
pMuteFairWt,           // subtree (equal height) replacement Mutation weight
pMuteFairLimit,        // max size of mutation tree created by MuteFair
pAnnealMuteWt,                 // Mut. Annealing weight
pAnnealCrossWt,         // Crossover Annealing weight
pCopyWt,         // Copy weight on generate
pDirCrossWt,         // Directed Crossover weight (given Crossover)
pSeeds,              // Number of seeds in initial pop
pMuteSeedWt,         // Chance of mutating seeds (except 1st)
pSelectMethod,         // can be tournament, ranked, proportional
pGenerational,         // 0 => steady state
pThreads,              // 0 no separate threads, otherwise num extra threads
pTournSize,         // tournament size 2-10
pTournComp,         // If select or kill cant choose, compare remaining
                    // winners with upto PTournComp in rest of pop
pTournOff,          // num generates after which pTournSize is reduced to 1
pNicheSize,         // 0 => no niches otherwise in same niche if identical
pElitist,           // <>0 keep best. Only for pareto so far
pPopWidth,      // if !=0 treat pop as being torroid with this as width
pDemeWidth,     // if !=0 treat deme as being rectangle. nb area =pMateRadius
pMateRadius,         // mating radius.  0 for panmictic
pGaussRegion,         // 0 for flat selection in region, 1 for gaussian
pRepeatEval,         // repeat evaluation on copy?  For sampled evals
pKillTourn,         // number for anti-tournament <=10
pMaxAge,         // age to start killing off a good guy
pParsimony,

pFitnessCases,     // number of fitness cases
pStoreFit,         //1=>use storefitness LIST:2=>across runs memory charging
pCPU,              //problem dependant
pNeighbours,       // number of neighbours to be analysed by ALL
pPDGPinit,         // chance per 100 of init with introns
pSAANWt,           // SAAN  weight on PDGP crossover
pSIANWt,           // SIAN  weight on PDGP crossover
pSAINWt,           // SAIN  weight on PDGP crossover
pSIINWt,           // SIIN  weight on PDGP crossover
pSSAANWt,          // SSAAN weight on PDGP crossover
pSSIANWt,          // SSIAN weight on PDGP crossover
pSSAINWt,          // SSAIN weight on PDGP crossover
pSSIINWt,          // SSIIN weight on PDGP crossover
pPDGP,             // chrome width at each depth
PARAM_COUNT }; // total number of parameters

enum {
pStartTemp,        // for SA, 0=> same or better <0 better only
pEndTemp,          // for SA
pSamePen,          // fitness penalty for doing as parent
pAllMinFit,        // min fitness for program to be analysed by ALL
PARAM_DBL_COUNT }; // total number of double parameters

enum {
pDumpfile1,        // basic name of dumpfile, can include path
pDumpfile2,        // basic name of dumpfile, both 1 & 2 NULL => dont use
pZipCmd,           // command to run compress on dumpfile output
pPopfile,          // name of population output file, can include path
PARAM_STR_COUNT }; // total number of string parameters

#define PARAM_DBL_OFF PARAM_COUNT
#define PARAM_STR_OFF (PARAM_COUNT+PARAM_DBL_COUNT)
#define PARAM_TOTAL   (PARAM_COUNT+PARAM_DBL_COUNT+PARAM_STR_COUNT)

#define pname_size 20
#define pTraceStatic  1
#define pTraceDynamic 2
#define pTraceDStats  4
#define pTracePareto  8
#define pTracePStats 16
#define pTraceTest   32
#define pTraceDupStats   64
#define pTraceAllTests  128
#define pGenStats       256
#define pTraceLPStats   512
#define pTraceZStats   1024

class ChromeParams {    // parameter set for initializing a chrome
#ifdef PDGP
        enum {PDGPMAXHEIGHT = 100};
	int pdgpwidth[PDGPMAXHEIGHT];
	int pdgpindex[PDGPMAXHEIGHT+1];
        int net_limits[NUM_TREES+1];     // start of each separate network
	int pdgpheight;
void loadpdgp(char* string);
void savepdgp(ostream& out) const;
public:
int  netdepth(const int net) const;
//friend node* Chrome::pdgp(const int x, const int y) const;
friend class Chrome;
friend class Pop;
#endif /*PDGP*/
public:
        enum {param_str_size = 512};

	int params[PARAM_COUNT];                //Initialised in constructor
        double params_dbl[PARAM_DBL_COUNT];                   // ditto
        char params_str[PARAM_STR_COUNT][param_str_size];     // ditto
	char pnames[PARAM_COUNT+PARAM_DBL_COUNT+PARAM_STR_COUNT][pname_size];
	int varcount;
	int funccount;

	ChromeParams();
	virtual ~ChromeParams();
	virtual void Edit() {};                 // Put your own edit code here
	void Save( ostream& out = cout);
	void Load( istream& in);
	void Load( char * s);

};

// Carrier for fitness values.
// Feel free to subclass this if you want to track your fitness cases more closely (by adding data)
// or minimize instead of maximize (by redefining IsBetter)
// or combine multiple scores
class FitnessValue {
public:
	float fvalue;			// Always provide this for reporting reasons
#ifdef INTRON64
	int exon;
	int nsol;
#endif /*INTRON64*/
	virtual BOOL IsBetter(FitnessValue* fv) {return fvalue>fv->fvalue;};
	virtual BOOL IsBetter(float fv) {return fvalue>fv;};
	virtual void update_rank() {;};
	virtual void Copy(FitnessValue* fv) {memcpy(this,fv,sizeof(FitnessValue));};
	virtual ~FitnessValue() {;};
	operator float() {return fvalue;};
	BOOL operator >(FitnessValue& fv) {return IsBetter(&fv);};
	BOOL operator <(FitnessValue& fv) {return fv.IsBetter(this);};
#ifdef PARETO
	virtual BOOL ParetoBest() {cerr<<"opps default ParetoBest() called\n";}
#endif
	virtual void write(ostream& ofile = cout) { ofile<<fvalue; };
	virtual BOOL unique_best() const
	        {cout<<"FitnessValue::unique_best\n";return TRUE;}
#ifdef GEN_MEM
	virtual void Save( ostream& out = cout);
	virtual int Load( istream& in); //LOAD_OK or otherwise
#endif
};

#define EXPRLEN 300000000  // maximum length of expression (NOT short int NOW)
enum {LOAD_OK,LOAD_BADFILE,LOAD_BADFUNC,LOAD_TOOFEW,LOAD_TOOMANY,LOAD_BADCONST,LOAD_TOOLONG,BAD_TREE_NAME,BAD_TREE_FORMAT,LOAD_BADLINK,LOAD_FUNCTOODEEP};   // subload error values
enum {PRETTY_NONE,PRETTY_PARENS,PRETTY_NOPARENS};       // source listing options



// The Chrome is the carrier for a program expression
// It includes the eval, initialization, reading, mutation, crossover, reading and writing

class Chrome     {               // GP program structure
public:
//#ifndef PDGP
        int ip;                  // Current instruction pointer, EXPRLEN
#ifdef MULTREE
        int tree_starts [NUM_TREES+1];  // start of each tree,  EXPRLEN, last element gives total size
#endif
//#endif /*PDGP*/
	node *expr;                                     // actual code
	//float lastfitness;                              // fitness on last eval
	INT32 birth;                                             // gencount when generated
        Problem* probl;
//	CSelector* funcbag;                              // select functions by weight
        FitnessValue* nfitness;                       // pointer to an allocated object holding the type of fitness
	Function** funclist;            // array of pointers to functions and operations
	ChromeParams* params;
	retval* constlist;              // array of constants
#ifdef TRACE
        short int num_children;// = 0;
	struct parent { INT32 birth;// = -1;
	         int xpoint;// = -1; //EXPRLEN
#ifdef MULTREE
                 int xtree;// = -1;  //EXPRLEN
#endif
#ifdef PTHREADS
                 int xsize;// = 0;   //EXPRLEN
#endif
                 parent(): birth(-1), xpoint(-1)
#ifdef MULTREE
                           ,xtree(-1)
#endif
#ifdef PTHREADS
                           ,xsize(0)
#endif
                 {};
               };
	parent mum;
	parent dad;
#endif
#ifdef PTHREADS
	int thread = -1;
#endif

	Chrome(ChromeParams* p,Problem* pr,Function** f,retval* c,BOOL doinit=TRUE, istream* fil = NULL );     // returns a randomly initialized condition
#ifdef PDGP
        int   pdgpi(const int x, const int y) const;
inline  node* pdgp(const int x, const int y) const {return &expr[pdgpi(x,y)];};
inline  int   max_depth(const int net) const 
                       {return params->net_limits[net+1]-1;};
        void  setpdgp(const int x, const int y, const node* a);
        void  setpdgp(const int x, const int y, const int op);
        void  writepdgpnode(const node* n, ostream& ostr) const;
        void  writepdgp(const int tree, ostream& ostr) const;
        int   activepdgpnode(const int net) const;
        int   activepdgp(const int spp, const depth, int size, int tree[EXPRLEN]) const;
        int   pdgpnode(const int net) const;
        void  copypdgp(const int depth, const int y, const int net,
		       const Chrome* in, const int depth2, const int y2, 
		       const BOOL full);
inline  int   X(const int ptr, const int net) const {
                  for(int d=params->net_limits[net];d<=max_depth(net);d++)
		    if(params->pdgpindex[d+1]>ptr) return d;
		  }
inline  int   Y(const int ptr, const int net) const {
	          const int y = ptr - params->pdgpindex[X(ptr,net)];
#ifdef assert
//		  cout<<" Y("<<ptr<<","<<net<<")="<<X(ptr,net)<<","<<y<<endl;
		  assert(y>=0 && y<params->pdgpwidth[X(ptr,net)]);
#endif /*assert*/
		  return y;
                  }
        Chrome* CrossNet(Chrome* mate, int net );
	int loadpdgpnode(istream& istr,node* buf, const int d, const int net) const;
	int SubLoadpdgp(istream& istr,node** buf, const int net) const;// load a network
#else
	void reinit_trees(BOOL which_trees[NUM_TREES]);
#endif /*PDGP*/
	virtual Chrome* Copy();                 // copy yourself.  Extend for subclasses, and call base.
	virtual void Dup(Chrome* source);  // make a copy without re-initializing.  Extend for subclass, and call base
	virtual ~Chrome();
#ifdef PDGP
        void SubInit(const int d, const int y, const int net, CSelector* funcbag, const int depth, const BOOL isfull, const int mindepth );
	void initpdgp();
	evalnode* SetupEvalpdgp(evalnode* ip, const int depth, const node* spp) const;
#endif /*PDGP*/
#ifdef MULTREE
	void SubInit(CSelector* funcbag, int argstogo, int maxlen, 
	             int full, BOOL start, int tree);
#else
	void SubInit(CSelector* funcbag, int argstogo,int maxlen,int full=FALSE, BOOL start=FALSE);     // Initialize a subtree half full
#endif
#ifdef RAND_TREE
	void lable_tree(const gtree* tp, CSelector* const funcbag);
	BOOL rand_tree(const int length, const int tree, const gtree* tp=NULL);
#endif
	BOOL IsBetter(Chrome* chrome) {return nfitness->IsBetter(chrome->nfitness);};

#ifndef FASTEVAL
	reteval eval() {ip++;return funclist[FUNCNUM(expr[ip])]->eval(this);} ;
#endif
#ifdef MULTREE
	reteval evalAll(int) ; // eval the whole expression anew
#else
	reteval evalAll() ; // eval the whole expression anew
#endif
	void SetupEval(const int tree = -1, evalnode* = NULL);// expand expression for multiple evals
	void Traverse();              // Skips over next subtree
	void TraverseGlobal();
	int Depth(int start, int end);
	int DepthMax(int start, int end);
	int DepthMax(const int start);
	int DepthMax(); //of subtree rooted at ip

typedef struct {
int backip; //EXPRLEN
int size;   //EXPRLEN
short int depth;  //
short int argnum; //max_arity
} nodeinfo;
        void route(const nodeinfo n[], const int nip, int route[]);
        int compareroutes(const nodeinfo a[], const int aip, const nodeinfo b[], const int bip );
        int treeinfo(const int ip, const int root, const int depth,const int a,
                     nodeinfo subtreeinfo[]) const;
        void InsertTree(const int tree, const int thislen,
                        const int thissubptr, const int thissublen,
  const node* mateexpr, const int matesubptr, const int matesublen);
#ifdef MULTREE
        int ChooseCrossTree(Chrome* mate);
	Chrome* CrossTree(Chrome* mate, int tree, const BOOL fair,
			  const BOOL onepoint, const BOOL strict);
	Chrome* JoinTree(Chrome* mate, int tree);
	int GetIntNode(int tree);  // get an internal node for crossover
	int GetAnyNode(int tree);  // get any node for crossover
#else
	Chrome* CrossTree(Chrome* mate, const BOOL fair);
	int GetIntNode();                                       // get an internal node for crossover
	int GetAnyNode();                               // get any node for crossover
#endif
	void Mutate();                                  // mutate nodes with rate r
	void MutateC();
	void MutateShrink();
	void MutateSubTree(BOOL fair);
	int SubLen(int startat=0)
	 {ip=startat;Traverse();return ip-startat;};    // return length of expression at startat
	int Length() { 
#ifdef PDGP
	  return (params->params[pPDGP])? params->params[pMaxExpr]       :
		 tree_starts[NUM_TREES-1] + SubLen(tree_starts[NUM_TREES-1]);
#endif
	  assert(tree_starts[NUM_TREES]>0);
	  return tree_starts[NUM_TREES];
		     };    // return length of chrome

#ifdef MULTREE
	void write(const int pretty,ostream& ofile = cout, 
		   const int last_tree=NUM_TREES-1, const int limit=MAXINT)
#else
	void write(int pretty,ostream& ofile = cout)
         {ip=-1;SubWrite(ofile,pretty,0);};}      // write the expression
#endif
        ;
#ifdef TRACE
        void write_trace(ostream& ofile = cout);
#endif
	void SubWrite(ostream& ostr,const int pretty = 0, const int depth = 0, const int limit=MAXINT,
	      int* size=NULL, int* maxdepth=NULL);
#ifdef MULTREE
	int SubLoad(istream& istr,node* buf, int tree);// load an expression.
	int FindFunc(const char* funcname, const int tree) const;  // find index of a function
#else
	int SubLoad(istream& istr,node* buf);   // load an expression.  Return 0 on success or error val
	int FindFunc(char* funcname);                   // find index of a function, or -1 if not found
#endif
	virtual int Load(istream& istr,int issource); // load from a stream.  Return success
	BOOL same(Chrome* other) const; //Is code for two chromes the same
	int hashcode() const;           //hash code
        BOOL iscomment(const int c,istream& istr) const;
};
extern void DisplayMStats( ostream& out = cout);

const int pcount_size = 4;//fix later

class Problem {         // Sets up a particular GP problem.  GA independent.
					// contains function set, variables, and fitness function
#ifdef RAND_TREE
	rand_tree_tables* rand_tree[NUM_TREES];
#endif
        enum {max_arity = 4};
	int arityc[NUM_TREES][max_arity+1];// number of fncs of this arity
	int arityop[NUM_TREES][max_arity+1];// opcode of 1st fnc of this arity
protected:
	retval* constlist;              // array of constants in order
	retval* varlist;                // array of variables
public:
        int funcconstmin;
        int funcconstmax;
	Function** funclist;                    // primitive functions

public:
	int funccount;
	CSelector* funcbag[NUM_TREES];  // select functions by weight

	Problem();     // allocate the arrays
	virtual ~Problem();
	void AddF(Function* f)          // add a function to the problem
	 {funclist[funccount++]=f;};
	Function** getfuncs() {return funclist;};
	retval* getconsts() {return constlist;};
//	CSelector* getfuncbag() ;
        void addtofuncbag(int tree, Function* f);
inline	int Arityc(const int t, const int a) const {return arityc[t][a];};
#ifdef CCCPOLY_LIB
inline	int  RndOp(const int t, const int a) const {
     return (a>0)?       arityop[t][a]+rnd(arityc[t][a]) : //Function
       (funcconstmin==0)?                                  //no constants
        ((arityc[t][0]==1) ? arityop[t][0]:                //First terminal 'X'
          arityop[t][0]+rnd(arityc[t][0]) ):               //random variable
        ((funcconstmin==arityop[t][0] || rnd(2)!=0)?       //no variables 
	  funcconstmin+rnd(1+funcconstmax-funcconstmin):   //or choose constant
	  arityop[t][0]+rnd(funcconstmin-arityop[t][0]) ); //random variable
#else
inline	int  RndOp(const int t, const int a) const {
                                  return arityop[t][a]+rnd(arityc[t][a]);};
#endif /*POLY_LIB*/
#ifdef RAND_TREE
        void init_1rand_tree(const int size); //one size
        void init_rand_tree(int size_limit);
        int* rand_tree_lookup(const int length, const int t ) const;
#endif
	// User hooks
	virtual FitnessValue* GetFitnessObj(Chrome* c=NULL);	// Can change behavior of FitnessValue
	virtual float fitness(Chrome* chrome);//, evalnode* = NULL);  // The active ingredient.  You add this.
	virtual void WriteTreeName(int tree, ostream& ostr = cout ) // The active ingredient.  You add this.
	{ostr<<" T"<<tree<<" =";}; //default You overwrite this	
virtual BOOL TreeNameMatch(int tree, char* s ) // The active ingredient.  You add this.
	{return (s[0]=='T');}; //default You overwrite this	
	virtual Chrome* NewChrome(ChromeParams* p, istream* fil = NULL){
cout<<"Problem::NewChrome\n";
	       return new Chrome(p,this,getfuncs(),getconsts(), TRUE, fil);
	};// Can change behavior
	virtual Chrome* Bestof(const PtrChrome  list[], const int listsize,
	                        int* bestinlist = NULL, const int target = 0 );
	virtual Chrome* Worstof(const PtrChrome list[], const int listsize,
	                        const int timenow, int* worstinlist,
				const int target = 0 );	
	virtual int static_check(Chrome* chrome, int tree) {return 0;};
#ifdef MAXTREEDEPTH
	virtual void ProbLPStats(int opcode_counts[][NUM_TREES][MAXTREEDEPTH+1],
	                 int array[pcount_size] ){;};
#else
	virtual void ProbLPStats(int opcode_counts[][NUM_TREES],
	                 int array[pcount_size] ){;};
#endif
	virtual int ProbLTStat(int tree, int start, int length, node* code){
                               return 0;};
};

class Pop {                     // a population for a particular GA.  Problem independent
					// controls the GA - selection, crossover, etc
        void select_candidates (int target, int tourn_size,
	                        Chrome* candidates[], int slots[], BOOL best ) const;
#ifdef MAXTREEDEPTH
	void GatherPStats( int opcode_counts[][NUM_TREES][MAXTREEDEPTH+1],
			   int pcount[pcount_size],
			   int tcount[NUM_TREES]           )   const;
#else
	void GatherPStats( int opcode_counts[][NUM_TREES],
			   int pcount[pcount_size],
			   int tcount[NUM_TREES]           )   const;
#endif
public:
#ifdef MULTREE
        int num_crosses[NUM_TREES];
	int num_crosses_cleared_at;// = -1;
#endif
public:
	Chrome** pop;           // allocated array holding the actual chromes
#ifdef GENERATIONAL
	Chrome** newpop;        // array holding pointers to new chromes
	Chrome** oldpop;        // array holding pointers or zero
#endif
	ChromeParams* params;
	Problem* problem;
	UINT popsize;           // size of array
	INT32 gencount;          // number of individuals generated
	time_t start;           // starting time
	time_t elapsed;         // elapsed time for last go_until
	BOOL isAborted;

#ifndef PARETO
	FitnessValue* BestFitness;      // Best current fitness value
	int   BestMember;       // Index of best current member
#endif
	int initeval;           // Flag: Still evaluating initial population?
	UINT nexteval;          // Next initial pop member to evaluate

	Pop(Problem* prob,ChromeParams* par, const int size,istream* = NULL );       // set up a population for a particular problem
	Pop(Problem* prob,ChromeParams* par);   // stub for derived classes
	virtual ~Pop();
#ifndef PARETO
    Chrome* best();                         // return the current best chrome
#endif
	virtual Chrome* selectParent(const int target, int* out = NULL) const; // select one parent. Return pointer to selected parent
						 // target identifies selection region in pop

#ifdef PTHREADS
//	void thread_fitness(const void *t);       // calculate multiple fitnesses in separate pThread
	void generation_fitness(const int start); // calculate all remaining fitnesses for this generation
#endif
#ifdef MULTREE
	virtual Chrome* generate(int tree);     // generate a new chrome.  Return fitness
#else
	virtual Chrome* generate();     // generate a new chrome.  Return fitness	
#endif
							// Alternate steady state GA code can go here
						// generate until reaching time max_time, evaluating maxevals individuals, or reaching maxfitness
#ifdef MULTREE
	virtual Chrome* go_until(time_t max_time, int maxevals,
                                 float maxfitness, int tree = -1);
#ifndef PDGP
        void reinitialise(int change_trees[NUM_TREES]);
#endif /*PDGP*/
#else
	virtual Chrome* go_until(time_t max_time, int maxevals, float maxfitness);
#endif
	int GetTarget();                // get a target member to replace with anti-tournament
	double Temperature() const;
	BOOL accept(const FitnessValue* f1, const FitnessValue*) const; // new ok for annealing?
	int GetAnnealTarget();                  // get a target for annealing, with selection tournament
	virtual BOOL Aborted(){return isAborted;};
	virtual float Fitness(Chrome* chrome);//, evalnode* xip = NULL);   // Evaluate a member using the current problem
	void UpdateBest(const int slot, const Chrome* NewChrome);
	void InsertMember(int slot,Chrome* NewChrome,int nodelete = FALSE);                     //insert a chrome and update BestMember values
	void DisplayStats( ostream& out = cout) const;
	void DisplayDStats( ostream& out = cout) const;
	void DisplayPStats( ostream& out = cout) const;
	void DisplayLPStats( ostream& out = cout) const;
	void DisplayDupStats( ostream& out = cout) const;
	void ClearCStats() {memset(num_crosses,0,sizeof(num_crosses));
	                    num_crosses_cleared_at = gencount;};
	void DisplayCStats( ostream& out, int gentime) const;
	void DisplaySStats() const;
	void DisplayZStats(ostream& out = cout) const;
	 int Parent(const INT32 birth) const;
	void DisplaySolStats( ostream& out = cout) const;
	void Write( ostream& out = cout, BOOL binary = FALSE );
#ifdef MULTITEST
	void Retest();
	void Save();
	void Restore();
#endif /*MULTITEST*/
};

#endif          // #ifndef _CHROME_H
////////////////////////////////////////////////////////////

//may not be needed eg Borland
int strcmpi(char* s1, char* s2);         // only checks equality

void clear_rank();

class ghashchain;

class ghashtable {
public:
  ghashchain** table;
  int size;
  inline int safe_index(int code) const { return (code%size+size)%size; }
  ghashtable(int s); 
  ~ghashtable(); 
  inline ghashchain* chainhead(int code) const
	  { return table[safe_index(code)]; };
  inline void setchainhead(int code, ghashchain* value) 
	  { table[safe_index(code)] = value; };
  ghashchain* find(const int code,BOOL (*same)(void*),BOOL& found) const;
//friend class hashchain;
};

class ghashchain { 
public:
  ghashchain* chain; 
  ghashchain(ghashtable* t, int code) : chain(NULL) {
    t->setchainhead(code,this);
  }
  ghashchain(ghashchain* endofchain)  : chain(NULL) {
    endofchain->chain = this;
  }
};

class costats;

class stats {						//WBL
	int cases;// = 0;
	float min;// = FLT_MAX;
	float max;// = -FLT_MAX;
	int  nmax;// = -1
	float sum;// = 0.0;
	double sqr;// = 0.0;

double	Variance();

public:
	stats(): cases(0),
	         min(FLT_MAX),
	         max(-FLT_MAX),
	         nmax(-1),
	         sum(0.0),
	         sqr(0.0){};

       void 	include ( float in );
inline int	Cases() { return cases; };
inline float	Min()   { return min; };
inline float	Max()   { return max; };
inline float	Sum()   { return sum; };
inline float	Mean()  { return (cases>0)? sum / cases : -FLT_MAX; };
inline float	SD()    { return (cases>0)? sqrt(Variance()) : -FLT_MAX; };
       void	display ( char* name,              ostream & fout = cout,
                          stats* own_stats = NULL, stats* co_stats = NULL,
                          costats* covar = NULL );
friend class costats;
};//end class stats

class costats { //WBL
  stats x;
  stats y;
  double product;
public:
  costats(): product(0) {;}
void include(float xx, float yy){x.include(xx); y.include(yy); product+=xx*yy;}
double Covariance() {return product/x.Cases() - x.Mean()*y.Mean();}
double Correlation() {
  double d = sqrt(x.Variance() * y.Variance());
  if ( d > 0 )			//protect rounding error
    return Covariance() / d;
  else
    return FLT_MAX;
}
};//end class costats

void update_fitter(const int fit, const float fitter);
void update_neighbours(const int fit, const int vworse, const int worse, 
		                     const int same,   const int better);

void display_neighbour_stats();
