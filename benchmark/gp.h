// gp.h   Include file for gp.cc and problem specific file
// W. Langdon cs.ucl.ac.uk 23 Jan 1995 from queue2.h version 1.08

// version "$Revision: 1.67 $"

//Modifications (in reverse order)
//WBL 19 Feb 2000  Make rawfit a double

//WBL  6 Sep 1999  Add vhits, vfit, passes

//WBL  4 Jun 1999  Add rawfit

//WBL  7 Sep 1997  Add FitnessAllTests

//WBL 12 Apr 1997  Make ProbLPStats support MAXTREEDEPTH

//WBL 13 Jun 1996  Go back to south wales size schedule

//WBL 18 Oct 1995  Support random test selection

//WBL 15 Oct 1995  Add IsBranch

//WBL 28 Sep 1995  Ensure Memory_errors works without STORE_FIT

//WBL 31 Aug 1995  Get STORE_FIT working without stop_on_error

//WBL 20 Jun 1995  Add ProbLTStat

//WBL 14 May 1995  Add cpu_time_max

//WBL 12 May 1995  Add memory_errors (and reduce cpu_time, so it fits)

//WBL  7 May 1995  Pack data in myfitnessvalue tighter (nb changes fitness)

//WBL  6 May 1995  Pack data in myfitnessvalue tighter

//WBL  3 May 1995  Add STORE_FIT. 
//                 Remove update_rank call from Copy and its unneeded argument

//WBL  3 May 1995  BUGFIX remove nfitness from mychrome. Always use base class

//WBL 30 Apr 1995  Remove chains from myfitnessvalue
//                 Leave Chrome* argument to myfitnessvalue, Copy
//                 and GetFitnessValue for the time being

//WBL 26 Apr 1995  Death by paging noted on hodgkin with 5-6 meg in machine.
//                 I am attributeding this to paging caused by scanning elite
//                 chain when remove an elite individual. Get populations where//                 almost all the pop is in one or more chains.
//                   Add Bchain

//WBL 23 Apr 1995  Add Ntests_run

//WBL 22 Apr 1995  Add cache_init, cache_stats, cache_stats_init

//WBL 20 Apr 1995  Add IsLoop

//WBL 18 Apr 1995  Add AddF8

//WBL 17 Apr 1995  Add print_funcnum

//WBL  6 Apr 1995  Add ARG1_funcnum and ARG2_funcnum, FUNC_funcnum

//WBL 31 Mar 1995  Add AddF7

//WBL 23 Mar 1995  Add Chrome* argument to GetFitnessObj, myfitnessvalue
//                 and Copy

//WBL 14 Mar 1995  Add chain to myfitness_value

//WBL 11 Mar 1995  Add arg1_funcnum and arg2_funcnum

//WBL 26 Feb 1995  Update definition of select_hits

//WBL 18 Feb 1995  Add reporting indetical scores in select_hits

//WBL 11 Feb 1995  Use clear() to initialise myfitnessvalue

//WBL 23 Jan 1995  New file

#ifdef STORE_FIT
#ifdef TRACE_RUN
#include "trace.h"
#endif
#endif

//extern indexed_memory* indexed_memory_ptr; quick max frig

class myfitnessvalue : public FitnessValue {
public:
#ifdef GEN_MEM
       retval genetic_memory[GENETIC_MEM_SIZE];
       void load_genetic_memory();
       void read_genetic_memory();
#endif
#ifdef GRID_LIB
       char schedule [42]; //NL
#endif
#ifdef RAWFITSTAT
	double rawfit;
#endif /*RAWFITSTAT*/
//#ifdef PARETO
       scoretype hits [NUM_OPERATIONS+NUM_OTHERS];
       int passes [1000];
       float vfit;
       scoretype vhits;
//#endif
       unsigned char mem_used; //=0
#ifndef STORE_FIT
       unsigned short int Memory_Errors; //=0
#endif
       unsigned short int Ntests_run; //=0
       int test_seed; //=0
#ifdef MINC_ADF1
       int adf1; //=0
#endif
#ifdef TIMINGS
#ifdef ANT_LIB
       float cpu_time; //=0
#else
       unsigned int cpu_time:24; //=0
#endif
#endif
#ifdef TRACE_RUN
#ifdef STORE_FIT
       unsigned char score[MAX_BREAKS]; 
       signed char test_phase; //=-1;
       int ok (int tree) const;
       int nak(int tree) const;
       void unpack_hits(); 
#else /*STORE_FIT*/
       unsigned char  ok[NUM_TREES];
       unsigned char nak[NUM_TREES];
#endif
#endif

#ifdef STORE_FIT
       struct per_test {
       unsigned char passes;
       unsigned char mem_used; //=0
       unsigned char Ntests_run; //=0
       unsigned char last_test_offset;
       unsigned char last_oktest_offset;
       unsigned char memory_errors;
#ifdef TIMINGS
       unsigned int cpu_time:16;
       enum {cpu_time_max = 65535};
#endif
#ifdef TRACE_RUN
       trace_packed trace;
#endif
       per_test(): mem_used(0), Ntests_run(0) {};
       };
       per_test test_data[NUM_TEST_SEQUENCES];
#endif
       myfitnessvalue(Chrome* p): mem_used(0)
#ifndef STORE_FIT
                ,Memory_Errors(0)
#endif
                ,Ntests_run(0)
                ,test_seed(0)
#ifdef MINC_ADF1
                ,adf1(0)
#endif
#ifdef TIMINGS
                ,cpu_time(0)
#endif
#ifdef STORE_FIT
                ,test_phase(-1)
#endif
	       {clear();};
#ifdef PARETO
       ~myfitnessvalue();
#endif //PARETO
       void clear();
#ifdef PARETO
       BOOL IsBetter(FitnessValue* fv);
#else
#ifdef SAMEPENALTY
       BOOL IsBetter(FitnessValue* fv);
#endif
#endif
       void write(ostream& fout = cout );
#ifdef PARETO
       int dominates(myfitnessvalue* target);
       void update_rank(); //only works for best so far
       BOOL unique_best() const;
#endif
       void Copy(FitnessValue* fv) {
	       memcpy(this,fv,sizeof(myfitnessvalue)); };
       BOOL ParetoBest();
#ifdef STORE_FIT
       int Memory_errors() const;
#else
       inline int Memory_errors() const {return Memory_Errors;};
#endif
};//end class myfitnessvalue;

typedef myfitnessvalue* ptrmyfitnessvalue;

#ifdef GRID_LIB
scoretype pareto_fitness(const ptrmyfitnessvalue x,const int i);
#endif

class mychrome : public Chrome {
public:
       mychrome(ChromeParams* p,Problem* pr,Function** f,retval* c,
		istream* fil = NULL):
	       Chrome(p,pr,f,c,TRUE,fil) {};
};//end class mychrome

class gp : public Problem {
public:
	int i0_funcnum;//=0
	int for_funcnum;//=0;
	int down_funcnum;//=0;
	int while_funcnum;//=0;
	int forwhile_funcnum;//=0;

	int arg1_funcnum;//=0;
	int arg2_funcnum;//=0;
	int ARG1_funcnum;//=0;
	int ARG2_funcnum;//=0;
	int FUNC_funcnum;//=0;

	int print_funcnum;//=0;

	const char* ProbVersion;

	gp();             // Initialize primitives,
	void WriteTreeName(int tree, ostream& ostr = cout );
	                               // parameters, and data
	BOOL TreeNameMatch(int tree, char* s);
	float fitness(Chrome* chrome); // GP fitness function
#ifdef MULTITEST
	void FitnessAllTests(Chrome* chrome, ostream& fout);
#endif

	Chrome* NewChrome(ChromeParams* p, istream* fil = NULL ){
//cout<<"queue::NewChrome\n";//debug
	       return new mychrome(p,this,getfuncs(),getconsts(),fil);
	};

	FitnessValue* GetFitnessObj(Chrome* chrome = NULL) {
//cout<<"queue::GetFitnessObj\n";//debug
return new myfitnessvalue(chrome);};

#ifdef PARETO
	Chrome* Bestof( const PtrChrome  list[], const int listsize, 
			int* bestinlist = NULL, const int target = 0 );
	Chrome* Worstof(const PtrChrome list[], const int listsize, 
			const int timenow, int* worstinlist,
			const int target = 0);
#endif
        int  static_check(Chrome* chrome, int tree);
	void AddF1    (int tree, Function* f );
	void AddFbar1 (int exclude_tree, Function* f );
	void AddF2    (int tree1, int tree2, Function* f );
	void AddF3    (int tree1, int tree2, int tree3, Function* f );
	void AddF4    (int tree1, int tree2, int tree3,
		       int tree4, Function* f );
	void AddF5    (int tree1, int tree2, int tree3,
		       int tree4, int tree5, Function* f );
	void AddF6    (int tree1, int tree2, int tree3,
		       int tree4, int tree5, int tree6, Function* f );
	void AddF7    (int tree1, int tree2, int tree3,
		       int tree4, int tree5, int tree6, 
		       int tree7, Function* f );
	void AddF8    (int tree1, int tree2, int tree3,
		       int tree4, int tree5, int tree6, 
		       int tree7, int tree8,            Function* f );
	void AddF9    (int tree1, int tree2, int tree3,
		       int tree4, int tree5, int tree6, 
		       int tree7, int tree8, int tree9, Function* f );
	void AddF11   (int tree1,  int tree2,  int tree3,  int tree4, 
                       int tree5,  int tree6,  int tree7,  int tree8,
                       int tree9,  int tree10, int tree11, Function* f );
	void AddF13   (int tree1,  int tree2,  int tree3,  int tree4, 
                       int tree5,  int tree6,  int tree7,  int tree8,
                       int tree9,  int tree10, int tree11, int tree12,
                       int tree13, Function* f );
	void AddFmain (Function* f );
	void AddFall  (Function* f );
	void LoadTests( istream& in=cin);
	void write_stats( ostream& out=cout);
#ifdef MAXTREEDEPTH
	void ProbLPStats( int f[][NUM_TREES][MAXTREEDEPTH+1], int pcount[pcount_size]);
#else
	void ProbLPStats( int f[][NUM_TREES], int pcount[pcount_size]);
#endif
	int  ProbLTStat(int tree, int start, int length, node* expr);
};//end class gp

extern gp* ThisProblem;

void find_elite(Chrome* list[], const int listsize, int& output_size,
	       scoretype best[num_ellite_pareto],int count[num_ellite_pareto]);

#ifdef PARETO
Chrome* select_hits (const PtrChrome  list[], const int listsize,
                     const int target, int* bestinlist, const BOOL best = TRUE,
	             int* select = NULL, int* duplicates = NULL, 
		     int* select_size = NULL);
void select_hits_write_stats(ostream& fout = cout);
#endif

float display_run(Chrome* chrome, const int gencount = 0);

BOOL IsLoop(int funcnum);
int  IsBranch(int funcnum);

#ifdef ADF_CACHE
void cache_init();
void cache_stats_init();
void cache_stats();
#endif

void write_problem_params();
