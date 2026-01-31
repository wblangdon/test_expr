// W. B. Langdon cs.ucl.ac.uk
// This code is released for non-commercial use only
#ifdef AVX
#define Version "AVX Sextic polynomial. $Revision: 1.39 $"
#else
#define Version "Sextic polynomial. $Revision: 1.39 $"
#endif

//WBL 22 Sep 2019 for g++ 8.3.1
//WBL 21 Sep 2019 Make sure only relevant binEval is defined
//WBL  8 Jan 2019 join AVX and PTHREADS
//WBL  5 Jan 2019 Use chrome.expr directly so replacing xip and SetupEval
//WBL 16 Dec 2018 For PTHREADS pass EvalSP,evalstack via EVAL
//WBL 10 Dec 2018 Add AVX instructions
//WBL 26 Nov 2018 Replace rand_0to1() by test::test() code poly.cc r1.17
//WBL 23 Nov 2018 Rewrite bool.cc r1.38 as avx.cc, cf poly.cc r1.17, primitiv.cxx r1.18
//WBL  5 Aug 2016 Replace fit64 with raw value
//WBL 15 Jul 2016 Add fit64
//WBL 12 Jul 2016 Add exon
//WBL  6 Jun 2016 for ~/gp/eden_mux add INTRON64. For simplificity for INTRON64
//put all opcodes here and dispense with primitive.cxx

//Modifications (reverse order):
// WBL 30 Apr 1999 Conversion for SiliconGraphics cc -c -64 -mips4

// WBL  9 Jan 1999 

// WBL 24 Jun 1998  Allow problem_order to be up to 13

// WBL 11 Jun 1998  Add EPARITY

// WBL  9 Mar 1998  Add reading prim.dat

// WBL  3 Mar 1998  based upon r1.35 antprob.cxx and ant.cc r1.29

 
// GPQUICK

#include "pch.h"
#include "chrome.h"

//https://gcc.gnu.org/onlinedocs/gcc-4.8.2/gcc/Restricted-Pointers.html
//tell each thread which data structures to use
//for performance ensure read-write per thread data are in different cache lines
#ifdef PTHREADS
#undef IP
#undef OPDEF
#undef GETVAL
//#undef CHROMEP
//#if __SIZEOF_POINTER__ != 8
//check pointer size used to align struct threads 
//#endif
//typedef struct {/*Chrome* chrome;*/ node* ip; /*retval* sp;*/ char padd[CACHE_LINE_SIZE /*- __SIZEOF_POINTER__ -*/  - __SIZEOF_POINTER__];} type_thread_info;
//type_thread_info threads[npthreads] __attribute__ ((aligned (CACHE_LINE_SIZE)));
EVALFUNC* Evalfunc  = NULL;
retval*   Constlist = NULL;
//#define CHROMEP threads[my_id].chrome
#define IP     (*ip)
#define EvalSP (*sp)
#define evalstackwidth MAXTESTCASES
#define inc_EvalSP  *sp += evalstackwidth
#define dec_EvalSP  *sp -= evalstackwidth
#define dec2_EvalSP *sp -= 2*evalstackwidth
#define OPDEF(Op)    reteval Op(node** ip, retval** sp)
#define OPDEF2(Op,f) reteval Op(node** ip, retval** sp,f)
#undef EVAL
#define EVAL  (Evalfunc[(++IP)->op])(ip,sp)
#define GETVAL Constlist[IP->op]
#define EVAL2        Eval2(ip,sp)
#define BINEVAL(f) binEval(ip,sp,f)
#else /*PTHREADS*/
#define OPDEF2(Op,f) reteval Op(f)
#define EVAL2 Eval2()
#define BINEVAL(f) binEval(f)
#endif /*PTHREADS*/

#include "primitiv.h"

#ifdef FASTEVAL
	extern evalnode* IpGlobal;
//	extern Chrome* ChromeGlobal;
#endif

//not used but keeps loader happy
int num_sol_found = 0;//used to report solution to each test phase found

//**************************************************************************
// Define your Problem class

#include <assert.h>
#include "prob.h"
#include "gp.h"
#include "test.h"
#ifdef AVX
#include <immintrin.h>
//#include <immintrin.h>
//#include <zmmintrin.h>
#endif

//**************************************************************************
// static variables that make up the evaluation environment
// Build your environment here as statics
// where the evaluation functions can get to the data

int evalstacksize = 0;
retval* evalstack = NULL;
#ifdef PTHREADS
#if retval != float
error evalstacksize assumes we are pushing and popping floats
#endif
#if((MAXTESTCASES*evalstacksize) % (CACHE_LINE_SIZE/__SIZEOF_FLOAT__) != 0)
error for performance reasons my_stackoff should be multiple of cache line size
#endif
//#define my_stackoff (my_id*MAXTESTCASES*evalstacksize)
int nthreads = 0;
//int* evalSP = NULL;
#else
#define my_stackoff 0
int EvalSP = 0;
#endif
float setmaxfit(const ChromeParams* p) {
  assert(p);
  assert(p->params);
  evalstacksize = (p->params[pMaxDepth]!=0)? p->params[pMaxDepth] : 
                                             10 * sqrt(p->params[pMaxExpr]+100);
  cout<<"evalstacksize: "<<evalstacksize<<" "<<flush;
  const int n = p->params[pThreads];
  nthreads = (n<=0)? 1 : n;
  cout<<"nthreads: "<<nthreads<<endl;
  assert(nthreads<=npthreads);
  //https://www.gnu.org/software/libc/manual/html_node/Aligned-Memory-Blocks.html
  evalstack = (retval*)aligned_alloc(512/8,nthreads*MAXTESTCASES*evalstacksize*sizeof(retval)); // new retval[MAXTESTCASES*evalstacksize];
  //evalSP = new int[nthreads];
#ifdef AVX
  assert((MAXTESTCASES% 8) == 0); // 8 needed by ConstEval
  assert((MAXTESTCASES%16) == 0); //16 needed by XEval etc
#endif
#ifdef PTHREADS
//assert(sizeof(type_thread_info) == CACHE_LINE_SIZE);
#endif
  return 0; //keep going for ever
}
retval* dataX = (retval*)aligned_alloc(512/8,MAXTESTCASES*sizeof(retval));  //new retval[MAXTESTCASES];
retval* dataY = (retval*)aligned_alloc(512/8,MAXTESTCASES*sizeof(retval));  //new retval[MAXTESTCASES];

//**************************************************************************
// Data handling functions

int current_phase = 0; 
int tests_run = 0;
int tests_apparently_run = 0;
int total_cpu = 0;

float sextic(const float X) {
  return (X*X*(X-1)*(X-1)*(X+1)*(X+1)); //poly.fitness.cc r1.5 and aigp3.tex r1.74
}

//gawk -f poly.test.awk -v "train=48,0,101" /dev/null > gp.test
//from poly.cc r.17
test::test(Problem* probl, istream& istr) {
cout<<"Loaded "<<flush;
int data_limit = -1;
{for(int t=0;t<MAXTESTCASES;t++) {
  istr>>dataX[t];
  istr>>dataY[t];
  data_limit = t;
  if(istr.peek()==EOF) break;
}}
cout<<(data_limit+1)<<" tests\n";
{for(int t=0;t<=data_limit;t++) {
  cout<<dataX[t]<<"\t"<<dataY[t]<<endl;
  assert(abs(dataY[t]-sextic(dataX[t]))<0.0001);
}}
assert((data_limit+1)==MAXTESTCASES);
}//end test::test
//end poly.cc r.17

/*int Reset() {
  **use fixed "random" selection cf AiGP3
%The 50 test points we used in all our experiments we chosen uniformly
%at random from the range -1 and 1.
%No granularity was imposed.
%Again by chance none of the three special
%values of -1, 0 or 1 where included and
%no value was repeated.**
  const int save = gpquick_seed;
  gpquick_seed = 713215;
  tests_run = 0;
  for(int i=0;i<MAXTESTCASES;i++) {
    //Xvalues[i] = 1;//better for debugging
    Xvalues[i] = 2*rand_0to1()-1;
    dataY[i]   = sextic(Xvalues[i]);
    cout<<"Xvalues["<<i<<"] "<<Xvalues[i]<<", dataY["<<i<<"] "<<dataY[i]<<endl;
  }
  gpquick_seed = save;
}
*/
//////////////////////////  Problem specific functions
#ifdef INTRON64
#ifndef FASTEVAL
error1!
#endif
#ifndef SMC_GP
error2!
#endif
//unsigned char intron64[EXPRLEN];
//unsigned char fit64[EXPRLEN];
//retval        raw[EXPRLEN];
//INT32 chrome64_age;
#endif

/*OPDEF(ConstEval) {** keep compiler happy ConstEval place holder never called *
  assert(1==0);
}*/
OPDEF(ConstEval) {
  //if(EvalSP>=evalstacksize) cout<<"ConstEval: EvalSP "<<EvalSP<<" evalstacksize "<<evalstacksize<<endl;
  assert(EvalSP < &evalstack[nthreads*MAXTESTCASES*evalstacksize]);
  const retval val = GETVAL;
#ifdef AVX
  for(int i=0;i<MAXTESTCASES;i+=8) {
    //https://software.intel.com/en-us/node/524140
    //No corresponding Intel AVX instruction.
    _mm256_store_ps(&EvalSP[i], _mm256_set1_ps(val));
  }
#else
  for(int i=0;i<MAXTESTCASES;i++) EvalSP[i] = val;
  //for(int i=0;i<MAXTESTCASES;i++) evalstack[my_stackoff+EvalSP*MAXTESTCASES+i] = val;
#endif /*AVX*/
  inc_EvalSP; //++EvalSP;
}

OPDEF(XEval){ 
  //cout<<"XEval"<<IP-ExprGlobal<<" "<<flush;
  //if(EvalSP>=evalstacksize) cout<<"XEval: EvalSP "<<EvalSP<<" evalstacksize "<<evalstacksize<<endl;
  assert(EvalSP < &evalstack[nthreads*MAXTESTCASES*evalstacksize]);
#ifdef AVX
  for(int i=0;i<MAXTESTCASES;i+=16) {
    //https://software.intel.com/en-us/node/523488
    //https://software.intel.com/en-us/node/523585
    _mm512_store_ps(&EvalSP[i],_mm512_load_ps(&dataX[i]));
    assert(EvalSP[i]   ==dataX[i]);
    assert(EvalSP[i+15]==dataX[i+15]);
  }
#else
  memcpy(EvalSP,dataX,MAXTESTCASES*sizeof(retval));
  //memcpy(&evalstack[my_stackoff+EvalSP*MAXTESTCASES],dataX,MAXTESTCASES*sizeof(retval));
#endif /*AVX*/
  inc_EvalSP; //++EvalSP;
}


inline
OPDEF (Eval2) {//EVAL two subtrees but leave 2 answers on top of evalstack
  assert(EvalSP >= evalstack);  //assert(EvalSP>=0);
  const retval* e0 = EvalSP;
  EVAL;
  assert(EvalSP==e0+evalstackwidth);
  EVAL;
  assert(EvalSP==e0+2*evalstackwidth);
  assert(EvalSP >= &evalstack[2*evalstackwidth]);
  dec2_EvalSP; //EvalSP -= 2;
  assert(EvalSP >= evalstack);
}
#ifdef AVX
inline
__m512 doadd(const __m512 a, const __m512 b){return _mm512_add_ps(a,b);}
inline
__m512 dosub(const __m512 a, const __m512 b){return _mm512_sub_ps(a,b);}
inline
__m512 domul(const __m512 a, const __m512 b){return _mm512_mul_ps(a,b);}

#if __GNUC__ > 7
#if (!defined(__OPTIMIZE__)) || (!defined(_AVX512FINTRIN_H_INCLUDED))
//from GCC 8.3.1 avx512fintrin.h allow to compile for GDB etc
extern __inline __mmask16
__attribute__ ((__gnu_inline__, __always_inline__, __artificial__))
_mm512_cmpneq_ps_mask (__m512 __X, __m512 __Y)
{
  return (__mmask16) __builtin_ia32_cmpps512_mask ((__v16sf) __X,
						   (__v16sf) __Y, _CMP_NEQ_UQ,
						   (__mmask16) -1,
						   _MM_FROUND_CUR_DIRECTION);
}
#endif
#else /*__GNUC__ <=7 */
//IntrinsicsGuide.txt: missing from GCC 7.3.1
inline
__mmask16 _mm512_cmpneq_ps_mask (const __m512 a, const __m512 b) {
  __mmask16 mask = 0;
  __mmask16 m    = 1;
  for(int i=0;i<16;i++,m*=2) {
    if(a[i] != b[i]) mask = mask | m;
  }
  return mask;
}
#endif
inline				// "Protected" division
__m512 dodiv(const __m512 numerator, const __m512 denominator){
  //https://software.intel.com/en-us/node/523793
  //const __m512    zero = _mm512_set1_ps(0.0f);
  const __m512    zero = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  //https://software.intel.com/en-us/node/523469
  const __mmask16 mask = _mm512_cmpneq_ps_mask(denominator,zero);
  //const __mmask16 mask = _mm512_cmpneq_epi32_mask(denominator,zero);
  //https://software.intel.com/en-us/node/523793
  const __m512    val  = _mm512_maskz_div_ps(mask,numerator,denominator);
  const __m512    one  = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
  //https://software.intel.com/en-us/node/523805
  const __m512 ans = _mm512_mask_blend_ps(mask,one,val);
  return ans;
}
inline
OPDEF2(binEval,__m512 f(const __m512 a, const __m512 b)) {
  EVAL2;
  for(int i=0;i<MAXTESTCASES;i+=16) {
    const __m512 sp0 = _mm512_load_ps(&EvalSP[i]);
    const __m512 sp1 = _mm512_load_ps(&EvalSP[MAXTESTCASES+i]);
    const __m512 val = f(sp0,sp1);
    _mm512_store_ps(&EvalSP[i],val);
  }
  inc_EvalSP; //++EvalSP;
}
#endif /*AVX*/
#ifndef AVX
inline
OPDEF2(binEval,float f(const retval a, const retval b)) {
  EVAL2;
  for(int i=0;i<MAXTESTCASES;i++) {
    const retval val = f(EvalSP[i],EvalSP[evalstackwidth+i]);
    EvalSP[i] = val;
  }
  inc_EvalSP;
}
inline
retval doadd(const retval a, const retval b){return a+b;}
inline
retval dosub(const retval a, const retval b){return a-b;}
inline
retval domul(const retval a, const retval b){return a*b;}
				// "Protected" division
//Nov 2018 should be ok for float
inline
retval dodiv(const retval numerator, const retval denominator){
	if (denominator !=0)
		return numerator/denominator;	//retval!=float
	else
		return 1.0f;			//retval!=float
}
#endif /*not AVX*/
OPDEF(AddEval) {return BINEVAL(doadd);}
OPDEF(SubEval) {return BINEVAL(dosub);}
OPDEF(MulEval) {return BINEVAL(domul);}
OPDEF(DivEval) {return BINEVAL(dodiv);}





//**************************************************************************
///////////////////////// Read Terminal and Function sets form file

//taken form ant.cc r1.29

BOOL firstline = TRUE;
void getline(istream& in, char* buffer,const int size, const int delim='\n')
{
BOOL rval;
do {
rval = TRUE;
int n = 0;
char ch;
//cout<<"getline"<<flush;
for(; (ch = in.get())!=delim && !in.rdstate(); ) {
  //cout<<" "<<in.rdstate()<<flush;
  if(n<size-1) buffer[n++] = ch;
  else rval = FALSE;
}
buffer[n] = '\0';
//cout<<" "<<n<<flush;
if(firstline) {
  firstline = FALSE;
  cout<<buffer<<endl;
}
} while(buffer[0]=='%');

if(!rval) cout<<"ERROR line `"<<buffer<<"' exceeds "<<size<<endl;
}//end getline

BOOL inline layout(const char c) { return (c==' '||c=='\t');}

//**************************************************************************
///////////////////////// PROBLEM definition

void gp::AddFall (Function* f )
{
AddF(f);
for ( int i = 0; i < NUM_TREES; i++ )
	addtofuncbag(i, f);
}


const int ConstEvalWt=2;
const int DataEvalWt =ConstEvalWt*250; //AiGP3 sextic should have 250 constants
#define ExtendConstlist(x) assert(x<CONSTARRAYSIZE)

gp::gp() : Problem()
{
	ProbVersion = Version;
	assert ( NUM_TREES >= NUM_OPERATIONS );
	AddFall(new ConstFunc(0));  // Required for GPQUICK, but zero weight

  ifstream inf("prim.dat");
  if (!inf) {
    cout<<"ERROR cannot open prim.dat\n";
    exit(0);
  }

  char buffer[1023+1]; //printing chars + \0
  char copy[1023+1];
  //from sym.cc
  BOOL constgiven = FALSE;

  do {
    getline(inf,buffer,sizeof(buffer));
    if(inf.rdstate()) break; //exit loop on ios::eof etc

    char* ibuff = &buffer[0];
    do {
      for(;(*ibuff!='\0')&&isspace(*ibuff);ibuff++); //skip layout
      char* c = &copy[0];
      for(;(*ibuff!='\0')&&!isspace(*ibuff);ibuff++) *c++=*ibuff; *c++='\0';

if(strcmp(copy,""   )==0) break;
else if(isdigit(copy[0]) || copy[0]=='-' || copy[0]== '.') {
  float constant;
  if(sscanf(copy,"%e",&constant)<=0) {
    cout<<"ERROR "<<copy<<" is not a number\n"; exit(0);
  }
  constgiven = TRUE;
  ExtendConstlist(funccount);
  if(funcconstmin==0) 
    funcconstmin       = funccount; //first constant
  else {
    if(funcconstmax != (funccount-1)) {
      cout<<"ERROR constants not grouped together\n"; exit(0); }
    if(constlist[funcconstmax] > constant) {
      cout<<"ERROR constants must be in ascending order "
	  <<"constlist["<<funcconstmax<<"]="
	  <<constlist[funcconstmax]<<" "<<constant<<endl; exit(0);}
  }
  constlist[funccount] = constant;
  funcconstmax         = funccount;
//cout<<"constlist["<<funccount<<"] "<<constlist[funccount]<<endl;
  AddFall(new Function(0,0,ConstEvalWt,ConstEval, copy));
}
else if(strcmp(copy,"X"   )==0) AddFall(new Function(0,0,DataEvalWt,XEval,"X"));
else if(strcmp(copy,"ADD" )==0) AddFall(new AddFunc(100));
else if(strcmp(copy,"SUB" )==0) AddFall(new SubFunc(100));
else if(strcmp(copy,"MUL" )==0) AddFall(new MulFunc(100));
else if(strcmp(copy,"DIV" )==0) AddFall(new DivFunc(100));
else {cout<<"unknown function "<<copy<<endl;}
} while(TRUE); //end do
} while(TRUE); //end do

	cout<<"\n% OPs ";
	cout<<"("<<funccount<<")";
	{int j;
	for (j=1;j<funccount;j++) {//ignore ConstFunc(0)
	  const int argnum = funclist[j]->argnum;
	  const char* name = funclist[j]->getname();

	  if(j<=20||j+20>funccount) {
	    cout<<name<<" "<<flush;
	    if(j%10==0) cout<<endl;
	  }

	  //if(argnum==0 && *name=='D' && isdigit(name[1])) 
	  //  num_inputs++;
	}
	if(j%10!=1) cout<<endl;
	}

	// If you need to load problem specific data for Evaluation
	// Or load an environment (a map, a scenario)
	// Do it here

	cout<<endl;

	//signal(SIGFPE,handle_SIGFPE);
  //end taken from sym.cc


	// If you need to load problem specific data for Evaluation
	// Or load an environment (a map, a scenario)
	// Do it here

	Evalfunc = new EVALFUNC[funccount];
	for(int i=0;i<funccount;i++){Evalfunc[i] = funclist[i]->evalfunc;}
	Constlist = constlist;
}
/*from sym.cc
BOOL SIGFPE_raised = 0;

void handle_SIGFPE(int x) {
if(SIGFPE_raised == 0) cout<<"handle_SIGFPE("<<x<<") "<<SIGFPE_raised<<endl;
else SIGFPE_raised = 1;
}
//end taken from sym.cc */

//from poly.cc r1.17, keep linker happy
//#include "timing.cc"
extern Chrome* time_chrome;
Chrome* time_chrome;
//end poly.cc

/***********************************************************************/
//inline
retval* evalall(Chrome* __restrict__ chrome, const int my_id/*, const** evalnode* __restrict__ xip*/) {
  assert(my_id>=0 && my_id<nthreads);
  //CHROMEP = chrome;
  node* expr = chrome->expr; //xip;
  node* IP = &expr;
  assert(Evalfunc[3]);
  assert(Constlist[255]<=1.0f);
  assert(IP);
  /*const*/ retval* start = &evalstack[my_id*evalstacksize*evalstackwidth];
  retval* stack = start;
  retval* EvalSP = &stack;
  --IP;
  EVAL;
  assert(EvalSP == start+evalstackwidth); //assert(EvalSP==1);
  //dec_EvalSP;
  return start;
}
#define testout cout
inline
float stub_fitness(Chrome* __restrict__ chrome, const int my_id/*, const** evalnode* __restrict__ xip*/)  // fitness function.
{
#include "avx.fitness.cc"
}
float gp::fitness(Chrome*  __restrict__ chrome/*, **const** evalnode* __restrict__ xip*/)  // fitness function.
{
  const int thread = chrome->thread;
  assert(thread < nthreads);
  const int       my_id = (thread<0)? 0 : thread;
  ///*const*/ evalnode* my_ip = (thread<0)? ExprGlobal : xip;
  const float fit = stub_fitness(chrome,my_id);//,my_ip);

#ifdef INTRON64
#if NUM_TREES != 1
Error here and DisplayZStats assumes one tree
#endif
  const int size = chrome->SubLen(0);
  int exon = 0;
  int nsol = 0;

  for(int i=0;i<size;i++) {
    if(intron64[i]==0) exon++;
    //if(fit64[i]==64)   nsol++;
    if(raw[i]==mux64[0]) nsol++;
  }
  chrome->nfitness->exon = exon;
  chrome->nfitness->nsol = nsol;
#endif /*INTRON64*/
  return fit;
}//end gp::fitness

/* not implemented **
float display_run0(Chrome* chrome, const int gencount) {
  evalnode* xip = IpGlobal;
  assert(xip);
#define DISPLAY
#include "avx.fitness.cc"
#undef DISPLAY
}
**/
float display_run(Chrome* chrome, const int gencount) {
  cout << "\ndisplay_run not implemented\n";
  return -1;
/* not implemented **

  cout << "\nExecuting "<< chrome << flush;
  chrome->write(PRETTY_NONE); cout<<endl;

  chrome->SetupEval();

  display_run0(chrome, gencount);

  //last leaf did not increment IP but needs to be included
  //const int treesize = (IP-ExprGlobal)+1;
  //for(int ip=0;ip<treesize; ip++) printf("%s",intron64[ip]? "i" : "_");
  //printf("#%d intron64\n",treesize);
**/
}
/***********************************************************************/

void test::write_stats(ostream& out) const
{
out<<"Ran " << tests_run << " trees. " << endl;
}//end write_stats

void write_problem_params()
{
//cout<<store_limit<<"["<<Array_bot<<".."<<Array_top<<"] ";
//cout<<mem_penalty_bot<<" "
  cout << "max fitness "<<max_fitness;
#ifdef TIMINGS
//  cout<<" "<<cpu_penalty_bot<<;
#endif
  cout<<endl;
}


// Define your FitnessValue class

void myfitnessvalue::write(ostream& fout)
{	fout << fvalue;
#ifdef MINC_ADF
	fout << " adf1ver " << adf1;
#endif
#ifdef TIMINGS
	fout << " cpu " << cpu_time;
#endif
//	fout << " Same " << Ntests_run;
#ifdef TRACE_RUN
	fout << " (";

	int notrun[NUM_TREES];
	missed(this, notrun);

	for (i=0; i < NUM_TREES; i++)
		{ char buf[80]; 
		  sprintf(buf,"%4d", ok(i)); fout << buf;
		  {if(Tree_OK(this,i))
			   fout<<"K"; else fout<<"-";}
		  sprintf(buf,"%-4d",nak(i)); fout << buf;
		  sprintf(buf,"%-4d",notrun[i]); fout << buf;		  
	        }
	fout << ")";
#endif
}//end myfitnessvalue::write

void cache_init()
{
}

void cache_stats_init()
{
}

void cache_stats()
{
}

//**************************************************************************
///////////////////////// PROBLEM definition

void gp::LoadTests(istream& in)
{ 
}

void gp::WriteTreeName(int tree, ostream& ostr)
{ostr<< "avx";
}//end queue::WriteTreeName()

BOOL gp::TreeNameMatch(int tree, char* s)
{return 0; //any string will do
}//end gp::TreeNameMatch()

int static_check_loop_check = 0;
int gp::static_check(Chrome* chrome, int tree)
// Perform static analysis of new chromsome. Just changed in tree
{
return 0; //Uncomment to disable static_check//ie ok
}//end queue::static_check


void gp::write_stats(ostream& out)
{
/*
#ifndef GP_RUN
cout<<"Frequency of TRUE\ntimes count\n";
{for(int i=0; i<=pow(2,problem_order); i++)
  cout<<i<<" ; "<<trues[i]<<endl;}
cout<<"Frequency of rules\nrule count\n";
{for(int i=0; i<=pow(2,problem_order); i++)
  cout<<i<<" : "<<frequency[i]<<endl;}
#endif /*GP_RUN*/
}

retval random_value(int& seed) //seed must be > 0
{
//return retval(F_tan_spread*tan(float(intrnd(seed) % 3141)/1000.0));
return retval(intrnd(seed)); //default -- not used
}

void myfitnessvalue::clear()
{
//NB clearing whole of myfitnessclass cause gcc produced code to crash

//In theory dont need to do this most of the time since on copy or
//crossover will soon be overwritten with data from parent

//cout<<"Sizeof myfitnessvalue "<<sizeof(*this)<<endl;

}

int mincpu_found; //lowest cpu used by any solution reported so far


int end_gens()
{
	//we never move on to next test phase
	return 0;
}//end_gens


int Chrome::ChooseCrossTree(Chrome* second_parent) {
//not actually in use
return rnd(NUM_TREES);
} //end Chrome::ChooseCrossTree


#ifdef MAXTREEDEPTH
void gp::ProbLPStats(int f[][NUM_TREES][MAXTREEDEPTH+1], int pcount[pcount_size] ) 
#else
void gp::ProbLPStats(int f[][NUM_TREES], int pcount[pcount_size] ) 
#endif /*MAXTREEDEPTH*/
{
pcount[0] += 1;//not used
pcount[1] += 1;//not used
pcount[2] += 1;//not used
pcount[3] += 1;//not used

}//end ProbLPStats

int gp::ProbLTStat(int tree, int start, int length, node* expr )
{
return -1; //ignore
}//end ProbLTStat

