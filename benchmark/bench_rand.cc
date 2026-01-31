// W. Langdon cs.ucl.ac.uk 9 Sep 2019 how fast in AVX interpretter?
//based on vienna_rna/gp/sse512/test_md.c r1.40
#define main_version "$Revision: 1.1 $"

//Modifications (in reverse order)
//WBL  23 Nov 2019 add timing for random_tree, remove eval etc
//WBL  9 Sep 2019

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include </usr/include/sys/param.h>
#include "tick.h"

#include <string.h>

/*stuff to sub GPQUICK, cf gp.cc r1.128*/

#include "pch.h"
#include "chrome.h"
#include "prob.h"
#include "gp.h"
#include "test.h"
#include "primitiv.h"

#define Ntests 1
node* test_expr[Ntests];
//retval test_ans[Ntests][48];
/*end GPQUICK*/

//retval* evalall1(node* __restrict__  expr);
//retval* evalall1(node* expr);
retval* evalall(Chrome* __restrict__ chrome, const int my_id);

gp* ThisProblem;
ChromeParams* ThisParams;
float max_fitness = 0;
Chrome* ThisChrome;

/*
void stub_fitness_(const int testid, const int verbose,
		      int* ok) {
  assert(testid>=0 && testid<Ntests);
//const retval* sp = evalall1(test_expr[testid]);
//const retval* sp = evalall1(ThisChrome->expr);
  const retval* sp = evalall(ThisChrome,0);
  retval save = 5555.55;
  for(int t=0;t<48;t++) {
    const retval answer = sp[t];
    //const retval error = test_ans[testid][t] - answer;

    //if(fabs(error) <= FLT_EPSILON) *ok=1;
    if(t==0 && answer==0 || answer==1.0f) {
      *ok=1;
      save = answer;
    }
    else if(answer!=save){
      if(verbose)
	printf("failed test %d[%d] evalall1 returned %g expected %g\n",
	       testid,t,answer,save);
//	       testid,t,answer,test_ans[testid][t]);
    }
    tests_run++;
  }
}

long long int do_test(const int testid, const int verbose,
		      int* ok) {
  const long long int t0 = getTickCount();
  stub_fitness_(testid,verbose, ok);  // fitness function.
  const long long int t1 = getTickCount();
  const long long int fitness_time = t1 - t0;
  return fitness_time;
}

/*ignore gcc:
  warning: passing argument 4 of ‘qsort’ from incompatible pointer type*/
int cmp(const long long int* a, const long long int* b){
  return (*a < *b)? -1 : (*a == *b)? 0 : 1;
//int cmp(const void* A, const void* B){
//  const long long int a = *A;
//  const long long int b = *B;
//  return (a < b)? -1 : (a == b)? 0 : 1;
}

/*
#define ntries_ 0
long long int do_testQ(const int testid, const int ntries, const int verbose,
		       int* ok, long long int times[ntries_]) {
  int i;
  for(i=0;i<ntries;i++) {
    times[i] = do_test(testid,verbose, ok);
    if(*ok==0) { //stop on first failure
      if(verbose) printf("do_test(%d,%d,*ok) failed ",testid,ntries,verbose);
      return times[i];
    }
  }
  //sort times, take 1st quartile
  qsort(times,ntries,sizeof(long long int),cmp);
  const int quartile = ntries/4;
  if(verbose) {
  printf("times[%d]%d ",ntries,quartile);
  for(i=0;i<ntries;i++) {
    if(i==quartile) printf(">");
    printf("%llu",times[i]);
    if(i==quartile) printf("<");
    printf(" ");
  }}
  return times[quartile];
}

int fit=0;
long long int total_time  = 0;
long long int total_noise = 0;

long long int absdiff(const long long int a, const long long int b) {
  return (a>=b)? a-b : b-a;
}
void do_testP(const int testid, const int ntries, const int verbose) {
  long long int times1[ntries];
  int ok1;
  ok1 = 0;
  const long long int q1 = do_testQ(testid, ntries, verbose, &ok1,times1);
  if(ok1) fit++;
  if(verbose) printf("\n");
  total_time  += q1;
  if(!ok1) return; //avoid using times1 since much of it will not have been set
  const int quartile = ntries/4;
  const int low  = (quartile)?          quartile - 1 : 0;
  const int high = (quartile<ntries-1)? quartile + 1 : ntries - 1;
 //dont divide by two to avoid dealing with fractions 
  const long long int noise = times1[high]-times1[low];
  total_noise += noise;
  if(verbose) printf("ntries %d noise %lld = times1[%d]-times1[%d] %lld %lld\n",
		     ntries,noise,high,low,times1[high],times1[low]);
  if(verbose) printf("\ntest %d\n",testid);
}
*/
const int max_arity = 4;
int random_tree(const int size, const int arity_count[], node *expr, Problem const *probl);
//extern int evalstacksize;

BOOL rand_tree(const int length, const int chrometree)
{
  Problem* probl = ThisProblem;
  //cout<<"rand_tree "<<length<<" "<<flush;
  //general case takes too long go straight to binary tree answer
  //const int* tac = probl->rand_tree_lookup(length, chrometree);
  //if(tac==NULL) return FALSE;
  const int tac[max_arity] = {0,length/2,0,0};
  //cout<<"tac[";
  //for(int i=1; i<=max_arity; i++) cout<<tac[i-1]<<", ";
  //cout<<"]\n";

  const long long int t0 = getTickCount();
  const int max_depth = random_tree(length,tac,ThisChrome->expr,probl);
  const long long int t1 = getTickCount();
  const long long int delta = t1 - t0;
  cout<<"random_tree "<<length<<" depth "<< max_depth << " took " << delta << " ns\n";

  /*
  if(max_depth > evalstacksize){
    cout<<"ERROR random_tree Max_Depth:"<<max_depth
        <<" exceeds evalstacksize:"<<evalstacksize<<"\n";
    return FALSE;
  }
  */
  return TRUE;
}//end rand_tree;
/*end rand_tree.cc r1.22 */


int ntests;
void do_tests(const int ntries, const int verbose){
/*program with just one leaf
  test_expr[0] = new node[1]; //chrome->expr;
  test_expr[0][0].op   = 255; ///*ConstEval; // Evalfunc[255];
  for(int t=0;t<48;t++) test_ans[0][t] = 0.997; //tail prim.dat Constlist[255];
  even for DIV 0 answer could be 0 or 1 for(int t=0;t<48;t++) test_ans[0][t] =
*/

  const int size = ThisParams->params[pMaxExpr];
  //test_expr[0] = new node[size]; //chrome->expr;
  ThisChrome->ip = 0;
  const BOOL ok = /*ThisChrome->*/rand_tree(size,0);
  if(!ok) return; //cout<<"rand_tree("<<size<<","<<0<<") failed"<<endl; exit(1);}
  /**/
  ThisChrome->tree_starts[0] = 0;
  ThisChrome->tree_starts[1] = size;
  ThisChrome->nfitness->fvalue = 0;
  ThisChrome->write_trace();
  ThisChrome->write(PRETTY_NONE,cout);
  cout<<endl;
  /**/
  /*
  assert(0<Ntests);
  ntests = 0;
  do_testP(0,ntries,verbose);
  ntests++;
  */
}//end do_tests

int main(int argc, char * argv[]) {
  /*from gp.cc r1.128*/
  cout<<"bench_avx" << flush;
  time_t start, end;
  time (&start);
  cout << " ("<<main_version<<") "<<flush;
  cout << "Node " << flush;
  char buffer[MAXHOSTNAMELEN+1];
  gethostname(buffer,sizeof(buffer));
  cout << buffer << endl;
  ThisProblem=new gp();
  cout <<ThisProblem->ProbVersion<<endl;

  ifstream gpParams ("gp.ini");
  if ( gpParams ) { 
    ThisParams = new ChromeParams();//cc.64
    ThisParams->Load(gpParams);
    gpParams.close();
  } else {
    cerr << "Could not open gp.ini\n";
    exit(0);
  }

  for(int i=1; i<argc;i++) {
    ThisParams->Load(argv[i]); }

  gpquick_seed = ThisParams->params[pPopSeed];
  //max_fitness = setmaxfit(ThisParams);
  cout<<"seed "<<gpquick_seed<<" ";
  cout<<"pMaxExpr "<<ThisParams->params[pMaxExpr]<<" "<<flush;

  if(ThisParams->params[pMaxExpr]%2==0) {
    cout<<" len must be odd for binary trees\n";
    return 1;
  }
  //general case takes too long go straight to binary tree answer
  //cout << "Creating random tables up to tree size " << max<<flush;
  //ThisProblem->init_rand_tree(max);
  //cout << " ... done"<<endl;
  /*end gp.cc r1.128*/
  ThisChrome = new Chrome(ThisParams,ThisProblem,
			  ThisProblem->getfuncs(),ThisProblem->getconsts(),
			  false, NULL);

  const int ntries     = 11;//(argc>1)? atoi(argv[1]) : 11;
  const int verbose    = 1;//(argc>2)? atoi(argv[2]) : 0;

  //printf("start tests %d %d\n",ntries,verbose);
  do_tests(ntries,verbose);
  /*fprintf(stdout,"ran %lld tests, passed %d duration %llu %g seconds noise %lld\n",
	  ntests,fit,
	  total_time, 1.0e-9*(double)total_time,
	  total_noise);
  */
  time (&end);
  cout << "\nGP took " << end - start << " secs. Finished at ";
  cout << ctime(&end);
  cout << endl;
  return 0;
}
