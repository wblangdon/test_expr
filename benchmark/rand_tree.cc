// rand_tree.cc   File for uniform random tree generator, with GPQUICK-2.1
// See ETL-TR-95-35 14 Nov 1995 Hitoshi Iba Random tree generator for GP
// W.B.Langdon cs.bham.ac.uk 9 August 1997
#define main_version "$Revision: 1.42 $"

//Modifications (in reverse order)
//WBL 19 Sep 19  Avoid using gtree,stack and lable_tree
//WBL 18 Sep 19  rename as RAND_TREE2_FASTER, rewrite onedom (make faster)
//WBL 15 Sep 19  Add gtree::Maxdepth()
//WBL 14 Sep 19  for g++ 7.3.1
//WBL 13 Sep 19  Add RAND_TREE_FASTER

//WBL  6 Aug 99  BUGFIX

//WBL ...        cc MIPS4 compiler not calling constructors for static data 
//               structues. do explicitly

//WBL 20 May 98  BUGFIX lookup() following email from Christian Igel
//               <igel@neuroinformatik.ruhr-uni-bochum.de>

//WBL  7 Jul 98  Use RndOp

//WBL 12 Jun 98  Allow caller to specify tree shape to rand_tree

//WBL 11 Mar 98  Add code to create tables only one tree size

//WBL 10 Mar 98  Speed up by, eg using static structure for gtree and stack

//WBL 21 Jan 98  Use double precision (should take above 1e07 trees created)
//               remove code to remove trees less than FLT_EPSILON
//               Potentially changes GPQUICK_SEED, previous verison r1.4

//WBL 29 Dec 97  Optionally include number of different terminals and functions
//               of each arity when calculating tree shape probabilties
//               tree_count replaced by version in ntrees.cc r1.1
//               Also ensure scope of for loop variables is limited

//WBL 12 Aug 97  Use log of factorial, add remove i[4] loop

//WBL 11 Aug 97  Extend range of size before overflow
//               BUGFIX limit (efficiency saving only), 

//WBL  1 Aug 97  File created from rand_tree-test.cc

#include <assert.h>
#include "pch.h"
#include "chrome.h"

#define DEBUG true

//dont change!
#define max_arity 4

double* log_fact;
int log_fact_size = 0;
void new_log_fact(const int x) {
  log_fact = new double[x+1];
  log_fact_size = x+1;
  double f = 0;
  log_fact[0] = 0;
  for(int i=1;i<log_fact_size;i++) log_fact[i] = (f += log(double(i)));
}
inline
double log_factorial(const int x) {
  assert (x<log_fact_size);
  return log_fact[x];
}

double tree_count(const int size, int arity_count[], const int arity[]) {

  //actual count of distinct tree shapes with a particular choice of
  //node arities is given by Cn(s) (8). However we are interested in
  //the proportion of the total rather than the absolute number. This
  //espacially important since Cn(s) overflows.

  // value returned is log(actual number of trees)
  // arity_count[0] set to number of terminals

  int residue = size;
  {for(int i=1;i<=max_arity;i++) residue -= arity_count[i];}
  arity_count[0] = residue;

  double trees = log_factorial(size)-log(size);
  {for(int i=0;i<=max_arity;i++) trees -= log_factorial(arity_count[i]);}

  //now include number of ways of labeling tree with given terminal+functions
  {for(int i=0;i<=max_arity;i++) 
    if(arity[i] != 0)  trees += arity_count[i] * log(arity[i]);}

  return trees;
}

int limit(const int size, const int arity[], const int index) {
  return (arity[index]==0)? 0 : size / index;
}

class gtree {
//public:
  int arity;
#ifdef DEBUG
  int max_depth = 1;
#endif
  gtree* child[max_arity];
public:
  gtree() {};
//use static structures for speed
//   inline gtree(const int a, gtree* c[]) {
//     arity = a;
//     for(int i=0; i<max_arity;i++)
//       child[i] = (i<a)? c[i] : NULL;
//   }
//   inline ~gtree() {
//     for(int i=0; i<arity;i++)
//       delete child[i];
//   }
  inline int Arity() const               { return arity; };
#ifdef DEBUG
  inline void setarity(const int c)      { assert(c>=0&&c<=max_arity);
					   arity = c; };
  inline gtree* Child(const int c) const { assert(c>=0&&c<arity);
					   return child[c]; };
  inline void setchild(const int c, /*const cc.64*/ gtree* t) { 
    assert(c>=0&&c<arity);
    child[c] = t;
    const int d = t->max_depth+1;
    if(d > max_depth) max_depth = d;
    /*if(max_depth>11207)
      cout<<"setchild("<<c<<",t "<<t->arity<<":"<<t->max_depth<<") "
    **
      cout<<"setchild "<<c<<" "<<t->arity<<" "<<t->max_depth<<" "
	  <<max_depth<<endl;
    /**/
  };
  inline int Maxdepth() const {return max_depth;};
#else
  inline void setarity(const int c)      { arity = c; };
  inline gtree* Child(const int c) const { return child[c]; };
  inline void setchild(const int c, const gtree* t) { child[c] = t; };
#endif
friend ostream& operator << (ostream& os, const gtree& t);
};//end class gtree

ostream& operator << (ostream& os, const gtree& t) {
  os<<"("<<t.arity;
  for(int i=0; i<t.arity;i++)
    os<<*t.child[i];
  os<<")"<<flush;
  return os;
}

gtree* fulllist    = NULL;
int fulllist_depth = 0;
int fulllist_size  = 0;

gtree* fulltree(const int depth, const int arity, int& fip) {
  if(fulllist_depth<depth) {
    int size = 1; for(int i=0;i<depth;i++) size *=2; size--;
    if(fulllist!=NULL) delete[] fulllist;
    fulllist       = new gtree[size];
    fulllist_depth = depth;
    fulllist_size  = size;
  }
  assert(fip<fulllist_size);  //use static structure
  const int ip = fip++;
  fulllist[ip].setarity((depth>1)? arity:0); //function or terminal?
  for(int a=0; a<fulllist[ip].Arity(); a++) {
//  cout<<"Arity["<<ip<<"] "<<fulllist[ip].Arity()<<endl;
    fulllist[ip].setchild(a,fulltree(depth-1,arity,fip));
  }
  return &fulllist[ip];
}//end fulltree

class onedom {
  int offset;
  int size;
  int* data;
public:
  onedom(const int s, const int in[]) {
    //cf Gambler's ruin and lattice paths Fig 5.12 p269 segdewick:1996:aa
    //An Introduction to the Analysis of Algorithms, Sedgewick and Flajolet
    size = s;
    data = in;
    //find offset such that data is 1-dominated
    //ie map last biggest knee to start
    offset = 0;
    int diag = 1;
    int x = 0;
    int y = 0;
    int i = 0;
    //cout<<"in["<<i<<"]"<<" "<<in[i]<<" "<<y<<endl;
    for(;i<size;i++) {
      if(in[i]==0) x++; else y++;
      const int d = x-y;
      //cout<<"in["<<i<<"]"<<" "<<in[i]<<" "<<x<<" "<<y<<" "<<d<<endl;
      if(diag>=d){ //take last
	diag = d;
	offset = i+1;
	//cout<<"i "<<i<<" diag "<<diag<<" offset "<<offset<<endl;
      }
    }
    //cout<<"offset "<<offset<<endl;
    if(!wellformed()) {cout<<"ERROR not wellformed "<<offset<<endl; exit(1);}
  }
  inline int item(const int index) const {return data[(index+offset)%size];}
  int wellformed() {
    //int x = 0;
    //int y = 0;
    int failed = FALSE;
    int c = 0;
    for(int i=0;i<size;i++) {
      c += 1 - item(i); 
      /*
      if(item(i)==0) x++; else y++;
      const int d = x-y;
      cout<<"item("<<i<<")"<<" "<<item(i)<<" "<<x<<" "<<y<<" "<<d<<endl;
      */
      if(c<=0) failed = TRUE; //return FALSE;
    }
    return (failed==FALSE && c==1);
  }
  //~onedom() { delete[] data;}
  int ip;  //NB step ip backwards
  int recurse(const int max){
    --ip;
    assert(ip>=0 && ip<size);
    const int a = item(ip);
    //cout<<"ip "<<ip<<" a "<<a<<" recurse("<<max<<")"<<endl;
    int mx = max+1;
    for(int i=0;i<a;i++) {
      const int m = recurse(max+1);
      if(m>mx) mx = m;
    }
    return mx;
  }
};//endclass onedom


class stack {
public:
  int sp;
  int size;
  gtree** data;
//gtree* data[EXPRLEN];
public:
//  stack() : sp(0), size(EXPRLEN) {};
//use static structures for speed
   stack(const int s) {
     data = new gtree*[s];
     size = s;
     sp = 0;
   }
#ifdef DEBUG
  ~stack()            {assert(sp==0); delete[] data;}
  gtree* pop()        {assert(sp>0);   return data[--sp];  }
  void push(gtree* t) {assert(sp<size);       data[sp++]=t;}
#else
  inline ~stack()            { delete[] data; }
  inline gtree* pop()        { return data[--sp];  }
  inline void push(gtree* t) { data[sp++]=t; }
#endif
};//endclass stack

gtree* tlist  = NULL;
int tlistsize = 0;
stack* sstack = NULL;

//int tp_max_depth = 0;
//gtree* 
int random_tree(const int size, const int arity_count[], node *expr, Problem const *probl) {
//Given numbers of each function by arity, create a random tree
int* dyck = new int[size];
#ifdef RAND_TREE2_FASTER
//Knuth shuffle https://en.wikipedia.org/wiki/Random_permutation
//for bench_avx.cc only want binary trees
assert(arity_count[0]==0);
assert(arity_count[2]==0);
assert(arity_count[3]==0);
{int i;
  //start with something close to desired random answer
  for(i=0; i<size; i++) dyck[i] = (i%2)? 2:0;
  for(i=0; i<size-1; i++) {
    const int j=i+rnd(size-i);
    //swap(dyck,i,j);
    const int k = dyck[i];
    dyck[i] = dyck[j];
    dyck[j] = k;
  }
}
#else /*RAND_TREE2_FASTER*/
{for(int i=0; i<size; i++) dyck[i] = -1;}

//cout<<" random_tree "<<size<<flush;
//for(i=0; i<max_arity; i++) cout<<" "<<arity_count[i]<<flush;

int u=0;
{for(int i=1; i<=max_arity; i++) {
  for(int n=0; n<arity_count[i-1]; n++) {
    int r = rnd(size-u);
    int j, k;
    for(j=-1, k=0;;k++) {
      if(dyck[k]<=0) j++;
      if(j==r) {
	dyck[k] = i;
	u++;
	//cout<<u<<" "<<r<<" dyck["<<k<<"]="<<i<<endl;
	break;
      }
    }
  }
}}
{for(int i=0; i<size; i++) {if(dyck[i] == -1) dyck[i]=0;}}
#endif /*RAND_TREE2_FASTER*/

//cout<<" Dyck ";
//{for(int i=0; i<size; i++) {
//  cout<<"x";
//  for(int j=0; j<dyck[i]; j++) cout<<"y";
//  cout<<dyck[i];
//}}
//cout<<flush<<endl;

//Before convert Dyck word to tree we need to convert it to
//1-dominated form. cf ETL-TR-95-35 

onedom* d = new onedom(size,dyck);
//onedom* d2 = new onedom(size,dyck); delete d2;

//cout<<" 1-dominated ";
//{for(int i=0; i<size; i++) {
//  cout<<"x";
//  for(int j=0; j<d->item(i); j++) cout<<"y";
//  cout<<d->item(i);
//}}
//cout<<flush<<endl;

#ifdef RAND_TREE2_FASTER
//do random bit sequence backwards as easier 
//but does not give same tree as forward recursive calls.
for(int i=size-1; i>=0; --i) {
  const int ip = size-1-i;
  const int a = d->item(i);
  SETNODE(expr[ip],probl->RndOp(0,a),0);
}

d->ip = size;
const int max_depth = d->recurse(0);
if(d->ip!=0) {cout<<"recurse gives non-zero ip "<<d->ip<<endl;exit(1);}
#else /*RAND_TREE2_FASTER*/

//stack* s = new stack(size);
if(tlistsize<size) {
  if(tlist !=NULL) delete[] tlist;
  if(sstack!=NULL) delete sstack;
  tlist     = new gtree[size];
  sstack    = new stack(size);
  tlistsize = size;
}
{for(int i=0; i<size; i++) {
  tlist[i].setarity(d->item(i)); //use static structure for speed
  for(int j=d->item(i); j>0; --j) tlist[i].setchild(j-1,sstack->pop());
  sstack->push(&tlist[i]);
//   gtree* t[max_arity];          original code
//   //cout<<"Popping "<<d->item(i)<<flush;
//   for(int j=d->item(i); j>0; --j) t[j-1] = s->pop();
//   gtree* n =  new gtree(d->item(i),t);
//   //cout<<" pushing ("<<n->arity<<")"<<flush;
//   s->push(n);
//   cout<<" i "<<i<<" stack "<<sstack->sp<<" deep"<<endl;
}}

//cout<<" tree "<<*final<<endl;
gtree* final = sstack->pop();
#endif /*RAND_TREE2_FASTER*/
#ifdef DEBUG
//tp_max_depth = final->Maxdepth();
cout<<"random_tree("<<size<<",arity_count["
    <<arity_count[0]<<","
    <<arity_count[1]<<","
    <<arity_count[2]<<","
    <<arity_count[3]<<"]) depth:"<<max_depth/*final->Maxdepth()*/<<endl;
#endif

//delete final;
//delete s;
delete d;
delete[] dyck;
//return final;
return max_depth;
}//end random_tree

class tree_chain {
public:
  double count;
  int* tree_arity_counts;
  tree_chain* next;
  tree_chain(const double c, int* const tac, tree_chain* const n):
    count(c), tree_arity_counts(tac), next(n) {};
};

class tree_table {
  int** tac;
  double* tac_index;
  int tac_size;
  double tac_total;
  void couttac(const int j) const {
    cout<<"tac_index["<<j<<"] "<<tac_index[j]<<"\t["<<flush;
    for(int k=0;k<max_arity-1;k++) cout<<tac[j][k]<<","<<flush;
    cout<<tac[j][max_arity-1]<<"]"<<endl;
  }
public:
  tree_table(const int size, const int arity[]) {
//all_arity
//Given tree size and Arities of functions
//return all possible combinations of arities 

  //max_arity = 4;
  double max = -FLT_MAX;
  tree_chain* head = NULL;
//cout<<"\nTrees of size "<<size<<flush;
//  for(int z=1;z<=4;z++) cout<<" "<<arity[z]; cout<<endl;
//cout<<"All possible combinations of arities for tree of size "<<size<<endl;
int i[5];
for(i[1]=0;i[1]<=limit((size-1),                   arity,1);i[1]++)
for(i[2]=0;i[2]<=limit((size-1)-i[1],              arity,2);i[2]++)
for(i[3]=0;i[3]<=limit((size-1)-i[1]-i[2]*2,       arity,3);i[3]++) {
  const int r3 = (size-1)-i[1]-i[2]*2-i[3]*3;
  i[4] = r3/4; 
  if(r3%4==0 && (arity[4]>0 || i[4]==0)) {
    //cout<<"("<<i[1]<<","<<i[2]<<","<<i[3]<<","<<i[4]<<") "<<flush;
    const double tree = tree_count(size,i,arity);
    //cout<<tree<<" "<<endl;
    if(tree>max) max = tree;
    int* savei = new int[max_arity];
    memcpy(savei,&i[1],sizeof(int)*max_arity);
    head = new tree_chain(tree,savei,head);
  }
}
//end all_arity

  //We rescale enties, by normalising we can convert back from logs
  //and discard many entries.  Set normalise count relative to best
  //inividual

    tree_chain* curr;
    double total = 0;
    for(curr=head; curr!=NULL; curr = curr->next) {
      curr->count = exp((curr->count)-max);
      total += curr->count;
    }

  //Normalise a second time so total is unity. We now discard small entries
    tree_chain* next;
    tree_chain* prev = NULL;
    tac_size = 0;
    double running = 0;
    for(curr=head; curr!=NULL; curr = next) {
      next = curr->next;
      const double normalise = curr->count / total;
/*    if(normalise < FLT_EPSILON) {//could be bigger to trim tables more
*	if(curr==head) head=next;
*	delete curr->tree_arity_counts;
*	delete curr;
*	if(prev!=NULL) prev->next = next;
*     }
*     else{ */
	prev = curr;
	tac_size++;
	running += normalise;
	curr->count = running;
/*    } */
    }
    if((running<.9999999||running>1.0000001) && tac_size>0)
      cout<<"ERROR "<<total<<" != "<<running;
    tac_total = running;

  //Now create tables and copy data into them
    tac       = new  int*[tac_size];
    tac_index = new double[tac_size];
    for(int j=0;j<tac_size;j++) {
      tac[j]       = head->tree_arity_counts;
      tac_index[j] = head->count;
      //couttac(j);
      next = head->next;
      delete head;
      head = next;
    }
    assert(head==NULL);
    //cout<<"\t"<<tac_size<<" table entries\n";
  }//end tree_table

  ~tree_table() {
    for(int j=0;j<tac_size;j++) delete[] tac[j];
    delete[] tac;
    delete[] tac_index;
  }

  void lookuperror(const char* s,const float f, const int i) {
    cout<<"ERROR "<<f<<s<<flush;
    couttac(i);
  }//end lookuperror

  int* lookup() {
    if(tac_size<=0) return NULL;
    const double f = drand_0to1()*tac_total;
    //based on time.h
    int bot = 0; //binary chop search
    int top = tac_size - 1;
    //cout<<"choping for "<<f<<" "<<flush;
    do {
      if(bot>top) lookuperror(" bot>top ",f,bot);
      const int i = (top+bot)/2;
      //cout<<tac_index[i]<<" "<<i<<" "<<bot<<" "<<top<<","<<flush;
      if     (i>0 && f < tac_index[i-1] ) top = i - 1;
      else if(       f > tac_index[i]   ) bot = i + 1;
      else {
	//cout<<" found "<<i<<" "<<tac_index[i]<<endl;
	return tac[i];
      }
    } while(TRUE);
  }//end lookup
};//end tree_table


class rand_tree_tables {
  tree_table** table;
  int lower_limit;
  int limit;
public:
  rand_tree_tables(const int maxlength, const int arity[]) {
    lower_limit = 1;
    limit = maxlength;
    table = new tree_table*[limit+1];
    for(int i=1;i<=limit;i++)
      table[i-1] = new tree_table(i,arity);
  }//end rand_tree_tables
  rand_tree_tables(const int arity[], const int len) { //one length only
    lower_limit = len;
    limit = len;
    table = new tree_table*[1];
    table[0] = new tree_table(len,arity);
  }//end rand_tree_tables

  ~rand_tree_tables() {
    for(int i=lower_limit;i<=limit;i++) delete table[i-lower_limit];
    delete[] table;
  }

  int* rand_tree_lookup(const int length) {
    return (length<lower_limit || length>limit)? NULL : 
                                          table[length-lower_limit]->lookup();
  }
};//end rand_tree_tables


//assumes ip set up  //based on SubInit
void Chrome::lable_tree(const gtree* tp, CSelector* const funcbag)
{
//assert(ip>=0&&ip<params->params[pMaxExpr]);
  const int a = tp->Arity();
  int i;

//  do {i=funcbag->roulette();} while (funclist[i]->argnum != a);

//  if (funclist[i]->varnum == 0)
      {SETNODE(expr[ip],probl->RndOp(0,a),0);} //faster!
//  else
//    {SETNODE(expr[ip],i,rnd(funclist[i]->varnum));}
  ip++;
  for(i=0;i<a;i++) lable_tree(tp->Child(i), funcbag);
}//end lable_tree

//assumes ip set up
BOOL Chrome::rand_tree(const int length, const int chrometree, const gtree* tp)
{
  //cout<<"rand_tree "<<length<<" "<<flush;
if(ip+length>params->params[pMaxExpr]) return FALSE; //no room in expr
if(tp!=NULL) 
  lable_tree(tp, probl->funcbag[chrometree]);
else {
  const int* tac = probl->rand_tree_lookup(length, chrometree);
  if(tac==NULL) return FALSE;
#ifndef RAND_TREE2_FASTER
  error needs care if want to use with new RAND_TREE2_FASTER
  WBL 19 Sep 2019

  gtree* tp = random_tree(length,tac);
  lable_tree(tp, probl->funcbag[chrometree]);
//  delete tp; now static structure for speedup
#endif /*RAND_TREE2_FASTER*/
}
  return TRUE;
}//end rand_tree;

//must be called after all functions have been added
//call this OR init_rand_tree(int size_limit) NOT BOTH
//create tables for one length only
void Problem::init_1rand_tree(const int size) {
  new_log_fact(size); 
  assert(NUM_TREES==1);
//for (int t=0;t<NUM_TREES;t++) { -- doesnt handle multiple trees!
    int arity[max_arity+1];
    memset(arity,0,sizeof(arity));
    for (int i=0;i<funcbag[0]->count();i++) //t
      {
	const FHANDLE f = funcbag[0]->GetHandle(i); //t
	const int a = funclist[f]->argnum;
	assert(a<=max_arity);
#ifdef RAND_TREE_ARITIES
	arity[a]++;
#else
	arity[a] = 1;
#endif
      } 
    rand_tree[0] = new rand_tree_tables(arity,size); //t
//  }
}

//must be called after all functions have been added
void Problem::init_rand_tree(int size_limit) {
  if(size_limit<0) size_limit = -size_limit;//use absolute value
  new_log_fact(size_limit);
  for (int t=0;t<NUM_TREES;t++) {
    //time_t now;
    //time (&now);
    //cout<<"Creating Random tree lookup table for tree "<<t<<" at ";
    //cout<<ctime(&now)<<flush;

    int arity[max_arity+1];
    memset(arity,0,sizeof(arity));
    for (int i=0;i<funcbag[t]->count();i++)
      {
	const FHANDLE f = funcbag[t]->GetHandle(i);
	const int a = funclist[f]->argnum;
	assert(a<=max_arity);
#ifdef RAND_TREE_ARITIES
	arity[a]++;
#else
	arity[a] = 1;
#endif
      } 
    //could look for common arity masks between trees, but simplier to
    //assume they are different and so no saving can be made. 
    rand_tree[t] = new rand_tree_tables(size_limit,arity);

    //time (&now);
    //cout<<" "<<size_limit<<" done "<<ctime(&now)<<endl;
  }
}//end init_rand_tree

//cannot be called before init_rand_tree
int* Problem::rand_tree_lookup(const int length, const int t) const
{
  return rand_tree[t]->rand_tree_lookup(length);
}
