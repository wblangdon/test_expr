// poly.fitness.cc  INCLUDE file
// W.Langdon@cs.ucl.ac.uk $Revision: 1.15 $

//Modifications (reverse order):
//WBL  8 Jan 2019 join AVX and PTHREADS
//WBL 14 Dec 2018 add AVX
//WBL 23 Nov 2018 based on poly.fitness.cc r1.5

  float f = 0;
  int hits = 0;
//int EvalSP = 0;
  const retval* sp = evalall(chrome, my_id);//,xip);
#ifdef AVX_xxxx
for compatibility between with and without avx do not do final calculation in parallel
  for(int t=0;t<MAXTESTCASES;t+=16) {
    const __m512 answer = _mm512_load_ps(&evalstack[my_stackoff+t]);
    const __m512 y      = _mm512_load_ps(&dataY[t]);
    const __m512 error  = _mm512_sub_ps(y,answer);
    const __m512 abs    = _mm512_abs_ps(error);
    const  float ff     = _mm512_reduce_add_ps(abs);
    f -= ff;
    for(int i=0;i<16;i++) {if(abs[i] < 0.01) hits++;}
  }
#ifdef DISPLAY
    cout<<"X"<<" = "<<"answer"<<" "<<"error"<<" "<<"sqr"<<" "<<hits<<endl;
#endif
tests_run += MAXTESTCASES;
#else
  for(int t=0;t<MAXTESTCASES;t++) {
    const retval answer = sp[t]; //evalstack[my_stackoff+t];
    const retval error = dataY[t] - answer;
    if(error>0){
      f -= error; //negate error, so bigger fitness are still better
      if(error < 0.01) hits++;
    }
    else {
      f += error;
      if(error > -0.01) hits++;
    }
#ifdef DISPLAY
    const retval sqr   = error*error;
//    f -= sqr;
    cout<</*X<<" = "<<*/answer<<" "<<error<<" "<<sqr<<" "<<hits<<endl;
#else
    tests_run++;
#endif
  }
#endif /*AVX*/
#ifdef DISPLAY
  float ff = 0;
  int hhits = 0;
  float save = dataX[0];     //only use first of MAXTESTCASES
  const int resolution = 100;
  {for(int t=0; t<2*resolution+1;t++) { //-1.0 .. 1.0
    const float X=-1.0+t/100.0;
    dataX[0] = X;
    //EvalSP = 0;
    assert(0); //not debugged yet
    const retval* sp = evalall(chrome,my_id);//,xip);
    const retval answer = sp[0];//evalstack[my_stackoff+0];
    const retval error = sextic(X) - answer;
    const retval sqr   = error*error;
//    ff -= sqr;
    if(error>0){
      ff -= error; //negate error, so bigger fitness are still better
      if(error < 0.01) hhits++;
    }
    else {
      ff += error;
      if(error > -0.01) hhits++;
    }
    testout<<gencount<<" "<<X<<" "<<answer<<" ";
    if(error>1)       testout<<"1"<<flush;  //trancate to make nice for gnuplot
    else if(error<-1) testout<<"-1"<<flush;
    else              testout<<answer<<flush;
    testout<<" "<<error<<" "<<sqr<<" "<<hhits<<endl;
  }}
  dataX[0] = save;
  //testout<<gencount<<" "<<ff/(2.0*resolution+1.0)<<endl;
  testout<<endl;
#else
  ptrmyfitnessvalue(chrome->nfitness)->Ntests_run = hits;
  chrome->nfitness->fvalue = f/MAXTESTCASES;
#endif
//cout<<"chrome->nfitness->fvalue="<<chrome->nfitness->fvalue<<" f="<<f<<" hits="<<hits<<endl;
  return chrome->nfitness->fvalue;

//end poly.test

