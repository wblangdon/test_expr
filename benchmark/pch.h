// pch.h   Adapted from GPQUICK-2.1
// W. Langdon cs.ucl.ac.uk 5 May 1994 (21 Feb 95 make it inline)
// Use Park-Miller random numbers

// QPQUICK
// Standard header files for compilers supporting Pre-Compiled Headers

//WBL 19 May 2016 for gcc version 4.8.5 of C++
//Remove include stream.h strstream.h iomanip.h
//remove ".h" on iostream.h fstream.h
//Add using namespace std;

#ifndef _PCH_H
#define _PCH_H

// Pick your compiler.  DOS/ANSI by default, now defined for UNIX
#define UNIX

#include <stdlib.h>
#ifdef UNIX
//#include <stream.h>
#include <strstream>
#else
#include <conio.h>
#include <strstrea.h>
#endif

#include <ctype.h>
#include <iostream>
//#include <iomanip.h>
#include <fstream>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <values.h>

using namespace std;

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define BOOL int

//////////  Random number functions cause portability problems.  Resolve here.
// Define	rand_0to1()  (random float from 0 to 1)
//			rnd(x) (random integer from 0 to x-1)
//			rndomize() (initialize with a time based seed) 			
// park-miller.cc   Park-Miller random number generator
// W. Langdon cs.ucl.ac.uk 5 May 1994

inline int intrnd (int& seed) // 1<=seed<=m
{
#ifdef LONG_GE46BITS
int const a    = 16807;      //ie 7**5
int const m    = 2147483647; //ie 2**31-1
	seed = (long(seed * a))%m;
	return seed;
#else
double const a    = 16807;      //ie 7**5
double const m    = 2147483647; //ie 2**31-1

	double temp = seed * a;
	seed = (int) (temp - m * floor ( temp / m ));
	return seed;
#endif
}

#ifdef UNIX
//////////////////// UNIX stuff
#define INT32 int
inline int kbhit() {return FALSE;}  // remove this MS DOS function
//May well work with Turbo C as well
	extern int gpquick_seed;

        inline float rand_0to1()
                { return (float)(intrnd(gpquick_seed)) / 2147483647.0; }

        inline double drand_0to1()
                { return (double)(intrnd(gpquick_seed)) / 2147483647.0; }

	inline int rnd(int __num)
#ifdef BIG_RND_BUGFIX
                { int temp;
		  do { //discard new randon number if % might cause bias
		    temp = intrnd(gpquick_seed);
		  } while((temp/__num)*__num>(2147483647-__num));
		  return (temp % __num);
		}
#else
		{ return (intrnd(gpquick_seed) % __num); }
#endif /*BIG_RND_BUGFIX*/
	inline void rndomize(void)
		{ do {
			gpquick_seed = (unsigned) time(NULL);
		     } while (gpquick_seed <= 0);
                }
#else
//////////////////// DOS/ANSI stuff
#define INT32 long
#define rand_0to1() (((float)rand())/RAND_MAX)
#ifndef __TURBOC__
	inline int rnd(int __num)
		 { return(int)(((long)rand()*__num)/(RAND_MAX+1)); }
	inline void rndomize(void) { srand((unsigned) time(NULL)); }


#else
//////////////////// Borland C++ stuff
#define rnd(x) random(x)
#define rndomize() randomize()
#endif

#endif

#endif
