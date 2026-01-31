############################################################################
#
#       Project		evaluate float function using MODIFIED gpquick-2.1,
# 			using Intel AVX-512 C/Unix/etc $Revision: 1.1 $ 
# 			Needs newish GCC
#	    
#       Author		W.B.Langdon
#
#	Created		9 Sep 2019, based upon avx.make Revision: 1.8
# 									   
#	Modifications:  
#
#	Module:         bench_avx.make
#
#	Description:          
#
#	Notes		ensure pch.h contain #define UNIX
#			#define MULTREE
#
############################################################################

NAME 	= bench_rand
OBJS 	= bench_rand.o	avx.o stats.o	chrome.o   rand_tree.o	selector.o   tick.o
LIBS 	= -lm
CC 	= g++
COMPILE = $(CC) -c $(CFLAGS)
CFLAGS 	= -O3 -DNDEBUG
#CFLAGS 	= -g
#CFLAGS = -g3 -Wall -check -pg
CFLAGS 	+= -pthread #for simplicity
CFLAGS 	+= -fpermissive #convert g++ errors to warning
#CFLAGS 	+= -fmax-errors=1
#pchuckle scl enable devtoolset-8 bash; tcsh; make -f bench_avx.make V=1
#try /opt/rh/devtoolset-7/root/usr/bin/g++ (GCC 7.3.1) for avx512f
#SSE     = -mavx512f #-DAVX
#SSE     = -march=skylake-avx512
#SSE     += -mvzeroupper
#SSE	+= -DAVX
#SSE	= -march=native -mtune=native
LINK 	= $(CC) -o 
LDFLAGS = -pthread -lpthread
#LDFLAGS += -g


$(NAME):$(OBJS) 
	$(LINK) $(NAME) $(LDFLAGS) $(OBJS) $(LIBS)

#C++ source files

bench_rand.o:	bench_rand.cc tick.h pch.h chrome.h prob.h avx.h gp.h test.h
	$(COMPILE) bench_rand.cc

#no AVX so can run on eden
avx.o:	avx.cc pch.h chrome.h selector.h avx.h primitiv.h test.h avx.fitness.cc
	$(COMPILE) $(SSE) avx.cc

stats.o:	stats.cc pch.h chrome.h selector.h avx.h gp.h
	$(COMPILE) stats.cc

chrome.o:	chrome.cxx pch.h chrome.h selector.h avx.h
	$(COMPILE) chrome.cxx

selector.o:	selector.cxx pch.h selector.h
	$(COMPILE) selector.cxx

rand_tree.o:	rand_tree.cc pch.h chrome.h selector.h avx.h
	$(COMPILE) rand_tree.cc -DRAND_TREE2_FASTER

tick.o:	tick.c tick.h
	$(COMPILE) tick.c

