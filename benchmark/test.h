// test.h   Include file for test.cc etc
// W. Langdon cs.ucl.ac.uk 30 Jan 1995

// version "$Revision: 1.32 $"

//Modifications (in reverse order)
//WBL 18 Oct 1995  Support rand_test

//WBL 17 Oct 1995  Support hurdle

//WBL 24 Sep 1995  STORE_FIT: always using 6 args on update_score

//WBL 27 Jun 1995  Support testing in multiple phases

//WBL 10 Jun 1995  Support gcc 2.5.7 by moving dynamic_data to test.run.cc

//WBL 27 May 1995  Suport Op_Ok

//WBL 22 May 1995  Remove lessthanperfect

//WBL 15 Feb 1995  Add update_passed, set_passed

//WBL 11 Feb 1995  Add stop_mask, display_tests(), max_passes

//WBL 20 Jan 1995  New file

#include <limits.h>


//statistical class data kept separate so can be updated by const functions
//and loaded from dumpfile

extern int current_phase;
extern int tests_run;
extern int tests_apparently_run;
extern int total_cpu;

class test {
//list data dependant (shouldnt be here but...) Used by OK and Unpack_hits
	scoretype Op_Max[NUM_TEST_PHASES][NUM_OPERATIONS]; 
        unsigned char Op_Ok[NUM_OPERATIONS][NUM_TEST_SEQUENCES];
	struct test_struct {
		int op; //0..255
		int res; //0..9 
		retval arg1;//Also print buffer iff prtlist
		retval arg2;
		retval* print;
	};
#ifndef GRID_LIB
        int stop_mask;// = 0;
#endif
        int rand_test_min_;// =0;
        int rand_test_max_fails_;
	int test_size;
	int max_passes[NUM_TEST_PHASES];
	test_struct* tests;

	retval random_data[10];
	retval unique_data[10];
//const   retval reserved_value_min = -505059;
//const   retval reserved_value_max = -505050;
enum {reserved_value_min = -505059,reserved_value_max = -505050};
	int seed;// = 1;
//	int tests_run = 0;

void	load_random();
inline retval  Data(retval data, retval dynamic_data[10+1]) const;
retval  read_data(istream&, int test_report_num);
inline retval  get_reserved_value(istream&, int test_report_num);
BOOL    reserved_value(retval data) const;
void    write_data(retval data, ostream& fout) const;
void	write_test(int test, Problem* probl, ostream& fout = cout) const;

enum { scoreq =     255};
enum { scorep =     254};
enum { break_test = 253};
enum { end_test =   252};
enum { start_test = 251};
enum { scorneq =    250};
enum { scoreq_ =    249};
enum { scorneq_ =   248};
enum { hurdle =     247};

void    set_max_fitness() const; //problem dependant
int     set_passed(BOOL passed[], Chrome* chrome, int first_test) const;

public:
#ifdef GRID_LIB
const int stop_mask;
#endif
int     num_phases;
int     sequence_start[NUM_TEST_SEQUENCES+1]; //+1 for end of last seq
int     phase_end_seq[NUM_TEST_PHASES];
int*    break_index;   //0 => not a break instruction
int     test_seed;
inline int rand_test_min() const {return rand_test_min_;};// = 0;
inline int rand_test_max_fails() const {return rand_test_max_fails_;};

	test(Problem* probl, istream& istr = cin);
	~test();
void	write(Problem* probl, ostream& fout = cout) const;
int     run_tests(Chrome* chrome, int test_sequence, int first_test) const;
int     display_tests(Chrome* chrome, int sequence, int first_test) const;
void    write_stats(ostream& fout = cout) const;
int     num_tests_run() const;
#ifdef TEST_STATS
float   sequence_fitness(int sequence, int last_test_run = 0) const;
void    sequence_stats(int sequence, int passes);
void    sequence_passed_stats(int s, int last_test_passed) const;
#endif
//problem specific

void    update_passed(BOOL passed[], int first_test, int last_test_plus_one)
	const;

inline void 
	update_score(Chrome* chrome, int test_sequence,
		     int first_test, int last_test_plus_one, scoretype score,
		     BOOL passed) const;
void sub_sequence(int start, int end, scoretype score, scoretype hits[]) const;
void ok_sub_sequences(int start, int end, scoretype hits[]) const;
#ifdef STORE_FIT

void update_seq_score(Chrome* chrome, int test_sequence, int seqpasses, 
		      int sequence_tests) const;
float   final_score (Chrome* chrome) const;

void not_run_seq(void* fitptr, int sequence, int count[NUM_OPERATIONS]) const;
void not_run(void* fitptr, int count[NUM_OPERATIONS]) const;
void write_not_run(void* fitptr, int sequence, ostream& os ) const;
BOOL Operation_OK(void* fitptr, int operation) const;//reads myfitnessvalue
void unpack_hits(Chrome*) const;//manipulates myfitnessvalue
#else
void update_seq_score(Chrome* chrome, int sequence, int passes)const;
float   final_score (Chrome* chrome, int total_passes, int tests)const;
#endif //STORE_FIT
};
