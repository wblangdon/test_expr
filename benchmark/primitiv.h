// primitiv.h   From GPQUICK-2.1

// W. Langdon cs.ucl.ac.uk
// For questions or upgrades contact:
// Andy Singleton, Creation Mechanics Inc.
// PO Box 248, Peterborough, NH 03458
// Internet: p00396@psilink.com
// Compuserve 73313,757
// Phone: (603) 563-7757


//Modification (reverse order)
//  3 Mar 1998   Restore Not, And, Or
// 30 Nov 1995   Add SUBR and DIVR
// 26 Oct 1995   Restore Mul, Div, Abs
// 25 Jan 1995   Restore Iflte

// 13 May 1994   Comment out unused classes (since they mess up gcc
//               without optimisation switch!)


//**************************************************************************
//  Standard Functions
//	includes ConstFunc (numeric constants)
//		ADD,SUB,MUL,DIV,SINE,ABS,ADD4,PROG4,BOOL,NOT,AND,OR,IF,IFLTE
//
// They are packaged as objects for easy initialization
//
// Evaluation functions are separated from FUNCTION objects for speed reasons.
// Each evaluation function goes with a FUNCTION object.
//
// Use the OPDEF,EVAL,IP and CURCHROME defines
// and the code will recompile for different eval methods
// Put problem specific functions in the <PROBLEM>.CPP file
// **********************************************************************

// ConstFunc
// REQUIRED function for GPQuick (required as the first function in a problem)
// Return a numeric constant from -128 to +127
OPDEF(ConstEval);
class ConstFunc : public Function {
public:
	ConstFunc(int w=300) {strcpy(name,"NUMBER");argnum = 0;varnum=256;weight=w;
	evalfunc=ConstEval;};
};

OPDEF(AddEval);	// Add two arguments
class AddFunc : public Function {
public:
	AddFunc(int w=100) {strcpy(name,"ADD");argnum=2;varnum=0;weight=w;evalfunc=AddEval;};
};

OPDEF(SubEval);	// Subtract two arguments
class SubFunc : public Function {
public:
	SubFunc(int w=100) {strcpy(name,"SUB");argnum=2;varnum=0;weight=w;evalfunc=SubEval;};
};

OPDEF(SubREval);	// Subtract two arguments (order reversed)
class SubRFunc : public Function {
public:
	SubRFunc(int w=100) {strcpy(name,"SUBR");argnum=2;varnum=0;weight=w;evalfunc=SubREval;};
};

OPDEF(MulEval);	// Multiply two arguments
class MulFunc : public Function {
public:
	MulFunc(int w=100) {strcpy(name,"MUL");argnum=2;varnum=0;weight=w;evalfunc=MulEval;};
};

OPDEF(DivEval);	// "Protected" division
class DivFunc : public Function {
public:
	DivFunc(int w=100) {strcpy(name,"DIV");argnum=2;varnum=0;weight=w;evalfunc=DivEval;};
};

OPDEF(DivREval);	// "Protected" division (order reversed)
class DivRFunc : public Function {
public:
	DivRFunc(int w=100) {strcpy(name,"DIVR");argnum=2;varnum=0;weight=w;evalfunc=DivREval;};
};

/***********************************************************************
OPDEF(SineEval);	// Return the Sine of one argument
class SineFunc : public Function {
public:
	SineFunc(int w=100) {strcpy(name,"SINE");argnum=1;varnum=0;weight=w;evalfunc=SineEval;};
};
**********************************************************************/

OPDEF(AbsEval);	// Return the absolute value
class AbsFunc : public Function {
public:
	AbsFunc(int w=100) {strcpy(name,"ABS");argnum=1;varnum=0;weight=w;evalfunc=AbsEval;};
};

/***********************************************************************
OPDEF(Prog4Eval);		// Do 4 in a row.  Return the last
class Prog4Func : public Function {
public:
	Prog4Func(int w=100) {strcpy(name,"PROG4");argnum=4;varnum=0;weight=w;evalfunc=Prog4Eval;};
};

OPDEF(Add4Eval);		// Add four arguments.  Useful for neural type behavior
class Add4Func : public Function {
public:
	Add4Func(int w=100) {strcpy(name,"ADD4");argnum=4;varnum=0;weight=w;evalfunc=Add4Eval;};
};
**********************************************************************/

OPDEF(IfEval);		// If(condition>0, dothis, otherwise dothat)
class IfFunc : public Function {
public:
	IfFunc(int w = 100) {strcpy(name,"IF");argnum=3;varnum=0;weight=w;evalfunc=IfEval;};
};


OPDEF(IflteEval);	//IfLTE(condition1islessthan,condition2,dothis,dothat)
					// This is the Koza conditional
class IflteFunc : public Function {
public:
	IflteFunc(int w = 100) {strcpy(name,"IFLTE");argnum=4;varnum=0;weight=w;evalfunc=IflteEval;};
};

/**********************************************************************
// **********************************************************************
// Boolean pack  BOOL, AND, OR, NOT
// for these functions, 0 or less is FALSE, greater than 0 is TRUE

OPDEF(BoolEval);	// Convert number to boolean value 0 or 1 - BOOL(arg)
class BoolFunc : public Function {
public:
	BoolFunc(int w = 100) {strcpy(name,"BOOL");argnum=1;varnum=0;weight=w;evalfunc=BoolEval;};
};
***********************************************************************/

OPDEF(NotEval);		// NOT(arg)
class NotFunc : public Function {
public:
	NotFunc(int w = 100) {strcpy(name,"NOT");argnum=1;varnum=0;weight=w;evalfunc=NotEval;};
};

OPDEF(AndEval);		// AND(arg1,arg2)
class AndFunc : public Function {
public:
	AndFunc(int w = 100) {strcpy(name,"AND");argnum=2;varnum=0;weight=w;evalfunc=AndEval;};
};


OPDEF(OrEval);		// OR(arg1,arg2)
class OrFunc : public Function {
public:
	OrFunc(int w = 100) {strcpy(name,"OR");argnum=2;varnum=0;weight=w;evalfunc=OrEval;};
};
