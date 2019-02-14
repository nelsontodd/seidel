#ifndef IterativeNL_Included
#define IterativeNL_Included

#include "matrix.h"
#include "iterativeLA.h" //For state definition

state fixedpt(void g(const Vector&,Vector&), Vector& x,
	      double tol, int& iter);
state newton(void f(const Vector&,Vector&), void df(const Vector&,Matrix&), Vector& x, double tol, int& iter);

#endif
