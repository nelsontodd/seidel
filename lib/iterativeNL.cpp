
///////////////////////////////////////////////////////////////////////////////
// FIXED POINT SUCCESSIVE SUBSTITUTION METHOD
//
// state fixedpt(void g(const Vector&,Vector&), Vector& x,
//               double tol, int& maxIter);
// Inputs:
//   g              The name of the function for which a fixed point is sought.
//   x              The initial guess at the solution.
//   tolerance      The convergence tolerance (must be > 0).
//   iter           The maximum number of iterations that can be taken.
// Outputs:
//   x              The solution.
//   iter           The number of iterations used.
// Return:
//   state          An error status code.
//     SUCCESS      Sucessful termination.
//     WONT_STOP    Error: Exceeded maximum number of iterations.
//     BAD_DATA     Error: The tolerance is <= 0
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
using namespace std;
#include "matrix.h"
#include "gaussElim.h"
#include "iterativeNL.h"

#define MONITOR 
static int prec;

state fixedpt(void g(const Vector&,Vector&), Vector& x,
	      double tol, int& iter) {
  if(iter < 1) iter = 1;
  if(tol <= 0) return BAD_DATA;

  int maxIter = iter;
  Vector y(x);
  double error;
  for(iter = 1; iter <= maxIter; iter++) {
		//Vector y is effectively x_0, and x is x_1. So y acts like a temp variable, and for error.
		g(y,x); // x = g(y) is the new guess      
    y -= x;      // y is now the change to x, this is used to calculate the error
    error = maxNorm(y);
#ifdef MONITOR
    cout << "Iter " << iter << ": x= " << x << ", err = " << error << endl;
#endif    
    if (error <= tol) return SUCCESS;
    y = x;
  }
  return WONT_STOP;
}

state newton(void f(const Vector&,Vector&), void df(const Vector&,Matrix&), Vector& x,
		double tol, int& iter) {
  if(iter < 1) iter = 1;
  if(tol <= 0) return BAD_DATA;

  int i;
	double error;
  int maxIter = iter;
  Vector y(x);
	Vector s(x);
  Permutation p(2);
  Matrix J(2,2);
#ifdef MONITOR
    cout.precision(12);
		cout << "Iter " << 0 << ": x= " << x << ", err = " << error << endl;
#endif    
	for(iter = 1; iter <= maxIter; iter++) {
		f(x,y);
		df(x,J);
		for(i = 0; i < y.n(); i++) {
      y(i) = -1*y(i);
		}
		s = y;
		ge_state solved = solve(J,p,s);
		for(i = 0; i < x.n(); i++) {
      x(i) = x(i) + s(i);
		}
    error = maxNorm(y);
#ifdef MONITOR
    cout.precision(12);
		cout << "Iter " << iter << ": x= " << x << ", err = " << error << endl;
#endif    
    if (error <= tol) return SUCCESS;
	}
  return WONT_STOP;
}
