#ifndef IterativeLA_Included
#define IterativeLA_Included

#include "matrix.h"

enum state {SUCCESS, WONT_STOP, BAD_DIAGONAL, BAD_DATA};

state jacobi(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol);

state gauss_seidel(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol);
state SOR(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol, double w);
#endif
