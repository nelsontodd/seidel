#ifndef IterativeLA_Included
#define IterativeLA_Included

#include "matrix.h"

enum state {SUCCESS, WONT_STOP, BAD_DIAGONAL, BAD_DATA, NOT_SYMMETRIC};

state jacobi(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol);

state gauss_seidel(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol);
state SOR(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol, double w);

state cgd(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol, string preconditioner);
state cgd_no_conditioning(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol);

#endif
