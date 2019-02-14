#include <iostream>
#include <cmath>
using namespace std;

#include "iterativeLA.h"
#include "matrix.h"

#define MONITOR

state jacobi(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol) {

  cout << "Using method JACOBI" << endl;

  // CHECK DATA

  int n = A.n(0);
  if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA; //size of vector and matrix dont match
  if(tol <= 0) return BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n; i++) {
    if(A(i,i) == 0) return BAD_DIAGONAL;
  }

  // APPLY JACOBI

  Vector xOld(x);

  for(int iter=0; iter<maxIter; iter++) {

    // Get new x
    for(int i=0; i<n; i++) {
      double sum = 0;
      for(int j=0; j<n; j++) {
	if(j==i) continue;
	sum += A(i,j)*xOld(j);
      }
      x(i) = ( -sum + b(i) ) / A(i,i);
    }

    // Check error tolerance
    xOld -= x;
    double l2error = l2norm(xOld) / (l2norm(x)+1e-16);
#ifdef MONITOR
    cout << "Iter " << iter+1 << ", l2-error " << l2error << endl;
#endif    
    if( l2error <= tol) {
      maxIter = iter+1;
      return SUCCESS;
    }
    xOld = x;
  }

  return WONT_STOP;
}



state gauss_seidel(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol) {

  cout << "Using method GAUSS SEIDEL" << endl;

  // CHECK DATA

  int n = A.n(0);
  if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA; //size of vector and matrix dont match
  if(tol <= 0) return BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n; i++) {
    if(A(i,i) == 0) return BAD_DIAGONAL;
  }

  // APPLY GAUSS-SEIDEL

  Vector xOld(x);

  for(int iter=0; iter<maxIter; iter++) {

    // Get new x
    for(int i=0; i<n; i++) {
      double sum = 0;
      for(int j=0; j<n; j++) {
	if(j==i) continue;
	sum += sum + A(i,j)*xOld(j);
      }
      x(i) = ( -sum + b(i) ) / A(i,i);
    }

    // Check error tolerance
    xOld -= x;
    double l2error = l2norm(xOld) / (l2norm(x)+1e-16);
#ifdef MONITOR
    cout << "Iter " << iter+1 << ", l2-error " << l2error << endl;
#endif    
    if( l2error <= tol) {
      maxIter = iter+1;
      return SUCCESS;
    }
    xOld = x;
  }

  return WONT_STOP;
}

state SOR(const Matrix& A, const Vector& b, Vector& x,
	     int& maxIter, double tol, double w) {
  cout << "Using method SOR" << endl;

  // CHECK DATA

  int n = A.n(0);
  if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA; //size of vector and matrix dont match
  if(tol <= 0) return BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n; i++) {
    if(A(i,i) == 0) return BAD_DIAGONAL;
  }

  // APPLY SOR

  Vector xOld(x);

  for(int iter=0; iter<maxIter; iter++) {

    // Get new x
    for(int i=0; i<n; i++) {
      double sum = 0;
      for(int j=0; j<n; j++) {
	if(j==i) continue;
	sum += sum + A(i,j)*xOld(j);
      }
      x(i) = (1 - w)*x(i) + w*( -sum + b(i) ) / A(i,i);
    }

    // Check error tolerance
    xOld -= x;
    double l2error = l2norm(xOld) / (l2norm(x)+1e-16);
#ifdef MONITOR
    cout << "Iter " << iter+1 << ", l2-error " << l2error << endl;
#endif    
    if( l2error <= tol) {
      maxIter = iter+1;
      return SUCCESS;
    }
    xOld = x;
  }

  return WONT_STOP;
}


state cgd(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol, string preconditioner="Jacobi") {
  // CHECK DATA

  int n = A.n(0);

  if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
  if(tol <= 0) return BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n-1; i++) {
    for(int j=i+1; j<n; j++) {
      if(A(i,j) != A(j,i) ) return NOT_SYMMETRIC;
    }
  }  

  // INITIALIZE CONJUGATE GRADIENTS

  Matrix Minv(n, n);
	if (preconditioner == "jacobi") {
		for(int i=0; i<n; i++) {
			Minv(i,i) = 1.0/A(i,i);
		}  
	}
	if (preconditioner == "gauss-seidel") {
		for(int i=0; i<n; i++) {
			Minv(i,i) = 1.0/A(i,i);
		}  
	}
  cout << "Matrix " << Minv << "\n" << endl;
  // Set initial residual r = b - Ax
  Vector r(n);
  cout << "Vector: " << r << "\n" << endl;

  matVecMult(A,x,r);
  r-=b; r*=(-1);

  // Set initial search direction d = r
  Vector d(n); // Creates d and sets d = M^-1*r
  matVecMult(Minv,r,d);
  Vector z(d); // Creates z and sets z = d
  double alpha = scDot(r,z);

  double tolSq = tol*tol;

  // CONJUGATE GRADIENT LOOP

  for(int iter=0; iter<maxIter; iter++) {

    if(scDot(d,d) <= tolSq) {
      maxIter = iter+1;
      return SUCCESS;
    }

    // Set u = Ad
    Vector u(n);

    matVecMult(A,d,u);

    // Update x = x + td and r = r - tu
    double t =  alpha / scDot(d,u);

    // Get new search direction d = r + s*d;
    for(int i=0; i<n; i++) {
      x(i) += t*d(i);
      r(i) -= t*u(i);
    }
		matVecMult(Minv, r, z);
    double beta = scDot(r,z);

#ifdef MONITOR
    cout << "Iter " << iter+1 << ", approx. solution: " 
         << x << ", res l2-error " << sqrt(beta) << endl;
#endif    

    if(beta <= tolSq) {
      maxIter = iter+1;
      return SUCCESS;
    }

    double s = beta / alpha;

    for(int i=0; i<n; i++) {
      d(i) = z(i) + s*d(i);
    }

    alpha = beta;
  }

  return WONT_STOP;
}

state cgd_no_conditioning(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol) {
  // CHECK DATA

  int n = A.n(0);

  if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
  if(tol <= 0) return BAD_DATA;
  if(maxIter <= 0) maxIter = 1;

  for(int i=0; i<n-1; i++) {
    for(int j=i+1; j<n; j++) {
      if(A(i,j) != A(j,i) ) return NOT_SYMMETRIC;
    }
  }  

  // INITIALIZE CONJUGATE GRADIENTS

  // Set initial residual r = b - Ax
  Vector r(n);

  matVecMult(A,x,r);
  r-=b; r*=(-1);

  double alpha = scDot(r,r);

  // Set initial search direction d = r
  Vector d(r); // Creates d and sets d = r

  double tolSq = tol*tol;

  // CONJUGATE GRADIENT LOOP

  for(int iter=0; iter<maxIter; iter++) {

    if(scDot(d,d) <= tolSq) {
      maxIter = iter+1;
      return SUCCESS;
    }

    // Set u = Ad
    Vector u(n);

    matVecMult(A,d,u);

    // Update x = x + td and r = r - tu
    double t = alpha / scDot(d,u);

    for(int i=0; i<n; i++) {
      x(i) += t*d(i);
      r(i) -= t*u(i);
    }

    // Get new search direction d = r + s*d;
    double beta = scDot(r,r);

#ifdef MONITOR
    cout << "Iter " << iter+1 << ", approx. solution: " 
         << x << ", res l2-error " << sqrt(beta) << endl;
#endif    

    if(beta <= tolSq) {
      maxIter = iter+1;
      return SUCCESS;
    }

    double s = beta / alpha;

    for(int i=0; i<n; i++) {
      d(i) = r(i) + s*d(i);
    }

    alpha = beta;
  }

  return WONT_STOP;
}
