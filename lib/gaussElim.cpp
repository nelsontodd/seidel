#include "gaussElim.h"
#include <cmath>

//#define TEST
//#define LARGE_TEST

static void swapRows(Matrix& a, int i, int j) {
  // E_i <--> E_j
  for(int k=0; k<a.n(1); k++) {
    double aa = a(i,k);
    a(i,k) = a(j,k);
    a(j,k) = aa;
  }
}

static void swap(Vector& v, int i, int j) {
  double vv = v(i);
  v(i) = v(j);
  v(j) = vv;
}

static void rowReplacement(Matrix& a, int i, int j) {
  // E_j --> E_j - (a(j,i)/a(i,i))E_i
  double f = a(j,i)/a(i,i);
  a(j,i) = f;
  for(int k=i+1; k<a.n(1); k++) {
    a(j,k) -= f*a(i,k);
  }
}

ge_state luFactorize(Matrix& a, Permutation& p) {
  int n = a.n(0);
  int i,j,k;

  if(a.n(1) != n || p.n() != n) return BADDATA; // error (wrong sizes)

  // Set up permutation
  p.identity();

  // Determine scale factors for scaled partial pivoting
  Vector s(n);
  for(int i=0; i<n; i++) {
    s(i) = fabs(a(i,0));
    for(int j=1; j<n; j++) {
      if( s(i) < fabs(a(i,j)) ) s(i) = fabs(a(i,j));
    }
  }

  // Loop on Columns
  for(j=0; j<n; j++) {

    // Get nonzero pivot (use scaled partial pivoting)
    double pivot = fabs(a(j,j))/s(j);
    i = j;
    for(k=j+1; k<n; k++) {
      double q = fabs(a(k,j))/s(k);
      if(q > pivot) {
	pivot = q;
        i = k;
      }
    }
    if(pivot == 0) return SINGULAR;
    if(i != j) {
      swapRows(a,i,j);
      p.swap(i,j);
      swap(s,i,j);
    }

    // Loop on rows
    for(i=j+1; i<n; i++) rowReplacement(a,j,i);
  }
  return GE_SUCCESS;
}

ge_state luSolve(const Matrix& a, const Permutation& p, Vector& x) {
  int n = a.n(0);
  int i,j,k;

  if(a.n(1) != n || p.n() != n || x.n() != n) return BADDATA;

  // Apply permutation to x
  p.permute(x);

  // FORWARD SUBSTITUTION

  // Loop on columns and rows of a
  for(j=0; j<n; j++) {
    for(i=j+1; i<n; i++) {
      x(i) -= a(i,j)*x(j);
    }
  }

  // BACKWARD SUBSTITUTION

  for(i=n-1; i>=0; i--) {
    for(j=i+1; j<n; j++) {
      x(i) -= a(i,j)*x(j);
    }
    x(i) /= a(i,i);
  }

  return GE_SUCCESS;
}

ge_state solve(Matrix& a, Permutation& p, Vector& x) {
  ge_state s = luFactorize(a,p);
  if(s != GE_SUCCESS) return s;
  return luSolve(a,p,x);
}
