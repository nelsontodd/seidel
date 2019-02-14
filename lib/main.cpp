#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "iterativeLA.h"
#include "iterativeNL.h"
#include "gaussElim.h"

void Fn(const Vector& x, Vector& y) {
  y(0) = 6*pow(x(0),3) + x(0)*x(1) - 3*pow(x(1),3) - 4;
  y(1) = pow(x(0),2) - 18*x(0)*pow(x(1),2) + 16*pow(x(1),3) + 1;
}

void Fprime(const Vector& x, Matrix& J) {
  J(0,0) = 18*pow(x(0),2) + x(1);
  J(0,1) = x(0) - 9*pow(x(1),2);
  J(1,0) = 2*x(0) - 18*pow(x(1),2);
  J(1,1) = -36*x(0)*x(1) + 48*pow(x(1),2);
}

int main() {
  
  state method_state;
  string method_name;
  int n;
  ifstream nfile("config/size.cfg");
  nfile >> n;
  ifstream Afile("config/matrix.cfg");
  Matrix A(n,n);
  Afile >> A;
  ifstream bfile("config/vector.cfg");
  Vector b(n);
  bfile >> b;

  Vector x(n);
  x = 0;
  ifstream itertolfile("config/maxiterandtolerance.cfg");
  int maxIter;
  double tolerance;
  itertolfile >> maxIter >> tolerance;

  if (std::getenv("METHOD"))
     method_name = std::getenv("METHOD");
  if (method_name == "GAUSS-SEIDEL")
    method_state = gauss_seidel(A,b,x,maxIter,tolerance);
  else if (method_name == "SOR")
    method_state = SOR(A,b,x,maxIter,tolerance,1.2);
  else if (method_name == "CGD") {
    method_state = cgd(A,b,x,maxIter,tolerance, "jacobi");
    cout << "Method state: " << method_name << endl;
  }
  else if (method_name == "CGD-NO-PRECONDITION") {
    method_state = cgd_no_conditioning(A,b,x,maxIter,tolerance);
    cout << "Method state: " << method_name << endl;
  }
  else if (method_name == "FIXED-POINT") {
    x = b;
    method_state = fixedpt(Fn, x, tolerance, maxIter);
    cout << "Method state: " << method_name << endl;
  }
  else if (method_name == "NEWTON") {
    x = b;
    method_state = newton(Fn, Fprime, x, tolerance, maxIter);
    cout << "Method state: " << method_name << endl;
  }
  else {
    method_name = "JACOBI";
    method_state = jacobi(A,b,x,maxIter,tolerance);
    cout << "WARN: METHOD VAR MUST BE IN: [JACOBI, GAUSS-SEIDEL, SOR, CGD, NEWTON, FIXED-POINT]" << endl;
  }


  switch(method_state) {
  case WONT_STOP:
    cout << "ERROR: Exceeded maximum number of iterations." << endl;
    return 1;
  case BAD_DIAGONAL:
    cout << "ERROR: A diagonal entry of A was 0." << endl;
    return 1;
  case NOT_SYMMETRIC:
    cout << "ERROR: A was not symmetric." << endl;
    return 1;
  default:
    cout << "ERROR: Unspecified." << endl;
    return 1;
  case SUCCESS:
    cout << "The solution is:" << endl;
    cout << x << endl;

    if (method_name == "JACOBI" || method_name == "GAUSS-SEIDEL" || method_name == "SOR" || method_name == "CGD") {
    Vector y(n);
    matVecMult(A,x,y);
    y -= b;
    cout << "The number of iterations is: " << maxIter << endl;
    cout << "The max-norm of residual is: " << maxNorm(y) << endl;
    cout << "The residual is: " << endl;
    cout << y << endl;
    }
    return 0;
  }
}

