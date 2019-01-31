#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

#include "iterativeLA.h"

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
  if (method_name == "GAUSS SEIDEL")
    method_state = gauss_seidel(A,b,x,maxIter,tolerance);
  else if (method_name == "SOR")
    method_state = SOR(A,b,x,maxIter,tolerance,1.2);
  else {
    method_name = "JACOBI";
    method_state = jacobi(A,b,x,maxIter,tolerance);
    cout << "WARN: METHOD VAR MUST BE IN: [JACOBI, GAUSS SEIDEL, SOR]" << endl;
  }


  switch(method_state) {
  case WONT_STOP:
    cout << "ERROR: Exceeded maximum number of iterations." << endl;
    return 1;
  case BAD_DIAGONAL:
    cout << "ERROR: A diagonal entry of A was 0." << endl;
    return 1;
  default:
    cout << "ERROR: Unspecified." << endl;
    return 1;
  case SUCCESS:
    cout << "The solution is:" << endl;
    cout << x << endl;

    Vector y(n);
    matVecMult(A,x,y);
    y -= b;
    cout << "The number of iterations is: " << maxIter << endl;
    cout << "The max-norm of residual is: " << maxNorm(y) << endl;
    cout << "The residual is: " << endl;
    cout << y << endl;
    return 0;
  }
}

