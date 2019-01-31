#include <iostream>
#include <fstream>
using namespace std;

#include "iterativeLA.h"

int read_input_file(string fpath) {
	ifstream reader(fpath) ;
	int i;
	char val;
	for ( i = 0 ; ! reader.eof() ; i++ ){
	reader.get(val);
	cout << val;
	}
	return 0;
}

int main() {
  
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
  cout << A;


  state s = jacobi(A,b,x,maxIter,tolerance);

  switch(s) {
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

