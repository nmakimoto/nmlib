// Usage of numerical differentiation code
// - Include "diff.h" and use namespace "nmlib". Nothing to link.
// - Features: gradient, Jacobian


#include <iostream>
#include "diff.h"  // gradient(), jacobian()
using namespace nmlib;


// Sample functions to be differentiated
inline double f11(double        x){ return x*x;          }  // R   --> R
inline double fn1(const Matrix& x){ return inner(x,x);   }  // R^n --> R
inline Matrix fnn(const Matrix& x){ return inner(x,x)*x; }  // R^n --> R^n


int main(void){
  bool romberg=true;  // if true, use higher order formula
  int j=2;  // axis for d/dXj

  double x0=2, dx=1.e-6, df;
  df=gradient(f11,x0,dx,romberg);  std::cout<<"df/dx:\t"<<x0<<'\t'<<df<<'\n';

  Matrix X0=Matrix(1,2,3), dX=Matrix(1,1,1)*1.e-6, dF, J;
  df=gradient(fn1,X0,j,dX(j),romberg);  std::cout<<"df/dXj:\t"       <<X0<<'\t'<<j<<'\t'<<df<<'\n';
  dF=gradient(fn1,X0,  dX,   romberg);  std::cout<<"(df/dXj)_j:\t"   <<X0<<'\t'<<dF<<'\n';
  J =jacobian(fnn,X0,j,dX(j),romberg);  std::cout<<"(dFi/dXj)_i:\t"  <<X0<<'\t'<<j<<'\t'<<J<<'\n';
  J =jacobian(fnn,X0,  dX,   romberg);  std::cout<<"(dFi/dXj)_ij:\t" <<X0<<'\t'<<J <<'\n';
  
  return 0;
}
