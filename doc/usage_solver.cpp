// Usage of equation solvers
// - Include "solver.h" and use namespace "nmlib". Nothing to link.
// - Features: non-linear equation solvers (Newton, bisection, multi-dim Newton)


#include <iostream>
#include "solver.h"  // solve*()
using namespace nmlib;


// Sample functions to be solved
inline double f(double        x){  return x*x;  }  // R-->R
inline Matrix F(const Matrix& x){  return Matrix(inner(x,x), x(1)/x(0), x(2)/x(0));  }  // R^3-->R^3


int main(void){
  // Solve f(x)=y0 by Newton
  double x0=2, y0=f(x0);  // answer
  double x1=3, dx=1.e-6;  // initial value, step-size of numerical differentiation
  double x=solve(f,y0,x1,dx);
  std::cout<<"Newton:\n"<<x0<<' '<<x<<'\n'<<y0<<' '<<f(x)<<'\n';

  // Solve f(x)=y0 by bisection
  double xmin=0, xmax=10;  // initial interval
  double xb=solve_bisect(f,y0,xmin,xmax);
  std::cout<<"Bisection:\n"<<x0<<' '<<xb<<'\n'<<y0<<' '<<f(xb)<<'\n';

  // Solve F(X)=Y0 by multi-dimensional Newton
  Matrix X0=Matrix(1,2,3), Y0=F(X0);
  Matrix X1=Matrix(2,3,4), dX=Matrix(1,1,1)*1.e-6;
  Matrix X=solve(F,Y0,X1,dX);
  std::cout<<"Newton(multi-dimensional):\n"<<X0<<'\n'<<X<<'\n'<<Y0<<'\n'<<F(X)<<'\n';
  
  return 0;
}
