// Non-linear equation solvers


#ifndef SOLVER_H
#define SOLVER_H


#include <cmath>
#include "diff.h"
#include "matrix.h"
namespace nmlib{


/******** Prototype ********/

template<class Func,class T> T solve       (const Func& f, const T& y0, const T& x0, const T& dx, const T& tol=1.e-8, int iter=100);
template<class Func,class T> T solve_bisect(const Func& f, const T& y0, const T& x1, const T& x2, const T& tol=1.e-8);
template<class Func,class T> matrix<T> solve(const Func& f, const matrix<T>& y0, const matrix<T>& x0, const matrix<T>& dx, const T& tol=1.e-8, int iter=100);


/******** Implementation ********/

// Solve f(x)=y0 near x0 by Newton method (scalar version)
template<class Func,class T> T solve(const Func& f, const T& y0, const T& x0, const T& dx, const T& tol, int iter){
  T x,y;
  x=x0;
  for(int i=0; i<iter; i++){
    y=f(x);
    x-=(y-y0)/gradient(f,x,dx,true);  // numerical differentiation
    if(std::abs(y-y0)<tol) break;
  }
  return x;
}

// Solve f(x)=y0 near x0 by Newton method (vecor version)
template<class Func,class T> matrix<T> solve(const Func& f, const matrix<T>& y0, const matrix<T>& x0, const matrix<T>& dx, const T& tol, int iter){
  matrix<T> x,y;
  x=x0;
  for(int i=0; i<iter; i++){
    y=f(x);
    x-=inv(jacobian(f,x,dx,true))*(y-y0);  // numerical differentiation
    if(norm(y-y0)<tol) break;
  }
  return x;
}

// Solve f(x)=y0 in [x1,x2] by bisection
template<class Func,class T> T solve_bisect(const Func& f, const T& y0, const T& x10, const T& x20, const T& tol){
  T x,y,x1(x10),x2(x20),y1(f(x10)),y2(f(x20)), eps=std::numeric_limits<T>::epsilon();
  if(y2<y1){ std::swap(x1,x2); std::swap(y1,y2); }
  while(1){
    x=(x1+x2)/T(2);
    y=f(x);
    if(std::abs(x2-x1)<tol || !((x1<x && x<x2) || (x2<x && x<x1))) break;  // x=x1=x2 (or nan or inf...)
    else if(y<y0){ x1=x; y1=y; }
    else if(y0<y){ x2=x; y2=y; }
    else break;  // y=y0 (or nan or inf...)
    if( std::abs(x2-x1) < (std::abs(x2)+std::abs(x1))*eps ) break;  // workaround for infinite loop
  }
  return x;
}


}  //namespace nmlib
#endif //SOLVER_H
