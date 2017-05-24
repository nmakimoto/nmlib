// Numerical differentiation


#ifndef DIFF_H
#define DIFF_H


#include "matrix.h"
namespace nmlib{


/******** Prototype ********/

// Gradient of scalar-valued function f: R^n --> R
template<class F,class T>        T  gradient(const F& f, const        T & x,        const T & dx , bool romberg=false);  // df/dx
template<class F,class T>        T  gradient(const F& f, const matrix<T>& x, int j, const T & dxj, bool romberg=false);  // df/dxj
template<class F,class T> matrix<T> gradient(const F& f, const matrix<T>& x, const matrix<T>& dx , bool romberg=false);  // (df/dxj)_j
template<class F,class T> matrix<T> gradient(const F& f, const matrix<T>& x,        const T & dt , bool romberg=false);  // (df/dxj)_j

// Jacobian of vector-valued function f: R^n --> R^m
template<class Func,class T> matrix<T> jacobian(const Func& f, const matrix<T>& x, int j, const T & dxj, bool romberg=false);  // (dfi/dxj)_i
template<class Func,class T> matrix<T> jacobian(const Func& f, const matrix<T>& x, const matrix<T>& dx , bool romberg=false);  // (dfi/dxj)_ij
template<class Func,class T> matrix<T> jacobian(const Func& f, const matrix<T>& x,        const T & dt , bool romberg=false);  // (dfi/dxj)_ij

// Hessian of scalar-valued function f: R^n --> R
template<class Func,class T> matrix<T> hessian(const Func& f, const matrix<T>& x, const matrix<T>& dx);  // (d^2f/dxidxj)_ij


/******** Implementation ********/

// df/dx
template<class Func,class T> T gradient(const Func& f, const T& x, const T& dx, bool romberg){
  if(romberg) return (T(4)*gradient(f,x,dx/T(2),false)-gradient(f,x,dx,false))/T(3);  // Romberg acceleration - error=O(dx^4)
  return (f(x+dx)-f(x-dx))/(T(2)*dx);  // difference quotient - error=O(dx^2)
}

// df/dxj
template<class Func,class T> T gradient(const Func& f, const matrix<T>& x, int j, const T& dxj, bool romberg){
  if(romberg) return (T(4)*gradient(f,x,j,dxj/T(2),false)-gradient(f,x,j,dxj,false))/T(3);
  matrix<T> x1(x), x2(x);
  x1(j)-=dxj;
  x2(j)+=dxj;
  return (f(x2)-f(x1))/(T(2)*dxj);
}

// (df/dxj)_j
template<class Func,class T> matrix<T> gradient(const Func& f, const matrix<T>& x, const matrix<T>& dx, bool romberg){
  if(romberg) return (T(4)*gradient(f,x,dx/T(2),false)-gradient(f,x,dx,false))/T(3);
  int r=x.nrow();
  matrix<T> dy(r);
  for(int j=0; j<r; j++) dy(j)=gradient(f,x,j,dx(j),romberg);
  return dy;
}
template<class Func,class T> matrix<T> gradient(const Func& f, const matrix<T>& x,        const T & dt, bool romberg){
  matrix<T> dx(x.dim());
  for(size_t i=0; i<dx.dim(); i++) dx(i)=dt;  // dx=(dt ... dt)^T
  return gradient(f,x,dx,romberg);
}

// (dfi/dxj)_i
template<class Func,class T> matrix<T> jacobian(const Func& f, const matrix<T>& x, int j, const T& dxj, bool romberg){
  if(romberg) return (T(4)*jacobian(f,x,j,dxj/T(2),false)-jacobian(f,x,j,dxj,false))/T(3);
  matrix<T> x1(x),x2(x);
  x1(j)=x(j)-dxj;
  x2(j)=x(j)+dxj;
  return (f(x2)-f(x1))/(T(2)*dxj);
}

// (dfi/dxj)_ij
template<class Func,class T> matrix<T> jacobian(const Func& f, const matrix<T>& x, const matrix<T>& dx, bool romberg){
  if(romberg) return (T(4)*jacobian(f,x,dx/T(2),false)-jacobian(f,x,dx,false))/T(3);
  int r=f(x).nrow(), c=x.nrow();
  matrix<T> dy(r,c);
  for(int j=0; j<c; j++) setvec(dy,j,jacobian(f,x,j,dx(j),romberg));
  return dy;
}
template<class Func,class T> matrix<T> jacobian(const Func& f, const matrix<T>& x,        const T & dt, bool romberg){
  matrix<T> dx(x.dim());
  for(size_t i=0; i<dx.dim(); i++) dx(i)=dt;  // dx=(dt ... dt)^T
  return jacobian(f,x,dx,romberg);
}

// (d^2f/dxidxj)_ij - note that too small dx may cause loss of digits
template<class Func,class T> matrix<T> hessian(const Func& f, const matrix<T>& x, const matrix<T>& dx){
  const size_t n=x.dim();
  matrix<T> hf(n,n);
  matrix<T> df=gradient(f,x,dx,true);
  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      if( j<i ){ hf(i,j)=hf(j,i); continue; }
      matrix<T> x1=x;
      if( j==i ){
	x1(i)=x(i)+dx(i);  hf(i,i)+=f(x1);
	x1(i)=x(i)-dx(i);  hf(i,i)+=f(x1);
	hf(i,i)-=T(2)*f(x);
	hf(i,i)/=dx(i)*dx(i);
      }
      else{
	x1(i)=x(i)+dx(i);  x1(j)=x(j)+dx(j);  hf(i,j)+=f(x1);
	x1(i)=x(i)+dx(i);  x1(j)=x(j)-dx(j);  hf(i,j)-=f(x1);
	x1(i)=x(i)-dx(i);  x1(j)=x(j)+dx(j);  hf(i,j)-=f(x1);
	x1(i)=x(i)-dx(i);  x1(j)=x(j)-dx(j);  hf(i,j)+=f(x1);
	hf(i,j)/=dx(i)*dx(j)*T(4);
      }
    }
  }
  return hf;
}


}  //namespace nmlib
#endif //DIFF_H
