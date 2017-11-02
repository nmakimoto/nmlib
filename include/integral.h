// Numerical integration

#ifndef INTEGRAL_H
#define INTEGRAL_H


#include <cmath>
#include <limits>
#include <stdexcept>
#include "random.h"      // Low-discrepancy sequence
#include "polynomial.h"  // Legendre
#include "matrix.h"
namespace nmlib{



/******** Utility I/F ********/

// Integral of f(x) on [a,b] by Simpson's rule
template<class Func> double integral(const Func& f, double a, double b, size_t n);


// Integral of f(x) on [a,b] by Gaussian quadrature
Matrix weight_gq(size_t n);  // nodes and weights
template<class Func> double integral_gq(const Func& f, double a, double b, size_t n, const Matrix& tw=weight_gq(5));


// Integral of f(x) on (a,b), [a,+inf), or (-inf,+inf) by DE (double exponential) formula
// typically used when f is singular at the boundary or slowly decreasing near infinity
template<class Func> double integral_de(const Func& f, double a, double b, int n, double dt);
template<class Func> double integral_de(const Func& f, double a, /*+inf,*/ int n, double dt);
template<class Func> double integral_de(const Func& f, /*-inf,*/ /*+inf,*/ int n, double dt);


// Multidim integral of f(x) on the region {x: a<x<b} and {x: a<x<b and g(x)=true} by pseudo Monte Carlo
template<class Func>
double integral_mc(const Func& f, const Matrix& a, const Matrix& b, const Lds::Uvec& base, size_t hist);
template<class Func, class Supp>
double integral_mc(const Func& f, const Matrix& a, const Matrix& b, const Lds::Uvec& base, size_t hist, const Supp& g);


/******** Implementation ********/

// Integral of f(x) on [a,b] by Simpson's rule
template<class Func> double integral(const Func& f, double a, double b, size_t n){
  double h=(b-a)/n, a1=a+h/2, sum=0;
  sum=f(a)+f(b);
  for(size_t k=1; k<n; k++) sum+=f(a +h*k)*2;
  for(size_t k=0; k<n; k++) sum+=f(a1+h*k)*4;
  return sum*h/6;
}


// Integral of f(x) on [a,b] by Gaussian quadrature
template<class Func>
double integral_gq(const Func& f, double a, double b, size_t n, const Matrix& tw){
  double sum=0, dx=(b-a)/(n*2);
  for(size_t i=0; i<n; i++){
    double x=a+(i*2+1)*dx;
    for(size_t k=0; k<tw.nrow(); k++){
      double t=tw(k,0), w=tw(k,1);
      sum += f(x+t*dx) * w;
    }
  }
  return sum*dx;
}


// Nodes and weights of Gaussian quadrature
inline Matrix weight_gq(size_t n){
  Polynomial p0=legendre<double>(n), q0=diff(p0), p=p0, q=q0;
  Matrix tw(n,2);

  double x=1, eps=std::numeric_limits<double>::epsilon();
  for(size_t k=0; k<n; k++){
    while(1){
      double dx=p(x)/q(x);
      if( x<=x-dx ) break;
      x=x-dx;
      if( dx<2*eps ) break;  // workaround for infinite loop
    }
    x-=p0(x)/q0(x);  // k-th largest zero of Legendre found by Newton

    tw(k,0)=x;  // Gauss node
    tw(k,1)=2/((1-x*x)*q0(x)*q0(x));  // weight
    p=p/Polynomial({-x,1});
    q=diff(p);
  }
  return tw;
}


// Integral of f(x) on [a,b] by DE formula
template<class Func> double integral_de(const Func& f, double a, double b, int n, double dt){
  double sum=0;
  for(int k=-n; k<=n; k++){
    double t=k*dt, c1=cosh(t)*M_PI/2, s1=sinh(t)*M_PI/2, c2=cosh(s1), x, y;
    x = tanh(s1)*(b-a)/2 + (b+a)/2;
    y = f(x)*c1/c2/c2;
    if( std::isfinite(x) && std::isfinite(y) ) sum += y;
  }
  return sum*dt*(b-a)/2;
}


// Integral of f(x) on [a,+inf) by DE formula
template<class Func> double integral_de(const Func& f, double a, int n, double dt){
  double sum=0;
  for(int k=-n; k<=n; k++){
    double t=k*dt, c1=cosh(t)*M_PI/2, s1=sinh(t)*M_PI/2, e2=exp(s1), x, y;
    x = exp(s1)+a;
    y = f(x)*c1*e2;
    if( std::isfinite(x) && std::isfinite(y) ) sum += y;
  }
  return sum*dt;
}


// Integral of f(x) on (-inf,+inf) by DE formula
template<class Func> double integral_de(const Func& f, int n, double dt){
  double sum=0;
  for(int k=-n; k<=n; k++){
    double t=k*dt, c1=cosh(t)*M_PI/2, s1=sinh(t)*M_PI/2, c2=cosh(s1), x, y;
    x = sinh(s1);
    y = f(x)*c1*c2;
    if( std::isfinite(x) && std::isfinite(y) ) sum += y;
  }
  return sum*dt;
}


// Multple integral of f(x) on {x: a<x<b} by pseudo Monte Carlo
template<class Func>
double integral_mc(const Func& f, const Matrix& a, const Matrix& b, const Lds::Uvec& base, size_t hist){
  auto g = [](const Matrix& /*x*/){ return true; };  // trivial support function
  return integral_mc(f, a, b, base, hist, g);
}


// Multiple integral of f(x) on {a<x<b and g(x)=true} by pseudo Monte Carlo
template<class Func, class Supp>
double integral_mc(const Func& f, const Matrix& a, const Matrix& b, const Lds::Uvec& base, size_t hist, const Supp& g){
  const size_t d=base.size();
  if( !(a.dim()==d && b.dim()==d && a.dim()==d) )
    throw std::domain_error("integral_mc(): invalid dimensions");

  double sum=0;
  size_t cnt=0;
  Lds lds;

  for(size_t it=0; it<hist; it++){
    Matrix x=lds(base);
    for(size_t k=0; k<d; k++) x(k)=(b(k)-a(k))*x(k)+a(k);
    if( !g(x) ) continue;
    sum+=f(x);
    cnt++;
  }

  sum/=hist;
  for(size_t k=0; k<d; k++) sum*=fabs(b(k)-a(k));
  return sum;
}


}  //namespace nmlib
#endif //INTEGRAL_H
