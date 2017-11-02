// Statistics utils


#ifndef STAT_H
#define STAT_H


#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include "matrix.h"
namespace nmlib{


/******** Prototype ********/

double pdf_normal(double x);     // probability density of N(0,1)
double cdf_normal(double x);     // cummulative probability of N(0,1)
double pvalue_normal(double p);  // p-value of N(0,1)

Matrix average    (const Matrix& xx);  // note: xx(i,j) is the i-th component of j-th sample
Matrix variance   (const Matrix& xx);
Matrix covariance (const Matrix& xx, const Matrix& yy);
double correlation(const Matrix& xx, const Matrix& yy);  // correlation of 1D sequence
Matrix fivenum    (const Matrix& xx);  // 0%(min),25%,50%(median),75%,100%(max) percentile
Matrix principal  (const Matrix& xx);  // eigenvectors of Var[X] with decreasing eigenvalues
Matrix regression (const Matrix& xx, const Matrix& yy, bool add_b=true);  // coefficent A of Y=AX or [A|b] of Y=AX+b


/******** Implementation ********/

inline double pdf_normal(double x){  return exp(-x*x/2)/sqrt(2*M_PI);  }
inline double cdf_normal(double x){  return erfc(-x/sqrt(2))/2;  }
inline double pvalue_normal(double p){
  if( !std::isfinite(p) ) return p;  // +-inf,nan
  if( !(0<p && p<1) ) return std::numeric_limits<double>::infinity()*(p>0 ? +1 : -1);  // not in (0,1)
  if( 0.5+1.e-12<p ) return -pvalue_normal(1-p);  // avoid cancellation and ensure monotonic convergence

  double x=0, eps=std::numeric_limits<double>::epsilon();
  while( cdf_normal(x-1)>p ) x-=1;
  while(1){
    double y=cdf_normal(x)-p, dy=pdf_normal(x), dx=y/dy;
    if( y==0 || dy==0 || !std::isfinite(dx) ) break;
    if( x<=x-dx ) break;
    x=x-dx;
    if( dx < 2*(1+std::abs(x))*eps ) break;  // workaround for infinite loop
  }
  return x;
}


inline Matrix average(const Matrix& xx){
  int n=xx.ncol();
  Matrix wt(n);
  wt.fill(1);
  return xx*wt*(1.0/n);
}
inline Matrix variance(const Matrix& xx){
  int n=xx.ncol();
  Matrix exx=xx*tp(xx)*(1.0/n), ex=average(xx);
  return exx-ex*tp(ex);
}
inline Matrix covariance(const Matrix& xx, const Matrix& yy){
  if(xx.ncol()!=yy.ncol()) throw std::domain_error("covariance(): x and y must have same number of samples");
  int n=xx.ncol();
  Matrix ex=average(xx), ey=average(yy), exy=xx*tp(yy)*(1.0/n);
  return exy-ex*tp(ey);
}
inline double correlation(const Matrix& xx, const Matrix& yy){
  if(!(xx.nrow()==1 && yy.nrow()==1)) throw std::domain_error("correlation(): x and y must be one dimensional");
  return covariance(xx,yy)(0) / sqrt(variance(xx)(0)*variance(yy)(0));
}


inline Matrix fivenum(const Matrix& xx){
  int n=xx.ncol(), d=xx.nrow();
  Matrix buf(n), ret(d,5);
  for(int i=0; i<d; i++){
    for(int j=0; j<n; j++) buf(j)=xx(i,j);
    std::sort(&buf(0),&buf(0)+n);
    ret(i,0) = buf(0);
    ret(i,1) = buf(n/4);
    ret(i,2) = (n%2 ? buf(n/2) : (buf(n/2)+buf(n/2-1))/2);
    ret(i,3) = buf(n*3/4);
    ret(i,4) = buf(n-1);
  }
  return ret;
}


inline Matrix principal(const Matrix& xx){
  Matrix v,u,d;
  v=variance(xx);
  u=eigen(v);  // eigenvectors of Var[X]
  d=tp(u)*v*u;

  size_t n=u.ncol();
  sort_columns_by_value(u,-getdiag(d));  // decreasing order of eigenvalue
  if( det(u)<0 ) setvec(u,n-1,-getvec(u,n-1));  // det=+1
  return u;
}


inline Matrix regression(const Matrix& xx, const Matrix& yy, bool add_b){
  if(add_b){
    Matrix one;
    one.resize(1,xx.ncol()).fill(1);
    return regression(vcat(xx,one),yy,false);
  }
  return yy*tp(xx)*inv(xx*tp(xx));
}


} //namespace nmlib
#endif //STAT_H
