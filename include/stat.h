// Statistics utils


#ifndef STAT_H
#define STAT_H


#include <cmath>
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
  if(p>0.5 && 0.5>1-p) return -pvalue_normal(1-p);
  double x=0,x0=1;
  while(x<x0){
    x0=x;
    x-=(cdf_normal(x)-p)/pdf_normal(x);
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

  // sort eivenvectors in decreasing order of eigenvalue
  size_t n=u.ncol();
  std::vector<std::pair<double,size_t> > d2j(n);
  for(size_t j=0; j<n; j++) d2j[j]=std::pair<double,size_t>(-d(j,j),j);
  std::sort(d2j.begin(),d2j.end());
  Matrix u1(n,n);
  for(size_t j=0; j<n; j++) setvec(u1,j,getvec(u,d2j[j].second));

  // even permutation (det=+1)
  bool odd=false;
  for(size_t j=0; j<n; j++)
    for(size_t k=j+1; k<n; k++)
      if(d2j[j].second>d2j[k].second) odd=!odd;
  if(odd) setvec(u1,n-1,-getvec(u1,n-1));

  return u1;
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