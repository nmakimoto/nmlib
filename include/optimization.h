// Optimization


#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


#include <vector>
#include <cmath>
#include "diff.h"
#include "matrix.h"
namespace nmlib {


/******** Prototype ********/

// Nelder-Mead optimization (a.k.a amoeba or downhill simplex) - find argmin_x f(x) near x0
//   x0,lambda: initial solution, initial simplex size
//   iter,tol: max iteration, tolerance(final simplex size)
template<class Func> Matrix opt_amoeba(const Func& f, const Matrix& x0, double lambda, size_t iter, double tol);
template<class Func> Matrix opt_amoeba(const Func& f, std::vector<Matrix>& xx, std::vector<double>& ff, size_t iter, double tol);

// Levenberg-Marquardt optimization - find argmin_x ||f(x)||^2  near x0
template<class Func> Matrix opt_lma(const Func& f, const Matrix& x0, const Matrix& dx, double lambda=1.0, double nu=4.0);

// Nonlinear least square curvefit using LMA - find argmin_th \sum_k |yy[k]-f(xx[k],th)|^2 near th0
template<class Func> Matrix curvefit(const Func& f, const Matrix& xx, const Matrix& yy, const Matrix& th0, const Matrix& dth);


/******** Implementation ********/

// Nelder-Mead optimization (simple version)
template<class Func> Matrix opt_amoeba(const Func& f, const Matrix& x0, double lambda, size_t iter, double tol){
  const size_t d=x0.dim(), n=d+1;
  std::vector<Matrix> xx(n,x0);
  std::vector<double> ff(n);
  for(size_t k=0; k<d; k++) xx[k+1](k)+=lambda;  // initial simplex = { x0, x0+lambda*e_i... }
  for(size_t k=0; k<n; k++) ff[k]=f(xx[k]);
  return opt_amoeba(f,xx,ff,iter,tol);
}

// Nelder-Mead optimization (low level version - vertices xx and values ff will be overwritten)
template<class Func> Matrix opt_amoeba(const Func& f, std::vector<Matrix>& xx, std::vector<double>& ff, size_t iter, double tol){
  if( xx.size()<2 || xx.size()!=ff.size() ) throw std::runtime_error("opt_amoeba1(): invalid input");
  const double alpha=1, gamma=2, rho=0.5, sigma=0.5;
  const size_t n=xx.size(), d=xx[0].dim();
  std::vector<int> ord(n);

  while(iter--){
    for(size_t k=0; k<n; k++) ord[k]=k;
    std::sort( ord.begin(), ord.end(),  [&](int i,int j){return ff[i]<ff[j];}  );  // sort ord[] by value ff[]
    const size_t kb=ord[0], kw=ord[n-1], kw2=ord[n-2];  // best/worst/2nd worst

    // centroid
    Matrix xo(d);
    for(size_t k=0; k<n; k++)
      if(k!=kw) xo+=xx[k];
    xo /= double(n-1);

    // check convergence by simplex diameter
    double r=0;
    for(auto x1: xx) r=std::max(r,norm(x1-xo));
    if( r*2<tol ) break;

    // reflection
    Matrix xr = xo+alpha*(xo-xx[kw]);
    double fr = f(xr);
    if( ff[kb]<=fr && fr<ff[kw2] ){ xx[kw]=xr; ff[kw]=fr; continue; }

    // expansion
    if( fr<ff[kb] ){
      Matrix xe = xo + gamma*(xo-xx[kw]);
      double fe = f(xe);
      if( fe<fr ){ xx[kw]=xe; ff[kw]=fe; continue; }
      else       { xx[kw]=xr; ff[kw]=fr; continue; }
    }

    // contraction
    Matrix xc = xo + rho*(xo-xx[kw]);
    double fc = f(xc);
    if( fc<ff[kw] ){ xx[kw]=xc; ff[kw]=fc; continue; }

    // shrink
    for(size_t k=0; k<n; k++){
      if( k==kb ) continue;
      xx[k] = xx[kb] + sigma*(xx[k]-xx[kb]);
      ff[k] = f(xx[k]);
    }
  }
  return xx[ord[0]];
}

// Levenberg-Marquardt optimization - find argmin_x ||f(x)||^2 near x0
template<class Func> Matrix opt_lma(const Func& f, const Matrix& x0, const Matrix& dx, double lambda, double nu){
  Matrix x=x0, y=f(x), jf=jacobian(f,x,dx), hf=tp(jf)*jf;
  while(true){
    Matrix x1,y1,hf1;
    hf1 = hf;
    for(size_t k=0; k<hf1.nrow(); k++) hf1(k,k)*=1+lambda;
    x1 = x - inv(hf1)*tp(jf)*y;
    y1 = f(x1);
    if( std::isnan(norm(x1)) || std::isinf(norm(x1)) ) break;
    if( norm(x1-x)==0 ) break;  // converged
    if( norm(y1)<norm(y) ){ x=x1; y=y1; jf=jacobian(f,x,dx); hf=tp(jf)*jf; if(1+lambda!=1) lambda/=nu; continue; }
    lambda*=nu;
  }
  return x;
}

// Nonlinear least square curvefit using LMA - find argmin_th \sum_k |yy[k]-f(xx[k],th)|^2 near th0
template<class Func> Matrix curvefit(const Func& f, const Matrix& xx, const Matrix& yy, const Matrix& th0, const Matrix& dth){
  class Residual{
    public:
    Residual(const Func& f0, const Matrix& xx0, const Matrix& yy0): f(f0), xx(xx0), yy(yy0) {
      if( xx.dim()!=yy.dim() ) throw std::runtime_error("Resudial: invalid data dimensions");
    }
    Matrix operator()(const Matrix& th) const{
      const size_t n=yy.dim();
      Matrix r(n);
      for(size_t k=0; k<n; k++) r(k)=f(xx(k),th)-yy(k);
      return r;
    }
    const Func& f;
    const Matrix& xx,yy;
  };
  Residual r(f,xx,yy);
  return opt_lma(r,th0,dth);
}


}  //namespace nmlib
#endif //OPTIMIZATION_H
