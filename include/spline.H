// Spline class template library


#ifndef SPLINE_H
#define SPLINE_H


#include <vector>
#include <map>
#include <algorithm>
namespace nmlib{


/******** Class I/F ********/

template<class T> class spline;
typedef spline<double> Spline;


template<class T> class spline{
public:
  spline  (void);
  spline  (const std::map<T,T>& x2y);
  void set(const std::map<T,T>& x2y);        // set control points (xi,yi)
  void get(      std::map<T,T>& x2y) const;  // get control points (xi,yi)

  T operator()(T x, int deg=3) const;  // interpolation at x (degree=0,1,3)

private:
  std::vector<T> xx,yy,aa,bb,cc,dd;  // control points and spline coefficients
  void init_abcd(void);              // calculation of spline coefficients
};


/******** Implementation ********/

template<class T> spline<T>::spline(void){}
template<class T> spline<T>::spline(const std::map<T,T>& x2y){ set(x2y); }

template<class T> void spline<T>::set(const std::map<T,T>& x2y){
  xx=yy=aa=bb=cc=dd=std::vector<T>();
  for(typename std::map<T,T>::const_iterator p=x2y.begin(); p!=x2y.end(); p++){
    xx.push_back(p->first);
    yy.push_back(p->second);
  }
  if(x2y.size()>=3) init_abcd();
}
template<class T> void spline<T>::get(std::map<T,T>& x2y) const{
  x2y=std::map<T,T>();
  for(int i=0; i<xx.size(); i++) x2y[xx[i]]=yy[i];
}

// Interpolation at x (degree=0,1,3)
template<class T> T spline<T>::operator()(T x, int deg) const{
  int n=xx.size();
  if(n<=1) return (n==1 ? yy[0] : 0);
  if(n==2) deg=std::min(deg,1);

  int k=std::upper_bound(xx.begin(),xx.end(),x)-xx.begin()-1;  // x in [xx[k],xx[k+1])
  k=std::max(std::min(k,n-2),0);  // extrapolation
  x-=xx[k];

  if(deg<1) return (x*2<xx[k+1]-xx[k] ? yy[k] : yy[k+1]);  // nearest neighbour
  if(deg<3) return yy[k]+(yy[k+1]-yy[k])*x/(xx[k+1]-xx[k]);  // linear
  return ((aa[k]*x+bb[k])*x+cc[k])*x+dd[k];  // cubic spline
}

// Calculation of spline coefficients
// Determined by the conditions that the spline is:
//   - a piecewise-polynomial Sk(x)=ak(x-xk)^3+bk(x-xk)^2+ck(x-xk)+dk (k=0..n-1)
//   - pass through given control points (x0,y0)...(xn,yn)
//   - continuous up to 2nd derivative
//   - "natural", i.e., 2nd derivative is zero at endpoinnts x0 and xn
template<class T> void spline<T>::init_abcd(void){
  int n=xx.size()-1;  // number of intervals
  std::vector<T> vv(n+1),pp(n),qq(n),rr(n);

  // solve for b[1]..b[n-1] - known to be a tridiagonal matrix equation
  for(int j=1; j<n; j++){
    T h0=xx[j]-xx[j-1], h1=xx[j+1]-xx[j], g0=yy[j]-yy[j-1], g1=yy[j+1]-yy[j];
    pp[j] = 2*(h1+h0);  // diagonal A(j,j)
    qq[j] = h1;  // lower subdiagonal A(j+1,j)
    rr[j] = h0;  // upper subdiagonal A(j-1,j)
    vv[j] = 3*(g1/h1-g0/h0);  // initial value of r.h.s
  }
  for(int j=1; j<n-1; j++){ T s=qq[j]/pp[j];  vv[j+1]-=s*vv[j];  qq[j]=0;  pp[j+1]-=s*rr[j+1]; }  // sweep lower side
  for(int j=n-1; j>1; j--){ T s=rr[j]/pp[j];  vv[j-1]-=s*vv[j];  rr[j]=0; }  // sweep upper side
  for(int j=1; j<n  ; j++){ vv[j]/=pp[j];  pp[j]=1; }  // normalize diagonal part
  vv[0]=vv[n]=0;  // natural spline boundary

  // calculate a,c,d; then store results
  aa=bb=cc=dd=std::vector<T>(n);
  for(int k=0; k<n; k++){
    T b0=vv[k], b1=vv[k+1], h=xx[k+1]-xx[k], g=yy[k+1]-yy[k], y0=yy[k];
    aa[k] = (b1-b0)/(3*h);
    bb[k] = b0;
    cc[k] = g/h-(b1+2*b0)*h/3;
    dd[k] = y0;
  }
}


}  //namespace nmlib
#endif //SPLINE_H
