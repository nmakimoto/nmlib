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
  T grad(T x) const;  // d/dx
  T integral(T xa, T xb, int deg=3) const;   // definite integral over [xa,xb]

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

  int k0=std::upper_bound(xx.begin(),xx.end(),x)-xx.begin()-1;  // x in [xx[k0],xx[k0+1])
  int k1=std::max(std::min(k0,n-2),0);  // valid interval [xx[k1],xx[k1+1]]
  if(deg<1)   return (x-xx[k1]<xx[k1+1]-x ? yy[k1] : yy[k1+1]);  // nearest neighbour
  if(k0== -1) return grad(xx[k1])*(x-xx[k1])+yy[k1];   // linear extrapolation (left)
  if(k0==n-1) return grad(xx[k0])*(x-xx[k0])+yy[k0];  // linear extrapolation (right)
  if(deg<3)   return yy[k1]+(yy[k1+1]-yy[k1])*(x-xx[k1])/(xx[k1+1]-xx[k1]);  // linear
  x-=xx[k1];
  return ((aa[k1]*x+bb[k1])*x+cc[k1])*x+dd[k1];  // cubic spline
}

template<class T> T spline<T>::grad(T x) const{
  int n=xx.size();
  int k=std::upper_bound(xx.begin(),xx.end(),x)-xx.begin()-1;  // x in [xx[k-1],xx[k])
  if(n<3) return (n<2 ? 0 : (yy[1]-yy[0])/(xx[1]-xx[0]));
  if(k  <0){ k=0;   x=xx[0]; }    // extrapolation (left)
  if(n-2<k){ k=n-2; x=xx[n-1]; }  // extrapolation (right)
  x-=xx[k];
  return (3*aa[k]*x+2*bb[k])*x+cc[k];
}

template<class T> T spline<T>::integral(T xa, T xb, int deg) const{
  int n=xx.size();
  T ya=operator()(xa), yb=operator()(xb);
  if(n<3)   return (ya+yb)*(xb-xa)/2;
  if(xb<xa) return -integral(xb,xa,deg);

  auto f3 = [&](T x,int k){ x-=xx[k]; return (((aa[k]/4*x+bb[k]/3)*x+cc[k]/2)*x+dd[k])*x; };
  auto f1 = [&](T x,int k){ return (yy[k]+operator()(x,1))*(x-xx[k])/2; };  // trapezoid
  auto fe = [&](T x,T x0){ return (operator()(x0) + operator()(x)) * (x-x0)/2; };  // linear extrapolation

  int ka = std::upper_bound(xx.begin(),xx.end(),xa)-xx.begin();
  int kb = std::upper_bound(xx.begin(),xx.end(),xb)-xx.begin();

  T sum=0;
  sum += (kb==0 ? fe(xb,xx[0]) : kb==n ? fe(xb,xx[n-1]) : deg<3 ? f1(xb,kb-1) : f3(xb,kb-1));
  sum -= (ka==0 ? fe(xa,xx[0]) : ka==n ? fe(xa,xx[n-1]) : deg<3 ? f1(xa,ka-1) : f3(xa,ka-1));

  for(int k=ka; k<kb; k++ )
    if(0<k && k<n)
      sum += (deg<3 ? f1(xx[k],k-1) : f3(xx[k],k-1));

  return sum;
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
