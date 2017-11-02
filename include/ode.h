// ODE solvers


#ifndef ODE_H
#define ODE_H


#include <map>
#include <cmath>
#include "matrix.h"
namespace nmlib{


/******** Module I/F ********/

// Solvers of ODE initial value problem dy/dt=f(t,y) on [t0,t1], y(t0)=y0
// where value type Y of y is Typically Matrix or double. Variations are:
//   - 4th order explicit Runge-Kutta (classical RK4)
//   - 4th order implicit Runge-Kutta (Gauss-Legendre)
//   - 4th order adaptive Runge-Kutta (Merson)
template<class Func,class Y> std::map<double,Y> solve_ode_rk4 (const Func& f, const Y& y0, double t0, double t1, int n);
template<class Func,class Y> std::map<double,Y> solve_ode_rk4i(const Func& f, const Y& y0, double t0, double t1, int n);
template<class Func,class Y> std::map<double,Y> solve_ode_rk4a(const Func& f, const Y& y0, double t0, double t1, int n, double tol);


/******** Implementation ********/

// private utils for implicit/adaptive RK methods
namespace nmode{
  template<class T> T norm(             T   x){ return std::abs(x); }
  template<class T> T norm(const matrix<T>& x){ return nmlib::norm(x); }
}


template<class Func,class Y>
std::map<double,Y> solve_ode_rk4(const Func& f, const Y& y0, double t0, double t1, int n){
  std::map<double,Y> t2y={{t0,y0},};
  double t=t0, h=(t1-t0)/n;
  Y      y=y0;
  while(n--){
    Y k1,k2,k3,k4;
    k1 = h*f(t,y);
    k2 = h*f(t+h/2, y+k1*0.5);
    k3 = h*f(t+h/2, y+k2*0.5);
    k4 = h*f(t+h  , y+k3    );
    t=t1-n*h;  //t+=h;
    y+=(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    t2y[t]=y;
  }
  return t2y;
}


template<class Func,class Y>
std::map<double,Y> solve_ode_rk4i(const Func& f, const Y& y0, double t0, double t1, int n){
  const double
    rt3=sqrt(3),
    a1 =(3-rt3)/6,   a2= (3+rt3)/6,
    b11=1.0/4,       b12=1.0/4-rt3/6,
    b21=1.0/4+rt3/6, b22=1.0/4,
    c1 =1.0/2,       c2 =1.0/2;
    //c21=(1+rt3)/2,   c22=(1-rt3)/2;

  std::map<double,Y> t2y={{t0,y0},};
  double t=t0, h=(t1-t0)/n;
  Y      y=y0;

  while(n--){
    Y k1,k2;
    k1=h*f(t+a1*h,y);
    k2=h*f(t+a2*h,y);
    double err=hypot(nmode::norm(k1),nmode::norm(k2));

    while(1){
      Y dk1,dk2;  // solve dk=0 by newton - may sometimes fail...
      dk1 = k1 - h*f(t+a1*h, y+b11*k1+b12*k2);
      dk2 = k2 - h*f(t+a2*h, y+b21*k1+b22*k2);
      double err2 = hypot(nmode::norm(dk1),nmode::norm(dk2));
      if( err<=err2 ) break;
      err=err2;
      k1-=dk1;
      k2-=dk2;
    }

    y+=c1*k1+c2*k2;
    t=t1-n*h;  //t+=h;
    t2y[t]=y;
  }
  return t2y;
}


template<class Func,class Y>
std::map<double,Y> solve_ode_rk4a(const Func& f, const Y& y0, double t0, double t1, int n, double tol){
  std::map<double,Y> t2y={{t0,y0},};
  double t=t0, h=(t1-t0)/n;
  Y      y=y0;

  while(n>0){
    Y k1,k2,k3,k4,k5,dy1,dy2;
    k1 = h*f(t,y);
    k2 = h*f(t+h/3, y + k1/3.);
    k3 = h*f(t+h/3, y + (k1 + k2)/6.);
    k4 = h*f(t+h/2, y + (k1 + 3.*k3)/8.);
    k5 = h*f(t+h  , y + (k1 - 3.*k3 + 4.*k4)/2.);
    dy1=(k1        +4.*k4 +k5)/6.;
    dy2=(k1 -3.*k3 +4.*k4    )/2.;

    double err=nmode::norm(dy1-dy2)/5;
    if( err>tol ){ n*=2; h/=2; continue; }
    n--;
    t=t1-n*h;  //t+=h;
    y+=dy1;
    t2y[t]=y;
    if( err<tol/64 && n%2==0 ){ n/=2; h*=2; }
  }
  return t2y;
}


} //namespace nmlib
#endif //ODE_H
