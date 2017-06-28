// Unit test (ODE solver)


#include <gtest/gtest.h>

#include "ode.h"
using namespace nmlib;


// equation: y'=y
static double f1(double /*t*/, double y){
  return y;
}
// squation: y(t) = y(0) exp(t)
static double g1(double t, double y0){
  return y0*exp(t);
}


// equation: y'' = -w^2 y + A cos(w0t)
static Matrix f2(double t, const Matrix& y){
  const double w=2, w0=1, a=3;
  Matrix dy(2);
  dy(0)=y(1);  // velocity
  dy(1)=-w*w*y(0)+a*cos(w0*t);  // acceleration
  return dy;
}
// solution: y(t) = B cos(w0t) + C1 cos(wt) + C2 sin(wt)  where  B=A/(w^2-w0^2), C1=y(0)-B, C2=y'(0)/w
static Matrix g2(double t, const Matrix& y0){
  const double w=2, w0=1, a=3;
  const double b=a/(w*w-w0*w0), c1=y0(0)-b, c2=y0(1)/w;
  Matrix y(2);
  y(0) =    b*cos(w0*t) +   c1*cos(w*t) +   c2*sin(w*t);
  y(1) =-w0*b*sin(w0*t) - w*c1*sin(w*t) + w*c2*cos(w*t);
  return y;
}


TEST(ode,rk4_1d){
  double t0=0, t1=1, y0=1, tol=1.e-6, ya,yb;
  int n=100;
  std::map<double,double> t2y;

  t2y = solve_ode_rk4(f1,y0,t0,t1,n);  // classical
  ya = t2y.rbegin()->second;
  yb = g1(t1,y0);
  EXPECT_NEAR(ya, yb, tol);

  t2y = solve_ode_rk4i(f1,y0,t0,t1,n);  // implicit
  ya = t2y.rbegin()->second;
  yb = g1(t1,y0);
  EXPECT_NEAR(ya, yb, tol);

  n=10;
  t2y = solve_ode_rk4a(f1,y0,t0,t1,n,tol);  // adaptive
  ya = t2y.rbegin()->second;
  yb = g1(t1,y0);
  EXPECT_NEAR(ya, yb, tol*t2y.size());
}


TEST(ode,rk4_vec){
  double t0=0, t1=10, tol=1.e-6;
  Matrix y0={1.2,3.4}, ya, yb;
  int n=1000;
  std::map<double,Matrix> t2y;

  t2y = solve_ode_rk4(f2,y0,t0,t1,n);  // classical
  ya = t2y.rbegin()->second;
  yb = g2(t1,y0);
  EXPECT_LE(norm(yb-ya), tol);

  t2y = solve_ode_rk4i(f2,y0,t0,t1,n);  // implicit
  ya = t2y.rbegin()->second;
  yb = g2(t1,y0);
  EXPECT_LE(norm(yb-ya), tol);

  n=10;
  t2y = solve_ode_rk4a(f2,y0,t0,t1,n,tol);  // adaptive
  ya = t2y.rbegin()->second;
  yb = g2(t1,y0);
  EXPECT_LE(norm(yb-ya), tol*t2y.size());
}


TEST(ode,rk4_order){
  double t0=0, t1=10;
  Matrix y0={1.2,3.4}, ya, yb, yc;
  int n=1000;
  std::map<double,Matrix> t2y;

  t2y = solve_ode_rk4(f2,y0,t0,t1,n);
  ya = t2y.rbegin()->second;
  t2y = solve_ode_rk4(f2,y0,t0,t1,n*2);
  yb = t2y.rbegin()->second;
  yc = g2(t1,y0);
  EXPECT_NEAR(norm(ya-yc)/norm(yb-yc), 16, 0.1);

  t2y = solve_ode_rk4i(f2,y0,t0,t1,n);
  ya = t2y.rbegin()->second;
  t2y = solve_ode_rk4i(f2,y0,t0,t1,n*2);
  yb = t2y.rbegin()->second;
  yc = g2(t1,y0);
  EXPECT_NEAR(norm(ya-yc)/norm(yb-yc), 16, 0.1);

  n=10;
  for(double tol=1.e-5; tol>1.e-10; tol/=10){
    t2y = solve_ode_rk4a(f2,y0,t0,t1,n,tol);
    ya = t2y.rbegin()->second;
    yb = g2(t1,y0);
    EXPECT_LE( norm(ya-yb), tol*t2y.size());
  }
}
