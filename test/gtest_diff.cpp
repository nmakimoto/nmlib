// Unit Test (numerical differentiation)


#include <gtest/gtest.h>

#include <cmath>
#include "diff.h"
#include "matrix.h"
using namespace nmlib;


static double f_11(double x){ return cos(x); }
static double f_n1(const Matrix& x){ return cos(x(0))*cos(2*x(1)); }
static Matrix f_nm(const Matrix& x){ Matrix y(3); y(0)=cos(x(0)-x(1)); y(1)=sin(x(0)-x(1)); y(2)=sin(x(0)); return y; }


TEST(diff,gradient11){
  for(int i=-10; i<=10; i++){
    double x=i*0.9, dx=1.e-3;
    EXPECT_NEAR(gradient(f_11,x,dx,false), -sin(x), 10.e-6);  // df/dx
    EXPECT_NEAR(gradient(f_11,x,dx,true ), -sin(x), 10.e-12);  // higher order version
  }
}


TEST(diff,gradientN1){
  for(int i=-10; i<=10; i++){
    Matrix x(2), dx(2), dy(2);
    x(0)=i*0.9;
    x(1)=1.2;
    dx(0)=1.e-3;
    dx(1)=1.e-3;
    dy(0)=-sin(x(0))*cos(2*x(1));
    dy(1)=-cos(x(0))*sin(2*x(1))*2;

    for(size_t j=0; j<x.nrow(); j++){
      EXPECT_NEAR(gradient(f_n1,x,j,dx(j),false), dy(j), 10.e-6);  // df/dxj
      EXPECT_NEAR(gradient(f_n1,x,j,dx(j),true ), dy(j), 10.e-12);
    }
    EXPECT_NEAR(norm(gradient(f_n1,x,dx,false)-dy), 0, 10.e-6);  // (df/dxj)_j
    EXPECT_NEAR(norm(gradient(f_n1,x,dx,true )-dy), 0, 10.e-12);
    EXPECT_NEAR(norm(gradient(f_n1,x,dx(0),false)-dy), 0, 10.e-6);
    EXPECT_NEAR(norm(gradient(f_n1,x,dx(0),true)-dy), 0, 10.e-12);
  }
}


TEST(diff,jacobian){
  for(int i=-10; i<=10; i++){
    Matrix x(2), dx(2), dy(3,2);
    x(0)=i*0.9;
    x(1)=1.2;
    dx(0)=1.e-3;
    dx(1)=1.e-3;
    dy(0,0)=-sin(x(0)-x(1));  dy(0,1)=+sin(x(0)-x(1));
    dy(1,0)=+cos(x(0)-x(1));  dy(1,1)=-cos(x(0)-x(1));
    dy(2,0)=+cos(x(0));       dy(2,1)=0;

    for(size_t j=0; j<x.nrow(); j++){
      EXPECT_NEAR(norm(jacobian(f_nm,x,j,dx(j),false)-getvec(dy,j)), 0, 10.e-6);  // (dfi/dxj)_i
      EXPECT_NEAR(norm(jacobian(f_nm,x,j,dx(j),true )-getvec(dy,j)), 0, 10.e-12);
    }
    EXPECT_NEAR(norm(jacobian(f_nm,x,dx,false)-dy), 0, 10.e-6);  // (dfi/dxj)_ij
    EXPECT_NEAR(norm(jacobian(f_nm,x,dx,true )-dy), 0, 10.e-12);
    EXPECT_NEAR(norm(jacobian(f_nm,x,dx(0),false)-dy), 0, 10.e-6);
    EXPECT_NEAR(norm(jacobian(f_nm,x,dx(0),true )-dy), 0, 10.e-12);
  }
}


TEST(diff,hessian){
  Matrix x({1,2}), dx({1.e-3,1.e-3}), h=hessian(f_n1,x,dx);
  Matrix h1(2,2);
  h1(0,0) = -cos(x(0))*cos(2*x(1));
  h1(1,1) = -cos(x(0))*cos(2*x(1))*4;
  h1(0,1) = h1(1,0) = +sin(x(0))*sin(2*x(1))*2;
  EXPECT_NEAR(norm(h-h1), 0, 1.e-5*norm(h));
  EXPECT_NEAR(norm(h-tp(h)), 0, 1.e-5*norm(h));
}
