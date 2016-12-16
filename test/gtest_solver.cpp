// Unit Test (nonlinear equation solver)


#include <gtest/gtest.h>
#include <iostream>

#include <cmath>
#include "solver.h"
#include "matrix.h"
using namespace nmlib;


inline double f_11(double x){
  return cos(x);
}
inline Matrix f_nn(const Matrix& x){
  Matrix y(2);
  y(0)=cos(x(0)*0.7+x(1)*0.3);
  y(1)=sin(x(0)*x(0)*0.3-x(1)*x(1)*0.4);
  return y;
}


TEST(solver,newton1){
  double x0=1.7, y0=f_11(x0), dx=1.e-3, x1=2.0, x, tol=1.e-8;
  x=solve(f_11,y0,x1,dx,tol);
  EXPECT_NEAR(x, x0, tol);
}


TEST(solver,newtonN){
  Matrix x0(2),x1(2),x(2),dx(2),y0(2);
  double tol=1.e-6;

  for(size_t i=0; i<x0.dim(); i++){
    x0(i)=i+1;
    dx(i)=1.e-6;
    x1(i)=x0(i)+0.1;
  }
  y0=f_nn(x0);

  x=solve(f_nn,y0,x1,dx,tol);
  EXPECT_NEAR(norm(y0-f_nn(x)),0,tol);
  EXPECT_NEAR(norm(x-x0),0,tol);
}


TEST(solver,bisect){
  double x0=1.7, y0=f_11(x0), x1=x0-0.3, x2=x0+0.8, x, tol=1.e-12;

  x=solve_bisect(f_11,y0,x1,x2,tol);
  EXPECT_TRUE(x0-tol<x and x<x0+tol);
  EXPECT_NEAR(std::abs(x0-x), 0, tol);

  x=solve_bisect(f_11,f_11(x0+2*tol),x1,x2,tol);
  EXPECT_FALSE(x0-tol<x and x<x0+tol);

  x=solve_bisect(f_11,f_11(x0-2*tol),x1,x2,tol);
  EXPECT_FALSE(x0-tol<x and x<x0+tol);

  x=solve_bisect(f_11,f_11(x0),x1,x2,-1.0);
  EXPECT_TRUE(x0-tol<x and x<x0+tol);
  EXPECT_NEAR(std::abs(x0-x), 0, tol);
  EXPECT_DOUBLE_EQ(x0,x);
}
