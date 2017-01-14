// Unit Test (spline)


#include <gtest/gtest.h>
#include <iostream>

#include <map>
#include <cmath>
#include <cstdlib>
#include "spline.h"
using namespace nmlib;


inline double urand(void){
  return (random()+0.5)/RAND_MAX;
}
inline double urand(double x0,double x1){
  return (x1-x0)*urand()+x0;
}
inline std::map<double,double> sample_data(void){
  std::map<double,double> x2y;
  double x,y;
  for(int i=-10; i<=10; i++){
    x=i*0.1*M_PI + urand(-1,1)*0.1;
    y=sin(x) + urand(-1,1)*0.1;
    x2y[x]=y;
  }
  return x2y;
}


TEST(spline,lessthan3pt){
  double x, x1=4.4, x2=6.6, a=3.5, b=-4.7;
  std::map<double,double> x2y;

  Spline f0;       // 0pt (zero)
  x2y[x1]=a*x1+b;
  Spline f1(x2y);  // 1pt (constant)
  x2y[x2]=a*x2+b;
  Spline f2(x2y);  // 2pt (linear)

  for(int i=-100; i<100; i++){
    x=0.11*i;

    EXPECT_DOUBLE_EQ(f0(x  ),0);
    EXPECT_DOUBLE_EQ(f0(x,3),0);
    EXPECT_DOUBLE_EQ(f0(x,1),0);
    EXPECT_DOUBLE_EQ(f0(x,0),0);

    EXPECT_DOUBLE_EQ(f1(x  ),a*x1+b);
    EXPECT_DOUBLE_EQ(f1(x,3),a*x1+b);
    EXPECT_DOUBLE_EQ(f1(x,1),a*x1+b);
    EXPECT_DOUBLE_EQ(f1(x,0),a*x1+b);

    EXPECT_NEAR(f2(x,3),a*x+b,1.e-12);
    EXPECT_NEAR(f2(x,1),a*x+b,1.e-12);
    if(x<(x1+x2)/2-1.e-12) EXPECT_DOUBLE_EQ(f2(x,0), a*x1+b);
    if(x>(x1+x2)/2+1.e-12) EXPECT_DOUBLE_EQ(f2(x,0), a*x2+b);
  }
  EXPECT_DOUBLE_EQ(f1(x2y.begin()->first), x2y.begin()->second);
  for(std::map<double,double>::const_iterator i=x2y.begin(); i!=x2y.end(); i++)
    EXPECT_NEAR(f2(i->first),i->second,1.e-12);
}


TEST(spline,gradient){
  std::map<double,double> x2y = sample_data();
  Spline f(x2y);
  double x0=x2y.begin()->first, xn=x2y.rbegin()->first;
  for(int i=-5; i<=15; i++){
    double x=x0+i*(xn-x0)/10;
    EXPECT_NEAR((f(x+1.e-6)-f(x-1.e-6))/2.e-6, f.grad(x), 1.e-6);  // f'(xk)
  }
}


TEST(spline,boundary){  // should be continuous up to 2nd derivative
  std::map<double,double> x2y = sample_data();
  Spline f(x2y);
  double x,y,dx=1.e-6;
  for(std::map<double,double>::const_iterator i0=x2y.begin(); i0!=x2y.end(); i0++){
    x=i0->first;
    y=i0->second;

    EXPECT_NEAR(f(x),y,1.e-12);  // f(xk)=yk
    EXPECT_NEAR(f.grad(x-dx), f.grad(x+dx), 1.e-4);  // f'(xk-0)=f'(xk+0)
    EXPECT_NEAR((f.grad(x)-f.grad(x-dx))/dx, (f.grad(x+dx)-f.grad(x))/dx, 1.e-4);  // f"(xk-0)=f"(xk+0)
  }
  x=x2y.begin()->first;
  EXPECT_NEAR((f.grad(x+dx)-f.grad(x-dx))/(2*dx), 0, 1.e-4);  // f"(x0)=0
  x=x2y.rbegin()->first;
  EXPECT_NEAR((f.grad(x+dx)-f.grad(x-dx))/(2*dx), 0, 1.e-4);  // f"(xn)=0
}


TEST(spline,extrapolation){  // should be linear
  std::map<double,double> x2y = sample_data();
  Spline f(x2y);
  double x0=x2y.begin()->first, xn=x2y.rbegin()->first;
  for(double dx=0; dx<10.5; dx+=1.0){
    EXPECT_NEAR(f(x0-dx), f(x0)-f.grad(x0)*dx, 1.e-4);
    EXPECT_NEAR(f(xn+dx), f(xn)+f.grad(xn)*dx, 1.e-4);
  }
}
