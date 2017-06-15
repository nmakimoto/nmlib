// Unit Test (optimization)


#include <gtest/gtest.h>
#include <iostream>

#include "optimization.h"
#include "matrix.h"
#include "random.h"
using namespace nmlib;


static double f(const Matrix& x){
  Matrix x0=Matrix({1.2,2.3,3.4});
  return inner(x-x0,x-x0)+4.5;
}

static Matrix g(const Matrix& x){
  Matrix x0=Matrix({1.2,2.3,3.4});
  Matrix y=x-x0;
  y(0)+= y(0)*y(0)*y(0);
  y(1) = exp(y(1))-1;
  y(2)+= y(0)*y(1);
  return y;
}

static double h(double x, const Matrix& th){
  return (x-th(0))*sqrt(fabs(th(1)));
}


TEST(optimization,amoeba){
  Matrix x, x0({-3,-2,-1}), x1({1.2,2.3,3.4});
  x = opt_amoeba(f,x0,0.1,500,1.e-4);
  EXPECT_NEAR( norm(x-x1), 0, 1.e-4);
  EXPECT_NEAR( f(x)-f(x1), 0, 1.e-8);
}


TEST(optimization,lma){
  Matrix x, x0({-3,-2,-1}), x1({1.2,2.3,3.4}), dx({1.e-4,1.e-4,1.e-4});
  x = opt_lma(g, x0, dx);
  EXPECT_NEAR( norm(x-x1), 0, 1.e-4);
  EXPECT_NEAR( norm(g(x)-g(x1)), 0, 1.e-8);
}


TEST(optimization,curvefit){
  size_t n=100;
  Matrix th({1.1,2.2}), th0({1,1}), dth({1.e-6,1.e-6});
  Matrix xx(n), yy(n);

  Rng rng(12345);
  for(size_t k=0; k<n; k++){
    xx(k)=rng.n();
    yy(k)=h(xx(k),th) + rng.n()*0.1;
  }

  Matrix th1=curvefit(h,xx,yy,th0,dth);

  for(size_t k=0; k<th1.dim(); k++){
    Matrix th2=th1;
    double err1, err2;

    err1=0;  for(size_t i=0; i<xx.dim(); i++){ double dy=h(xx(i),th1)-yy(i); err1+=dy*dy; }

    th2(k)=th1(k)+1.e-8;
    err2=0;  for(size_t i=0; i<xx.dim(); i++){ double dy=h(xx(i),th2)-yy(i); err2+=dy*dy; }
    EXPECT_LT(err1,err2);

    th2(k)=th1(k)-1.e-8;
    err2=0;  for(size_t i=0; i<xx.dim(); i++){ double dy=h(xx(i),th2)-yy(i); err2+=dy*dy; }
    EXPECT_LT(err1,err2);
  }
}
