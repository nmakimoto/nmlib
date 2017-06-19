// Unit test (numerical integration)


#include <gtest/gtest.h>

#include "integral.h"
using namespace nmlib;


static double f1(double x){ return 1/(1+x*x); }  // atan(x)'
static double f2(double x){ return (1/sqrt(x)+1/sqrt(1-x))/2; }  // sqrt(x)'-sqrt(1-x)'
static double f3(const Matrix& x){ return exp(-inner(x,x)/2)/(2*M_PI); }
static bool   f4(const Matrix& x){ return norm(x)<1; }  // support


TEST(integral,simpson){
  EXPECT_NEAR(integral(exp,0,+1,10), exp(1)-1, 1.e-6);
  EXPECT_NEAR(integral(f1,-1,+1,10), M_PI/2  , 1.e-6);
}


TEST(integral,gauss){
  EXPECT_NEAR(integral_gq(exp,0, 1,10,weight_gq(5)), exp(1)-1,1.e-12);
  EXPECT_NEAR(integral_gq(f1, 0, 1,10,weight_gq(5)), M_PI/4,   1.e-6);
  EXPECT_NEAR(integral_gq(f1,-1,+1,10,weight_gq(5)), M_PI/2,   1.e-6);
}


TEST(integral,deformula){
  EXPECT_NEAR(integral_de(f1, 1,    100,0.1), M_PI/4, 1.e-6);  // [a,+inf)
  EXPECT_NEAR(integral_de(f1,       100,0.1), M_PI  , 1.e-6);  // (-inf,+inf)
  EXPECT_NEAR(integral_de(f2, 0, 1, 100,0.1), 2     , 1.e-6);  // (singular,singular)
}


TEST(integral,montecarlo){
  EXPECT_NEAR(integral_mc(f3,{-1,-1},{+1,+1},{2,3},10000,f4), 1-exp(-0.5), 1.e-2);  // -exp(-r^2/2)|[0,1]
  EXPECT_NEAR(integral_mc(f4,{-1,-1},{+1,+1},{2,3},10000)   , M_PI       , 1.e-2);
}
