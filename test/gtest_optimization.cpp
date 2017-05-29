// Unit Test (optimization)


#include <gtest/gtest.h>
#include <iostream>

#include "optimization.h"
#include "matrix.h"
using namespace nmlib;


inline double f(const Matrix& x){
  Matrix x0=Matrix({1.2,2.3,3.4});
  return inner(x-x0,x-x0)+4.5;
}


TEST(optimization,amoeba){
  Matrix x, x0({-3,-2,-1}), x1({1.2,2.3,3.4});
  x = opt_amoeba(f,x0,0.1,500,1.e-4);
  EXPECT_NEAR( norm(x-x1), 0, 1.e-4);
  EXPECT_NEAR( f(x)-f(x1), 0, 1.e-8);
}
