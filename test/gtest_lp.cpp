// Unit Test (LP solver)


#include <gtest/gtest.h>

#include <cmath>
#include "lp.h"
#include "matrix.h"
#include "io.h"
using namespace nmlib;


static bool nonnegative(const Matrix& v, double tol=0){
  for(size_t i=0; i<v.dim(); i++)
    if(v(i)+tol<0) return false;
  return true;
}


TEST(lp,optimal){
  Matrix a,b,c,x,y;
  LP lp;
  int ret;

  str2any("2 3   3 2 1   2 5 3", a);
  str2any("2 1   10 15", b);
  str2any("3 1   2 3 4", c);

  // orig: optimal
  lp.init(a,b,c);
  ret=lp.solve();
  x  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(x,1.e-8));
  EXPECT_TRUE(nonnegative(b-a*x,1.e-8));

  // dual: optimal
  lp.init(-tp(a),-c,-b);
  ret=lp.solve();
  y  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(y));
  EXPECT_TRUE(nonnegative(tp(a)*y-c,1.e-8));

  EXPECT_NEAR(inner(c,x), inner(b,y), 1.e-8);  // max cx = min by
  EXPECT_NEAR(inner(c,x), 20, 1.e-8);  // max cx = 20
}


TEST(lp,unbounded){
  Matrix a,b,c,x,y;
  LP lp;
  int ret;

  str2any("3 2   -2 -1   -1 -1   -1 -2", a);
  str2any("3 1   -8 -6 -8", b);
  str2any("2 1   4 3", c);

  // orig: unbounded
  lp.init(a,b,c);
  ret=lp.solve();
  x  =lp.vertex();
  EXPECT_EQ(ret,2);
  EXPECT_TRUE(nonnegative(x,1.e-8));
  EXPECT_TRUE(nonnegative(b-a*x,1.e-8));

  // dual: infeasible
  lp.init(-tp(a),-c,-b);
  ret=lp.solve();
  y  =lp.vertex();
  EXPECT_EQ(ret,-1);
  EXPECT_FALSE(nonnegative(tp(a)*y-c,1.e-8));
}


TEST(lp,infeasible){
  Matrix a,b,c,x,y;
  LP lp;
  int ret;

  str2any("2 2   1 -2   -3 6", a);
  str2any("2 1   -2 -1", b);
  str2any("2 1   1 1", c);

  // orig: infeasible
  lp.init(a,b,c);
  ret=lp.solve();
  x  =lp.vertex();
  EXPECT_EQ(ret,-1);
  EXPECT_FALSE(nonnegative(b-a*x,1.e-8));

  // dual: infeasible
  lp.init(-tp(a),-c,-b);
  ret=lp.solve();
  y  =lp.vertex();
  EXPECT_EQ(ret,-1);
  EXPECT_FALSE(nonnegative(tp(a)*y-c));
}


TEST(lp,degenerate){
  Matrix a,b,c,x,y;
  LP lp;
  int ret;

  // feasible region is a thin triangle near 45-deg line
  // x=(1/2,1/2) is optimal, max cx=1
  int n=101;
  a=Matrix(n,2);
  b=Matrix(n);
  c=Matrix(2);
  for(int i=0; i<n/2; i++){
    a(i,0)=   -(n-i);  a(i,1)=+1;  // -kx + y <= 0 (k>1)
    a(i,0)=1.0/(n-i);  a(i,1)=-1;  // x/k - y <= 0 (k>1)
  }
  a(n/2,0)=a(n/2,1)=1;  // x+y<=1
  b(n/2)=1;  // b_i=0 for i!=n/2 (degenerate)
  c(0)=c(1)=1;

  lp.init(a,b,c);
  ret=lp.solve();
  x  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(x,1.e-8));
  EXPECT_TRUE(nonnegative(b-a*x,1.e-8));

  lp.init(-tp(a),-c,-b);
  ret=lp.solve();
  y  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(y,1.e-8));
  EXPECT_TRUE(nonnegative(tp(a)*y-c,1.e-8));

  EXPECT_NEAR(inner(c,x), inner(b,y), 1.e-8);  // max cx = min by
  EXPECT_NEAR(inner(c,x), 1, 1.e-8);  // max cx = 1
}


TEST(lp,cyclic){
  Matrix a,b,c,x,y;
  LP lp;
  int ret;

  // example where cycling may occur with a usual pivot rule
  str2any("3 4   1 -11 -5 18   1 -3 -1 2   1 0 0 0", a);
  str2any("3 1   0 0 1", b);
  str2any("4 1   10 -57 -9 -24", c);

  // orig: optimal
  lp.init(a,b,c);
  ret=lp.solve();
  x  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(x,1.e-8));
  EXPECT_TRUE(nonnegative(b-a*x,1.e-8));

  // dual: optimal
  lp.init(-tp(a),-c,-b);
  ret=lp.solve();
  y  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(y));
  EXPECT_TRUE(nonnegative(tp(a)*y-c,1.e-8));

  EXPECT_NEAR(inner(c,x), inner(b,y), 1.e-8);  // max cx = min by
}


TEST(lp,large){
  Matrix a,b,c,x,y;
  LP lp;
  int ret;

  // feasible region tightly bounds a disk |x-2c|<=|c|
  // x=3c is one of the optimal solutions and the optimal value is 3|c|^2
  int n=1000;
  a=Matrix(n,2);
  b=Matrix(n);
  c=Matrix(2);
  for(int i=0; i<n; i++){
    a(i,0)=cos(2*i*M_PI/n);
    a(i,1)=sin(2*i*M_PI/n);
  }
  //for(int i=0; i<n; i++) b(i)=1;  // center=0, radius=1
  for(int j=0; j<2; j++) c(j)=1;
  for(int i=0; i<n; i++) b(i)=inner(tp(getsub(a,i,0,1,2)),c)*2+norm(c);  // center=2c, radius=|c|

  // orig: optimal
  lp.init(a,b,c);
  ret=lp.solve();
  x  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(x,1.e-8));
  EXPECT_TRUE(nonnegative(b-a*x,1.e-8));

  // dual: optimal
  lp.init(-tp(a),-c,-b);
  ret=lp.solve();
  y  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(y,1.e-8));
  EXPECT_TRUE(nonnegative(tp(a)*y-c,1.e-8));

  EXPECT_NEAR(inner(c,x), inner(b,y), 1.e-8);  // max cx = min by
  EXPECT_NEAR(inner(c,x), 3*inner(c,c), 1.e-8);  // x=3c is optimal

  // minimize
  lp.init(a,b,-c);
  ret=lp.solve();
  x  =lp.vertex();
  EXPECT_EQ(ret,1);
  EXPECT_TRUE(nonnegative(x,1.e-8));
  EXPECT_TRUE(nonnegative(b-a*x,1.e-8));
  EXPECT_NEAR(inner(c,x), 1*inner(c,c), 1.e-8);  // x=c is optimal
}
