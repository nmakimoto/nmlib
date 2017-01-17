// Unit Test (random number generator)


#include <gtest/gtest.h>

#include <map>
#include "random.h"
using namespace nmlib;


TEST(random,uniform){
  int n=1000000, k=10;
  double sgm1 = sqrt((1.0/k)*(1-1.0/k));  // stddev of binomial dist B(1/k)
  std::map<int,int> count;

  init_rand(12345ULL);
  for(int i=0; i<n; i++) count[int(urand()*k)]++;
  for(int i=0; i<k; i++) EXPECT_NEAR(count[i]*1.0/n, 1.0/k, 3*sgm1/sqrt(n));
}


TEST(random,white){
  int n=1000000;
  double s01=0,x0=0,x1=0;

  init_rand(12345ULL);
  for(int i=0; i<n; i++){
    x0=x1;
    x1=urand()-0.5;
    s01+=x0*x1;
  }
  s01/=n;  // autocorrelation E X(t)X(t+1)

  double sgm1=1.0/12;  // stddev of x0*x1 = (stddev of x0)^2
  EXPECT_NEAR(fabs(s01), 0, 3*sgm1/sqrt(n));
}


TEST(random,lowdiscrepancy){
  int n=10000, k=10;

  // van der Corput sequence
  for(int base: {2,3,5,7,11,13}){
    std::map<int,int> count;
    for(int seq=1; seq<=n; seq++) count[int(corput(seq,base)*k)]++;
    for(int j=0; j<k; j++) EXPECT_NEAR(count[j], n/k, log(n));  // uniform
  }

  // Halton sequence
  {
    std::vector<uint_fast64_t> bb={2,3,5,7,11,13};
    const int dim=bb.size();
    Matrix s1(dim), s2(dim,dim);
    for(int seq=1; seq<=n; seq++){
      Matrix x=halton(seq, bb);
      s1+=x;
      s2+=x*tp(x);
    }
    s1/=double(n);
    s2=s2/double(n)-s1*tp(s1);
    for(int i=0; i<dim; i++){
      EXPECT_NEAR(s1(i), 0.5, log(n)/n);  // unbiased
      for(int j=0; j<dim; j++)
        if(i!=j) EXPECT_NEAR(s2(i,j), 0.0, log(n)/n);  // uncorrelated
    }
  }
}
