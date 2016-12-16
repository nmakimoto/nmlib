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
  int n=10000, k=10, bb[]={2,3,5,7,11,13}, dim=sizeof(bb)/sizeof(bb[0]);

  // van der Corput sequence
  for(int i=0; i<dim; i++){
    int base=bb[i];
    std::map<int,int> count;
    for(int seq=1; seq<=n; seq++) count[int(corput(seq,base)*k)]++;
    for(int j=0; j<k; j++) EXPECT_NEAR(count[j], n/k, log(n)/log(base)*5);
  }

  // Halton sequence
  {
    std::map<int,int> cnt[dim];
    for(int seq=1; seq<=n; seq++){
      Matrix x=halton(seq,dim);
      for(int d=0; d<dim; d++) cnt[d][int(x(d)*k)]++;
    }
    for(int d=0; d<dim; d++)
      for(int j=0; j<k; j++) EXPECT_NEAR(cnt[d][j], n/k, log(n)/log(bb[d])*5);  // uniform

    double cov=0;
    for(int seq=1; seq<=n; seq++){
      Matrix x=halton(seq,2);
      cov+=x(0)*x(1);
    }
    cov=cov/n-0.5*0.5;
    EXPECT_NEAR(cov,0,5.0/n);  // no correlation

    EXPECT_THROW(halton(0,7), std::runtime_error);  // max dimension
  }
}
