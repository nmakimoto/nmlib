// Unit Test (random number generator)


#include <gtest/gtest.h>

#include <map>
#include "random.H"
using namespace nmlib;


TEST(random,uniformity){
  int n=1000000, k=10;
  double sgm1 = sqrt((1.0/k)*(1-1.0/k));  // stddev of binomial dist B(1/k)
  std::map<int,int> count;

  init_rand(12345ULL);
  for(int i=0; i<n; i++) count[int(urand()*k)]++;
  for(int i=0; i<k; i++) EXPECT_NEAR(count[i]*1.0/n, 1.0/k, 3*sgm1/sqrt(n));
}


TEST(random,whiteness){
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
  int n=1000000, k=10, bb[]={2,3,5,7,11};

  for(int i=0; i<4; i++){
    int base=bb[i];
    std::map<int,int> count;
    for(int seq=0; seq<n; seq++) count[int(lds(seq,base)*k)]++;
    for(int i=0; i<k; i++) EXPECT_NEAR(count[i], n/k, log(n)/log(base));
  }
}
