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

  // van der Corput sequence
  for(int i=0; i<4; i++){
    int base=bb[i];
    Corput c(base);
    std::map<int,int> count;
    for(int seq=0; seq<n; seq++) count[int(c()*k)]++;
    for(int j=0; j<k; j++) EXPECT_NEAR(count[j], n/k, log(n)/log(base));
  }

  // Halton sequence
  {
    Halton h(3);
    std::map<int,int> cnt0,cnt1,cnt2;
    for(int seq=0; seq<n; seq++){
      Matrix x=h();
      cnt0[int(x(0)*k)]++;
      cnt1[int(x(1)*k)]++;
      cnt2[int(x(2)*k)]++;
    }
    for(int j=0; j<k; j++) EXPECT_NEAR(cnt0[j], n/k, log(n)/log(2));
    for(int j=0; j<k; j++) EXPECT_NEAR(cnt1[j], n/k, log(n)/log(3));
    for(int j=0; j<k; j++) EXPECT_NEAR(cnt2[j], n/k, log(n)/log(5));

    EXPECT_THROW(Halton(6), std::runtime_error);
  }
}
