// Unit Test (random number generator)


#include <gtest/gtest.h>

#include <map>
#include "random.h"
using namespace nmlib;


TEST(random,init){
  Rng rng;
  std::map<double,int> xx1,xx2,xx3;

  for(int i=0; i<10; i++){ rng.init();    xx1[rng()]=1; }  // random seed
  for(int i=0; i<10; i++){ rng.init(i);   xx2[rng()]=1; }
  for(int i=0; i<10; i++){ rng.init(123); xx3[rng()]=1; }
  EXPECT_EQ(xx1.size(), 10U);
  EXPECT_EQ(xx2.size(), 10U);
  EXPECT_EQ(xx3.size(),  1U);
}


TEST(random,uniform){
  int n=1000000, k=10;
  double sgm1 = sqrt((1.0/k)*(1-1.0/k));  // stddev of binomial dist B(1/k)
  std::map<int,int> count;

  Rng rng(12345ULL);
  for(int i=0; i<n; i++) count[int(rng()*k)]++;
  for(int i=0; i<k; i++) EXPECT_NEAR(count[i]*1.0/n, 1.0/k, 3*sgm1/sqrt(n));
}


TEST(random,white){
  int n=1000000;
  double s01=0,x0=0,x1=0;

  Rng rng(12345ULL);
  for(int i=0; i<n; i++){
    x0=x1;
    x1=rng()-0.5;
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


TEST(random,rejection){
  int n=10000, div=10;
  std::vector<int> cnt(div);


  Rng rng(12345);
  cnt=std::vector<int>(div);
  for(int i=0; i<n; i++){
    double x=rng.rand_pdf(sqrt,1,0,1);
    int j=int(x*div);
    cnt[j]++;
  }
  for(int j=0; j<div; j++){
    double dx=1.0/div, x=j*dx, prob=pow(x+dx,1.5)-pow(x,1.5);
    EXPECT_NEAR(prob*n, cnt[j], sqrt(cnt[j])*3);
  }

  Lds lds(12345);
  cnt=std::vector<int>(div);
  for(int i=0; i<n; i++){
    double x=lds.rand_pdf(sqrt,1,0,1, {2,3});
    int j=int(x*div);
    cnt[j]++;
  }
  for(int j=0; j<div; j++){
    double dx=1.0/div, x=j*dx, prob=pow(x+dx,1.5)-pow(x,1.5);
    EXPECT_NEAR(prob*n, cnt[j], 5);
  }
}


TEST(random,converter){
  int n=10000;
  Lds lds;
  std::vector<int> cnt;

  // Box-Muller
  cnt=std::vector<int>(8);
  double x1=3.186394e-01, x2=6.744898e-01, x3=1.150349e+00;  // p=5/8, 6/8, 7/8
  for(int i=0; i<n; i++){
    Matrix xu=lds({2,3});
    double x=box_muller(xu);
    int seg=0;
    if(x<0){ seg+=0x4; x=-x; }
    if     (x>x3) seg+=0x3;
    else if(x>x2) seg+=0x2;
    else if(x>x1) seg+=0x1;
    cnt[seg]++;
  }
  for(auto c: cnt) EXPECT_NEAR(c, n/8, 10);

  // I^2 --> S^2
  cnt=std::vector<int>(8);
  for(int i=0; i<n; i++){
    Matrix x=conv_i2s2(lds({2,3}));
    int seg=0;
    if(x(0)>0) seg+=0x1;
    if(x(1)>0) seg+=0x2;
    if(x(2)>0) seg+=0x4;
    cnt[seg]++;
  }
  for(auto c: cnt) EXPECT_NEAR(c, n/8, 5);
}
