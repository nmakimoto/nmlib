// Unit test (stat)


#include <gtest/gtest.h>

#include <cstdlib>  // random
#include "stat.h"
using namespace nmlib;


// Sample data generators
static double urand1(void){
  return (random()%RAND_MAX+0.5)/RAND_MAX;
}
static Matrix mrand(size_t r, size_t c){
  Matrix m(r,c);
  for(size_t k=0; k<m.dim(); k++) m(k)=urand1()*2-1;
  return m;
}


// Normal distribution
TEST(stat,normal){
  double xx[] = { 0.0,    0.5,    1.0,    1.5,    2.0    };
  double pp[] = { 0.5000, 0.6915, 0.8413, 0.9332, 0.9772 };
  for(int i=0; i<5; i++) EXPECT_NEAR(cdf_normal( xx[i]),   pp[i], 1.e-4);
  for(int i=0; i<5; i++) EXPECT_NEAR(cdf_normal(-xx[i]), 1-pp[i], 1.e-4);

  double dx=1.e-6;
  for(double x=-3.0;  x<3.01;   x+=0.1  ) EXPECT_NEAR(pdf_normal(x), (cdf_normal(x+dx)-cdf_normal(x-dx))/(2*dx), 1.e-6);
  for(double x=-3.0;  x<3.01;   x+=0.1  ) EXPECT_NEAR(pvalue_normal(cdf_normal(x)), x, 1.e-12);
  for(double p=0.001; p<0.9991; p+=0.001) EXPECT_NEAR(cdf_normal(pvalue_normal(p)), p, 1.e-12);

  // boundary
  for(double p: {0.5-1.e-10, 0.5-0.9e-10, 0.5, 0.5+0.9e-10, 0.5+1.e-10, 0.5+1.1e-10})
    EXPECT_NEAR(cdf_normal(pvalue_normal(p)), p, 1.e-12);
  for(double p: {+1.e-10, 1-1.e-10}) EXPECT_TRUE(std::isfinite(pvalue_normal(p)));
  for(double p: {-1.e-10, 1+1.e-10}) EXPECT_TRUE(std::isinf   (pvalue_normal(p)));
}


// Average and variance
TEST(stat,avgvar){
  Matrix data=mrand(3,100);
  Matrix avg=average(data), cov=variance(data);
  Matrix avg2(3),cov2(3,3);

  for(int j=0; j<100; j++) avg2+=getvec(data,j)*0.01;
  for(int j=0; j<100; j++) cov2+=(getvec(data,j)-avg2)*tp(getvec(data,j)-avg2)*0.01;
  EXPECT_NEAR(norm(avg2-avg), 0, 1.e-8);
  EXPECT_NEAR(norm(cov2-cov), 0, 1.e-8);
}


// Fivenum - 0%(min), 25%, 50%(median), 75%, 100%(max)
TEST(stat,fivenum){
  Matrix data(3,100);
  for(int i=0; i<3; i++){
    for(int j=0; j<100; j++) data(i,j)=i*1000+(99-j);
    for(int k=0; k<100; k++)
      std::swap(data(i,int(urand1()*100)%100), data(i,int(urand1()*100)%100));
  }

  Matrix y;
  y=fivenum(data);
  EXPECT_TRUE(y.nrow()==3 && y.ncol()==5);
  for(int i=0; i<3; i++){
    EXPECT_NEAR(y(i,0), i*1000+ 0, 1.e-8);
    EXPECT_NEAR(y(i,1), i*1000+25, 1.e-8);
    EXPECT_NEAR(y(i,2), i*1000+49.5, 1.e-8);
    EXPECT_NEAR(y(i,3), i*1000+75, 1.e-8);
    EXPECT_NEAR(y(i,4), i*1000+99, 1.e-8);
  }
}


// Linear regression
TEST(stat,regression){
  int nx=5, ny=3, k=10;
  Matrix a,b,ab,xx,yy,xx1,one,c1,c2;

  a=mrand(ny,nx);
  b=mrand(ny,1);
  ab=hcat(a,b);
  one.resize(1,k).fill(1);
  xx=mrand(nx,k);
  xx1=vcat(xx,one);

  yy=ab*xx1;  // inhomogeneous
  c1=regression(xx,yy,true);
  for(int i=0; i<ny; i++){
    for(int j=0; j<nx+1; j++){
      c2=c1;
      EXPECT_NE(c1(i,j),0);
      c2(i,j)=c1(i,j)+1.e-8;  EXPECT_LT(norm(yy-c1*xx1), norm(yy-c2*xx1));
      c2(i,j)=c1(i,j)-1.e-8;  EXPECT_LT(norm(yy-c1*xx1), norm(yy-c2*xx1));
    }
  }

  yy=a*xx;  // homogeneous
  c1=regression(xx,yy,false);
  for(int i=0; i<ny; i++){
    for(int j=0; j<nx; j++){
      c2=c1;
      EXPECT_NE(c1(i,j),0);
      c2(i,j)=c1(i,j)+1.e-8;  EXPECT_LT(norm(yy-c1*xx), norm(yy-c2*xx));
      c2(i,j)=c1(i,j)-1.e-8;  EXPECT_LT(norm(yy-c1*xx), norm(yy-c2*xx));
    }
  }
}


// PCA
TEST(stat,pca){
  for(int k=0; k<10; k++){
    Matrix data,u,v,d;
    data=mrand(3,1000);
    u=principal(data);
    v=variance(data);
    d=tp(u)*v*u;

    EXPECT_NEAR(norm(tp(u)*u-1.0), 0, 1.e-8);  // U^T U=1
    EXPECT_NEAR(inner(outer(getvec(u,0),getvec(u,1)),getvec(u,2)), 1, 1.e-8);  // detU=+1
    for(size_t k=0; k+1<d.ncol(); k++) EXPECT_GE(d(k,k), d(k+1,k+1));  // sorted by eigenvalues
    setdiag(d,Matrix(3));
    EXPECT_NEAR(norm(d), 0, 1.e-8);  // U^T V U=diag
  }
}
