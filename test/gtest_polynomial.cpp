// Unit Test (polynomial class)


#include <gtest/gtest.h>
#include <stdexcept>

#include <map>
#include <cmath>
#include "polynomial.h"
#include "io.h"
#include "solver.h"  // for Legendre zeros
using namespace nmlib;


// Abbreviations
#define FOR1(P,I)  for(int I=P.deg(); I>=0; I--)
#define EEQ(A,B)   EXPECT_EQ(A,B)
#define ELE(A,B)   EXPECT_LE(A,B)
#define EDQ(A,B)   EXPECT_DOUBLE_EQ(A,B)
#define ENR(A,B,C) EXPECT_NEAR(A,B,C)
#define ETH(A,B)   EXPECT_THROW(A,B)


Polynomial P,Q,X,Z;


void init_polynomial_testdata(void){
  P=Q=X=Z=Polynomial();
  for(int i=0; i<=6; i++) P.c(i)=(i+1)* 1;
  for(int i=0; i<=3; i++) Q.c(i)=(i+1)*10;
  X.c(1)=1;
  Z=Polynomial();
}


// Access to coefficients
TEST(polynomial, coefficient){
  Polynomial p;

  // default (empty)
  EEQ(p.deg(),-1);  // deg=-1
  EDQ(p(10),0);

  // auto resize
  p.c(5)=5;
  EEQ(p.deg(),5);
  for(int i=0; i<=5; i++) EEQ(p.c  (i),(i==5 ? 5 : 0));
  ETH(p.c(-1),std::domain_error);

  // set coefficients
  for(int i=0; i<=p.deg(); i++) p.c(i)=i+1;
  EEQ(p.deg(),5);
  EDQ(p(10),654321);
  p.c(6)=7;  // modify p...
  EDQ(p(10),7654321);  // existing coefficients remain unchanged

  // copy
  Polynomial q=p;
  EEQ(q.deg(),p.deg());
  FOR1(p,i) EDQ(q.c(i),p.c(i));
  FOR1(p,i) p.c(i)++;  // modify p...
  FOR1(p,i) EDQ(p.c(i),q.c(i)+1);  // q remains unchanged

  // c() const version
  const Polynomial r=p;
  EEQ(r.deg(),p.deg());  FOR1(p,i) EDQ(r.c(i),p.c(i));
  EEQ(r.c(10),0);
  EEQ(r.c(-1),0);
}


// Arithmetic operations (incremental)
TEST(polynomial,arithmetic_incr){
  init_polynomial_testdata();
  Polynomial p(P),q(Q),r,s;
  double x=1.23;

  // scalar operation
  r=p; r+=x;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)+(i==0 ? x : 0));
  r=p; r-=x;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)-(i==0 ? x : 0));
  r=p; r*=x;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)*x);
  r=p; r/=x;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)/x);

  // sum/difference
  r=p; r+=q;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)*(i<=3 ? 11 : 1));
  r=q; r+=p;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)*(i<=3 ? 11 : 1));
  r=p; r-=q;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)*(i<=3 ? -9 : 1));
  r=q; r-=p;  EEQ(r.deg(),p.deg()); for(int i=0; i<=6; i++) EDQ(r.c(i),(i+1)*(i<=3 ? +9 :-1));

  // multiplivstion
  r=p; r*=q;
  EEQ(r.deg(),6+3);
  EDQ(r.c(0),p.c(0)*q.c(0));
  EDQ(r.c(1),p.c(1)*q.c(0)+p.c(0)*q.c(1));
  EDQ(r.c(5),p.c(5)*q.c(0)+p.c(4)*q.c(1)+p.c(3)*q.c(2)+p.c(2)*q.c(3));
  EDQ(r.c(9),p.c(6)*q.c(3));

  // division
  r=q;  r/=p;  EEQ(r.deg(),-1);
  s=q;  s%=p;  EEQ(s.deg(),q.deg());  FOR1(q,i) EDQ(s.c(i),q.c(i));
  r=p;  r/=q;  EEQ(r.deg(),p.deg()-q.deg());
  s=p;  s%=q;  ELE(s.deg(),q.deg()-1);
  r*=q; r+=s;  EEQ(r.deg(),p.deg());  FOR1(p,i) EDQ(r.c(i),p.c(i));  // p=(p/q)*q+(p%q)

  // shift
  r=p;  r<<=3;  EEQ(r.deg(),p.deg()+3);  FOR1(r,i) EDQ(r.c(i),(i<3 ? 0 : p.c(i-3)));
  r=p;  r>>=3;  EEQ(r.deg(),p.deg()-3);  FOR1(r,i) EDQ(r.c(i),p.c(i+3));
  r=p;  r>>=10;  EEQ(r.deg(),-1);
  r=p;  ETH(r<<=-1, std::domain_error);
  r=p;  ETH(r>>=-1, std::domain_error);

  // operation with itself
  r=p; r=-r;  EEQ(r.deg(),p.deg());  FOR1(r,i) EDQ(r.c(i),-p.c(i));
  r=p; r+=r;  EEQ(r.deg(),p.deg());  FOR1(r,i) EDQ(r.c(i),p.c(i)*2); 
  r=p; r-=r;  EEQ(r.deg(),p.deg());  FOR1(r,i) EDQ(r.c(i),0); 
  r=p; r*=r;  EEQ(r.deg(),p.deg()*2);
  r=p; r/=r;  EEQ(r.deg(),0);  EDQ(r.c(0),1); 
  r=p; r%=r;  ELE(r.deg(),p.deg()-1);  FOR1(r,i) EDQ(r.c(i),0);

  // operations with zero
  r=p;  s=Z;  r*=s;  EEQ(r.deg(),-1);
  r=Z;  s=p;  r*=s;  EEQ(r.deg(),-1);
  r=s=Z;      r*=s;  EEQ(r.deg(),-1);

  r=Z;  s=p;  r/=s;  EEQ(r.deg(),-1);
  r=p;  s=Z;  ETH(r/=s, std::domain_error);
  r=s=Z;      ETH(r/=s, std::domain_error);

  r=Z;  s=p;  r%=p;  EEQ(r.deg(),-1);
  r=p; s=Z;   ETH(r%=s, std::domain_error);
  r=s=Z;      ETH(r%=s, std::domain_error);
}


// Arithmetic operations
TEST(polynomial, arithmetic){
  init_polynomial_testdata();
  Polynomial p(P),q(Q),r,s;
  double x=1.23;

  r=-p; s=p; s*=-1.0; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));

  r=p+x; s=p; s+=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p-x; s=p; s-=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p*x; s=p; s*=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p/x; s=p; s/=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));

  r=x+p; s=p; s+=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=x-p; s=p; s*=-1.0; s+=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=x*p; s=p; s*=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  //r=x/p; s=p; s/=x; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));

  r=p+q; s=p; s+=q; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p-q; s=p; s-=q; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p*q; s=p; s*=q; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p/q; s=p; s/=q; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p%q; s=p; s%=q; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));

  r=p<<3; s=p; s<<=3; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
  r=p>>3; s=p; s>>=3; EEQ(r.deg(),s.deg());  FOR1(r,i) EDQ(r.c(i),s.c(i));
}


// Set/get
TEST(polynomial,io){
  init_polynomial_testdata();
  Polynomial p(P),q(Q);
  std::vector<double> v;
  ELE(q.deg(),p.deg());
  p.get(v);
  q.set(v);
  EEQ(p.deg(),q.deg());
  FOR1(p,i) ENR(q.c(i),p.c(i),1.e-5);
}


// Utilities
TEST(polynomial,util){
  init_polynomial_testdata();
  Polynomial p(P),q(Q), r;
  for(int i=20; i>=0; i--) r.c(i)=21-i;
  q=shrink(r, 10.5);  EEQ(q.deg(),10);  for(int i=q.deg(); i>=0; i--) EDQ(q.c(i),21-i);
  q=mirror(p);     for(double x=-1; x<1.01; x+=0.1) ENR(q(x),p(-x)      ,1.e-8);
  q=shift(p,0.7);  for(double x=-1; x<1.01; x+=0.1) ENR(q(x),p(x+0.7)   ,1.e-8);
  q=scale(p,0.7);  for(double x=-1; x<1.01; x+=0.1) ENR(q(x),p(x*0.7)   ,1.e-8);
  q=power(p,3);    for(double x=-1; x<1.01; x+=0.1) ENR(q(x),pow(p(x),3),1.e-8);
  q=diff(p);   EEQ(q.deg(),p.deg()-1);  FOR1(q,i) EDQ(q.c(i),p.c(i+1)*(i+1));
  q=integ(p);  EEQ(q.deg(),p.deg()+1);  FOR1(p,i) EDQ(q.c(i+1),p.c(i)/(i+1));
}


// Overall: Gaussian quadrature
TEST(polynomial,overall){
  init_polynomial_testdata();
  Polynomial p=legendre<double>(5), q=p, r=X, s;

  // nodes x_i and weights w_i...
  std::map<double,double> x2w;
  for(int i=0; i<p.deg(); i++){
    double x,w;
    x=solve(q,0.0,1.0,1.e-6);  // x_i = zero of Legendre
    q/=(r-x);
    s=(1.0-r*r)*power(diff(p),2)/2.0;
    w=1/s(x);  // w_i = 2/[(1-x^2) P'(x_i)^2]
    x2w[x]=w;

  }

  // integral of (n+1) x^n  on [-1,+1] - exact up to n=deg(P)*2-1
  for(int n=0; n<p.deg()*2; n++){
    s = (n+1.0) * (Z+1.0)<<n;  // (n+1) x^n
    double theoretical=integ(s)(1)-integ(s)(-1);
    double numerical=0;
    for(std::map<double,double>::const_iterator i=x2w.begin(); i!=x2w.end(); i++)
      numerical += s(i->first) * i->second;
    ENR(theoretical, numerical, 1.e-12);
  }
}
