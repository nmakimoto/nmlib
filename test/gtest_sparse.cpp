// Unit test (sparse matrix)


#include <gtest/gtest.h>

#include <cstdlib>  // rand
#include "sparse.h"
using namespace nmlib;


// Random sample generators
static double urand1(void    ){
  return (random()%RAND_MAX+0.5)/RAND_MAX;
}
static Matrix random_vector(size_t n){
  Matrix x(n);
  for(size_t i=0; i<n; i++) x(i)=2*urand1()-1;
  return x;
}
static Sparse random_sparse(size_t n, size_t w, size_t k){
  Sparse a(n,n);
  for(size_t i=0; i<n; i++)
    for(size_t j=0; j<n; j++)
      if( j-w<=i && i<=j+w ) a(i,j)=urand1()*2-1;  // w=band width
  for( ; k>0; k--)
    a(size_t(urand1()*n), size_t(urand1()*n)) = urand1()*2-1;  // random position
  return a;
}


// Access to components
TEST(sparse,component){
  size_t n=10, w=2;
  Sparse a=random_sparse(n,w,n*4)+10.0;
  Matrix ma=dense(a);
  EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
  EXPECT_NEAR(norm(a), norm(ma), 1.e-8);

  for(size_t k=0; k<n*4; k++){
    size_t i=size_t(urand1()*n), j=size_t(urand1()*n);
    double c=urand1()*2-1;
    a(i,j)=ma(i,j)=c;
  }
  EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
}


// Arithmetic operations
TEST(sparse,arithmetic){
  size_t n=10, w=2;
  Sparse a, a2;
  Matrix ma, ma2, b;

  a =random_sparse(n,w,n*4)+10.0;  ma =dense(a );
  a2=random_sparse(n,w,n*6)+10.0;  ma2=dense(a2);
  b =random_vector(n);
  double t=2;

  a+=t;   ma+=t;    EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
  a-=t;   ma-=t;    EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
  a*=t;   ma*=t;    EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
  a/=t;   ma/=t;    EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
  a+=a2;  ma+=ma2;  EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
  a-=a2;  ma-=ma2;  EXPECT_TRUE(norm(dense(a)-ma)<1.e-8);
  EXPECT_TRUE(norm(dense(-a)-(-ma))<1.e-8);
  EXPECT_TRUE(norm(dense(a+t)-(ma+t))<1.e-8);
  EXPECT_TRUE(norm(dense(a-t)-(ma-t))<1.e-8);
  EXPECT_TRUE(norm(dense(a*t)-(ma*t))<1.e-8);
  EXPECT_TRUE(norm(dense(a/t)-(ma/t))<1.e-8);
  EXPECT_TRUE(norm(dense(t+a)-(t+ma))<1.e-8);
  EXPECT_TRUE(norm(dense(t-a)-(t-ma))<1.e-8);
  EXPECT_TRUE(norm(dense(t*a)-(t*ma))<1.e-8);
  EXPECT_TRUE(norm(dense(a+a2)-(ma+ma2))<1.e-8);
  EXPECT_TRUE(norm(dense(a-a2)-(ma-ma2))<1.e-8);
  EXPECT_TRUE(norm(dense(a*a2)-(ma*ma2))<1.e-8);
  EXPECT_TRUE(norm((a*b)-(ma*b))<1.e-8);
  EXPECT_TRUE(norm(dense(tp(a))-tp(ma))<1.e-8);
  EXPECT_TRUE(norm(tpab(a,b)-tp(ma)*b)<1.e-8);
}


// Preconditioner (ILU0)
TEST(sparse,preconditioner){
  size_t n=10, w=2;
  Sparse a;
  Matrix ma,ml,mu;
  SparseConf sc;

  a=random_sparse(n,w,0)+10.0;  // band matrix
  ma=dense(a);

  sc.init_ilu0(a);  // ILU0 info is stored in sc.lu and retrieved below
  ml=mu=Matrix(n,n);
  for(size_t i=0; i<n; i++)
    for(size_t j=0; j<n; j++)
      if(i<=j) mu(i,j)=sc.lu(i,j);
      else     ml(i,j)=sc.lu(i,j);
  ml+=1.0;

  // ILU0 related operations
  Matrix v=random_vector(n);
  EXPECT_NEAR(norm( lux  (sc.lu,v)-       ml*mu  *v), 0, 1.e-8);
  EXPECT_NEAR(norm(sluxb (sc.lu,v)-inv(   ml*mu) *v), 0, 1.e-8);
  EXPECT_NEAR(norm(slutxb(sc.lu,v)-inv(tp(ml*mu))*v), 0, 1.e-8);

  // ILU0 result
  EXPECT_NEAR(norm(ma-ml*mu)/norm(ma), 0, 1.e-8);  // |A-LU|/|A|=0
  EXPECT_NEAR(norm(ma*inv(ml*mu)-1.0)/sqrt(n*1.0), 0, 1.e-8);  // |A/(LU)-I|/|I|=0
}


// Solver (CG/BCG/PBCG)
TEST(sparse,solver){
  size_t n=100, w=2;
  Sparse a;
  Matrix b,x;

  a=random_sparse(n,w,4*n)+10.0;
  b=random_vector(n);

  x=solve_bcg (a,b);  EXPECT_NEAR(norm(a*x-b)/norm(b), 0, 1.e-8);
  x=solve_pbcg(a,b);  EXPECT_NEAR(norm(a*x-b)/norm(b), 0, 1.e-8);
  a+=tp(a);
  x=solve_cg  (a,b);  EXPECT_NEAR(norm(a*x-b)/norm(b), 0, 1.e-8);
}
