// Unit test (matrix decomposition)


#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include "matrix.h"
#include "matrix_decomp.h"
using namespace std;
using namespace nmlib;


typedef double R;
typedef std::complex<R> C;
typedef matrix<R> MatR;
typedef matrix<C> MatC;


// Random sample
template<class T> T urand1(T /*c*/=0){
  return (random()+0.5)/RAND_MAX;
}
template<class T> std::complex<T> urand1(std::complex<T> /*c*/=0){
  return std::complex<T>(urand1<T>(),urand1<T>());
}
template<class T> matrix<T> urandm(size_t r, size_t c, T /*c*/=0){
  matrix<T> m(r,c);
  for(size_t k=0; k<m.dim(); k++) m(k)=urand1<T>(T(0));
  return m;
}
template<class T> matrix<T> urandm_sym(size_t r, size_t c, T /*c*/=0){
  matrix<T> m=urandm<T>(r,c);
  return (m+tp(m))/T(2);
}


template<class T> void verify_lu(const matrix<T>& m, const matrix<T>& l, const matrix<T>& u, const matrix<T>& p0, bool pivote){
  size_t n=m.nrow();
  matrix<T> p=(pivote ? p0 : matrix<T>(n,n)+T(1));
  EXPECT_NEAR(norm(l*u-p*m), 0, 1.e-8);
  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<i; j++){
      EXPECT_DOUBLE_EQ(std::abs(u(i,j)), 0);
      EXPECT_DOUBLE_EQ(std::abs(l(j,i)), 0);
    }
    EXPECT_DOUBLE_EQ(std::abs(l(i,i)-T(1)), 0);
  }
  EXPECT_NEAR(norm(tp(p)*p-T(1)), 0, 1.e-8);
  for(size_t i=0; i<n; i++)
    for(size_t j=0; j<n; j++)
      EXPECT_DOUBLE_EQ(std::abs(p(i,j)*(p(i,j)-T(1))), 0);
}


template<class T> void verify_qr(const matrix<T>& m, const matrix<T>& q, const matrix<T>& r){
  size_t n=m.nrow();
  EXPECT_NEAR(norm(q*r-m), 0, 1.e-8);
  EXPECT_NEAR(norm(tp(q)*q-T(1)), 0, 1.e-8);
  for(size_t i=0; i<n; i++)
    for(size_t j=0; j<i; j++)
      EXPECT_NEAR(std::abs(r(i,j)), 0, 1.e-8);
}


template<class T> void verify_svd(const matrix<T>& m, const matrix<T>& u, const matrix<T>& d, const matrix<T>& v){
  EXPECT_NEAR(norm(u*d*tp(v)-m), 0, 1.e-8);
  EXPECT_NEAR(norm(tp(u)*u-T(1)), 0, 1.e-8);
  EXPECT_NEAR(norm(tp(v)*v-T(1)), 0, 1.e-8);
  for(size_t i=0; i<d.nrow(); i++)
  for(size_t j=0; j<d.ncol(); j++)
    if(i!=j){  EXPECT_NEAR(std::abs(d(i,j)), 0, 1.e-8);  }
}


template<class T> void verify_eigen(const matrix<T>& m, const matrix<T>& u, const matrix<T>& d){
  size_t n=m.nrow();
  EXPECT_NEAR(norm(tp(u)*m*u-d), 0, 1.e-8);
  EXPECT_NEAR(norm(tp(u)*u-T(1)), 0, 1.e-8);
  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      if(i!=j){  EXPECT_NEAR(std::abs(d(i,j)), 0, 1.e-8);  }
      if(i==j){  EXPECT_NEAR(std::imag(d(i,j)), 0, 1.e-8);  }
    }
  }
}


TEST(matrix_decomp,lu){
  size_t n=10;

  matrix<R> mr0,mr,lr,ur,pr;
  mr=mr0=urandm<R>(n,n);
  lu_decomp (lr,ur,mr);           verify_lu(mr,lr,ur,pr,false);
  lup_decomp(lr,ur,pr,mr,false);  verify_lu(mr,lr,ur,pr,false);
  lup_decomp(lr,ur,pr,mr,true);   verify_lu(mr,lr,ur,pr,true);

  matrix<C> mc0,mc,lc,uc,pc;
  mc=mc0=urandm<C>(n,n);
  lu_decomp (lc,uc,mc);           verify_lu(mc,lc,uc,pc,false);
  lup_decomp(lc,uc,pc,mc,false);  verify_lu(mc,lc,uc,pc,false);
  lup_decomp(lc,uc,pc,mc,true);   verify_lu(mc,lc,uc,pc,true);

  EXPECT_DOUBLE_EQ(norm(mr-mr0)+norm(mc-mc0),0);
}


TEST(matrix_decomp,qr){
  size_t n=10;

  matrix<R> mr0,mr,qr,rr;
  mr=mr0=urandm<R>(n,n);
  qr_decomp(qr,rr,mr);  verify_qr(mr,qr,rr);

  matrix<C> mc0,mc,qc,rc;
  mc=mc0=urandm<C>(n,n);
  qr_decomp(qc,rc,mc);  verify_qr(mc,qc,rc);

  EXPECT_DOUBLE_EQ(norm(mr-mr0)+norm(mc-mc0),0);
}


TEST(matrix_decomp,svd){
  int n=10, n1=7, n2=12;

  matrix<R> mr0, mr,ur,dr,vr;
  mr=mr0=urandm<R>(n ,n );  svd_decomp(ur,dr,vr,mr);  verify_svd(mr,ur,dr,vr);
  mr=mr0=urandm<R>(n1,n2);  svd_decomp(ur,dr,vr,mr);  verify_svd(mr,ur,dr,vr);
  mr=mr0=urandm<R>(n2,n1);  svd_decomp(ur,dr,vr,mr);  verify_svd(mr,ur,dr,vr);

  matrix<C> mc0, mc,uc,dc,vc;
  mc=mc0=urandm<C>(n ,n );  svd_decomp(uc,dc,vc,mc);  verify_svd(mc,uc,dc,vc);
  mc=mc0=urandm<C>(n1,n2);  svd_decomp(uc,dc,vc,mc);  verify_svd(mc,uc,dc,vc);
  mc=mc0=urandm<C>(n2,n1);  svd_decomp(uc,dc,vc,mc);  verify_svd(mc,uc,dc,vc);

  EXPECT_DOUBLE_EQ(norm(mr-mr0)+norm(mc-mc0),0);
}


TEST(matrix_decomp,eigen_jacobi){
  int n=10;

  matrix<R> mr0,mr,ur,dr;
  mr=mr0=urandm_sym<R>(n,n);
  ur=eigen(mr);  dr=tp(ur)*mr*ur;  verify_eigen(mr,ur,dr);

  matrix<C> mc0,mc,uc,dc;
  mc=mc0=urandm_sym<C>(n,n);
  uc=eigen(mc);  dc=tp(uc)*mc*uc;  verify_eigen(mc,uc,dc);

  EXPECT_DOUBLE_EQ(norm(mr-mr0)+norm(mc-mc0),0);
}


TEST(matrix_decomp,eigen_qr){
  int n=10;

  matrix<R> mr0,mr,ur,dr;
  mr=mr0=urandm_sym<R>(n,n);
  eigen_qr(ur,dr,mr,true);               verify_eigen(mr,ur,dr);
  eigen_qr(ur,dr,mr,false,1.e-8,10000);  verify_eigen(mr,ur,dr);  // convergence may be slow

  matrix<C> mc0,mc,uc,dc;
  mc=mc0=urandm_sym<C>(n,n);
  eigen_qr(uc,dc,mc,true);               verify_eigen(mc,uc,dc);
  eigen_qr(uc,dc,mc,false,1.e-8,10000);  verify_eigen(mc,uc,dc);  // convergence may be slow

  EXPECT_DOUBLE_EQ(norm(mr-mr0)+norm(mc-mc0),0);
}
