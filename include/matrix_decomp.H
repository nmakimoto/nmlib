// Matrix decomposition tools - LU/LUP, QR, tridiagonalization


#ifndef MATRIX_DECOMP_H
#define MATRIX_DECOMP_H


#include "matrix.H"
namespace nmlib{


  /******** Prototype ********/

  // LU and LUP decomposition of A=(P)LU
  template<class T> void lu_decomp   (matrix<T>& l, matrix<T>& u,               const matrix<T>& a);
  template<class T> void lup_decomp  (matrix<T>& l, matrix<T>& u, matrix<T>& p, const matrix<T>& a, bool pivote=true);

  // QR decompoosition of A=QR
  template<class T> void qr_decomp   (matrix<T>& q, matrix<T>& r,               const matrix<T>& a);

  // Tridiagonalization D=U'AU of symmetric A
  template<class T> void tridiag     (matrix<T>& u, matrix<T>& d,               const matrix<T>& a, double tol=0);


  /******** Implementation ********/

  // LU decomposition of A (A=LU, L=unit lower triangular, U=upper triangluar)
  template<class T> void lu_decomp(matrix<T>& l, matrix<T>& u, const matrix<T>& a){
    matrix<T> p_dummy;
    lup_decomp(l,u,p_dummy,a,false);
  }


  // LUP decomposition with pivoting (PA=LU, P=pivoting of rows)
  template<class T> void lup_decomp(matrix<T>& l, matrix<T>& u, matrix<T>& p, const matrix<T>& a, bool pivote){
    size_t n=a.nrow();
    std::vector<size_t> perm(n);  // P(i,j)=1 if j=perm[i], 0 otherwise

    l=u=p=matrix<T>(n,n);
    for(size_t k=0; k<n; k++) perm[k]=k;

    for(size_t k=0; k<n; k++){
      // Pivoting (will not affect the already established portion of L and U)
      if(pivote){
	for(size_t k1=k; k1<n; k1++)
	  if(std::abs(a(perm[k],k))<std::abs(a(perm[k1],k))) std::swap(perm[k],perm[k1]);
      }
      p(k,perm[k])=1;

      // Determine upper left portion of L and U
      l(k,k)=T(1);
      for(size_t j=0; j<k; j++){  // L(k,j=0..k-1)
	T s=a(perm[k],j);  // (PA)(k,j)
	for(size_t i=0; i<j; i++) s-=l(k,i)*u(i,j);
	l(k,j)=s/u(j,j);  // (handle division by zero...)
      }
      for(size_t i=0; i<=k; i++){  // U(i=0..k,k)
	T s=a(perm[i],k);  // (PA)(i,k)
	for(size_t j=0; j<i; j++) s-=l(i,j)*u(j,k);
	u(i,k)=s;
      }
    }
  }


  // QR decomposition of M  (M=QR, Q:unitary, R:upper triangular)
  template<class T> void qr_decomp(matrix<T>& q, matrix<T>& r, const matrix<T>& a){
    size_t n=a.nrow();
    r=a;
    for(size_t k=0; k<n-1; k++){
      // v = Householder reflection direction
      matrix<T> b=getsub(r,k,k,n-k,1), v=b;
      v(0)+=norm(v)*sgn(v(0));  // v/|v| <--> ek
      v/=T(norm(v));  // (handle division by zero...)

      // reflection - h=1-2vv^T; q=qh; r=hr;
      matrix<T> qv=getsub(q,0,k,n,n-k)*v, vr=tp(v)*getsub(r,k,0,n-k,n);
      for(size_t i=0; i<n; i++){
	for(size_t j=0; j<n-k; j++){
	  q(i,k+j)-=T(2)*qv(i)*cnj(v(j));  // q=qh
	  r(k+j,i)-=T(2)*v(j)*vr(i);  // r=hr
	}
      }
    }

    // normalizization - r(i,i)>=0
    for(size_t i=0; i<n; i++){
      T s=sgn(r(i,i));
      for(size_t j=0; j<n; j++){
	r(i,j)*=cnj(s);
	q(j,i)*=s;
      }
    }
  }


  // Tridiagonalization D=U'AU of symmetric A (U:orthogonal, D:tridiagonal)
  template<class T> void tridiag(matrix<T>& u, matrix<T>& d, const matrix<T>& a, double tol){
    size_t n=a.nrow();
    u=matrix<T>(n,n)+T(1);
    d=a;

    for(size_t k=1; k<n-1; k++){
      // reflection vector v s.t. b:=lower (n-k) of (k-1)-th row <--> |b|ek
      matrix<T> v=getsub(d,k,k-1,n-k,1);
      v(0)+=norm(v)*sgn(v(0));
      if(std::abs(v(0))<=tol) continue;  // eigen
      v/=T(norm(v));

      // h=1-2vv^T; d=uh; d=hdh;
      matrix<T> uv=getsub(u,0,k,n,n-k)*v;
      for(size_t i=0; i<n; i++)
	for(size_t j=0; j<n-k; j++) u(i,k+j)-=T(2)*uv(i)*cnj(v(j));  // u=uh
      matrix<T> dv=getsub(d,0,k,n,n-k)*v;
      for(size_t i=0; i<n; i++)
	for(size_t j=0; j<n-k; j++) d(i,k+j)-=T(2)*dv(i)*cnj(v(j));  // d=dh
      matrix<T> vd=tp(v)*getsub(d,k,0,n-k,n);
      for(size_t i=0; i<n; i++)
	for(size_t j=0; j<n-k; j++) d(k+j,i)-=T(2)*v(j)*vd(i);  // d=hd
    }
  }


}  //namespace nmlib
#endif //MATRIX_DECOMP_H
