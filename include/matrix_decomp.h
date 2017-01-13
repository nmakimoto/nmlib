// Matrix decomposition tools - LU/LUP, QR, tridiagonalization, eigensolver


#ifndef MATRIX_DECOMP_H
#define MATRIX_DECOMP_H


#include "matrix.h"
namespace nmlib{


  /******** Prototype ********/

  // LU and LUP decomposition of A=(P)LU
  template<class T> void lu_decomp   (matrix<T>& l, matrix<T>& u,               const matrix<T>& a);
  template<class T> void lup_decomp  (matrix<T>& l, matrix<T>& u, matrix<T>& p, const matrix<T>& a, bool pivote=true);

  // QR decompoosition of A=QR
  template<class T> void qr_decomp   (matrix<T>& q, matrix<T>& r,               const matrix<T>& a);

  // Tridiagonalization D=U'AU of symmetric A
  template<class T> void tridiag     (matrix<T>& u, matrix<T>& d,               const matrix<T>& a, double tol=0);

  // Eigen solver D=U'AU of symmetric/tridiagonal A (QR algorithm)
  template<class T> void eigen_qr     (matrix<T>& u, matrix<T>& d,              const matrix<T>& a, bool shift=false, double tol=1.e-8, int iter=1000);
  template<class T> void eigen_tridiag(matrix<T>& u, matrix<T>& d,              const matrix<T>& a, bool shift=false, double tol=1.e-8, int iter=1000);

  // SVD(Singular Value Decomposition) of A=UDV^T
  template<class T> void svd_decomp(matrix<T>& u, matrix<T>& d, matrix<T>& v, const matrix<T>& a);


  /******** Implementation ********/

  // LU decomposition of A (A=LU, L=unit lower triangular, U=upper triangluar)
  template<class T> void lu_decomp(matrix<T>& l, matrix<T>& u, const matrix<T>& a){
    matrix<T> p_dummy;
    lup_decomp(l,u,p_dummy,a,false);
  }


  // LUP decomposition with pivoting (PA=LU, P=pivoting of rows)
  template<class T> void lup_decomp(matrix<T>& l, matrix<T>& u, matrix<T>& p, const matrix<T>& a, bool pivote){
    if(a.nrow()!=a.ncol()) throw std::domain_error("lup_decomp(): A is not square");
    size_t n=a.nrow();

    std::vector<size_t> perm(n);  // P(i,j)=1 if j=perm[i], 0 otherwise
    l=u=p=matrix<T>(n,n);
    for(size_t i=0; i<n; i++) perm[i]=i;

    for(size_t i=0; i<n; i++){
      // Pivote rows below i (will not affect the already established portion of L and U)
      if(pivote){
	for(size_t i1=i; i1<n; i1++)
	  if(std::abs(a(perm[i],i))<std::abs(a(perm[i1],i))) std::swap(perm[i],perm[i1]);
      }
      p(i,perm[i])=1;

      // Determine i-th row of L and U
      l(i,i)=T(1);
      for(size_t j=0; j<i; j++){  // L(i,j=0..i-1)
	T s=a(perm[i],j);  // (PA)(i,j)
	for(size_t k=0; k<j; k++) s-=l(i,k)*u(k,j);
	l(i,j)=s/u(j,j);  // (handle division by zero...)
      }
      for(size_t j=i; j<n; j++){  // U(i,j=i..n-1)
	T s=a(perm[i],j);  // (PA)(i,j)
	for(size_t k=0; k<i; k++) s-=l(i,k)*u(k,j);
	u(i,j)=s;
      }
    }
  }


  // QR decomposition of M  (M=QR, Q:unitary, R:upper triangular)
  template<class T> void qr_decomp(matrix<T>& q, matrix<T>& r, const matrix<T>& a){
    if(a.nrow()!=a.ncol()) throw std::domain_error("qr_decomp(): A is not square");
    size_t n=a.nrow();

    q=matrix<T>(n,n)+T(1);
    r=a;

    for(size_t k=0; k<n-1; k++){
      // v = Householder reflection direction
      matrix<T> b=getsub(r,k,k,n-k,1), v=b;
      v(0)+=norm(v)*nm_sign(v(0));  // v/|v| <--> ek
      v/=T(norm(v));  // (handle division by zero...)

      // reflection - h=1-2vv^T; q=qh; r=hr;
      matrix<T> qv=getsub(q,0,k,n,n-k)*v, vr=tp(v)*getsub(r,k,0,n-k,n);
      for(size_t i=0; i<n; i++){
	for(size_t j=0; j<n-k; j++){
	  q(i,k+j)-=T(2)*qv(i)*std::conj(v(j));  // q=qh
	  r(k+j,i)-=T(2)*v(j)*vr(i);  // r=hr
	}
      }
    }

    // normalizization - r(i,i)>=0
    for(size_t i=0; i<n; i++){
      T s=nm_sign(r(i,i));
      for(size_t j=0; j<n; j++){
	r(i,j)*=std::conj(s);
	q(j,i)*=s;
      }
    }
  }


  // Tridiagonalization D=U'AU of symmetric A (U:orthogonal, D:tridiagonal)
  template<class T> void tridiag(matrix<T>& u, matrix<T>& d, const matrix<T>& a, double tol){
    if(a.nrow()!=a.ncol()) throw std::domain_error("tridiag(): A is not square");
    size_t n=a.nrow();

    u=matrix<T>(n,n)+T(1);
    d=a;

    for(size_t k=1; k<n-1; k++){
      // reflection vector v s.t. b:=lower (n-k) of (k-1)-th row <--> |b|ek
      matrix<T> v=getsub(d,k,k-1,n-k,1);
      v(0)+=norm(v)*nm_sign(v(0));
      if(std::abs(v(0))<=tol) continue;  // eigen
      v/=T(norm(v));

      // h=1-2vv^T; d=uh; d=hdh;
      matrix<T> uv=getsub(u,0,k,n,n-k)*v;
      for(size_t i=0; i<n; i++)
	for(size_t j=0; j<n-k; j++) u(i,k+j)-=T(2)*uv(i)*std::conj(v(j));  // u=uh
      matrix<T> dv=getsub(d,0,k,n,n-k)*v;
      for(size_t i=0; i<n; i++)
	for(size_t j=0; j<n-k; j++) d(i,k+j)-=T(2)*dv(i)*std::conj(v(j));  // d=dh
      matrix<T> vd=tp(v)*getsub(d,k,0,n-k,n);
      for(size_t i=0; i<n; i++)
	for(size_t j=0; j<n-k; j++) d(k+j,i)-=T(2)*v(j)*vd(i);  // d=hd
    }
  }


  // Eigen solver D=U'AU of symmetric A (U:orthogonal, D:diagonal) (QR algorithm)
  template<class T> void eigen_qr(matrix<T>& u, matrix<T>& d, const matrix<T>& a, bool shift, double tol, int iter){
    matrix<T> u1,d1;
    tridiag(u1,d1,a,tol);
    eigen_tridiag(u,d,d1,shift,tol,iter);
    u=u1*u;
  }


  // Eigen solver D=U'AU of tridiagonal A (U:orthogonal, D:diagonal) (QR algorithm)
  template<class T> void eigen_tridiag(matrix<T>& u, matrix<T>& d, const matrix<T>& a, bool shift, double tol, int iter){
    if(a.nrow()!=a.ncol()) throw std::domain_error("eigen_tridiag(): A is not square");
    size_t n=a.nrow(), n1=n;

    u=matrix<T>(n,n)+T(1);
    d=a;

    for(int it=0; it<iter; it++){
      // shift
      T lam=0;
      if(shift){
	while(n1>1 && std::norm(d(n1-1,n1-2))+std::norm(d(n1-2,n1-1))<std::norm(tol)/n) n1--;
	lam=d(n1-1,n1-1);
      }
      d-=lam;

      // Givens rotation G by th=atan2(dji,dii)
      matrix<T> cc(n),ss(n);
      for(size_t i=0; i<n1-1; i++){
	size_t j=i+1;
	T c,s,th;  // c=|dii|/r, s= dji/r (dii/|dii|)~
	th=atan2(std::abs(d(j,i)),std::abs(d(i,i)));
	c=cc(i)=cos(th);
	s=ss(i)=sin(th) * nm_sign(std::conj(d(i,i))) * nm_sign(d(j,i));  // choose sign to avoid underflow and to keep G unitary
	for(size_t k=0; k<n1; k++){ T p=d(i,k),q=d(j,k); d(i,k)=c*p+std::conj(s)*q; d(j,k)=-s*p+c*q; }  // D=G*D
      }
      for(size_t i=0; i<n1-1; i++){
	size_t j=i+1;
	T c=cc(i),s=ss(i),p,q;
	for(size_t k=0; k<n; k++){ p=d(k,i); q=d(k,j); d(k,i)=c*p+s*q; d(k,j)=-std::conj(s)*p+c*q; }  // D=D*G^T  G=| c s~|
	for(size_t k=0; k<n; k++){ p=u(k,i); q=u(k,j); u(k,i)=c*p+s*q; u(k,j)=-std::conj(s)*p+c*q; }  // U=U*G^T    |-s c |
      }
      d+=lam;

      T err=0; for(size_t k=1; k<n1; k++) err+=std::norm(d(k-1,k))+std::norm(d(k,k-1));
      if( std::abs(err)<std::norm(tol)/n ) break;  // (note: not monotonically decreasing)
    }
  }


  // SVD(Singular Value Decomposition) of A=UDV^T (U,V:orthogonal, D:diagonal)
  template<class T> void svd_decomp(matrix<T>& u, matrix<T>& d, matrix<T>& v, const matrix<T>& a){
    if(a.nrow()<a.ncol()){ svd_decomp(v,d,u,tp(a)); return; }
    v=eigen(tp(a)*a);
    u=orth(a*v);
    d=tp(u)*a*v;
  }


}  //namespace nmlib
#endif //MATRIX_DECOMP_H
