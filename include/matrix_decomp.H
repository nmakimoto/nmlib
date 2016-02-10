// Matrix decomposition tools - LU/LUP


#ifndef MATRIX_DECOMP_H
#define MATRIX_DECOMP_H


#include "matrix.H"
namespace nmlib{


  /******** Prototype ********/

  // LU and LUP decomposition of A=(P)LU
  template<class T> void lu_decomp   (matrix<T>& l, matrix<T>& u,               const matrix<T>& a);
  template<class T> void lup_decomp  (matrix<T>& l, matrix<T>& u, matrix<T>& p, const matrix<T>& a, bool pivote=true);


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


}  //namespace nmlib
#endif //MATRIX_DECOMP_H
