// Sparse matrix and solvers (CG/BiCG/PBCG)


#include <iostream>
#include <map>
#include <stdexcept>
#include <sstream>
#include "sparse.H"
#include "matrix.H"
//#include "io.H"
namespace nmlib{


// Abbreviations
#define FOR_K(M)    for(size_t K=0; K<M.ncol(); K++)
#define FOR_IJ(M)   for(size_t I=0; I<M.nrow(); I++)  for(CI C=M.columns(I).begin(); C!=M.columns(I).end(); C++)
#define FOR_IJ2(M)  for(size_t I=0; I<M.nrow(); I++)  for(IT C=M.columns(I).begin(); C!=M.columns(I).end(); C++)
#define J     C->first
#define M_IJ  C->second

// Abbreviations
typedef Sparse S;
typedef Matrix M;
typedef double T;
typedef std::map<size_t,T> Row;
typedef Row::const_iterator  CI;
typedef Row::iterator        IT;


/******************************** Range checks ********************************/

// Range check for each operation
void chk(size_t r1, size_t c1, size_t r2, size_t c2, char op){
  if( (op=='@' and 0<=r1 and r1<r2 and 0<=c1 and c1<=c2) or
      (op=='+' and r1==r2 and c1==c2)  or
      (op=='-' and r1==r2 and c1==c2)  or
      (op=='*' and c1==r2)  or
      (op=='/' and r1==c1 and c1==r2 and c2==1) ) return;
  std::stringstream msg;
  msg<<"Sparse: ("<<r1<<","<<c1<<") "<< op << " ("<<r2<<","<<c2<<")";
  throw std::domain_error(msg.str());
}
void chk(const S& a, const S& b, char op){ chk(a.nrow(),a.ncol(),b.nrow(),b.ncol(),op); }
void chk(const S& a, const M& b, char op){ chk(a.nrow(),a.ncol(),b.nrow(),b.ncol(),op); }


/*********************************** Basic operations ***********************************/

// Class methods
S::Sparse(size_t r, size_t c): row(r),col(c) { if(r<0 or c<0) throw std::domain_error("Sparse(r,c): size<0"); val=std::vector<Row>(r); }
size_t  S::nrow(void) const{ return row; }
size_t  S::ncol(void) const{ return col; }
T&   S::operator()(size_t i, size_t j)       { chk(i,j,row,col,'@'); return val[i][j]; }
T    S::operator()(size_t i, size_t j) const { chk(i,j,row,col,'@'); CI c=val[i].find(j); return (c==val[i].end() ? 0 : c->second); }

void SparseConf::init_ilu0(const S& a){ lu=a; ilu0(lu); }


// Incremental operations
S& operator+=(S& a, T s){ FOR_K(a) a(K,K)+=s; return a; }
S& operator-=(S& a, T s){ a+=-s; return a; }
S& operator*=(S& a, T s){ FOR_IJ2(a) M_IJ*=s; return a; }
S& operator/=(S& a, T s){ a*=1/s; return a; }
S& operator+=(S& a, const S& b){ chk(a,b,'+'); FOR_IJ(b) a(I,J)+=M_IJ; return a; }
S& operator-=(S& a, const S& b){ chk(a,b,'-'); FOR_IJ(b) a(I,J)-=M_IJ; return a; }


// Basic operations
S operator-(const S& a)     { S b(a); b*=-1.0; return b; }
S operator+(const S& a, T s){ S b(a); b+=s; return b; }
S operator-(const S& a, T s){ S b(a); b-=s; return b; }
S operator*(const S& a, T s){ S b(a); b*=s; return b; }
S operator/(const S& a, T s){ S b(a); b/=s; return b; }
S operator+(T s, const S& a){ S b(a); b+=s; return b; }
S operator-(T s, const S& a){ S b(a); b*=-1.0; b+=s; return b; }
S operator*(T s, const S& a){ S b(a); return (b*=s); }
//S operator/(T s, const S& a);
S operator+(const S& a, const S& b){ S c(a); c+=b; return c; }
S operator-(const S& a, const S& b){ S c(a); c-=b; return c; }
S operator*(const S& a, const S& b){ chk(a,b,'*'); S c(a.nrow(),b.ncol()); FOR_K(c) FOR_IJ(a) c(I,K)+=M_IJ*b(J,K); return c; }
M operator*(const S& a, const M& b){ chk(a,b,'*'); M c(a.nrow(),b.ncol()); FOR_K(c) FOR_IJ(a) c(I,K)+=M_IJ*b(J,K); return c; }


// Utilities
S tp   (const S& a){ S b(a.ncol(),a.nrow()); FOR_IJ(a) b(J,I)=M_IJ; return b; }
M dense(const S& a){ M b(a.nrow(),a.ncol()); FOR_IJ(a) b(I,J)=M_IJ; return b; }
M tpab (const S& a, const M& b){ chk(a.ncol(),a.nrow(),b.nrow(),b.ncol(),'*'); M c(b.nrow(),b.ncol()); FOR_K(c) FOR_IJ(a) c(J,K)+=M_IJ*b(I,K); return c; }


/******************************** Sparse solvers ********************************/

// Solve Ax=b (CG - Conjugate Gradient for symmetric positive definite A)
M solve_cg(const S& a, const M& b, const SparseConf& cf){
  chk(a,b,'/');

  M x,r,p,q;
  T t,rr;

  x=M(a.ncol());  if(cf.x0.dim()>0) x=cf.x0;
  p=r=b-a*x;
  for(int j=0; j<cf.loop; j++){
    if(cf.verb) std::cerr<<"CG["<<j<<"]\t"<<norm(r)/norm(b)<<"\t"<<norm(a*x-b)<<"\n";
    if(norm(r)<=norm(b)*cf.tol) break;  // convergence
    rr=inner(r,r);
    q=a*p;
    t=rr/inner(p,q);
    x+=t*p;
    r-=t*q;  // r=b-Ax
    t=inner(r,r)/rr;
    p=r+t*p;
  }
  return x;
}


// Solve Ax=b (BiCG - Bi-Conjugate Gradient)
M solve_bcg(const S& a, const M& b, const SparseConf& cf){
  chk(a,b,'/');

  M x,r1,r2,p1,p2,q1,q2;
  T t,rr;

  x=M(a.ncol());  if(cf.x0.dim()>0) x=cf.x0;
  p1=r1=b-a*x;
  p2=r2=r1;
  for(int j=0; j<cf.loop; j++){
    if(cf.verb) std::cerr<<"BCG["<<j<<"]\t"<<norm(r1)/norm(b)<<"\t"<<norm(a*x-b)<<"\n";
    if(norm(r1)<=norm(b)*cf.tol) break;
    rr=inner(r1,r2);
    q1=a*p1;
    q2=tpab(a,p2);  // A^T*p2
    t =rr/inner(q1,p2);
    x +=t*p1;
    r1-=t*q1;  // r1=b-Ax
    r2-=t*q2;
    t =inner(r1,r2)/rr;
    p1=r1+t*p1;
    p2=r2+t*p2;
  }
  return x;
}


// Solve Ax=b (PBCG - Preconditioned BiCG)
M solve_pbcg(const S& a, const M& b, const SparseConf& cf){
  chk(a,b,'/');

  M x,r1,r2,p1,p2,q1,q2;
  T t,rr;

  Sparse lu0;  if(cf.lu.nrow()==0){ lu0=a; ilu0(lu0); }
  const Sparse& lu=(cf.lu.nrow()>0 ? cf.lu : lu0);

  x=M(a.ncol());  if(cf.x0.dim()>0) x=cf.x0;
  p1=r1=sluxb(lu,b-a*x);  //lu.prod_inv  (b-a*x);
  p2=r2=slutxb(lu,b-a*x);  //lu.prod_tpinv(b-a*x);

  for(int j=0; j<cf.loop; j++){
    if(cf.verb) std::cerr<<"PBCG["<<j<<"]\t"<<norm(a*x-b)/norm(b)<<"\t"<<norm(a*x-b)<<"\t"<<"\n";
    if(norm(a*x-b)<=norm(b)*cf.tol) break;
    rr=inner(r1,r2);
    q1=sluxb(lu,a*p1);  //lu.prod_inv(a*p1);                 // (C^-1 A) * p1
    q2=tpab(a,slutxb(lu,p2));  //tpab(a,lu.prod_tpinv(p2));  // (C^-1 A)^T * p2
    t =rr/inner(q1,p2);
    x +=t*p1;
    r1-=t*q1;  // r1=C^-1 (b-Ax)
    r2-=t*q2;
    t =inner(r1,r2)/rr;
    p1=r1+t*p1;
    p2=r2+t*p2;
  }
  return x;
}


/******************************* Pre-conditioning *******************************/

// b=LUx
Matrix lux(const Sparse& lu, const Matrix&  x){
  // y=Ux
  Matrix y(x.nrow());
  for(size_t i=0; i<lu.nrow(); i++)
    for(CI c1=lu.columns(i).find(i),c2=lu.columns(i).end(); c1!=c2; c1++)
      y(i) += c1->second * x(c1->first);

  // z=Ly=LUx
  Matrix z(y);  // Lii=1
  for(size_t i=0; i<lu.nrow(); i++)
    for(CI c1=lu.columns(i).begin(),c2=lu.columns(i).find(i); c1!=c2; c1++)
      z(i) +=  c1->second * y(c1->first);

  return z;
}


// (LU)^-1 b (solve LUx=b)
Matrix sluxb(const Sparse& lu, const Matrix& b){
  Matrix x(b);

  // x --> L^-1 x
  for(size_t i=0; i<lu.nrow(); i++){
    CI c1=lu.columns(i).begin(), c2=lu.columns(i).find(i);
    for(; c1!=c2; c1++) x(i) -= c1->second * x(c1->first);
    //x(i)/=1;
  }

  // x --> U^-1 x
  for(size_t i=lu.nrow()-1; /*i>=0*/; i--){
    CI c1=lu.columns(i).find(i), c2=lu.columns(i).end();
    for(c1++; c1!=c2; c1++) x(i) -= c1->second * x(c1->first);
    x(i)/=lu(i,i);  // = (bi - \sum_{j>i} Uij xj) / Uii
    if(i==0) break;
  }

  return x;
}


// (LU)^-T b (solve LU^T x=b)
Matrix slutxb(const Sparse& lu, const Matrix& b){
  Matrix x(b);

  // x --> U^-T x
  for(size_t i=0; i<lu.ncol(); i++){
    CI c1=lu.columns(i).find(i), c2=lu.columns(i).end();
    x(i) /= c1->second; ;  // xi /= Uii
    for(c1++; c1!=c2; c1++) x(c1->first) -= c1->second * x(i);  // xj -= Uji*xi
  }

  // x --> L^-T x
  for(size_t i=lu.ncol()-1; /*i>=0*/; i--){
    CI c1=lu.columns(i).begin(), c2=lu.columns(i).find(i);
    //x(i)/=1;  // xi /== Lii
    for( ; c1!=c2; c1++) x(c1->first) -= c1->second * x(i);  // xj -= Lji*xi
    if(i==0) break;
  }

  return x;
}


// ILU(0) - incomplete LU decomposition
// - A typical preconditioner without fill-in
// - Decomposition is complete when A is a band matrix
void ilu0(Sparse& a){
  for(size_t i=0; i<a.nrow(); i++){
    std::map<size_t,double>& ai=a.columns(i);
    IT aii=ai.find(i);  // Uii
    if( aii==ai.end() ) throw std::runtime_error("ILU(0): division by Aii=0");

    for(IT aik=ai.begin(); aik!=aii; aik++){
      size_t k=aik->first;
      Row& ak=a.columns(k);
      aik->second /= a(k,k);  // Lik

      IT aij=aik;
      for(aij++; aij!=ai.end(); aij++){
	size_t j=aij->first;
	IT akj=ak.find(j);
	if( akj!=ak.end() ) aij->second -= aik->second * akj->second;  // Uij
      }
    }
  }
}


}  //namespace nmlib
