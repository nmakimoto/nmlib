// Sparse matrix and solvers (CG/BiCG/PBCG)


#ifndef SPARSE_H
#define SPARSE_H


#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include "matrix.h"

namespace nmlib{


/******** Class I/F ********/

// Sparse matrices
class Sparse{
public:
  explicit Sparse(size_t r=0, size_t c=0);

  size_t nrow(void) const;  // number of rows
  size_t ncol(void) const;  // number of columns
  double  operator()(size_t i, size_t j) const;  // get s=A(i,j)
  double& operator()(size_t i, size_t j);        // set A(i,j)=s

  const std::map<size_t,double>& columns(size_t i) const{ return val[i]; }  // i-th row (map j=>Aij)
        std::map<size_t,double>& columns(size_t i)      { return val[i]; }  // i-th row (non-const version)

private:
  size_t row, col;  // dimentions
  std::vector<std::map<size_t,double> > val;  // componets val[i][j]  CRS(Condensed Row Storage)
};


// Configurations for sparse solvers
class SparseConf{
public:
  SparseConf(void): tol(1.e-8), loop(100), x0(), lu(), verb(false) {}
  void init_ilu0(const Sparse& a);  // initialize preconditioner(ILU0)

  double  tol;   // tolerance |Ax-b|/|b|
  int     loop;  // max iteration
  Matrix  x0;    // initial value
  Sparse  lu;    // preconditioner
  bool    verb;  // verbosity (show progress to std::cerr)
};


// Configuration for sparse eigensolver
struct SparseEigenConf{
  SparseEigenConf(void): x0(), mu0(0), tol(1.e-8), loop(100), warmup(0), pbcg(false), verb(false) {}

  Matrix x0;    // initial guess of an eigenvector
  double mu0;   // initial guess of an eigenvalue
  double tol;   // tolerance for convergence check
  int    loop;  // max iteration
  int    warmup;// number of iteration to magnify eigens near mu0
  bool   pbcg;  // use PBCG with ILU(0) preconditioner
  bool   verb;  // verbosity
};


/******** Utility I/F ********/

// Incremental operations
Sparse& operator+=(Sparse& a, double s);   // A=A+sI
Sparse& operator-=(Sparse& a, double s);   // A=A-sI
Sparse& operator*=(Sparse& a, double s);   // A=A*s
Sparse& operator/=(Sparse& a, double s);   // A=A/s
Sparse& operator+=(Sparse& a, const Sparse& b);  // A=A+B
Sparse& operator-=(Sparse& a, const Sparse& b);  // A=A-B
//Sparse& operator*=(Sparse& a, const Sparse& b);
//Sparse& operator/=(Sparse& a, const Sparse& b);

// Scalar operations
Sparse operator-(const Sparse& a);            // -A (unary)
Sparse operator+(const Sparse& a, double s);  // A+sI
Sparse operator-(const Sparse& a, double s);  // A-sI
Sparse operator*(const Sparse& a, double s);  // A*s
Sparse operator/(const Sparse& a, double s);  // A/s
Sparse operator+(double s, const Sparse& a);  // sI+A
Sparse operator-(double s, const Sparse& a);  // sI-A
Sparse operator*(double s, const Sparse& a);  // s*A
//Sparse operator/(double s, const Sparse& a);  // s*inv(A)

// Sparse-Sparse operations
Sparse operator+(const Sparse& a, const Sparse& b);  // A+B
Sparse operator-(const Sparse& a, const Sparse& b);  // A-B
Sparse operator*(const Sparse& a, const Sparse& b);  // A*B
Matrix operator*(const Sparse& a, const Matrix& b);  // A*b (dense Matrix version)

// Utilities
double norm (const Sparse& a);  // |A|
Sparse tp   (const Sparse& a);  // A^T
Matrix dense(const Sparse& a);  // conversion to dens matrix - use with care
Matrix tpab (const Sparse& a, const Matrix& b);  // A^T * b

// I/O
std::istream& operator>>(std::istream& s,       Sparse& a);  // text I/O
std::ostream& operator<<(std::ostream& s, const Sparse& a);  // text I/O
std::istream& read      (std::istream& s,       Sparse& a);  // binary I/O
std::ostream& write     (std::ostream& s, const Sparse& a);  // binary I/O

// Sparse solvers
Matrix solve_cg  (const Sparse& a, const Matrix& b, const SparseConf& cf=SparseConf());  // CG (for symmetric positive definit A)
Matrix solve_bcg (const Sparse& a, const Matrix& b, const SparseConf& cf=SparseConf());  // BiCG
Matrix solve_pbcg(const Sparse& a, const Matrix& b, const SparseConf& cf=SparseConf());  // Preconditioned BiCG

// Preconditioners
void   ilu0(Sparse& a);  // incomplete LU decomposition ILU(0) (overwrites A)
Matrix lux   (const Sparse& lu, const Matrix& x);  // LUx
Matrix sluxb (const Sparse& lu, const Matrix& b);  // (LU)^-1 b (solve LUx=b)
Matrix slutxb(const Sparse& lu, const Matrix& b);  // (LU)^-T b (solve LU^T x=b)

// Sparse eigensolver
Matrix eigenvec(const Sparse& a0, const SparseEigenConf& cf1=SparseEigenConf(), const SparseConf& cf20=SparseConf());


/******** Implementation ********/

// Class methods
inline Sparse::Sparse(size_t r, size_t c): row(r),col(c),val(r) {}
inline size_t  Sparse::nrow(void) const{ return row; }
inline size_t  Sparse::ncol(void) const{ return col; }
inline double& Sparse::operator()(size_t i, size_t j){
  if( !(i<nrow() && j<ncol()) ) throw std::domain_error("Sparse::operator(): out of range");
  return val[i][j];
}
inline double  Sparse::operator()(size_t i, size_t j) const {
  if( !(i<nrow() && j<ncol()) ) throw std::domain_error("Sparse::operator() const: out of range");
  const auto& c=val[i].find(j);
  return (c==val[i].end() ? 0 : c->second);
}

inline void SparseConf::init_ilu0(const Sparse& a){ lu=a; ilu0(lu); }


// Incremental operations
inline Sparse& operator+=(Sparse& a, double s){
  for(size_t k=0; k<a.ncol(); k++) a(k,k)+=s;
  return a;
}
inline Sparse& operator-=(Sparse& a, double s){
  a+=-s;
  return a;
}
inline Sparse& operator*=(Sparse& a, double s){
  for(size_t i=0; i<a.nrow(); i++)
    for(auto& c: a.columns(i))
      c.second*=s;
  return a;
}
inline Sparse& operator/=(Sparse& a, double s){
  a*=1/s;
  return a;
}
inline Sparse& operator+=(Sparse& a, const Sparse& b){
  if( !(a.nrow()==b.nrow() && a.ncol()==b.ncol()) )
    throw std::domain_error("operator +=(Sparse,Sparse): invalid dimensions");
  for(size_t i=0; i<b.nrow(); i++)
    for(const auto& c: b.columns(i))
      a(i,c.first)+=c.second;
  return a;
}
inline Sparse& operator-=(Sparse& a, const Sparse& b){
  if( !(a.nrow()==b.nrow() && a.ncol()==b.ncol()) )
    throw std::domain_error("operator -=(Sparse,Sparse): invalid dimensions");
  for(size_t i=0; i<b.nrow(); i++)
    for(const auto& c: b.columns(i))
      a(i,c.first)-=c.second;
  return a;
}


// Basic operations
inline Sparse operator-(const Sparse& a)          { Sparse b(a); b*=-1.0; return b; }
inline Sparse operator+(const Sparse& a, double s){ Sparse b(a); b+=s; return b; }
inline Sparse operator-(const Sparse& a, double s){ Sparse b(a); b-=s; return b; }
inline Sparse operator*(const Sparse& a, double s){ Sparse b(a); b*=s; return b; }
inline Sparse operator/(const Sparse& a, double s){ Sparse b(a); b/=s; return b; }
inline Sparse operator+(double s, const Sparse& a){ Sparse b(a); b+=s; return b; }
inline Sparse operator-(double s, const Sparse& a){ Sparse b(a); b*=-1.0; b+=s; return b; }
inline Sparse operator*(double s, const Sparse& a){ Sparse b(a); return (b*=s); }
//inline Sparse operator/(double s, const Sparse& a);
inline Sparse operator+(const Sparse& a, const Sparse& b){ Sparse c(a); c+=b; return c; }
inline Sparse operator-(const Sparse& a, const Sparse& b){ Sparse c(a); c-=b; return c; }
inline Sparse operator*(const Sparse& a, const Sparse& b){
  if( a.ncol()!=b.nrow() )
    throw std::domain_error("operator *(Sparse,Sparse): invalid dimensions");
  Sparse ret(a.nrow(),b.ncol());
  for(size_t j=0; j<ret.ncol(); j++)
    for(size_t i=0; i<a.nrow(); i++)
      for(const auto& c: a.columns(i))
	ret(i,j)+=c.second*b(c.first,j);
  return ret;
}
inline Matrix operator*(const Sparse& a, const Matrix& b){
  if( a.ncol()!=b.nrow() )
    throw std::domain_error("operator *(Sparse,Matrix): invalid dimensions");
  Matrix ret(a.nrow(),b.ncol());
  for(size_t j=0; j<ret.ncol(); j++)
    for(size_t i=0; i<a.nrow(); i++)
      for(const auto& c: a.columns(i))
	ret(i,j)+=c.second*b(c.first,j);
  return ret;
}


// Utilities
inline double norm (const Sparse& a){
  double s=0;
  for(size_t i=0; i<a.nrow(); i++){
    for(const auto& c: a.columns(i)){
      double t=c.second;
      s+=t*t;
    }
  }
  return sqrt(s);
}
inline Sparse tp   (const Sparse& a){
  Sparse b(a.ncol(),a.nrow());
  for(size_t i=0; i<a.nrow(); i++)
    for(const auto& c: a.columns(i))
      b(c.first,i)=c.second;
  return b;
}
inline Matrix dense(const Sparse& a){
  Matrix b(a.nrow(),a.ncol());
  for(size_t i=0; i<a.nrow(); i++)
    for(const auto& c: a.columns(i))
      b(i,c.first)=c.second;
  return b;
}
inline Matrix tpab (const Sparse& a, const Matrix& b){
  if( a.nrow()!=b.nrow() )
    throw std::domain_error("tpab(Sparse,Matrix): invalid dimensions");
  Matrix ret(b.nrow(),b.ncol());
  for(size_t j=0; j<ret.ncol(); j++)
    for(size_t i=0; i<a.nrow(); i++)
      for(const auto& c: a.columns(i))
	ret(c.first,j)+=c.second*b(i,j);
  return ret;
}


// Solve Ax=b (CG - Conjugate Gradient for symmetric positive definite A)
inline Matrix solve_cg(const Sparse& a, const Matrix& b, const SparseConf& cf){
  if( !(a.nrow()==a.ncol() && a.ncol()==b.nrow() && b.ncol()==1) )
    throw std::domain_error("solve_cg(Sparse,Matrix): invalid dimensions");

  Matrix x,r,p,q;
  double t,rr;

  x=(cf.x0.dim()>0 ? cf.x0 : Matrix(a.ncol()));
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
inline Matrix solve_bcg(const Sparse& a, const Matrix& b, const SparseConf& cf){
  if( !(a.nrow()==a.ncol() && a.ncol()==b.nrow() && b.ncol()==1) )
    throw std::domain_error("solve_bcg(Sparse,Matrix): invalid dimensions");

  Matrix x,r1,r2,p1,p2,q1,q2;
  double t,rr;

  x=Matrix(a.ncol());  if(cf.x0.dim()>0) x=cf.x0;
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
inline Matrix solve_pbcg(const Sparse& a, const Matrix& b, const SparseConf& cf){
  if( !(a.nrow()==a.ncol() && a.ncol()==b.nrow() && b.ncol()==1) )
    throw std::domain_error("solve_pbcg(Sparse,Matrix): invalid dimensions");

  Matrix x,r1,r2,p1,p2,q1,q2;
  double t,rr;

  Sparse lu0;  if(cf.lu.nrow()==0){ lu0=a; ilu0(lu0); }
  const Sparse& lu=(cf.lu.nrow()>0 ? cf.lu : lu0);

  x=Matrix(a.ncol());  if(cf.x0.dim()>0) x=cf.x0;
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


// ILU(0) - incomplete LU decomposition
// - A typical preconditioner without fill-in
// - Decomposition is complete when A is a band matrix
inline void ilu0(Sparse& a){
  if( a.nrow()!=a.ncol() ) throw std::domain_error("ilu0(Sparse): invalid dimensions");

  for(size_t i=0; i<a.nrow(); i++){
    std::map<size_t,double>& ai=a.columns(i);
    auto aii=ai.find(i);  // Uii
    if( aii==ai.end() ) throw std::runtime_error("ilu0(Sparse): division by Aii=0");

    for(auto aik=ai.begin(); aik!=aii; aik++){
      size_t k=aik->first;
      auto& ak=a.columns(k);
      aik->second /= a(k,k);  // Lik

      auto aij=aik;
      for(aij++; aij!=ai.end(); aij++){
	size_t j=aij->first;
	auto akj=ak.find(j);
	if( akj!=ak.end() ) aij->second -= aik->second * akj->second;  // Uij
      }
    }
  }
}


// b=LUx
inline Matrix lux(const Sparse& lu, const Matrix&  x){
  if( !(lu.nrow()==lu.ncol() && lu.ncol()==x.nrow() && x.ncol()==1) )
    throw std::domain_error("lux(Sparse,Matrix): invalid dimensions");

  // y=Ux
  Matrix y(x.nrow());
  for(size_t i=0; i<lu.nrow(); i++)
    for(auto c1=lu.columns(i).find(i),c2=lu.columns(i).end(); c1!=c2; c1++)
      y(i) += c1->second * x(c1->first);

  // z=Ly=LUx
  Matrix z(y);  // Lii=1
  for(size_t i=0; i<lu.nrow(); i++)
    for(auto c1=lu.columns(i).begin(),c2=lu.columns(i).find(i); c1!=c2; c1++)
      z(i) +=  c1->second * y(c1->first);

  return z;
}


// (LU)^-1 b (solve LUx=b)
inline Matrix sluxb(const Sparse& lu, const Matrix& b){
  if( !(lu.nrow()==lu.ncol() && lu.ncol()==b.nrow() && b.ncol()==1) )
    throw std::domain_error("sluxb(Sparse,Matrix): invalid dimensions");

  Matrix x(b);

  // x --> L^-1 x
  for(size_t i=0; i<lu.nrow(); i++){
    auto c1=lu.columns(i).begin(), c2=lu.columns(i).find(i);
    for(; c1!=c2; c1++) x(i) -= c1->second * x(c1->first);
    //x(i)/=1;
  }

  // x --> U^-1 x
  for(size_t i=lu.nrow()-1; /*i>=0*/; i--){
    auto c1=lu.columns(i).find(i), c2=lu.columns(i).end();
    for(c1++; c1!=c2; c1++) x(i) -= c1->second * x(c1->first);
    x(i)/=lu(i,i);  // = (bi - \sum_{j>i} Uij xj) / Uii
    if(i==0) break;
  }

  return x;
}


// (LU)^-T b (solve LU^T x=b)
inline Matrix slutxb(const Sparse& lu, const Matrix& b){
  if( !(lu.nrow()==lu.ncol() && lu.ncol()==b.nrow() && b.ncol()==1) )
    throw std::domain_error("slutxb(Sparse,Matrix): invalid dimensions");

  Matrix x(b);

  // x --> U^-T x
  for(size_t i=0; i<lu.ncol(); i++){
    auto c1=lu.columns(i).find(i), c2=lu.columns(i).end();
    x(i) /= c1->second; ;  // xi /= Uii
    for(c1++; c1!=c2; c1++) x(c1->first) -= c1->second * x(i);  // xj -= Uji*xi
  }

  // x --> L^-T x
  for(size_t i=lu.ncol()-1; /*i>=0*/; i--){
    auto c1=lu.columns(i).begin(), c2=lu.columns(i).find(i);
    //x(i)/=1;  // xi /== Lii
    for( ; c1!=c2; c1++) x(c1->first) -= c1->second * x(i);  // xj -= Lji*xi
    if(i==0) break;
  }

  return x;
}


// Sparse eigensolver by shifted inverse iteration
inline Matrix eigenvec(const Sparse& a0, const SparseEigenConf& cf1, const SparseConf& cf20){
  Sparse amu=a0;
  SparseConf cf2=cf20;
  Matrix x, x_old;
  double mu=cf1.mu0, err=0, err_old;

  x=(cf1.x0.dim() ? cf1.x0 : cf2.x0.dim() ? cf2.x0 : Matrix(a0.nrow()).fill(1));
  x/=norm(x);

  for(int i=0; i<cf1.loop; i++){
    x_old=x;
    err_old=norm(a0*x-mu*x);

    for(size_t i=0; i<a0.nrow(); i++) amu(i,i)=a0(i,i)-mu;
    cf2.x0=x;
    if( cf1.pbcg ){ cf2.init_ilu0(amu); x=solve_pbcg(amu,x,cf2); }
    else x=solve_cg(amu,x,cf2);
    x/=norm(x);
    if( cf1.warmup<=i ) mu=inner(a0*x,x);  // use mu0 at an early stage to magnify eigens near mu0
    err=norm(a0*x-mu*x);

    if( cf1.warmup<i && err_old<cf1.tol && !(err<err_old) ) return x_old;
    if( cf1.verb ) std::cerr << "EIGEN["<<i<<"]\terr="<<err<<"\tmu="<<mu<<"\tcos="<<inner(x,x_old)<<"\n";
    if( !std::isfinite(err) ) break;
  }

  if( cf1.verb ) std::cerr<<"EIGEN: perhaps poor convergence.\n";
  return (std::isfinite(err) ? x : x_old);
}


}  //namespace nmlib
#endif //SPARSE_H
