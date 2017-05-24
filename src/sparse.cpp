// Sparse matrix and solvers (CG/BiCG/PBCG)


#include <iostream>
#include <map>
#include <stdexcept>
#include <sstream>
#include "sparse.h"
#include "matrix.h"
//#include "io.h"
namespace nmlib{


/******************************** Range checks ********************************/

// Range check for each operation
void chk(size_t r1, size_t c1, size_t r2, size_t c2, char op){
  if( (op=='@' && r1<r2  && c1<=c2)  ||
      (op=='+' && r1==r2 && c1==c2)  ||
      (op=='-' && r1==r2 && c1==c2)  ||
      (op=='*' && c1==r2)  ||
      (op=='/' && r1==c1 && c1==r2 && c2==1) ) return;
  std::stringstream msg;
  msg<<"Sparse: ("<<r1<<","<<c1<<") "<< op << " ("<<r2<<","<<c2<<")";
  throw std::domain_error(msg.str());
}
void chk(const Sparse& a, const Sparse& b, char op){ chk(a.nrow(),a.ncol(),b.nrow(),b.ncol(),op); }
void chk(const Sparse& a, const Matrix& b, char op){ chk(a.nrow(),a.ncol(),b.nrow(),b.ncol(),op); }


/*********************************** Basic operations ***********************************/

// Class methods
Sparse::Sparse(size_t r, size_t c): row(r),col(c) { val=std::vector<std::map<size_t,double>>(r); }
size_t  Sparse::nrow(void) const{ return row; }
size_t  Sparse::ncol(void) const{ return col; }
double& Sparse::operator()(size_t i, size_t j)       { chk(i,j,row,col,'@'); return val[i][j]; }
double  Sparse::operator()(size_t i, size_t j) const { chk(i,j,row,col,'@'); auto c=val[i].find(j); return (c==val[i].end() ? 0 : c->second); }

void SparseConf::init_ilu0(const Sparse& a){ lu=a; ilu0(lu); }


// Incremental operations
Sparse& operator+=(Sparse& a, double s){ for(size_t k=0; k<a.ncol(); k++) a(k,k)+=s; return a; }
Sparse& operator-=(Sparse& a, double s){ a+=-s; return a; }
Sparse& operator*=(Sparse& a, double s){ for(size_t i=0; i<a.nrow(); i++) for(auto& c: a.columns(i)) c.second*=s; return a; }
Sparse& operator/=(Sparse& a, double s){ a*=1/s; return a; }
Sparse& operator+=(Sparse& a, const Sparse& b){ chk(a,b,'+'); for(size_t i=0; i<b.nrow(); i++) for(const auto& c: b.columns(i)) a(i,c.first)+=c.second; return a; }
Sparse& operator-=(Sparse& a, const Sparse& b){ chk(a,b,'-'); for(size_t i=0; i<b.nrow(); i++) for(const auto& c: b.columns(i)) a(i,c.first)-=c.second; return a; }


// Basic operations
Sparse operator-(const Sparse& a)          { Sparse b(a); b*=-1.0; return b; }
Sparse operator+(const Sparse& a, double s){ Sparse b(a); b+=s; return b; }
Sparse operator-(const Sparse& a, double s){ Sparse b(a); b-=s; return b; }
Sparse operator*(const Sparse& a, double s){ Sparse b(a); b*=s; return b; }
Sparse operator/(const Sparse& a, double s){ Sparse b(a); b/=s; return b; }
Sparse operator+(double s, const Sparse& a){ Sparse b(a); b+=s; return b; }
Sparse operator-(double s, const Sparse& a){ Sparse b(a); b*=-1.0; b+=s; return b; }
Sparse operator*(double s, const Sparse& a){ Sparse b(a); return (b*=s); }
//Sparse operator/(double s, const Sparse& a);
Sparse operator+(const Sparse& a, const Sparse& b){ Sparse c(a); c+=b; return c; }
Sparse operator-(const Sparse& a, const Sparse& b){ Sparse c(a); c-=b; return c; }
Sparse operator*(const Sparse& a, const Sparse& b){
  chk(a,b,'*');
  Sparse ret(a.nrow(),b.ncol());
  for(size_t j=0; j<ret.ncol(); j++)
    for(size_t i=0; i<a.nrow(); i++)
      for(const auto&  c: a.columns(i))  ret(i,j)+=c.second*b(c.first,j);
  return ret;
}
Matrix operator*(const Sparse& a, const Matrix& b){
  chk(a,b,'*');
  Matrix ret(a.nrow(),b.ncol());
  for(size_t j=0; j<ret.ncol(); j++)
    for(size_t i=0; i<a.nrow(); i++)
      for(const auto&  c: a.columns(i))  ret(i,j)+=c.second*b(c.first,j);
  return ret;
}


// Utilities
double norm (const Sparse& a){ double s=0; for(size_t i=0; i<a.nrow(); i++) for(const auto& c: a.columns(i)){ double t=c.second; s+=t*t; } return sqrt(s); }
Sparse tp   (const Sparse& a){ Sparse b(a.ncol(),a.nrow()); for(size_t i=0; i<a.nrow(); i++) for(const auto& c: a.columns(i)) b(c.first,i)=c.second; return b; }
Matrix dense(const Sparse& a){ Matrix b(a.nrow(),a.ncol()); for(size_t i=0; i<a.nrow(); i++) for(const auto& c: a.columns(i)) b(i,c.first)=c.second; return b; }
Matrix tpab (const Sparse& a, const Matrix& b){
  chk(a.ncol(),a.nrow(),b.nrow(),b.ncol(),'*');
  Matrix ret(b.nrow(),b.ncol());
  for(size_t j=0; j<ret.ncol(); j++)
    for(size_t i=0; i<a.nrow(); i++)
      for(const auto& c: a.columns(i))  ret(c.first,j)+=c.second*b(i,j);
  return ret;
}


/******************************** Sparse solvers ********************************/

// Solve Ax=b (CG - Conjugate Gradient for symmetric positive definite A)
Matrix solve_cg(const Sparse& a, const Matrix& b, const SparseConf& cf){
  chk(a,b,'/');

  Matrix x,r,p,q;
  double t,rr;

  x=Matrix(a.ncol());  if(cf.x0.dim()>0) x=cf.x0;
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
Matrix solve_bcg(const Sparse& a, const Matrix& b, const SparseConf& cf){
  chk(a,b,'/');

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
Matrix solve_pbcg(const Sparse& a, const Matrix& b, const SparseConf& cf){
  chk(a,b,'/');

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


/******************************* Pre-conditioning *******************************/

// b=LUx
Matrix lux(const Sparse& lu, const Matrix&  x){
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
Matrix sluxb(const Sparse& lu, const Matrix& b){
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
Matrix slutxb(const Sparse& lu, const Matrix& b){
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


// ILU(0) - incomplete LU decomposition
// - A typical preconditioner without fill-in
// - Decomposition is complete when A is a band matrix
void ilu0(Sparse& a){
  for(size_t i=0; i<a.nrow(); i++){
    std::map<size_t,double>& ai=a.columns(i);
    auto aii=ai.find(i);  // Uii
    if( aii==ai.end() ) throw std::runtime_error("ILU(0): division by Aii=0");

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


}  //namespace nmlib
