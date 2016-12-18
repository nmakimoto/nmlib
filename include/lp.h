// LP solver class library


#ifndef LP_H
#define LP_H


#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "matrix.h"
namespace nmlib {


/******** Class I/F ********/

class LP{
public:
  LP(void);
  LP(const Matrix& a, const Matrix& b, const Matrix& c);  // max cx s.t. Ax<=b and x>=0

  void init(const Matrix& a, const Matrix& b, const Matrix& c);
  int  solve      (int max_iter=-1, bool bland_rule=true);
  int  solve_stdlp(int max_iter=-1, bool bland_rule=true, bool aux_prob=false);
  int  solve_auxlp(int max_iter=-1, bool bland_rule=true);
  int  choose_axis(int& i0, int& j0, bool aux_prob, bool bland_rule);
  void pivot(int i0, int j0);

  Matrix vertex(void) const;  // basic solution
  int    status(void) const;

  Matrix tbl;
  matrix<int> idx;
};


/******** Utility I/F ********/

// stream I/O
std::istream& operator>>(std::istream& str,       LP& lp);
std::ostream& operator<<(std::ostream& str, const LP& lp);
std::ostream& read (std::ostream& str,       LP& lp);
std::ostream& write(std::ostream& str, const LP& lp);


/******** Implementation ********/

inline LP::LP(void){ init(Matrix(),Matrix(),Matrix()); }
inline LP::LP(const Matrix& a, const Matrix& b, const Matrix& c){ init(a,b,c); }


// initialize simplex tableau
//   |A   -1 b|  ...  constraints
//   |c^T  0 0|  ...  objectives (original)
//   |0   -1 0|  ...  objectives (auxiliary)
inline void LP::init(const Matrix& aa, const Matrix& bb, const Matrix& cc){
  int r=aa.nrow(), c=aa.ncol();
  if( !(aa.ncol()==cc.dim() && bb.dim()==aa.nrow()) )
    throw std::domain_error("LP::init(): bad dimensions");

  // conefficients of original problem
  tbl=Matrix(r+2,c+2);
  for(int i=0; i<r; i++)
    for(int j=0; j<c; j++) tbl(i,j)=aa(i,j);
  for(int i=0; i<r; i++) tbl(i,c+1)=bb(i);
  for(int j=0; j<c; j++) tbl(r,j)=cc(j);

  // coefficients of auxiliary problem
  for(int i=0; i<r; i++) tbl(i,c)=-1;
  tbl(r+1,c)=-1;

  // variable indices
  idx=matrix<int>(r+c+1);
  for(int j=0; j<c; j++) idx(j)=j;  // nonbasic vars
  idx(c)=-1;  // aux var
  for(int i=0; i<r; i++) idx(c+1+i)=c+i;  // basic vars
}


// solve LP: max cx s.t. Ax<=b, x>=0
// status: -1(infeasible), 0(feasible), +1(optimal), +2(unbounded)
inline int LP::solve(int iter, bool bland){
  if(solve_auxlp(iter,true)<0) return -1;  // infeasible
  return  solve_stdlp(iter,bland,false);  // feasible(0), optimal(1) or unbounded(2)
}


// solve standard LP: max cx s.t. Ax+s=b, x>=0, s>=0  where b>=0
inline int LP::solve_stdlp(int iter, bool bland, bool aux){
  int r=tbl.nrow()-2, c=tbl.ncol()-2;

  for(int i=0; i<r; i++)
    if(tbl(i,c+1)<0) return -1;  // infeasible (not b>=0)
  for(int i=0; i<r; i++)
    if(!aux && idx(c+1+i)<0) return -1;  // infeasible (aux var is basic)
  while(iter--){
    int i0,j0;
    choose_axis(i0,j0,aux,bland);
    if(j0<0) return 1;  // optimal
    if(i0<0) return 2;  // unbounded
    pivot(i0,j0);
  }
  return 0;  // just feasible
}


// solve auxiliary LP: max -t s.t. Ax-t+s=b and x>=0, s>=0, t>=0
inline int LP::solve_auxlp(int iter, bool bland){
  int r=tbl.nrow()-2, c=tbl.ncol()-2;
  if(idx(c)>=0) throw std::runtime_error("LP::solve_auxlp(): aux variable not found at c-th column");

  int i0=0;  // i0 := argmin_i b(i)
  for(int i=0; i<r; i++)
    if(tbl(i,c+1)<tbl(i0,c+1)) i0=i;
  if(tbl(i0,c+1)<0){
    pivot(i0,c);
    solve_stdlp(iter,bland,true);  // bland=true is recommended to make sure aux var of feasible LP becomes nonbasic
  }

  int j0=-1;  // j0 := column of aux var
  for(int j=0; j<c+1; j++)
    if(idx(j)<0) j0=j;
  if(j0<0) return -1;  // aux var is basic, i.e., infeasible (unless degenerate)

  for(int i=0; i<r+2; i++) tbl(i,j0)=0;  // forget about auxiliary problem
  return 0;
}


// choose axis for pivoting
inline int LP::choose_axis(int& i0, int &j0, bool aux, bool bland){
  int r=tbl.nrow()-2, c=tbl.ncol()-2, r1=r, c1=c+1;
  if(aux) r1++;

  j0=i0=-1;

  for(int j=0; j<c1; j++){
    if(tbl(r1,j)<=0) continue;
    else if(j0<0) j0=j;
    else if(bland && idx(j)<idx(j0)) j0=j;  // j0:=argmin_j idx(j) s.t. c(j)>0
    else if(tbl(r1,j0)<tbl(r1,j)) j0=j;  // j0:=argmax_i c(j)
  };
  if(j0<0) return 1;  // optimal

  for(int i=0; i<r; i++){
    if(tbl(i,j0)<=0) continue;
    else if(i0<0) i0=i;
    else if(bland && idx(c1+i)<idx(c1+i0) && tbl(i,c1)==0) i0=i;  // j0:=argmin_i idx(i) s.t. b(i)=0
    else if(tbl(i,c1)/tbl(i,j0)<tbl(i0,c1)/tbl(i0,j0)) i0=i;  // i0:=argmin_i b(i)/A(i,j0) where A(i,j0)>0
  }
  if(i0<0) return 2;  // unbounded

  return 0;
}


// pivote tableau on A(i0,j0)
inline void LP::pivot(int i0, int j0){
  int r=tbl.nrow()-2, c=tbl.ncol()-2;
  double t;
  
  for(int i=0; i<r+2; i++){
    if(i==i0) continue;
    t=tbl(i,j0)/tbl(i0,j0);
    for(int j=0; j<c+2; j++) tbl(i,j)-=t*tbl(i0,j);  // sweep (i,j0)
    tbl(i,j0)=-t;  // swap variables
  }

  t=1/tbl(i0,j0);
  for(int j=0; j<c+2; j++) tbl(i0,j)*=t;  // normalize i0-th row
  tbl(i0,j0)=t;  // swap variables

  std::swap(idx(j0),idx(c+1+i0));
}


// retrieve basic solution from current tableau
inline Matrix LP::vertex(void) const{
  int r=tbl.nrow()-2, c=tbl.ncol()-2;
  Matrix x(c);
  for(int i=0; i<r; i++){
    int j=idx(c+1+i);
    if(0<=j && j<c) x(j)=tbl(i,c+1);
  }
  return x;
}


// status of tableau: -1(not feasible), 0(feasible), +1(optimal), +2(unbounded)
inline int LP::status(void) const {
  int r=tbl.nrow()-2, c=tbl.ncol()-2, i, j;

  for(i=0; i<r; i++)
    if(tbl(i,c+1)<0) return -1;  // currently not feasible (not b>=0)

  for(i=0; i<r; i++)
    if(idx(c+1+i)<0) return -1;  // infeasible (aux var is basic)

  for(j=0; j<c+1; j++)
    if(0<tbl(r,j)) break;
  if(j==c+1) return 1;  // optimal (c<=0)

  for(j=0; j<c+1; j++){
    if(idx(j)<0 || tbl(r,j)<=0) continue;
    for(i=0; i<r; i++)
      if(0<tbl(i,j)) break;
    if(i==r) return 2;  // unbounded (c_j>0 and a_*j<=0)
  }

  return 0;  // currently just feasible, no further knowledge
}


// stream I/O
inline std::istream& operator>>(std::istream& str, LP& lp){
  str >> lp.tbl >> lp.idx;
  return str;
}
inline std::ostream& operator<<(std::ostream& str, const LP& lp){
  str << lp.tbl << std::endl << lp.idx << std::endl;
  return str;
}
inline std::istream& read(std::istream& str, LP& lp){
  read(str,lp.tbl);
  read(str,lp.idx);
  return str;
}
inline std::ostream& write(std::ostream& str, const LP& lp){
  write(str,lp.tbl);
  write(str,lp.idx);
  return str;
}


}  //namespace nmlib
#endif //LP_H
