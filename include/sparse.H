// Sparse matrix and solvers (CG/BiCG/PBCG)


#ifndef SPARSE_H
#define SPARSE_H


#include <iostream>
#include <vector>
#include <map>
#include "matrix.H"
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


}  //namespace nmlib
#endif //SPARSE_H
