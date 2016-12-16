// Usage of Sparse matrix library
// - Include sparse.H and matrix.H into your code. Use namespace nmlib.
// - Compile and link your code with sparse.C.
// - See the code below for more on usage.


#include <iostream>
#include "sparse.H"
#include "matrix.H"
using namespace nmlib;


int main(void){
  size_t n=10000;
  Sparse a(n,n);
  Matrix b(n),x;

  // example: 1d heat equation (A=laplacian, b=external source)
  a+=2.0;  // diagonal
  for(size_t k=1; k<n; k++) a(k-1,k)=a(k,k-1)=-1;  // subdiagonal
  for(size_t k=0; k<n; k++) b(k)=1/(1.0+k);

  // solve Ax=b either by...
  x = solve_cg  (a,b);  // CG (Conjugate Gradient) - only for symmetric positive definite A
  x = solve_bcg (a,b);  // BiCG (Bi-conjugate Gradient)
  x = solve_pbcg(a,b);  // PBCG (Preconditioned BiCG) - band matrices are solved without iteration

  // if convergense is not satisfactory, then try any of SparseConf options below.
  SparseConf cf;
  //cf.tol  = 1.e-8;  // tolerance |b-Ax|/|b|
  //cf.loop = 1000;   // max number of iteration
  //cf.verb = true;   // show progress on std::cerr
  //cf.x0   = ...;    // initial approximate solution
  //cf.lu   = ...;    // preconditioner for PBCG - set cf.lu=(L-1)+U for given ILU of A
  x = solve_cg  (a,b,cf);
  x = solve_bcg (a,b,cf);
  x = solve_pbcg(a,b,cf);

  std::cout << x << std::endl;
  return 0;
}
