// Usage of matrix class library
// - Just include "matrix.h" to use nmlib::matrix<T>. No linking is required.
// - Features: arithmetic operations, symmetric eigensolvers, stream I/O, etc.
// - See the code below for more on usage. Enjoy!


#include <iostream>
#include "matrix.h"
#include "io.h"  // str2any()
using namespace nmlib;


/******** Simple Example ********/

int main(void){
  Matrix A(2,2),b(2),x;
  A(0,0)=2; A(0,1)=0; b(0)=2;  // A=|2 0| b=|2|
  A(1,0)=0; A(1,1)=1; b(1)=1;  //   |0 1|   |1|
  x=inv(A)*b;  // solve Ax=b
  std::cout << x << std::endl;
}


/******** Cheat Sheet ********/

void cheatsheet(const Matrix& A, const Matrix& B){
  Matrix M,U,V;
  size_t r=3, c=3, i=0, j=0, k=0, n=-3;
  double x=1.1, y=2.2, z=3.3, t;
  std::vector<double> vec={1,2,3};

  // Initialization & I/O
  M=A;
  M=Matrix(r,c);        // zero matrix of size r by c
  V=Matrix(r);          // zero matrix of size r by 1 (column vector)
  V=Matrix(vec);        // column vector
  V={x,y,z};            // column vector
  M={{x,y,z},{x,y,z}};  // matrix consisting of multiple column vectors
  M=str2any<Matrix>("3  2  0.0  0.1  1.0  1.1  2.0  2.1");  // (see "io.h")
  std::string str=any2str(M);
  std::cin  >> M;
  std::cout << M;

  // Getting/Setting Components
  t=M(i,j);
  t=V(k);
  M(i,j)=t;
  V(k)=t;

  // Arithmetic Operations
  M=-A;
  M+=t;   M-=t;   M*=t;   M/=t;   // (note: M+t:=M+tI)
  M=A+t;  M=A-t;  M=A*t;  M=A/t;
  M=t+A;  M=t-A;  M=t*A;  M=t/A;
  M=A+B;  M=A-B;  M=A*B;

  // Other Useful Operations
  M=inv(A);    // A^-1 (inverse of A)
  M=pow(A,n);  // A^n  (n-th power of A)
  M=tp(A);     // A^T  (transpose of A)
  U=orth(A);   // U^T A=upper triangular  (Gram-Schmidt orthonormalization of A)
  U=eigen(A);  // U^T A U=diagonal  (an orthonormal frame of eigenvectors of A)
  M=hcat(A,B);    // horizontal concatenation  (r x cA+cB matrix)
  M=vcat(A,B);    // vertical concatenation    (rA+rB x c matrix)
  M.resize(3,4).fill(5.0);
  V=getvec(M,j);  // get j-th column vector of M
  setvec(M,j,V);  // set j-th column vector of M to V
  //...

  t=norm(A);     // |A|:=<A,A>^(-1/2)
  t=inner(A,B);  // <A,B>:=tr(A^T B)
  t=det(A);

  // 3D-Specific Operations
  Matrix v1={1,0,0}, v2={0,1,0}, v3;
  v3=outer(v1,v2);  // (0,0,1)^T ("cross" product)
  M=vec2rot(v3*M_PI);  // rotation around v3 by 180 degrees
  //...
}
