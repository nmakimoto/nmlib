// Usage of Polynomial class
// - Include "polynomial.h" and use namespace "nmlib". Nothing to compile/link.
// - Features: arithmetic and other operations


#include <iostream>
#include "polynomial.h"
using namespace nmlib;


int main(void){
  int n=2;
  double c=3, x=4, y=5, eps=1.e-10;

  // Accessing coefficients
  Polynomial P,Q;  // zero polynomials (even without constant term)
  P.c(n) = c;   // P(x) = 3 x^2
  n = P.deg();  // polynomial degree = 2
  c = P.c(n);   // coefficient of x^n
  y = P(x);     // value at x

  // Arithmetic operations
  Q = (P*P-c*P+c*c)*(P+c);
  Q = (Q%P);

  // Other useful operations
  Q = power(P,n);    // P(x)^n
  Q = mirror(P);     // P(-x)
  Q = shift(P,c);    // P(x+c)
  Q = scale(P,c);    // P(cx)
  Q = diff(P);       // dP(x)/dx
  Q = integ(P);      // \int P(x)dx
  Q = compose(P,Q);  // P(Q(x))
  Q = shrink(Q,eps); // terms with small coefficient are omitted

  // Miscellaneous
  P = legendre<double>(n);   // Legendre polynoimal Pn
  x=1;
  for(int i=0; i<10; i++) x-=P(x)/diff(P)(x);  // zero of Pn by Newton
  std::cout << "Pn=" << P(x) << "=0 at " << x << '\n';
  
  return 0;
}
