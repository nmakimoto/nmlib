// Usage of I/O utilities
// - Include "io.h" and use namespace "nmlib". No linking is required.
// - Features: string conversion to/from generic T, stream I/O of vector<T> and map<K,T>


#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "io.h"
using namespace nmlib;


// I/O test class - anything will do as long as it has >> and <<
struct Dummy{ float a,b; };
std::istream& operator>>(std::istream& str,       Dummy& x){ str>>x.a     >>x.b;      return str; }
std::ostream& operator<<(std::ostream& str, const Dummy& x){ str<<x.a<<' '<<x.b<<' '; return str; }


int main(void){
  Dummy x1;
  std::vector<Dummy> xx;
  std::map<int,Dummy> k2x;
  std::string s;
  int prec=4;  // number of digits in scientific notation, used in any2str()

  // class T I/O
  //std::cin >> x1;
  str2any("1.1 1.2",x1);
  s=any2str(x1,prec);
  std::cout << x1 << '\n'+s+'\n';

  // vector<T> I/O (format: N x1 ... xN)
  //std::cin >> xx;
  str2any("3  2.1 2.2 3.1 3.2 4.1 4.2", xx);
  s=any2str(xx,prec);
  std::cout << xx << '\n'+s+'\n';

  // map<K,T> I/O (format: N k1 x1 ... kN xN)
  //std::cin >> k2x;
  str2any("3 5 1.1 1.2 6 2.1 2.2 7 3.1 3.2", k2x);
  s=any2str(k2x,prec);
  std::cout << k2x<< '\n'+s+'\n';

  return 0;
}
