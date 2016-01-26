// Polynomials


#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H


#include <vector>
#include <cmath>
#include <stdexcept>


/******** class I/F ********/

template<class T> class polynomial;
typedef polynomial<double> Polynomial;


template<class T> class polynomial{
public:
  explicit polynomial(const std::vector<T>& cc=std::vector<T>());  // p(x):=\sum_k cc[k] x^k
  void set(const std::vector<T>& cc);
  void get(      std::vector<T>& cc) const;

  int deg(void)  const;              // degree
  T&  c(int k);                      // p.c(k)=c
  T   c(int k) const;                // c=p.c(k)
  T   operator()(const T& x) const;  // y=p(x)

private:
  std::vector<T> cc;  // p(x):=\sum_k cc[k] x^k
};


/******** utility I/F ********/

// incremental operations
template<class T> polynomial<T>& operator+=(polynomial<T>& p, const T& s);  // p+=s
template<class T> polynomial<T>& operator-=(polynomial<T>& p, const T& s);  // p-=s
template<class T> polynomial<T>& operator*=(polynomial<T>& p, const T& s);  // p*=s
template<class T> polynomial<T>& operator/=(polynomial<T>& p, const T& s);  // p/=s
template<class T> polynomial<T>& operator+=(polynomial<T>& p, const polynomial<T>& q);  // p+=q
template<class T> polynomial<T>& operator-=(polynomial<T>& p, const polynomial<T>& q);  // p-=q
template<class T> polynomial<T>& operator*=(polynomial<T>& p, const polynomial<T>& q);  // p*=q
template<class T> polynomial<T>& operator%=(polynomial<T>& p, const polynomial<T>& q);  // (p/q)*q+(p%q)=p, deg(p%q)<deg(q)
template<class T> polynomial<T>& operator%=(polynomial<T>& p, const polynomial<T>& q);
template<class T> polynomial<T>& operator<<=(polynomial<T>& p, int n);  // p(x)*x^n
template<class T> polynomial<T>& operator>>=(polynomial<T>& p, int n);  // p(x)/x^n

// polynomial-scalar operations
template<class T> polynomial<T> operator-(const polynomial<T>& p);  // unary -p
template<class T> polynomial<T> operator+(const polynomial<T>& p, const T& s);  // p+s
template<class T> polynomial<T> operator-(const polynomial<T>& p, const T& s);  // p-s
template<class T> polynomial<T> operator*(const polynomial<T>& p, const T& s);  // p*s
template<class T> polynomial<T> operator/(const polynomial<T>& p, const T& s);  // p/s
template<class T> polynomial<T> operator+(const T& s, const polynomial<T>& p);  // s+p
template<class T> polynomial<T> operator-(const T& s, const polynomial<T>& p);  // s-p
template<class T> polynomial<T> operator*(const T& s, const polynomial<T>& p);  // s*p
//template<class T> polynomial<T> operator/(const T& s, const polynomial<T>& p);  // s/p

// polynomial-polynomial operations
template<class T> polynomial<T> operator+ (const polynomial<T>& p, const polynomial<T>& q);  // p+q
template<class T> polynomial<T> operator- (const polynomial<T>& p, const polynomial<T>& q);  // p-q
template<class T> polynomial<T> operator* (const polynomial<T>& p, const polynomial<T>& q);  // p*q
template<class T> polynomial<T> operator/ (const polynomial<T>& p, const polynomial<T>& q);  // p/q
template<class T> polynomial<T> operator% (const polynomial<T>& p, const polynomial<T>& q);  // p%q
template<class T> polynomial<T> operator<<(const polynomial<T>& p, int n);  // p*x^n
template<class T> polynomial<T> operator>>(const polynomial<T>& p, int n);  // p/x^n

// polynomial manipulations
template<class T> polynomial<T> shrink (const polynomial<T>& p, const T& eps=T(0), const T& x=T(1));  // omit terms s.t. |c_i|<|eps*x^i|
template<class T> polynomial<T> compose(const polynomial<T>& p, const polynomial<T>& q);  // p(q(x))
template<class T> polynomial<T> power  (const polynomial<T>& p, int n);       // p(x)^n
template<class T> polynomial<T> mirror (const polynomial<T>& p);              // p(-x)
template<class T> polynomial<T> shift  (const polynomial<T>& p, const T& c);  // p(x+c)
template<class T> polynomial<T> scale  (const polynomial<T>& p, const T& c);  // p(cx)
template<class T> polynomial<T> diff   (const polynomial<T>& p);              // dp(x)/dx
template<class T> polynomial<T> integ  (const polynomial<T>& p);              // \int p(x)dx
template<class T> polynomial<T> legendre(int n);  // Legendre polynomial Pn(x)


/******** implementation ********/

template<class T>      polynomial<T>::polynomial(const std::vector<T>& cc0){ set(cc0); }
template<class T> void polynomial<T>::set(const std::vector<T>& cc0)       { cc=cc0; }
template<class T> void polynomial<T>::get(      std::vector<T>& cc1) const { cc1=cc; }
template<class T> int  polynomial<T>::deg(void ) const{ return cc.size()-1; }  // note: deg=-1 if cc is empty
template<class T> T&   polynomial<T>::c  (int k)      { if(k<0) throw std::domain_error("polynomial::c(k): k<0"); if(deg()<k) cc.resize(k+1,0); return cc[k]; }
template<class T> T    polynomial<T>::c  (int k) const{ return ((0<=k and k<=deg()) ? cc[k] : 0); }
template<class T> T    polynomial<T>::operator()(const T& x) const{ T y=0; for(int i=deg(); i>=0; i--) y=y*x+cc[i]; return y; }

template<class T> polynomial<T>& operator+=(polynomial<T>& p, const T& s){ p.c(0)+=s; return p; }
template<class T> polynomial<T>& operator-=(polynomial<T>& p, const T& s){ p.c(0)-=s; return p; }
template<class T> polynomial<T>& operator*=(polynomial<T>& p, const T& s){ for(int i=p.deg(); i>=0; i--) p.c(i)*=s; return p; }
template<class T> polynomial<T>& operator/=(polynomial<T>& p, const T& s){ for(int i=p.deg(); i>=0; i--) p.c(i)/=s; return p; }
template<class T> polynomial<T>& operator+=(polynomial<T>& p, const polynomial<T>& q){ for(int j=q.deg(); j>=0; j--) p.c(j)+=q.c(j); return p; }
template<class T> polynomial<T>& operator-=(polynomial<T>& p, const polynomial<T>& q){ for(int j=q.deg(); j>=0; j--) p.c(j)-=q.c(j); return p; }
template<class T> polynomial<T>& operator*=(polynomial<T>& p, const polynomial<T>& q){
  polynomial<T> r;
  for(int i=p.deg(); i>=0; i--)
    for(int j=q.deg(); j>=0; j--)
      r.c(i+j)+=p.c(i)*q.c(j);
  p=r;
  return p;
}
template<class T> polynomial<T>& operator/=(polynomial<T>& p, const polynomial<T>& q){
  if(q.deg()<0) throw std::domain_error("polynomial p/q: q=0");
  polynomial<T> r;
  for(int k=p.deg()-q.deg(); k>=0; k--){
    T t=p.c(k+q.deg())/q.c(q.deg());
    r.c(k)=t;
    for(int j=q.deg(); j>=0; j--) p.c(j+k)-=q.c(j)*t;
  }
  p=r;
  return p;
}
template<class T> polynomial<T>& operator%=(polynomial<T>& p, const polynomial<T>& q){
  if(q.deg()<0) throw std::domain_error("polynomial p%q: q=0");
  if(q.deg()>p.deg()) return p;
  polynomial<T> r=p-(p/q)*q;
  int n=p.deg()-q.deg()-1;
  p=polynomial<T>();
  for(int k=n; k>=0; k--) p.c(k)=r.c(k);
  return p;
}
template<class T> polynomial<T>& operator<<=(polynomial<T>& p, int n){
  if(n<0) throw std::domain_error("polynomial p<<n: n<0");
  polynomial<T> q;
  for(int i=p.deg(); i>=0; i--) q.c(i+n)=p.c(i);
  p=q;
  return p;
}
template<class T> polynomial<T>& operator>>=(polynomial<T>& p, int n){
  if(n<0) throw std::domain_error("polynomial p>>n: n<0");
  polynomial<T> q;
  for(int i=p.deg(); i>=n; i--) q.c(i-n)=p.c(i);
  p=q;
  return p;
}

template<class T> polynomial<T> operator-(const polynomial<T>& p)     { polynomial<T> q(p); q*=-T(1); return q; }
template<class T> polynomial<T> operator+(const polynomial<T>& p, const T& s){ polynomial<T> q(p); q+=s; return q; }
template<class T> polynomial<T> operator-(const polynomial<T>& p, const T& s){ polynomial<T> q(p); q-=s; return q; }
template<class T> polynomial<T> operator*(const polynomial<T>& p, const T& s){ polynomial<T> q(p); q*=s; return q; }
template<class T> polynomial<T> operator/(const polynomial<T>& p, const T& s){ polynomial<T> q(p); q/=s; return q; }
template<class T> polynomial<T> operator+(const T& s, const polynomial<T>& p){ polynomial<T> q(p); q+=s; return q; }
template<class T> polynomial<T> operator-(const T& s, const polynomial<T>& p){ polynomial<T> q(p); q*=-T(1); q+=s; return q; }
template<class T> polynomial<T> operator*(const T& s, const polynomial<T>& p){ polynomial<T> q(p); q*=s; return q; }
//template<class T> polynomial<T> operator/(T s, const polynomial<T>& p);

template<class T> polynomial<T> operator+(const polynomial<T>& p, const polynomial<T>& q){ polynomial<T> r(p); r+=q; return r; }
template<class T> polynomial<T> operator-(const polynomial<T>& p, const polynomial<T>& q){ polynomial<T> r(p); r-=q; return r; }
template<class T> polynomial<T> operator*(const polynomial<T>& p, const polynomial<T>& q){ polynomial<T> r(p); r*=q; return r; }
template<class T> polynomial<T> operator/(const polynomial<T>& p, const polynomial<T>& q){ polynomial<T> r(p); r/=q; return r; }
template<class T> polynomial<T> operator%(const polynomial<T>& p, const polynomial<T>& q){ polynomial<T> r(p); r%=q; return r; }
template<class T> polynomial<T> operator<<(const polynomial<T>& p, int n){ polynomial<T> q(p); q<<=n; return q; }
template<class T> polynomial<T> operator>>(const polynomial<T>& p, int n){ polynomial<T> q(p); q>>=n; return q; }

template<class T> polynomial<T> shrink (const polynomial<T>& p, const T& eps, const T& x){
  polynomial<T> q;
  T xn=T(1);
  for(int i=0; i<p.deg(); i++,xn*=x)
    if(std::abs(eps)<std::abs(p.c(i)*xn)) q.c(i)=p.c(i);
  return q;
}
template<class T> polynomial<T> compose(const polynomial<T>& p, const polynomial<T>& q){
  polynomial<T> r;
  for(int i=p.deg(); i>=0; i--) r=r*q+p.c(i);
  return r;
}
template<class T> polynomial<T> power  (const polynomial<T>& p0, int n){
  if(n<0) throw std::domain_error("polynomial power(p,n): n<0");
  polynomial<T> q,p;
  for(p=p0, q.c(0)=1; n>0; n/=2,p*=p) if(n%2) q*=p;
  return q;
}
template<class T> polynomial<T> shift  (const polynomial<T>& p, const T& c){ polynomial<T> q; q.c(1)=1; q.c(0)=c; return compose(p,q); }
template<class T> polynomial<T> scale  (const polynomial<T>& p, const T& c){ polynomial<T> q; q.c(1)=c;           return compose(p,q); }
template<class T> polynomial<T> mirror (const polynomial<T>& p)            { polynomial<T> q; q.c(1)=-1;          return compose(p,q); }
template<class T> polynomial<T> diff   (const polynomial<T>& p)            { polynomial<T> q; for(int i=p.deg(); i>=1; i--) q.c(i-1)=p.c(i)*i;     return q; }
template<class T> polynomial<T> integ  (const polynomial<T>& p)            { polynomial<T> q; for(int i=p.deg(); i>=0; i--) q.c(i+1)=p.c(i)/(i+1); return q; }
template<class T> polynomial<T> legendre(int n){
  if(n<0) throw std::domain_error("polynomial legendre(n): n<0");
  polynomial<T> p0,p1,p2;
  p0.c(0)=T(1);  if(n==0) return p0;
  p1.c(1)=T(1);  if(n==1) return p1;
  for(int k=2; k<=n; k++){
    p2=(T(2*k-1)*(p1<<1)-T(k-1)*p0)/T(k);  // Bonnet recursion
    p0=p1;
    p1=p2;
  }
  return p2;
}


#endif //POLYNOMIAL_H
