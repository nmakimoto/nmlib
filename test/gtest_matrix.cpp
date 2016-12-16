// Unit test (matrix)


#include <gtest/gtest.h>
#include <stdexcept>

#include <cstdlib>  // random
#include "matrix.h"
using namespace nmlib;


// Sample data generators
double urand1(void){
  return (random()%RAND_MAX+0.5)/RAND_MAX;
}
Matrix mrand(size_t r, size_t c){
  Matrix m(r,c);
  for(size_t k=0; k<m.dim(); k++) m(k)=urand1()*2-1;
  return m;
}
Matrix exp(const Matrix& m){
  Matrix s=m+1.0;
  Matrix t=m;
  for(int i=2; i<100; i++){
    t=t*m/double(i);
    s+=t;
  }
  return s;
}


// Initialization
TEST(matrix,init){
  Matrix m1,m2;

  // constructors...
  m1=Matrix();
  EXPECT_TRUE(m1.nrow()==0 and m1.ncol()==0 and m1.dim()==0);

  m1=Matrix(3);
  EXPECT_TRUE(m1.nrow()==3 and m1.ncol()==1 and m1.dim()==3);
  for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m1(k),0);

  m1=Matrix(3,4);
  EXPECT_TRUE(m1.nrow()==3 and m1.ncol()==4 and m1.dim()==12);
  for(size_t i=0; i<m1.nrow(); i++)  for(size_t j=0; j<m1.ncol(); j++)  EXPECT_DOUBLE_EQ(m1(i,j), 0);

  m1=Matrix(1.1, 2.2, 3.3);
  EXPECT_TRUE(m1.nrow()==3 and m1.ncol()==1 and m1.dim()==3);
  for(size_t k=0; k<m1.dim(); k++) EXPECT_NEAR(m1(k),1.1*(k+1),1.e-12);

  double aa[]={1.11,2.22,3.33,4.44,5.55};
  m1=Matrix(aa,aa+5);
  EXPECT_TRUE(m1.nrow()==5 and m1.ncol()==1 and m1.dim()==5);
  for(size_t i=0; i<m1.dim(); i++) EXPECT_NEAR(m1(i), (i+1)*1.11d, 1.e-8);

  std::vector<double> vv(aa,aa+5);
  m2=Matrix(vv);
  EXPECT_TRUE(m1.nrow()==5 and m1.ncol()==1 and m1.dim()==5);
  for(size_t i=0; i<m2.dim(); i++) EXPECT_NEAR(m2(i), (i+1)*1.11d, 1.e-8);

  // type converters...
  vv=std::vector<double>(m2);
  EXPECT_TRUE(vv.size()==5);
  for(size_t i=0; i<vv.size(); i++) EXPECT_NEAR(vv[i], (i+1)*1.11d, 1.e-8);

  matrix<std::complex<double> > mc = matrix<std::complex<double> >(m1);
  EXPECT_TRUE(mc.nrow()==5 and mc.ncol()==1 and mc.dim()==5);
  for(size_t i=0; i<vv.size(); i++) EXPECT_NEAR(mc(i).real(), (i+1)*1.11d, 1.e-8);
  for(size_t i=0; i<vv.size(); i++) EXPECT_DOUBLE_EQ(mc(i).imag(), 0);
}


// Component access
TEST(matrix,component){
  Matrix m1(3,4);
  for(size_t i=0; i<m1.nrow(); i++)  for(size_t j=0; j<m1.ncol(); j++)  m1(i,j)=i*10+j+0.1;
  for(size_t i=0; i<m1.nrow(); i++)  for(size_t j=0; j<m1.ncol(); j++)  EXPECT_DOUBLE_EQ(m1(i,j), i*10+j+0.1);
  for(size_t k=0; k<m1.dim(); k++)  m1(k)=k+0.1;
  for(size_t k=0; k<m1.dim(); k++)  EXPECT_DOUBLE_EQ(m1(k), k+0.1);
}


// Norm and inner product
TEST(matrix,norm){
  Matrix m1,m2;

  m1=Matrix(3,2);
  EXPECT_DOUBLE_EQ(norm(m1), 0);

  m1=Matrix(3,2);  for(size_t k=0; k<m1.dim(); k++) m1(k)=k+1;
  m2=Matrix(3,2);  for(size_t k=0; k<m2.dim(); k++) m2(k)=pow(10,k);
  EXPECT_DOUBLE_EQ(inner(m1,m2), 654321);
  EXPECT_DOUBLE_EQ(inner(m1,m1), 6*7*13/6);  // n(n+1)(2n+1)/6
  EXPECT_DOUBLE_EQ(norm(m1),sqrt(inner(m1,m1)));
}


// Arithmetic operations
TEST(matrix,arithmetic){
  Matrix m1,m2,m3;
  double t;

  t =-2;
  m1=Matrix(3,4);  for(size_t k=0; k<12; k++) m1(k)=(k+1.0)* 1;
  m2=Matrix(3,4);  for(size_t k=0; k<12; k++) m2(k)=(k+1.0)*10;

  // scalar multiple
  m3=m1; m3*=t;   for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)*2);
  m3=m1; m3/=t;   for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)/2);
  m3=-m1;         for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0));
  m3=m1/t;        for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)/2);
  m3=m1*t;        for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)*2);
  m3=t*m1;        for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)*2);

  // matrix +/- matrix
  m3=m1; m3+=m2;  for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k), (k+1.0)*11);
  m3=m1; m3-=m2;  for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)* 9);
  m3=m1+m2;       for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k), (k+1.0)*11);
  m3=m1-m2;       for(size_t k=0; k<m1.dim(); k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)* 9);

  // matrix * matrix
  m1=mrand(4,3);
  m2=mrand(3,5);
  m3=m1*m2;
  for(size_t i=0; i<m3.nrow(); i++)
    for(size_t j=0; j<m3.ncol(); j++)
      EXPECT_DOUBLE_EQ(m3(i,j), m1(i,0)*m2(0,j) + m1(i,1)*m2(1,j) + m1(i,2)*m2(2,j));

  t=-2;
  m1=Matrix(3,3);  for(size_t k=0; k<9; k++) m1(k)=(k+1.0)*1;

  // matrix +/- scalar matrix
  m3=m1; m3+=t;  for(size_t k=0; k<9; k++) EXPECT_DOUBLE_EQ(m3(k), (k+1.0)+(k/3==k%3 ? -2 : 0));
  m3=m1; m3-=t;  for(size_t k=0; k<9; k++) EXPECT_DOUBLE_EQ(m3(k), (k+1.0)+(k/3==k%3 ?  2 : 0));
  m3=m1+t;       for(size_t k=0; k<9; k++) EXPECT_DOUBLE_EQ(m3(k), (k+1.0)+(k/3==k%3 ? -2 : 0));
  m3=m1-t;       for(size_t k=0; k<9; k++) EXPECT_DOUBLE_EQ(m3(k), (k+1.0)+(k/3==k%3 ?  2 : 0));
  m3=t+m1;       for(size_t k=0; k<9; k++) EXPECT_DOUBLE_EQ(m3(k), (k+1.0)+(k/3==k%3 ? -2 : 0));
  m3=t-m1;       for(size_t k=0; k<9; k++) EXPECT_DOUBLE_EQ(m3(k),-(k+1.0)+(k/3==k%3 ? -2 : 0));
}


// Inverse and power
TEST(matrix,inverse){
  double t=-2;
  Matrix m1=mrand(3,3)*0.2+1.0;

  // inverse
  EXPECT_EQ(inv(m1).dim(), m1.dim());
  EXPECT_NEAR(norm(inv(m1)*m1-1.0)    , 0, 1.e-12);
  EXPECT_NEAR(norm((t/m1)-(t*inv(m1))), 0, 1.e-12);

  // power
  EXPECT_NEAR(norm(pow(m1, 5)-m1*m1*m1*m1*m1), 0, 1.e-12);
  EXPECT_NEAR(norm(pow(m1,-5)-inv(pow(m1,5))), 0, 1.e-12);
  EXPECT_NEAR(norm(pow(m1, 0)-1.0)           , 0, 1.e-12);
}


// Submatrix (set/get)
TEST(matrix,submat){
  Matrix m1,m2,m3;
  m1=mrand(3,4);
  m2=mrand(3,4);
  for(size_t j=0; j<m1.ncol(); j++) setvec(m2,j,getvec(m1,j));
  EXPECT_DOUBLE_EQ(norm(m1-m2), 0);

  m2=m1;
  m2.resize(5,6);
  for(size_t k=0;  k<12;   k++) EXPECT_DOUBLE_EQ(m2(k),m1(k));
  for(size_t k=12; k<30; k++) EXPECT_DOUBLE_EQ(m2(k),0);
  m2.resize(4,3);
  for(size_t k=0; k<12; k++) EXPECT_DOUBLE_EQ(m2(k),m1(k));
  m2.fill(11.1);
  for(size_t k=0; k<12; k++) EXPECT_DOUBLE_EQ(m2(k),11.1);

  m1=mrand(7,7);
  m2=getdiag(m1);
  for(size_t k=0; k<7; k++) EXPECT_DOUBLE_EQ(m2(k),m1(k,k));
  for(size_t k=0; k<7; k++) m1(k,k)=0;
  setdiag(m1,m2);
  for(size_t k=0; k<7; k++) EXPECT_DOUBLE_EQ(m2(k),m1(k,k));
}


// Concatenation (hcat/vcat)
TEST(matrix,concat){
  Matrix m1,m2,m3;

  m1=mrand(3,4);
  m2=mrand(3,5);
  m3=hcat(m1,m2);
  EXPECT_EQ(m3.nrow(),m1.nrow());
  EXPECT_EQ(m3.ncol(),m1.ncol()+m2.ncol());
  for(size_t i=0; i<m3.nrow(); i++)
    for(size_t j=0; j<m3.ncol(); j++)
      EXPECT_DOUBLE_EQ(m3(i,j), (j<m1.ncol() ? m1(i,j) : m2(i,j-m1.ncol())));

  m1=mrand(3,4);
  m2=mrand(2,4);
  m3=vcat(m1,m2);
  EXPECT_EQ(m3.nrow(),m1.nrow()+m2.nrow());
  EXPECT_EQ(m3.ncol(),m1.ncol());
  for(size_t i=0; i<m3.nrow(); i++)
    for(size_t j=0; j<m3.ncol(); j++)
      EXPECT_DOUBLE_EQ(m3(i,j), (i<m1.nrow() ? m1(i,j) : m2(i-m1.nrow(),j)));
}


// Transpose
TEST(matrix,transpose){
  Matrix m1,m2,m3;

  m1=mrand(4,3);
  m2=tp(m1);
  EXPECT_TRUE(m1.nrow()==m2.ncol() and m1.ncol()==m2.nrow());
  EXPECT_NEAR(norm(tp(tp(m1))-m1), 0, 1.e-12);  // M^T^T = M
  for(size_t i=0; i<m1.nrow(); i++)
    for(size_t j=0; j<m1.ncol(); j++)
      EXPECT_DOUBLE_EQ(m2(j,i),m1(i,j));
}


// Orthogonalization
TEST(matrix,orthgonal){
  Matrix m1,m2,m3;

  m1=mrand(4,3);
  m2=orth(m1);
  m3=tp(m2)*m1;  // should be upper triangular
  EXPECT_NEAR(norm(tp(m2)*m2 - 1.0), 0, 1.e-12);  // U^T U = I
  for(size_t i=0; i<m3.nrow(); i++) EXPECT_LE(0,m3(i,i));
  for(size_t i=0; i<m3.nrow(); i++) for(size_t j=0; j<i; j++) EXPECT_NEAR(m3(i,j),0,1.e-12);
}


// 3D operation
TEST(matrix,threedim){
  Matrix m,v0,v1,v2;
  double th;

  m=mrand(3,3);
  v0=getvec(m,0);
  v1=getvec(m,1);
  v2=getvec(m,2);
  th=urand1()*M_PI*2;

  // outer product
  EXPECT_NEAR(norm(outer(v0,v1)+outer(v0,v2)-outer(v0,v1+v2)), 0, 1.e-8);  // linear
  EXPECT_NEAR(norm(outer(v0,v1)+outer(v1,v0)), 0, 1.e-8);  // alternate
  EXPECT_NEAR(inner(outer(Matrix(1,0,0),Matrix(0,1,0)),Matrix(0,0,1)), 1, 1.e-8);  // normalized

  // v in R^3 <--> A in so(3)
  EXPECT_NEAR(norm(asym2vec(vec2asym(th*v0))-th*v0), 0, 1.e-8);  // V<-->A
  EXPECT_NEAR(norm(vec2asym(th*v0)+tp(vec2asym(th*v0))), 0, 1.e-8);  // A'=-A
  EXPECT_NEAR(norm(vec2asym(th*v0)*v0), 0, 1.e-8);  // AV=0
  EXPECT_NEAR(norm(vec2asym(th*v0)*v1-outer(th*v0,v1)), 0, 1.e-8);  // AX=VxX

  // v in R^3 <--> R in SO(3)  (rotation about v by angle |v|)
  m=orth(mrand(3,3)+3.0);
  v0=getvec(m,0);
  v1=getvec(m,1);
  v2=getvec(m,2);
  th=urand1()*M_PI*2 * 0.1;
  EXPECT_NEAR(norm(rot2vec(vec2rot(th*v0))-th*v0), 0, 1.e-8);  // V<-->R
  EXPECT_NEAR(norm(vec2rot(rot2vec(orth(1.0+m*1.e-6)))-orth(1.0+m*1.e-6)), 0, 1.e-12);  // avoid loss of significance near R=1
  EXPECT_NEAR(norm(rot2vec(Matrix(3,3)+1.0+1.e-15)), 0, 1.e-20);  // avoid NaN error near R=1
  EXPECT_NEAR(norm(rot2vec(vec2rot(Matrix(0,0,M_PI*0.99))) - Matrix(0,0,M_PI*0.99)), 0,  1.e-12);  // almost 180[deg] rotation
  EXPECT_NEAR(norm(tp(vec2rot(th*v0))*vec2rot(th*v0)-1.0), 0, 1.e-8);  // R'R=1
  EXPECT_NEAR(norm(vec2rot(th*v0)*v0-v0), 0, 1.e-8);  // RV=V
  EXPECT_NEAR(norm(vec2rot(th*v0)*v1-( cos(th)*v1+sin(th)*v2)), 0, 1.e-8);
  EXPECT_NEAR(norm(vec2rot(th*v0)*v2-(-sin(th)*v1+cos(th)*v2)), 0, 1.e-8);
  th=5.e-9;  // very small angle
  EXPECT_NEAR(norm(vec2rot(th*v0)*v0-v0), 0, 1.e-12);
  EXPECT_NEAR(norm(vec2rot(th*v0)*v1-( cos(th)*v1+sin(th)*v2)), 0, 1.e-12);
  EXPECT_NEAR(norm(vec2rot(th*v0)*v2-(-sin(th)*v1+cos(th)*v2)), 0, 1.e-12);
  EXPECT_NEAR(norm(vec2rot(th*v0)-orth(vec2rot(th*v0))), 0, 1.e-16);

  // rotation about coord axis
  th=urand1()*M_PI*2;
  for(int k=0; k<3; k++){
    m=rotabout(k,th);
    int i=(k+1)%3, j=(k+2)%3;
    v0=Matrix(3); v0(i)=1;
    v1=Matrix(3); v1(j)=1;
    v2=Matrix(3); v2(k)=1;
    EXPECT_NEAR(norm( cos(th)*v0+sin(th)*v1 - m*v0), 0, 1.e-8);
    EXPECT_NEAR(norm(-sin(th)*v0+cos(th)*v1 - m*v1), 0, 1.e-8);
    EXPECT_NEAR(norm(v2-m*v2), 0, 1.e-8);
  }
}


// Mixed calculations (exponential)
TEST(matrix,misc){
  Matrix m1,m2;

  m1=mrand(4,4);
  m2=m1-tp(m1);
  EXPECT_NEAR(norm(exp(m1)*exp(m1)-exp(2.0*m1)), 0, 1.e-12);
  EXPECT_NEAR(norm(exp(m1)*exp(-m1)-1.0)       , 0, 1.e-12);
  EXPECT_NEAR(norm(exp(m2)*tp(exp(m2))-1.0)    , 0, 1.e-12);
}


// Exceptions (size check)
TEST(matrix,exception){
  Matrix m1,m2;
  double t;

  m1=mrand(3,4);
  m2=mrand(5,6);
  t =urand1();

  EXPECT_THROW(m1+=t, std::domain_error);
  EXPECT_THROW(m1-=t, std::domain_error);
  EXPECT_THROW(m1+t , std::domain_error);
  EXPECT_THROW(m1-t , std::domain_error);
  EXPECT_THROW(t+m1 , std::domain_error);
  EXPECT_THROW(t-m1 , std::domain_error);
  EXPECT_THROW(t/m1 , std::domain_error);

  EXPECT_THROW(m1+=m2, std::domain_error);
  EXPECT_THROW(m1-=m2, std::domain_error);
  EXPECT_THROW(m1+m2 , std::domain_error);
  EXPECT_THROW(m1-m2 , std::domain_error);
  EXPECT_THROW(m1*m2 , std::domain_error);

  EXPECT_THROW(inner(m1,m2), std::domain_error);
  EXPECT_THROW(inv(m1)     , std::domain_error);
  EXPECT_THROW(pow(m1,3)   , std::domain_error);
  EXPECT_THROW(orth(m1)    , std::domain_error);

  EXPECT_THROW(getsub(m1,0,4,1,1), std::domain_error);
  EXPECT_THROW(getsub(m1,3,0,1,1), std::domain_error);
  EXPECT_THROW(getsub(m1,0,0,4,1), std::domain_error);
  EXPECT_THROW(getsub(m1,0,0,1,5), std::domain_error);
  EXPECT_THROW(setsub(m2,3,0,m1) , std::domain_error);
  EXPECT_THROW(setsub(m2,0,3,m1) , std::domain_error);
  EXPECT_THROW(getvec(m2,6)      , std::domain_error);
  EXPECT_THROW(setvec(m2,1,m1)   , std::domain_error);

  EXPECT_THROW(hcat(m1,m2), std::domain_error);
  EXPECT_THROW(vcat(m1,m2), std::domain_error);

  EXPECT_THROW(outer(Matrix(3,3),Matrix(3,1)), std::domain_error);
  EXPECT_THROW(vec2rot(Matrix(3,3)), std::domain_error);
  EXPECT_THROW(rot2vec(Matrix(3,1)), std::domain_error);
  EXPECT_THROW(vec2asym(Matrix(3,3)), std::domain_error);
  EXPECT_THROW(asym2vec(Matrix(3,1)), std::domain_error);
  EXPECT_THROW(rotabout(4,0.0), std::domain_error);
}


// Comlex matrix operations
TEST(matrix,complex){
  typedef std::complex<double> Cpx;
  typedef matrix<Cpx> MatC;
  MatC m,m1,u,d;

  m = MatC(mrand(3,3)) + MatC(mrand(3,3))*Cpx(0,1);

  // inner
  m1= MatC(mrand(3,3)) + MatC(mrand(3,3))*Cpx(0,1);
  Cpx c(urand1(),urand1());
  EXPECT_NEAR(std::abs(inner(c*m,m1)-nm_conj(c)*inner(m,m1)), 0, 1.e-8);
  EXPECT_NEAR(std::abs(inner(m,c*m1)-        c *inner(m,m1)), 0, 1.e-8);
  EXPECT_NEAR(std::abs(inner(m,m1) - 
		       (inner(real(m),real(m1)) + inner(imag(m),imag(m1)) +
			Cpx(0,1) * (-inner(imag(m),real(m1)) + inner(real(m),imag(m1))))),
	      0, 1.e-8);

  // tp
  m1=tp(m);
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<m.ncol(); j++)
      EXPECT_DOUBLE_EQ(std::norm(m(i,j)-nm_conj(m1(j,i))), 0);

  // inv
  EXPECT_NEAR(std::abs(norm(inv(m)*m-Cpx(1.0))), 0, 1.e-8);

  // orth
  u = orth(m);
  d = tp(u)*m;
  EXPECT_GT  (inner(m,m).real(), 0);
  EXPECT_NEAR(inner(m,m).imag(), 0, 1.e-8);
  EXPECT_NEAR(std::abs(norm(tp (u)*u-Cpx(1.0))), 0, 1.e-8);
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<i; j++)
      EXPECT_NEAR(std::abs(d(i,j)),0,1.e-8);

  // eigen
  m+=tp(m);
  u=eigen(m);
  d=tp(u)*m*u;
  EXPECT_NEAR(norm(tp(u)*u-Cpx(1)),0, 1.e-8);
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<m.ncol(); j++){
      if(i==j) EXPECT_NEAR(d(i,j).imag()   , 0, 1.e-8);
      else     EXPECT_NEAR(std::abs(d(i,j)), 0, 1.e-8);
    }
}
