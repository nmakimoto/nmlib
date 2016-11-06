// kinematics of 6-axis manipulator


#ifndef ROBOT_H
#define ROBOT_H


#include <iostream>
#include "matrix.H"
namespace nmlib{


/******** Class I/F ********/

// 6-axis manipulator
class Robot{
public:
  // kinematics
  Matrix fk(const Matrix& th, int k=6) const;  // forward kinematics (of k-th part)
  Matrix ik(const Matrix& hom, const Matrix& th_ini, double err=1.e-8, int iter=100) const;  // inverse kinematics
  Matrix jacobian(const Matrix& th, int k=6) const;  // jacobian of fk (of k-th part)

  // configuration
  void set_zero(const Matrix& th0);  // set offset of joint angles
  void set_base(const Matrix& hom);  // set position and attitude of the robot
  Matrix& link(int k);               // set parameters of k-th joint
  Matrix  link(int k) const;         // get parameters of k-th joint

private:
  Matrix  lk(int k, double th) const;  // homog trsf of k-th part relative to (k-1)-th when k-th joint angle=theta
  Matrix dlk(int k, double th) const;  // dlk/dth

  // homogeneous transformations between adjacent parts, connected by k-th joint, at zero angles
  //   - lk0[k] is the homog trsf from (k-1)-th part to k-th part when k-th joint angle=0
  //   - rotation axis of k-th joint is "Z"-axis of k-th part
  //   - for k=0, "(k-1)-th part" is the world
  //   - for k=6, "k-th part" is a fixed end effector
  Matrix lk0[6+1];
};


/******** Utility I/F ********/

// homogeneous transformation matrix <--> (position,rotation)
Matrix hom2pos(const Matrix& hom);
Matrix hom2rot(const Matrix& hom);
Matrix homtrsf(const Matrix& pos, const Matrix& rot);
Matrix inv_homtrsf(const Matrix& hom);
bool   is_homtrsf (const Matrix& hom, double err=1.e-8);

// stream I/O
std::istream& operator>>(std::istream& str,       Robot& robot);
std::ostream& operator<<(std::ostream& str, const Robot& robot);


/******** Implementation ********/

// forward kinematics (of k-th part)
inline Matrix Robot::fk(const Matrix& th, int k) const{
  Matrix h=Matrix(4,4)+1.0;
  for(int j=0; j<=k and j<6; j++) h=h*lk(j,th(j));
  if(k==6) h=h*lk(k,0);
  return h;
}


// inverse kinematics
inline Matrix Robot::ik(const Matrix& hom0, const Matrix& th_ini, double err, int iter) const{
  Matrix th=th_ini, p0=hom2pos(hom0), r0=hom2rot(hom0);
  while(iter--){
    Matrix t=fk(th), dp=hom2pos(t)-p0, dr=rot2vec(hom2rot(t)*tp(r0));
    if(norm(dp)+norm(dr) < err) return th;
    th -= inv(jacobian(th)) * vcat(dp,dr);  // newton
  }
  return th;  // poor convergence; should be verified in the user code
}


// jacobian of fk (of k-th part)
inline Matrix Robot::jacobian(const Matrix& th, int k) const{
  if(k==6) k--;
  Matrix jac(6,6), hj(4,4), hk=fk(th,k), rk=hom2rot(hk), jac1;
  hj+=1.0;
  for(int j=0; j<=k; j++){
    jac1 = hj*dlk(j,th(j));  // T_0 ... T_{j-1} dT_j
    hj   = hj* lk(j,th(j));  // T_0 ... T_j
    jac1 = jac1*inv_homtrsf(hj)*hk;  // T_0 ... T_{j-1} dT_j T_{j+1} ... T_k
    setsub(jac, 0,j, hom2pos(jac1));  // velocity = dpos/dth_j
    setsub(jac, 3,j, asym2vec(hom2rot(jac1) * tp(rk)));  // angular velocity = drot/dth_j rot^-1
  }
  return jac;
}


// configuration
inline void Robot::set_zero(const Matrix& th0){ for(int j=0; j<6; j++) link(j)=lk(j,th0(j)); }
inline void Robot::set_base(const Matrix& hom){ link(0)=hom; }
inline Matrix& Robot::link(int k)      { return lk0[k]; }
inline Matrix  Robot::link(int k) const{ return lk0[k]; }


// homogeneous transformation of k-th part relative to (k-1)-th when joint angle=theta
inline Matrix Robot::lk(int k, double th) const{
  if(k==6 or th==0) return lk0[k];
  Matrix m(4,4);
  m(0,0)=cos(th); m(0,1)=-sin(th);
  m(1,0)=sin(th); m(1,1)= cos(th);
  m(2,2)=m(3,3)=1;
  return lk0[k]*m;
}
// dlk/dth
inline Matrix Robot::dlk(int k, double th) const{
  Matrix m(4,4);
  m(0,0)=-sin(th); m(0,1)=-cos(th);
  m(1,0)= cos(th); m(1,1)=-sin(th);
  return lk0[k]*m;
  
}


// homogeneous transformation matrix T <--> (position P,rotation R)
inline Matrix hom2pos(const Matrix& t){ return getsub(t,0,3,3,1); }
inline Matrix hom2rot(const Matrix& t){ return getsub(t,0,0,3,3); }
inline Matrix homtrsf(const Matrix& p, const Matrix& r){
  Matrix t(4,4);
  setsub(t,0,0,r);  // T = |R P|  P:3-vector, R:orthogonal
  setsub(t,0,3,p);  //     |0 1|
  t(3,3)=1;
  return t;
}
// inverse of homg trsf
inline Matrix inv_homtrsf(const Matrix& t0){
  Matrix t(4,4), r0=hom2rot(t0), p0=hom2pos(t0);
  setsub(t,0,0, tp(r0));     // T^-1 = |R^T -R^T P|
  setsub(t,0,3,-tp(r0)*p0);  //        |0   1     |
  t(3,3)=1;
  return t;
}
// returns true if hom is actually homog trsf matrix
inline bool is_homtrsf(const Matrix& hom, double err){
  return norm(hom*inv_homtrsf(hom)-1.0)<err;
}


// stream I/O
inline std::istream& operator>>(std::istream& str,       Robot& robot){
  for(int i=0; i<6+1; i++) str>>robot.link(i);
  return str;
}
inline std::ostream& operator<<(std::ostream& str, const Robot& robot){
  for(int i=0; i<6+1; i++) str<<robot.link(i)<<"\n";
  return str;
}


}  //namespace nmlib
#endif //ROBOT_H
