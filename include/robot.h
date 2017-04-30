// kinematics of 6-axis manipulator


#ifndef ROBOT_H
#define ROBOT_H


#include <iostream>
#include "matrix.h"
namespace nmlib{


/******** Class I/F ********/

// 6-axis manipulator
class Robot{
public:
  Robot(void);

  // kinematics
  Matrix fk      (const Matrix& th, int k=7) const;  // forward kinematics (of k-th part)
  Matrix jacobian(const Matrix& th, int k=7) const;  // jacobian of fk (of k-th part)
  Matrix ik      (const Matrix& hom_dst, const Matrix& th_ini, double err=1.e-8, int iter=100) const;  // inverse kinematics

  // configuration
  Matrix  zero(void)  const;  // get offset of joint angles
  Matrix& zero(void);         // set ...
  Matrix  hom0(int k) const;  // get position and attitude of k-th part
  Matrix& hom0(int k);        // set ...

private:
  Matrix hm0[6+2];  // relative position and attitude of adjacent parts
  Matrix th0;       // offset of joint angles

  // Notes:
  // - For k=0..7, hm0[k] represents the homog. trsf. from (k-1)-st part to k-th part when (k-1)-st joint angle is 0.
  //   By convention, -1st and 7th parts represent the world and an end effector, respectively.
  // - For k=1..6, the Z axis of k-th part is bound to the rotation axis of (k-1)-st joint.
  //   Note that there are some freedom of choice of the origin, X axis, and Y axis of k-th part.
  // - If relevant geometry is given in the world coordinate, then initialization code may look like:
  //     robot.hom0(k) = inv_homtrsf(T[k-1])*T[k];
  //     robot.zero()  = robot.ik(T_ref,th_ini);
  // - IK may fail (e.g. near singularities). User should verify convergence.
};


/******** Utility I/F ********/

// stream I/O
std::istream& operator>>(std::istream& str,       Robot& robot);
std::ostream& operator<<(std::ostream& str, const Robot& robot);

// homogeneous transformation matrix <--> (position,rotation)
Matrix hom2pos(const Matrix& hom);
Matrix hom2rot(const Matrix& hom);
Matrix homtrsf(const Matrix& pos, const Matrix& rot);
Matrix inv_homtrsf(const Matrix& hom);
bool   is_homtrsf (const Matrix& hom, double err=1.e-8);


/******** Implementation ********/

inline Robot::Robot(void){
  // approximate model of Yaskawa Motoman ES280D-230
  Matrix eye=rotabout(2,0.0), x90=rotabout(0,M_PI/2), yzx=rotabout(1,M_PI/2)*rotabout(2,M_PI/2);
  hm0[0] = homtrsf(Matrix({   0,    0,  0}),eye);
  hm0[1] = homtrsf(Matrix({ 285,    0,650}),tp(x90));
  hm0[2] = homtrsf(Matrix({   0,-1150,  0}),eye);
  hm0[3] = homtrsf(Matrix({-250,-1015,  0}),x90);
  hm0[4] = homtrsf(Matrix({   0,    0,  0}),tp(x90));
  hm0[5] = homtrsf(Matrix({   0, -250,  0}),x90);
  hm0[6] = homtrsf(Matrix({   0,    0,200}),yzx);
  hm0[7] = homtrsf(Matrix({   0,    0,  0}),eye);

  th0 = Matrix(6);
}


// forward kinematics (of k-th part)
inline Matrix Robot::fk(const Matrix& th, int k) const{
  if(k<0 || 7<k) throw std::domain_error("Robot::fk: bad part number");
  Matrix h=hm0[0], rz=Matrix(4,4)+1.0;
  for(int j=0; j<6 && j<k ; j++){
    double thj=th(j)+th0(j);
    setsub(rz, 0,0, rotabout(2,thj));
    h=h*hm0[j+1]*rz;
  }
  if(k==7) h=h*hm0[7];
  return h;
}


// jacobian of fk (of k-th part)
inline Matrix Robot::jacobian(const Matrix& th, int k) const{
  if(k<0 || 7<k) throw std::domain_error("Robot::jacobian: bad part number");
  Matrix jac(6,6), hj=hm0[0], hk=fk(th,k), rk=hom2rot(hk), jac1, rz(4,4), drz(4,4);
  for(int j=0; j<6 && j<k; j++){
    double thj=th(j)+th0(j);
    setsub( rz, 0,0, rotabout(2,thj));          rz(2,2)= rz(3,3)=1;
    setsub(drz, 0,0, rotabout(2,thj+M_PI/2));  drz(2,2)=drz(3,3)=0;
    jac1 = hj*hm0[j+1]*drz;           // T_0 ... T_j dT_{j+1}
    hj   = hj*hm0[j+1]* rz;           // T_0 ... T_j  T_{j+1}
    jac1 = jac1*inv_homtrsf(hj)*hk;   // T_0 ... T_j dT_{j+1} T_{j+2} ... T_k
    setsub(jac, 0,j, hom2pos(jac1));  // velocity = dX/dth
    setsub(jac, 3,j, asym2vec(hom2rot(jac1)*tp(rk)));  // angular velocity = dR/dth R^-1
  }
  return jac;
}


// inverse kinematics
inline Matrix Robot::ik(const Matrix& hom_dst, const Matrix& th_ini, double err, int iter) const{
  Matrix th=th_ini, p0=hom2pos(hom_dst), r0=hom2rot(hom_dst);
  while(iter--){
    Matrix t=fk(th), dp=hom2pos(t)-p0, dr=rot2vec(hom2rot(t)*tp(r0)), dth;
    if(norm(dp)+norm(dr) < err) return th;
    dth = inv(jacobian(th)) * vcat(dp,dr);  // newton
    dth = (norm(dth)<1 ? dth : norm(dth)>1 ? dth/norm(dth) : Matrix(6).fill(0.01));  // avoid singularity
    th -= dth;
  }
  return th;  // poor convergence; should be verified in the user code
}


// configuration
inline Matrix  Robot::zero(void)  const { return th0; }
inline Matrix& Robot::zero(void)        { return th0; }
inline Matrix  Robot::hom0(int k) const { return hm0[k]; }
inline Matrix& Robot::hom0(int k)       { return hm0[k]; }


// stream I/O
inline std::istream& operator>>(std::istream& str,       Robot& robot){
  for(int i=0; i<6+2; i++) str>>robot.hom0(i);
  str>>robot.zero();
  return str;
}
inline std::ostream& operator<<(std::ostream& str, const Robot& robot){
  for(int i=0; i<6+2; i++) str<<robot.hom0(i)<<"\n";
  str<<robot.zero()<<"\n";
  return str;
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


}  //namespace nmlib
#endif //ROBOT_H
