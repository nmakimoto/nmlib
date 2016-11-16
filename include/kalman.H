// Extended Kalman filter


#ifndef KALMAN_H
#define KALMAN_H


#include <functional>
#include "diff.H"
#include "matrix.H"
namespace nmlib{


/******** Class I/F ********/

// Extended Kalman filter
class Kalman{
public:
  template<class TM> void predict(double        t1, const TM& f, const Matrix& Q);  // predict phase
  template<class OM> void update (const Matrix& y1, const OM& h, const Matrix& R);  // update phase
  template<class TM> Matrix mean(double t1, const TM& f)                  const;  // state at time t1
  template<class TM> Matrix cov (double t1, const TM& f, const Matrix& Q) const;  // covariance at time t1

  double t;   // current time
  Matrix x;   // current state vector
  Matrix P;   // current covariance
  Matrix dx;  // pitch width for numerical differenciation
};


/******** Implementation ********/

// execute predict phase
// t1: new time, f(x,delta_t): transition model, Q: covariance of system noise
template<class TM> void Kalman::predict(double t1, const TM& f, const Matrix& Q){
  Matrix F=jacobian(std::bind(f,std::placeholders::_1,t1-t),x,dx);
  P=F*P*tp(F)+Q;
  x=f(x,t1-t);
  t=t1;
};

// execute update phase
// y1: new observation data, h(x): observation model, R: covariance of observation error
template<class OM> void Kalman::update(const Matrix& y1, const OM& h, const Matrix& R){
  Matrix H=jacobian(h,x,dx);
  Matrix K=P*tp(H)*inv(H*P*tp(H)+R);
  x+=K*(y1-h(x));
  P =(1.0-K*H)*P;
};

// calculate expected state vector at future time t1
template<class TM> Matrix Kalman::mean(double t1, const TM& f) const{
  return f(x,t1-t);
}

// calculate expected covariance at future time t1
template<class TM> Matrix Kalman::cov(double t1, const TM& f, const Matrix& Q) const{
  Matrix F=jacobian(std::bind(f,std::placeholders::_1,t1-t),x,dx);
  return P=F*P*tp(F)+Q;
}


} //namespace
#endif //KALMAN_H
