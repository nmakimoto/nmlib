// Chronograph (stopwatch)


#ifndef CHRONO_H
#define CHRONO_H


#include <chrono>
namespace nmlib{


class Chrono{
public:
  Chrono(void){ t0=now(); }
  void   reset(void)       { t0=now(); }
  double lap  (void) const { return std::chrono::duration_cast<std::chrono::microseconds>(now()-t0).count() * 1.e-6; }

private:
  std::chrono::system_clock::time_point t0;
  std::chrono::system_clock::time_point now(void) const { return std::chrono::system_clock::now(); }
};


}  //namespace nmlib
#endif  //CHRONO_H
