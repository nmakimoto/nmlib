// Chronograph (stopwatch) for Linux


#ifndef CHRONO_H
#define CHRONO_H


#include <sys/time.h>


class Chrono{
public:
  Chrono(void){ t0=now(); }
  double now  (void) const { struct timeval tv;  gettimeofday(&tv,0);  return tv.tv_sec+tv.tv_usec*(1.e-6); }
  double lap  (void) const { return now()-t0; }
  void   reset(void)       { t0=now(); }

private:
  double t0;
};


#endif  //CHRONO_H
