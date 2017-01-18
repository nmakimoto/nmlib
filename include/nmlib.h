// All-in-one header of nmlib


#ifndef NMLIB_H
#define NMLIB_H


#include "chrono.h"
#include "diff.h"
#include "fft.h"
#include "io.h"
#include "kalman.h"
#include "lp.h"
#include "matrix.h"
#include "matrix_decomp.h"
#include "polynomial.h"
#include "random.h"  // compile and link random.cpp to use non-class RNG functions
#include "robot.h"
#include "route.h"
#include "solver.h"
#include "sparse.h"  // compile and link sparse.cpp to use sparse solvers
#include "spline.h"
#include "stat.h"

#ifdef NMLIB_DISABLE_NAMESPACE
using namespace nmlib;
#endif


#endif //NMLIB_H
