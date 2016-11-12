// All-in-one header of nmlib


#ifndef NMLIB_H
#define NMLIB_H


#ifdef __linux__
#include "chrono.H"
#endif

#include "diff.H"
#include "fft.H"
#include "io.H"
#include "lp.H"
#include "matrix.H"
#include "matrix_decomp.H"
#include "polynomial.H"
#include "random.H"  // compile and link random.C to use pseudo RNGs
#include "robot.H"
#include "solver.H"
#include "sparse.H"  // compile and link sparse.C to use sparse solvers
#include "spline.H"
#include "stat.H"

#ifdef NMLIB_DISABLE_NAMESPACE
using namespace nmlib;
#endif


#endif //NMLIB_H
