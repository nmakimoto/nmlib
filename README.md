# nmlib
A tiny C++ math library


---


## Introduction


### What's this?

This is a C++ header library for mathematical calculations. Main features include:

  - Matrix calculations (basic operations, eigenvalues, etc.)
  - Nonlinear equation solvers and optimizers
  - ODE solvers
  - LP solvers
  - Sparse matrix solvers
  - Statistics (basic statistics, linear regression, PCA)
  - Numerical integration and differentiation
  - Spline interpolation
  - FFT
  - Random number generators and low-discrepancy sequences
  - Dijkstra shortest path algorithm
  - Kalman filter
  - Robot kinematics
  - Chronograph
  - Generic I/O utilities

and others. Design policy is "simple and easy." Enjoy!


### How to use

- C++11 compiler is required.
- Copy `include/*.h` to an include path and add the following lines in your code. All is ready.
    `#include "nmlib.h"`
    `using namespace nmlib;`
- The code should work on any OS, except that binary I/O may fail on Windows due to CR/LF conversion.
  Unit test environment: ubuntu 16.04 / g++ 5.4.0 with -std=c++11 option / googletest 1.7.0.
- Below is an example usage. See `include/*.h` for full functionality.


---


## matrix.h
This is a matrix template library, a cruicial part of nmlib. It is also used throughout the library.

- Basic operations (+-*, inverse, sub/super matrices, component access, etc.)
- Eigensolvers (symmetric and hermitian, Jacobi and QR)
- 3D geometric operations (cross product, rotation, etc.)
- Decompositions (LU/LUP, QR, SVD)

```c++
Matrix a={{3,4},{4,3}}, b={5,6}, x;
x=inv(a)*b;

Matrix u,d;
u=eigen(a);   // orthoginal matrix U consisting of unit eigenvectors of A
d=tp(u)*a*u;  // diagonal matrix D=U'AU of corresponding eigenvalues
```


## chrono.h
- Simple stopwatch

```c++
Chrono ch;
//...task...
double t=ch.lap();  // elapsed time[sec] since constructed or last reset
ch.reset();         // set lap() to zero for subsequent measurements
```


## fft.h
- FFT and IFFT (1D / 2D)

```
int n=1024;  // must be a power of 2
std::vector<std::complex<double>> x(n),y;
for(int i=0; i<n; i++) x[i]=...;
y= fft(x);
x=ifft(y);
```


## io.h
- Generic converters of type T to/from std::string
- Generic I/O utils of std::map<T> and std::vector<T> to/from std::iostream

```c++
typedef double T;  // any type with stream I/O operators (<< and >>), such as Matrix and double
T x;
std::string s="3.14";
x=str2any<T>(s);
s=any2str(x);
```


## integral.h
- Numerical integration by Simpson's rule / Gaussian quadrature / DE (double exponential) formula
- Pseudo Monte Carlo multiple integration

```c++
double a=0, b=1, h=0.1, y;
int ndiv=100, iter=1000000;

double f(double x);
y=integral   (f,a,b,ndiv);  // integral of f(x) over [a,b] by Simpson
y=integral_gq(f,a,b,ndiv);  // integral of f(x) over [a,b] by Gaussian quadrature
y=integral_de(f,ndiv,h);    // infinite integral of f(x) by DE formula

double g(const Matrix& x);
y=integral_mc(g,{a,a,a},{b,b,b},{2,3,5},iter);  // integral of g(x) over [a,b]^3 by pseudo MC ({2,3,5}: base of Halton)
```


## lp.h
- LP solver by 2-step simplex

```c++
Matrix a,b,c;
std::cin >> a >> b >> c;

LP lp(a,b,c);
LP::Status s=lp.solve();  // solve problem: max_x cx s.t. Ax<=b and x>=0
                          // status: Infeasible / Feasible (but non-optimal) / Optimal / Unbounded / Unknown
Matrix x=lp.vertex();     // retrieve current basic solution
```


## ode.h
- Generic solvers for ODE initial value problem: dy/dt=f(t,y) on [t0,t1], y(t0)=y0
- 4th order classical / impllicit / adaptive Runge-Kutta

```c++
typedef Matrix Y;                             // type of y (typically Matrix or double)
Y      f(double t, const Y& y);
Y      y0={1.0,2.0,3.0};
double t0=0.0, t1=3.0, tol=1.e-8;
int    ndiv=100;

std::map<double,Y> t2y;                       // solution holder (map from t to y(t))
t2y = solve_ode_rk4 (f,y0,t0,t1,ndiv);        // classical RK4
t2y = solve_ode_rk4i(f,y0,t0,t1,ndiv);        // implicit
t2y = solve_ode_rk4a(f,y0,t0,t1,ndiv,tol);    // adaptive (aka embedded)

//Spline g(t2y);  // can be used to interpolate scalar-valued solution (when Y=double)
//for(double t=0; t<3.01; t+=0.1) std::cout << t << '\\t' << g(t) << '\\t' << g.grad(t) << '\\n';
```


## optimization.h
- Nonlinear optimization by Nelder-Mead (aka amoeba algorithm)
- Nonlinear least squares by Levenberg-Marquardt (aka LMA)
- Nonlinear least square curvefit by LMA

```c++
Matrix th_opt, th_ini={...}, dth={...};  // optimal parameter \\theta, initial guess, step of numerical diff.

double f(const Matrix& th);
int    iter=1000;
double r0=5.0, tol=1.e-6;                // initial "amoeba" size, tolerance
th_opt=opt_amoeba(th_ini,r0,iter,tol);   // argmin_th f(th)

Matrix g(const Matrix& th);
th_opt=opt_lma(g,th_ini,dth);            // argmin_th ||g(th)||^2

double h(double x, const Matrix& th);
Matrix xx, yy;
std::cin >> xx >> yy;                    // sample data for 1D curvefit y=h(x|th)
th_opt=curvefit(h,xx,yy,th_ini,dth);     // argmin_th \\sum_k |y_k- h(x_k,th)|^2
```


## random.h
- Random number generator (uniform, gaussian, exponential, multi-dimensional)
- Low-discrepancy sequence (van der Corput, Halton)
- Rejection sampling from given PDF

```c++
Rng rng;
Lds lds;
auto pdf = [](double x){ return sqrt(1-x*x); };  // example probability density function
double pmax=1;                                   // and its upper bound
  
for(int i=0; i<100; i++){
  double x1,x2,x3,x4;
  x1=lds();    // LDS on (0,1)
  x2=rng();    // U(0,1) uniform
  x3=rng.n();  // N(0,1) gaussian
  x4=rng.rand_pdf(pdf, pmax,-1.0,1.0);  // sampling from given PDF
  //...
}
```


## robot.h
- Forward / Inverse kinematics of 6-axis manipulator
- Example robot model is built-in and used by the default constructor Robot()

```c++
Matrix angle0, homog0, angle, homog;  // joint angle vectors and homogeneous transformation matrices
angle0 = {0,0,0,0,0,0};
homog0 = homtrsf({1000,1000,1000},{{1,0,0},{0,1,0},{0,0,1}});

Robot robot;                          // built-in robot model (can be configured)
homog = robot.fk(angle0);             // calculate FK at angle0
angle = robot.ik(homog0,angle0);      // solve IK near angle0
```


## route.h
- Dijkstra shortest-path algorithm on oriented graph

```c++
Route route;

route.dist[5][7]=10;  // set distance of adjacent nodes as necessary
route.dist[7][5]=15;  // note that symmetry is not assumed
route.dist[7][2]=12;
//route.dist[i][j]=d_ij;
//...

auto path=route(5,2);
for(auto node: path) std::cout << node << '\\n';
```


## solver.h
- 1D/multidim nonlinear equation solver (Newton)
- 1D nonlinear equation solver (bisection)

```c++
typedef Matrix T;              // type of variables (typically Matrix or double)
T f(const T& x);
T y0=..., x0=..., dx=...;      // target, initial guess, step of numerical diff.
T x=solve(f,y0,x0,dx);         // solve f(x)=y0 near x0 by newton

//double  x1=-10, x2=+10;      // initial interval of bisection
//x=solve_bisect(f,y0,x1,x2);  // solve f(x)=y0 on [x1,x2] by bisection (when T is scalar)
```


## sparse.h
This is an experimental implementation of sparse matrix solvers.
- Krylov sparse solvers (CG, BiCG, PBCG)
- Symmetric eigensolver (shifted inverse iteration)

```c++
int n=10000;
Sparse a(n,n);
Matrix b(n),x;

a+=4.0;  // diagonal part
for(int k=1; k<n; k++) a(k-1,k)=a(k,k-1)=1;
for(int k=0; k<n; k++) b(k)=10;

SparseConf conf;
conf.init_ilu0(a);      // use ILU(0) preconditioner (optional)
x=solve_pbcg(a,b,conf)  // solve Ax=b by PBCG

SparseEigenConf conf2;
conf2.mu0=1.0;          // use initial guess mu0 of an eigenvalue (optional)
x=eigenvec(a,conf2);    // a unit eigenvector of A
double mu=inner(x,a*x);
```


## spline.h
- Natural spline / linear / nearest neighbour interpolation
- Gradient of spline

```c++
std::map<double,double> x2y={{1,1},{2,4},{3,3},{5,4}};
Spline f(x2y);
for(double x=0; x<5.01; x+=0.1) std::cout << x << '\\t' << f(x) << '\\t' << f.grad(x) << '\\n';
```


## stat.h
- Basic statistics (average, variance, major percentile points, etc.)
- Linear regression
- Principal component analysis (PCA)
- Density / cummulative probability / p-value of normal distribution

```c++
Matrix xx, yy;
std::cin >> xx >> yy;         // j-th column vector represents j-th sample
const size_t dim=xx.nrow(), ndata=xx.ncol();

Matrix reg, err, one;
reg=regression(xx,yy,true);  // regression coefficients with a constant term
one=Matrix(1,ndata).fill(1);
err=variance(reg*vcat(xx,one)-yy);

Matrix eig;
eig=principal(xx);           // principal components (sorted eigenvectors of Var(X))
```
