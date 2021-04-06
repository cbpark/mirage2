
/** \file numerics.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: Integration of ODEs by Runge Kutta, minimum finding and
                derivative calculation

   $Log: numerics.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.7  2005/07/15 15:10:47  allanach
   Added analytic dilog routine

   Revision 1.6  2005/06/16 13:57:04  allanach
   Added a Cauchy-distributed random number generator

   Revision 1.5  2005/05/13 16:07:27  allanach
   Edited precision due to double precision

   Revision 1.4  2005/04/13 14:53:55  allanach
   Corrected binning procedure so it does what you expect

   Revision 1.3  2005/04/12 10:44:20  allanach
   Added bin function to calculate the bin number of some data

   Revision 1.2  2005/04/11 14:06:36  allanach
   Added random number routines

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.7  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.6  2003/08/14 09:25:49  allanach
   Used new standard convention for included files

   Revision 1.5  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.4  2002/10/30 16:22:48  allanach
   Fixed bug in f function

   Revision 1.3  2002/10/01 11:52:16  allanach
   Bug-fixed bound-finding routines. MGUT always determined.

   Revision 1.4  2001/08/08 09:52:33  allanach
   Added dilogarithm function - could be speeded up....

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef NUMERICS_H
#define NUMERICS_H

#include "utils.h"
#include "mycomplex.h"
#include <iostream>
using std::cout;
using std::endl;
using std::flush;
#include "def.h"
#include "linalg.h"

/// A single step of Runge Kutta (5th order), input: 
/// y and dydx (derivative of y), x is independent variable. yout is value
/// after step. derivs is a user-supplied function
void rungeKuttaStep(const DoubleVector & y, const DoubleVector & dydx, 
	     double x, double h, DoubleVector & yout, DoubleVector & yerr, 
	     DoubleVector (*derivs)(double, const DoubleVector &));

/// organises the variable step-size for Runge-Kutta evolution
int odeStepper(DoubleVector & y, const DoubleVector & dydx, double *x, double
		htry, double eps, DoubleVector & yscal, double *hdid, 
		double *hnext,		
		DoubleVector (*derivs)(double, const DoubleVector &));

/// Organises integration of 1st order system of ODEs
int integrateOdes(DoubleVector & ystart, double x1, double x2, double eps,
		  double h1, double hmin, 
		  DoubleVector (*derivs)(double, const DoubleVector &),
		  int (*rkqs)
		  (DoubleVector & y, const DoubleVector & dydx, double *x,
		   double htry, double eps, DoubleVector & yscal, double
		   *hdid, double *hnext, 
		   DoubleVector (*derivs)(double, const DoubleVector &)));

/// func is user-supplied, h is an estimate of what step-size to start with
/// and err returns error flags
double calcDerivative(double (*func)(double), 
		     double x, double h, double *err);

/// f is user-defined function, minimum value returned in xmin
double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin);

void shft2(double & a, double & b, double & c); ///< a=b and b=c
/// a=b, b=c and c=d
void shft3(double & a, double & b, double & c, double & d); 

/// For calculation of PV functions
double integrandThreshbn(double x);
/// Returns real part of b function, less accurate than analytic expressions
double bIntegral(int n, double p, double m1, double m2, double mt);
DoubleVector dd(double x, const DoubleVector & y);

/// Passarino-Veltman function definition
double b0(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b1(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b22(double p,  double m1, double m2, double q);
/// Passarino-Veltman function definition
double c0(double m1, double m2, double m3);
/// Passarino-Veltman function definition
double d27(double m1, double m2, double m3, double m4);
/// Passarino-Veltman function definition
double d0(double m1, double m2, double m3, double m4);

// inlined PV functions
inline double a0(double m, double q) {
  if (m == 0.0) return 0.0;
  return sqr(m) * (1.0 - log(sqr(m / q)));
}

inline double ffn(double p, double m1, double m2, double q) {
  return a0(m1, q) - 2.0 * a0(m2, q) - 
    (2.0 * sqr(p) + 2.0 * sqr(m1) - sqr(m2)) * 
    b0(p, m1, m2, q);
}

inline double gfn(double p, double m1, double m2, double q) {
  return (sqr(p) - sqr(m1) - sqr(m2)) * b0(p, m1, m2, q) - a0(m1, q) 
    - a0(m2, q); 
}

inline double hfn(double p, double m1, double m2, double q) {
  return 4.0 * b22(p, m1, m2, q) + gfn(p, m1, m2, q);
}

inline double b22bar(double p, double m1, double m2, double q) {
  return b22(p, m1, m2, q) - 0.25 * a0(m1, q) - 0.25 * a0(m2, q);
}

inline double fB(const Complex & x) {
  return (log(1.0 - x) - x * log(1.0 - 1.0 / x) - 1.0).real();
}

double dilogarg(double t);
double dilog(double x);

double integrandThreshbnr(double x);
Complex fnfn(double x);

/// Gaussian deviated random number, mean 0 variance 1. Don't re-set idum once
/// you've initially set it. Initialise with a NEGATIVE integer
double gasdev(long & idum);
/// Normally distributed random number between 0 and 1. Don't re-set idum once
/// you've initially set it. Initialise with a NEGATIVE integer
double ran1(long & idum);
/// Cauchy distribution ie 1 / [ pi gamma (1 + x^2/gamma^2) ]. For a width,
/// you must multiply the x coming out by the width. 
double cauchyRan(long & idum);

/// Returns the number of a bin that the data is in: from 1 to numBins in the
/// range (bins other than this range are also possible - you must deal with
/// them outside the function...)
int bin(double data, double start, double end, int numBins);

/// Adds logs of two numbers in a more careful way that avoids underflow
double logOfSum(double a, double b);
#endif

