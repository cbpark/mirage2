
/** \file utils.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: A few handy bits and pieces - little mathematical functions 
                and the like

   $Log: utils.h,v $
   Revision 1.3  2005/11/09 14:12:25  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.10  2004/09/28 10:39:07  allanach
   Added Les Houches output function for output to a file

   Revision 1.9  2004/05/12 19:18:00  allanach
   Made the close function safer: it's no longer sensitive to zeros

   Revision 1.8  2004/04/07 16:19:09  allanach
   Added proper NaN tests

   Revision 1.7  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.6  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.5  2003/03/28 15:59:35  allanach
   Two-loop stuff relegated to seperate file

   Revision 1.4  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.3  2002/05/01 16:02:20  allanach
   Added checkTolerance subroutine

*/

#ifndef UTILS_H
#define UTILS_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <def.h>
#include <stdio.h>
using namespace std;

/// Standard theta function: 1 is a>0, 0 otherwise
int theta(double a);
/// Sets number of decimal places to display in scientific output
void outputCharacteristics(int);
/// square of a number
inline double sqr(double a) { return a * a; }
/// maximum of two numbers
inline double maximum(double a, double b) { return ((a > b) ? a : b); }
/// minimum of a and b
inline double minimum(double a, double b) { return ((a < b) ? a : b); }
/// minimum of a and b
inline int minimum(int a, int b) { return ((a < b) ? a : b); }
/// Finds fractional difference between |a| and |b|
double toleranceCheck(double sTin, double sTout);

/// checks if ABSOLUTE (or squared) values are closer than tol
bool close(double m1, double m2, double tol);

/// Returns |a| with sign of b in front
inline double sign(double a, double b) 
{ return ((b) >= 0.0 ? fabs(a) : -fabs(a)); }

/// gives the sign (+-1) of x
inline int sgn(double x)
{ return (x >= 0.0 ? 1 : -1); }

/// Outputs a space if greater than zero, a minus otherwise.
/// Useful for outputting negative numbers in rows
void printRow(double x);
void printRow(fstream & f, double x);

/// Returns true if f's a nan
bool testNan(double f);
#endif

