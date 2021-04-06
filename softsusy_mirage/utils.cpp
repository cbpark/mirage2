
/** \file utils.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: utils.cpp,v $
   Revision 1.3  2005/11/09 14:12:25  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.2  2005/07/15 15:17:04  allanach
   Checking one calculation by using an underflow rather than compare with zero

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.15  2004/09/28 10:39:07  allanach
   Added Les Houches output function for output to a file

   Revision 1.14  2004/05/12 19:18:00  allanach
   Made the close function safer: it's no longer sensitive to zeros

   Revision 1.13  2004/04/07 16:19:09  allanach
   Added proper NaN tests

   Revision 1.12  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.11  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.10  2003/03/28 15:59:35  allanach
   Two-loop stuff relegated to seperate file

   Revision 1.9  2003/03/07 11:32:48  allanach
   Added tested 2-loop alpha_t^2, alpha_s alpha_t, alpha_s alpha_b corrections
   to REWSB and bug-fixed the Higgs mass calculation

   Revision 1.8  2003/03/04 12:47:44  allanach
   Deleted call to f2c.h and included it directly

   Revision 1.7  2003/02/26 10:11:06  allanach
   Put two-loop Higgs routines from Slavich and co in

   Revision 1.6  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.3  2002/05/01 16:02:20  allanach
   Added checkTolerance subroutine

 */

#include "utils.h"

int theta(double a) {
  int temp = 0;
  if (a > 0.0) temp = 1;
  return temp;
}

// Just sets precision and format of outputs
void outputCharacteristics(int n) {
  cin.setf(ios::scientific, ios::floatfield);
  cin.precision(n);
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(n);
  cerr.setf(ios::scientific, ios::floatfield);
  cerr.precision(n);
}

// Finds fractional difference between |a| and |b|
double toleranceCheck(double a, double b) {
  double sTin = fabs(a), sTout = fabs(b);
  double maxx = maximum(sTin, sTout);

  const double underflow = 1.0e-20;

  if (maxx < underflow) return 0.0;
  return fabs(1.0 - minimum(sTin, sTout) / maxx);
}

// Outputs a space if greater than zero, a minus otherwise.
// Useful for outputting negative numbers in rows
void printRow(double x) {

  // make it return a character when you've worked out the equivalent of printf

  double underflow = 1.0e-120;
  if (fabs(x) < underflow) x = 0.0; // Traps -0.0
  if (x >= 0.0) cout << " " << x;
  else cout << x;
}

bool testNan(double f) {
  if (f >= 0.0 || f <= 0.0) return false;
  return true;
}

bool close(double m1, double m2, double tol) {
  return (fabs(maximum(fabs(m1), fabs(m2)) - fabs(minimum(fabs(m1), fabs(m2))))
	  <= tol * maximum(fabs(m1), fabs(m2)));
}

// Outputs a space if greater than zero, a minus otherwise.
// Useful for outputting negative numbers in rows
void printRow(fstream & f, double x) {
  // make it return a character when you've worked out the equivalent of printf
  double underflow = 1.0e-120;
  if (fabs(x) < underflow) x = 0.0; // Traps -0.0
  if (x >= 0.0) f << " " << x;
  else f << x;
}

