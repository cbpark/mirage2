
/** \file rge.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: rge.cpp,v $
   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.12  2003/12/18 16:32:05  allanach
   Cleaned up and corrected assumed accuracy of running

   Revision 1.11  2003/10/28 13:58:32  allanach
   Bug-fixed running to negative scales

   Revision 1.10  2003/10/27 15:40:58  allanach
   If you try to run to a negative number, it will evolve to the absolute
   magnitude of the scale

   Revision 1.9  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.8  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.7  2002/10/22 13:14:01  allanach
   Bug fixed

   Revision 1.6  2002/09/09 10:42:54  allanach
   TOLERANCE replaces EPS as being more user-friendly

   Revision 1.5  2002/09/04 13:54:54  allanach
   Gauge unification condition now imposed if mx is input as negative

   Revision 1.4  2002/09/03 14:30:46  allanach
   Got rid of superfluous extern statements

   Revision 1.3  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, TOLERANCE out of def.h to make it quicker to 
   compile once they are changed.

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info

*/

#include "rge.h"

static RGE * tempRge;

// runto/run functions return >0 if there's a problem with the running
int RGE::runto(double x2, double eps) {
  double tol;
  if (eps < 0.0) tol = TOLERANCE;
  else if (eps < EPSTOL) tol = EPSTOL;
  else tol = eps;
  double x1 = this -> displayMu();
  return run(x1, x2, tol);
}

int RGE::run(double x1, double x2, double eps) {
  double tol;

  if (eps < 0.0) tol = TOLERANCE;
  else if (eps < EPSTOL) tol = EPSTOL;
  else tol = eps;

  tempRge = this;
  DoubleVector y(this -> display());
  int err = callRK(x1, x2, y, allDerivs, tol);
  if (err == 0) {
    this -> set(y);
    tempRge -> setMu(x2);
  }
  return err;
}

// Runge-Kutte user defined routine: given log renorm scale x and the dependent
// variables y of an RGE, will calculate the derivitives dydx.
DoubleVector allDerivs(double x, const DoubleVector & y)
{
  //  cout << " x=" << exp(x) << " y=" << y; // DEBUG
  tempRge->setMu(exp(x));
  tempRge->set(y);
  return tempRge->beta();
}

//Does the actual calling of Runge Kutta: default precision is TOLERANCE defined in
//def.h
//Returns >0 if there's a problem with the running
int RGE::callRK(double x1, double x2, DoubleVector & v,
		DoubleVector (*derivs)(double, const DoubleVector &), 
		double eps) {
  double tol;
  if (eps < 0.0) tol = TOLERANCE;
  else if (eps < EPSTOL) tol = EPSTOL;
  else tol = eps;
  // x1 == x2 with high precision
  if (close(fabs(x1), fabs(x2), EPSTOL)) return 0;

  // RGE in terms of natural log of renormalisation scale
  double from = log(fabs(x1));
  double to = log(fabs(x2));

  double guess = (from - to) * 0.1; //first step size
  double hmin = (from - to) * tol * 1.0e-5; 

  int err =
    integrateOdes(v, from, to, tol, guess, hmin, derivs, odeStepper);
  
  setMu(x2);
  return err;
}
 

 
