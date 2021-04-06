
/** \file rge.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: RGE objects consisting of energy scale and parameters 
                and loops (order in perturbation theory) and thresholds
		(accuracy of calculation, handles decoupling ususally).

   $Log: rge.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.3  2005/06/09 10:14:04  allanach
   Deleted general boundary conditions -- they were confusing. Also,

   Revision 1.2  2005/05/30 14:22:14  allanach
   Fixed stau mixing in ISAWIG interface to 7.64

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.8  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.7  2003/12/18 16:32:05  allanach
   Cleaned up and corrected assumed accuracy of running

   Revision 1.6  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.5  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.4  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, EPS out of def.h to make it quicker to compile once
   they are changed.

   Revision 1.3  2002/04/13 15:51:01  allanach
   Starting to implement anomalous dimensions and beta functions

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef RGE_H
#define RGE_H

#include "def.h"
#include "utils.h"
class DoubleVector; 
#include "linalg.h"
#include "numerics.h"

/// Describes a set of parameters and RGEs in general. 
class RGE
{
private:
  double mu; ///< Renormalisation scale
  int numPars; ///< Number of parameters
  int loops; ///< To what order does the RG evolution run
  int thresholds; ///< Threshold flag: not always used

public:
  /// Default constructor fills data with zeroes
  RGE(): mu(0.0), numPars(0), loops(0), thresholds(0) {};
  
  /// Sets renormalisation scale to e
  void setMu(double e) { mu = e; };
  /// Default is to just set mu, but if you've got a derived class, 
  /// you may like to call several setMu's.
  virtual void setmu(double f) { setMu(f); };
  /// Set number of parameters in RGE object
  void setPars(int i) { numPars = i; };
  /// Set number of loops used
  void setLoops(int l) { loops = l; };
  /// Set level of threshold approximation used. It's typically used in the
  /// beta functions or boundary conditions.
  void setThresholds(int t) { thresholds = t; };

  /// Return renomalisation scale
  double displayMu() const { return mu; };
  /// Return number of loops
  int displayLoops() const { return loops; };
  /// Return number of parameters
  int howMany() const { return numPars; };
  /// Return level of threshold approximation
  int displayThresholds() const { return thresholds; };  
  /// Displays all RGE parameters in a double vector. Obligatory function for
  /// derived classes.
  virtual DoubleVector display() const = 0;
  /// Sets all parameters in an RGE object to the contents of s. It's an
  /// obligatory function for a derived class.
  virtual void set(const DoubleVector & s) = 0;
  /// Beta functions - all RGE objects must have them.
  virtual DoubleVector beta() const = 0;

  /// Runs parameters from scale "from" to "to, eps is fractional accuracy
  virtual int run(double from, double to, double eps= -1.0);
  /// Runs parameters to scale "to", eps is accuracy
  virtual int runto(double to, double eps = -1.0);
  
  /// Does the actual calling of Runge Kutta: default precision is TOLERANCE.
  /// Returns >0 if there's a problem with the running.
  /// IO parameters: x1=start renormalisation scale, x2=end renormalisation
  /// scale, v=boundary condition on parameters to be run on input, but gets
  /// changed to calculated final values on output, function to
  /// calculate derivatives are encoded in *derivs, eps is the numerical
  /// accuracy 
  int callRK(double x1, double x2, DoubleVector & v, 
	     DoubleVector (*derivs)(double, const DoubleVector &), 
	     double eps = -1.0);

};  

/// Runge-Kutte user defined routine: given log renorm scale x and the
/// dependent variables y of an RGE, will calculate the derivitives dydx.
DoubleVector allDerivs(double, const DoubleVector &);

#endif






