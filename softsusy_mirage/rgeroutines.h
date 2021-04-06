
/** \file rgeroutines.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: Does the calling of iteration routines etc, NOT for public 
                release!

*/

#ifndef RGEROUTINES_H
#define RGEROUTINES_H

#include <iostream>
#include <utils.h>
#include <softsusy.h>
#include <softpars.h>
#include <susy.h>
#include <lowe.h>
#include <linalg.h>
#include <def.h>
#include <string.h>
#include <rpvsusypars.h>
#include <rpvsoft.h>
using namespace std;

void gaussErrors();
void spsPoints(const QedQcd & oneset);
void scanPars2(double starttb, double endtb, double startM12, double endM12, 
               double m0, int numPoints, int sgnMu, int numRoutine);
void scanPars(double startM0, double endM0, double startM12, double endM12, 
	      double tanb, int numPoints, int sgnMu, int numRoutine);
double findChiSqSugra(double m0, double m12, double a0, double tanb, 
		      int sgnMu, int numRoutine);
double chiSqP1(const MssmSoftsusy & r);
double chiSqP2(const MssmSoftsusy & r);
double chiSqP3(const MssmSoftsusy & r);
double chiSqP4(const MssmSoftsusy & r);
double chiSqP5(const MssmSoftsusy & r);
double chiSqP6(const MssmSoftsusy & r);
double chiSqSPS1a(const MssmSoftsusy & r);
double mhqMax(const MssmSoftsusy & r);
double mllMax(const MssmSoftsusy & r);
double mlqMax(const MssmSoftsusy & r);
double edgeRatio(const MssmSoftsusy & r);
double mhqMin(const MssmSoftsusy & r);
double mllqMin(const MssmSoftsusy & r);
double mttMax(const MssmSoftsusy & r);
double mll4max(const MssmSoftsusy & r);
double mlqMaxFar(const MssmSoftsusy & r);
double mllqMax(const MssmSoftsusy & r);
double mllbMin(const MssmSoftsusy & r);

void scaledParsSugra();

MssmSoftsusy pointDomination(double tanb, double
				   mgut, int sgnMu, int accuracy,
			       DoubleVector & ft, double m32,
				   double epsilon = 0.0);
void contDomination(double, int, int);
void m12Scan();
void translateNonUniGauginos(DoubleVector & pars, double m0, double m1, 
			     double m2, double m3, double a0);
void nonUniGauginos(MssmSoftsusy & m, const DoubleVector & inputParameters);
void contUniversalM0m12(int sgnMu, int accuracy);
void pointUniversalM0m12(double m0, double m12, double tanb, double a0,
			    double mgut, int sgnMu, int accuracy);
MssmSoftsusy doPointUniversalM0m12(double m0, double m12, double tanb,
					double a0, double mgut, int sgnMu,
					int accuracy, DoubleVector &);
double getChisqTbtau(const QedQcd &, double);
void cycleStandard(const QedQcd & one, const MssmSusy & run, double, int);
inline void translateSugra(DoubleVector & pars, double m0, double m12, double
			    a0);
void iterateBound(int maxTries, double tol, double & oldBound, 
		  double m0, 
		  double m12, double a0, double tanb, double &mgut, int sgnMu,
		  int & err, RpvCouplings couplingType, int i, int j, int k,
		  const QedQcd & oneset, double neutrinoBound, 
		  RpvSoftsusy & kw); 
void tachyonBound(int maxTries, double tol, double & oldStart, double & oldEnd,
		  bool & startIsOk, bool & endIsOk, double m0, 
		  double m12, double a0, double tanb, double &mgut, int sgnMu,
		  int & err, RpvCouplings couplingType, int i, int j, int k,
		  const QedQcd & oneset, double neutrinoBound, 
		  RpvSoftsusy & kw);

void unifyBounds(RpvCouplings couplingType, int i, int j, int k, double m0, 
		 double m12, double a0, double tanb, double & mgut, 
		 int sgnMu, double & tachBound, double & neutBound, 
		 double & tachBoundMZ, double & neutBoundMZ, 
		 const QedQcd & oneset, RpvSoftsusy & kw);

RpvSoftsusy rpvPoint(RpvCouplings couplingType, int i, int j, int k, 
		       double tanb, 
		       double lambda, double m12, int sgnMu, 
		     const QedQcd & oneset, double & mgut);

void rpvScan(RpvCouplings couplingType, int i, int j, int k, double m12Start, 
		 double m12End, double lStart, double lEnd, double tanb, 
	     int sgnMu, const QedQcd & oneset);

/** Main program for RPV SUGRA paper. Does: 
       - no-scale scans:
         -# R-parity conserving
         -# lambda_231 
         -# lambda''_323 
	 -# lambda'_333
       - Scanning bounds from tachyons/neutrinos:
         -# lambda'_333
	 -# lambda_231
       - total bounds on all L-violating operators
         -# all possibilities for mixing
*/
void organiseBounds();

// as rpvPoint, but including m0 and a0 as arguments
RpvSoftsusy rpvPointII(RpvCouplings couplingType, int i, int j, int k, 
		       double tanb, double lambda, double m12, double m0, 
		       double a0, int sgnMu, const QedQcd & oneset, 
		       double & mgut);

#endif
