
/** 
   Project:     SOFTSUSY 2.0
   File:        main.cpp
   Author:      Ben Allanach 
   Manual:      B.C. Allanach,hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   Description: main calling program example: performs a scan of tan beta 
   (starting at mSUGRA point SPS1a) and prints out Higgs masses as a result in
   the format
       <tan beta>      <mh0>     <mA0>    <mH0>     <mH+/->
*/

#include <iostream>
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"
#include <math.h>

using std::cout;
using std::endl;

/// global variable declaration
/// no quark mixing (dominant third family approx), and no verbose output
int MIXING = -1, PRINTOUT = 0;
/// fractional accuracy required
double TOLERANCE = 1.0e-3;
/// decay constant of muon
double GMU = 1.16637e-5; 
/// there are two possible conventions: if QEWSB > MZ, its value is assumed
/// in GeV and used as a constant MSUSY. Otherwise, it MULTIPLIES the usual 
/// MSUSY value, of root(mstop1 mstop2)
double QEWSB = 1.0; 
/// Do we include 2-loop RGEs of *all* scalar masses and A-terms, or only the
/// scalar mass Higgs parameters? (Other quantities all 2-loop anyway): the
/// default in SOFTSUSY 2.x is to include all 2-loop terms
bool INCLUDE_2_LOOP_SCALAR_CORRECTIONS = true;
/// number of loops used to calculate Higgs mass and tadpoles. They should be
/// identical for a consistent calculation
int numHiggsMassLoops = 2, numRewsbLoops = 2;
/// end of global variable declaration

int main() {
 /// Sets format of output: 6 decimal places
  outputCharacteristics(6);
  int numPoints = 1;
  double qMax = 0.;
	
  cerr << "SOFTSUSY" << VERSION << " test program, Ben Allanach 2002\n";
  cerr << "If you use SOFTSUSY, please refer to B.C. Allanach, \n";
  cerr << " Comput. Phys. Commun. 143 (2002) 305, hep-ph/0104145\n\n";

  /// Parameters used: Mirage Mediation parameters
  double alphac = 1. ,M0=500., am=0.5, ah=0.0, cm=0.5, ch=0.0 ,mGutGuess=2.0e16,
         tanb = 10.0;
  int sgnMu = 1;
  
  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
///  double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
  double alphasMZ = 0.1187, mtop = 172.4, mbmb = 4.2,mtau=1.777;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMass(mBottom, mbmb);
  oneset.setMass(mTau,mtau);
  oneset.toMz();      ///< Runs SM fermion masses to MZ

  
    MssmSoftsusy r; 
    DoubleVector pars(6);
    pars(1) = alphac; pars(2) = M0; pars(3) = am;
    pars(4) = ah    ; pars(5) = cm; pars(6) = ch;
    bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
    bool altEwsb = true;
//    bool ewsbBCscal = true;   

    /// Calculate the spectrum
    /// For user_define BCS
    
    r.lowOrg(userDefinedBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni); 
    
//    r.lesHouchesAccordOutput(nonUniversal, pars, sgnMu, tanb, qMax,
//		                            numPoints, mbmb, mtau, mgut, altEwsb);
    r.printLong();


    
/*    /// check the point in question is problem free: if so print the output
//    if (!r.displayProblem().test()) 
 //     cout << tanb << " " << r.displayPhys().mhiggs(1) << " " 
//	   << r.displayPhys().mhiggs(2) << " " 
	   << r.displayPhys().mhiggs(3) << " " 
	   << r.displayPhys().mhiggs(4) << endl;
    else
      /// print out what the problem(s) is(are)
      cout << tanb << " " << r.displayProblem() << endl;
  }*/
    
}


