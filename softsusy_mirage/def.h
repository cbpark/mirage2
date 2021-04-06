
/** \file def.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: Contains switches and default parameters
   */

#ifndef DEF_H
#define DEF_H

#include <cmath>
const char VERSION[] = "2.0.11";

/// uncomment if you want checking of vector/matrices bounds: slows code down
#define ARRAY_BOUNDS_CHECKING 

/// Make true if you want to include the 2-loop RGE corrections to scalar mass
/// squared parameters and trilinear terms: they slow it down by a factor of
/// 3. Note that gaugino and Higgs mass parameters are evolved to 2-loops by
/// default anyway.
extern bool INCLUDE_2_LOOP_SCALAR_CORRECTIONS;
/// Set to number of loops to use for calculation of Higgs mass 
/// (currently up to 2, the default)
extern int numHiggsMassLoops;
/// Set to number of loops to use for REWSB condition up to the default of 2
extern int numRewsbLoops;

const double EPSTOL = 1.0e-11; ///< underflow accuracy
const double PI = atan(1.0) * 4.0; ///< or 3.141592653589793 longhand;
const double root2 = sqrt(2.0);

extern double GMU;

/// particle data book 2004 central value. Is just used for intialisation etc
const double MW = 80.410; 
/// particle data book 2004 central value. Is just used for intialisation etc
const double MZ = 91.1876; 
/// variable for level of output and amount of quark: 0-3, higher numbers
/// giving more diagnostics. Set by user in file "massIn"
extern int PRINTOUT;
/// quark mixing flag: set by user in file "massIn":
/// 0=no quark mixing, 1=in up sector, 2=in down sector, -1=3rd family
/// approximation (all at MZ)
extern int MIXING; 
/// overall accuracy required: set by user in file "massIn"
extern double TOLERANCE;
/// SUSY breaking scale - if set by user
extern double QEWSB;

#endif
