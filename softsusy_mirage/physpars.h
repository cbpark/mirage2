
/** \file physpars.h
   - Project:     SOFTSUSY 
   - File:        physpars.h
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: Pole masses of sparticles and Higgs', mixing parameters

   $Log: physpars.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.2  2005/07/26 10:59:55  allanach
   Made comments more accurate

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.14  2004/03/19 21:03:04  allanach
   Added 1-loop tadpole part

   Revision 1.13  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.12  2003/08/12 15:18:28  allanach
   Bug-fixed problem flags - some were forgotten

   Revision 1.11  2003/07/22 09:16:04  allanach
   Added running MW, MZ to definition of drbar parameters

   Revision 1.10  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.9  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.8  2003/02/24 14:26:39  allanach
   Implementing DRbar parameters in loop corrections. Half way though:
   drbarpars now changes to MPZ notation, but need to get rid of pole masses
   in loop corrections (and re-calculate tadpoles).

   Revision 1.7  2003/02/21 17:59:36  allanach
   Added drbar parameter class and calculation, starting to move to DRbar
   parameters in the 1-loop corrections

   Revision 1.6  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.4  2002/10/22 13:12:10  allanach
   Introduced new problem flag for infra-red quasi fixed points

   Revision 1.3  2002/08/30 12:31:55  allanach
   Test RPV version that works!

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef PHYSPARSH
#define PHYSPARSH

#include <iostream>
using std::ostream;
using std::istream;
using std::endl;
#include "linalg.h"

/// Masses of the physical particles. We have made assumption that scalar
/// mixing is small except in third generation and Higgs sector
struct sPhysical
{
  /// in order: h^0, A^0, H^0, H^+-
  DoubleVector mhiggs;
  /// sneutrino masses
  DoubleVector msnu;
  /// chargino/neutralino masses: ordered
  DoubleVector mch, mneut;
  /// Gluino mass
  double mGluino;
  /// neutralino mixing
  DoubleMatrix mixNeut;
  /// chargino and third family mixing angles
  double thetaL, thetaR, thetat, thetab, thetatau;
  /// sparticle masses in order (i=L/R, family)
  DoubleMatrix mu, md, me;
  /// Higgs mixing angle (alpha)
  double thetaH;
  
  /// DRbar tadpoles evaluated at MSusy
  double t1OV1Ms, t2OV2Ms, t1OV1Ms1loop, t2OV2Ms1loop;
  
  /// Feynman rules
  DoubleMatrix aChi0ChicW, bChi0ChicW;
  
  sPhysical(); ///< Constructor: initialises with zeroes
  sPhysical(const sPhysical &); ///< Constructor copies another object
  /// Displays the object in situations where a const is required
  sPhysical displaysPhysical() const { return *this; };
  /// Sets whole contents to those of another object s
  void setsPhysical(const sPhysical &s) { *this = s; };
  /// Sets whole contents to those of another object s
  const sPhysical & operator = (const sPhysical & s);
  
  /// Displays contents in a C-style convention *a (starts with index zero)
  void display(double *a) const;
}; 

/// Formatted printout
ostream & operator <<(ostream &, const sPhysical &); 

/// Various boolean values to flag any problems with the parameter point
struct sProblem {
  bool badConvergence; ///< Nowhere near a decent solution
  bool irqfp; ///< Infra-red quasi fixed point breached
  bool noRhoConvergence; ///< Couldn't calculate electroweak rho parameter
  bool noConvergence; ///< Iteration did not converge: not always serious
  bool tachyon; ///< Tachyonic point
  bool muSqWrongSign; ///< mu^2 came out with wrong sign; no REWSB
  bool b; ///< b came out with wrong sign; no REWSB
  bool higgsUfb; ///< Higgs potential inconsistent with a good minimum
  bool nonperturbative; ///< Running went non-perturbative
  bool noMuConvergence; ///< mu couldn't be calculated 
  bool test() const 
  {return (irqfp || noConvergence || tachyon || muSqWrongSign ||
	   higgsUfb || nonperturbative || noRhoConvergence || 
	   noMuConvergence || b || badConvergence);}; ///< returns true if there's any problem
  /// Only returns true if there's a serious problem
  bool testSeriousProblem() const 
  {return (irqfp || tachyon || muSqWrongSign || higgsUfb || nonperturbative 
	   || noRhoConvergence || noMuConvergence || b || badConvergence);}; 

  inline sProblem(); ///< constructor full of false values
  /// Constructor that sets flags equal to those of s
  inline sProblem(const sProblem & s);
  /// Sets flags equal to those of s
  const sProblem & operator = (const sProblem &);
  
};
/// Formatted output, but won't print unflagged problems
ostream & operator <<(ostream &st, const sProblem & p);
/// Formatted input of physical parameters
istream & operator >>(istream & left, sPhysical &s);

/// DRbar values of masses and mixings in MSSM
struct drBarPars: public sPhysical { 
  double mz, mw;       /// Running electroweak gauge boson masses
  double mt, mb, mtau; /// Running top, bottom and tau mass
  /// BPMZ convention mixing matrices for neutralinos and charginos
  ComplexMatrix nBpmz, uBpmz, vBpmz; 
  /// positive definite masses for neutralinos and charginos
  DoubleVector mnBpmz, mchBpmz; 

  inline drBarPars(); ///< Initialises with zero values
  inline drBarPars(const drBarPars &); ///< Initialises with another object
  /// Sets contents equal to those of another object
  const drBarPars & operator = (const drBarPars &s);

  /// Returns mixing matrix o and neutralino masses mn in the MPZ convention
  /// (hep-ph/9606211), n is 4 by 4 and mneut is 1->4, ie
  /// Calculates mnBpmz, nBpmz
  void mpzNeutralinos(); 
  /// Calculates uBpmz, vBpmz, mchBpmz, ie
  /// Returns mixing matrices u,v and neutralino masses mneut in the MPZ
  /// convention (hep-ph/9606211),  u+v are (2,2) and mch is 1->2.
  void mpzCharginos();
};

ostream & operator <<(ostream &, const drBarPars &); 

// ---------------------- inline class members ------------------------
inline drBarPars::drBarPars()
  : sPhysical(), mz(0.0), mw(0.0), 
    mt(0.0), mb(0.0), mtau(0.0), nBpmz(4, 4), uBpmz(2, 2), 
    vBpmz(2, 2), mnBpmz(4), mchBpmz(2)
{}

inline drBarPars::drBarPars(const drBarPars &s)
  : sPhysical(s.displaysPhysical()), mz(s.mz), mw(s.mw),
    mt(s.mt), mb(s.mb), mtau(s.mtau), 
    nBpmz(s.nBpmz), uBpmz(s.uBpmz), vBpmz(s.vBpmz), mnBpmz(s.mnBpmz), 
    mchBpmz(s.mchBpmz)
{}

inline sPhysical::sPhysical()
  : mhiggs(4), msnu(3), mch(2), mneut(4), mGluino(0.0),
    mixNeut(4, 4), thetaL(0.0), thetaR(0.0), thetat(0.0), thetab(0.0),
    thetatau(0.0), mu(2, 3), md(2, 3), me(2, 3), thetaH(0.0), 
    t1OV1Ms(0.0), t2OV2Ms(0.0), t1OV1Ms1loop(0.), t2OV2Ms1loop(0.), 
    aChi0ChicW(4, 2), bChi0ChicW(4, 2)
{}

inline sPhysical::sPhysical(const sPhysical & s)
  : mhiggs(s.mhiggs), msnu(s.msnu), mch(s.mch), 
    mneut(s.mneut), mGluino(s.mGluino), mixNeut(s.mixNeut), thetaL(s.thetaL),
    thetaR(s.thetaR), thetat(s.thetat), thetab(s.thetab),
    thetatau(s.thetatau), mu(s.mu), md(s.md), me(s.me), thetaH(s.thetaH),
    t1OV1Ms(s.t1OV1Ms), t2OV2Ms(s.t2OV2Ms),    
    t1OV1Ms1loop(s.t1OV1Ms1loop), t2OV2Ms1loop(s.t2OV2Ms1loop),
    aChi0ChicW(s.aChi0ChicW), bChi0ChicW(s.bChi0ChicW)
{}

inline sProblem::sProblem()
  : badConvergence(false), 
    irqfp(false), noRhoConvergence(false), noConvergence(false),
    tachyon(false), muSqWrongSign(false), b(false), higgsUfb(false), 
    nonperturbative(false), noMuConvergence(false)
{}

inline sProblem::sProblem(const sProblem & s)
  : badConvergence(s.badConvergence), 
    irqfp(s.irqfp), noRhoConvergence(s.noRhoConvergence), 
    noConvergence(s.noConvergence),
    tachyon(s.tachyon), muSqWrongSign(s.muSqWrongSign), b(s.b),
    higgsUfb(s.higgsUfb), nonperturbative(s.nonperturbative), 
    noMuConvergence(s.noMuConvergence)
{}

#endif
