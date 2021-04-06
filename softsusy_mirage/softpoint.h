
/** \file softpoint.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: main calling program: command line interface

   $Log: softpoint.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.14  2004/03/23 00:16:24  allanach
   Added option to set MGUT=mEWSB

   Revision 1.13  2004/03/21 20:43:05  allanach
   Added alternative electroweak symmetry breaking conditions. Added possibility
   of EWSB=SUSY breaking boundary condition scale in SUSY Les Houches Accord.

   Revision 1.12  2004/03/19 21:27:25  allanach
   Dispensed with global variable storing 1-loop value of tadpoles. Intermediate
   step is just to include them in physpars structure.

   Revision 1.11  2004/03/19 20:41:34  allanach
   Corrected alternative EWSB conditions

   Revision 1.10  2004/03/17 13:05:23  allanach
   Added stau corrections to different REWSB condition

   Revision 1.9  2004/03/17 12:46:17  allanach
   Debugging different EWSB condition. 2-loop stuff commented better now,
   and getVev not used since it introduces the leading-log inaccuracy.

   Revision 1.8  2003/10/24 16:30:26  allanach
   Deleted old higgsVevMs variable, using running variable instead

   Revision 1.7  2003/08/19 14:26:22  allanach
   Changing lowOrg to be more sensible about gauge unification. Should now be
   called with POSITIVE mgut and a flag for gauge unification.

   Revision 1.6  2003/08/19 14:02:03  allanach
   Can now use "unified" to flag gauge unification

   Revision 1.5  2003/07/24 14:55:28  allanach
   Implemented les Houches input and output properly in the usual command-line
   interface

   Revision 1.4  2003/07/21 14:00:18  allanach
   MZ fully implemented as an input now. Kept MZ as the central PDG 2002 value,
   for defaults etc

   Revision 1.3  2003/05/22 12:45:41  allanach
   Fortran version of 2-loop Higgs corrections included with corrections:
   alpha_s (alpha_t + alpha_b) + alpha_t^2 + alpha_b^2 + alpha_t alpha_b

   Revision 1.2  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.1  2003/04/02 09:19:42  allanach
   Relegated MssmSoftsusy2 definition (and other headers) here

*/

#include <iostream>
#include <sstream>
#include <string>
using namespace std;
#include <mycomplex.h>
#include <def.h>
#include <linalg.h>
#include <lowe.h>
#include <rge.h>
#include <softsusy.h>
#include <softpars.h>
#include <susy.h>
#include <utils.h>
#include <numerics.h>
#include <twoloophiggs.h>
#include <cstring>

/// Does the user require gauge unification or not -- gaugeUnification changed
/// to be correct value
inline double mgutCheck(char * a, bool & gaugeUnification,
                        bool & ewsbBCscale) {
  gaugeUnification = false; ewsbBCscale = false;
  if (!strcmp(a, "?") || !strcmp(a,"unified")) {
    gaugeUnification = true;
    return 2.0e16;
  }
  if (!strcmp(a, "msusy")) {
    ewsbBCscale = true;
    return 1.0e3;
  }
  else return atof(a);
}

/// Incorrect input: gives advice on how to supply it
void errorCall();

/// Just like the usual MSSM SOFTSUSY class, but
/// enforces a different REWSB condition - given by mu(MS) and MA(pole)
/// conditions, instead of specifying mh1(MGUT) and mh2(MGUT)
class MssmSoftsusy2: public MssmSoftsusy {
private:
public:
  double mAcond, muCond; ///< user set conditions on mA and mu at M_SUSY
  /// Just sets the mu parameter to what the user wanted
  virtual void rewsbTreeLevel(int sgnMu) {
    setSusyMu(muCond);
  };
  /// "Inverts" the relation between SUSY breaking parameters and mA in order
  /// that mA(pole) is set by the user constraint: mH1^2 and mH2^2 are set in
  /// order to be consistent with mu and mA given.
  virtual void rewsb(int sgnMu, double mt) {

    setSusyMu(muCond);
    calcDrBarPars();
    double sinthDRbarMS = calcSinthdrbar();
    double tanb = displayTanb(), beta = atan(tanb);
    double mzRun = displayMzRun();
    doTadpoles(mt, sinthDRbarMS);

    double piaa = piAA(mAcond, displayMu());

    double mAsq =  sqr(displayDrBarPars().mhiggs(2));

    double gstrong = displayGaugeCoupling(3),
      rmtsq = sqr(displayDrBarPars().mt), scalesq = sqr(displayMu()),
      vev2 = sqr(displayHvev()), tbeta = displayTanb(),
      amu = -displaySusyMu(), mg = displayGaugino()(3);

    double p2s = 0., p2w = 0., p2b = 0., p2tau = 0.;
    if (numHiggsMassLoops > 1) {
      // two-loop Higgs corrections
      double sintau = sin(displayDrBarPars().thetatau),
        costau = cos(displayDrBarPars().thetatau);
      double msnusq = sqr(displayDrBarPars().msnu(3));
      double sxb = sin(displayDrBarPars().thetab),
        cxb = cos(displayDrBarPars().thetab);
      double msb1sq = sqr(displayDrBarPars().md(1, 3)),
        msb2sq = sqr(displayDrBarPars().md(2, 3));
      double mstau1sq = sqr(displayDrBarPars().me(1, 3)),
        mstau2sq = sqr(displayDrBarPars().me(2, 3));
      double cotbeta = 1.0 / tbeta;
      double rmbsq = sqr(displayDrBarPars().mb);
      double rmtausq = sqr(displayDrBarPars().mtau);
      double sxt = sin(displayDrBarPars().thetat),
        cxt = cos(displayDrBarPars().thetat);
      double mst1sq = sqr(displayDrBarPars().mu(1, 3)),
        mst2sq = sqr(displayDrBarPars().mu(2, 3));

      dszodd_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
              &tbeta, &vev2, &gstrong, &p2s);
      ddsodd_(&rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2,
              &p2w);
      dszodd_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
              &cotbeta, &vev2, &gstrong, &p2b);
      tausqodd_(&rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
                &costau, &scalesq, &amu, &tanb, &vev2, &p2tau);
    }

    double dMA = p2s + p2b + p2w + p2tau;

    double newMh1sq, newMh2sq;
    newMh1sq = sqr(sin(beta)) * (sqr(mAcond) + piaa + sqr(mzRun) - dMA)
      - (sqr(displaySusyMu()) + 0.5 * sqr(mzRun)) +
      displayTadpole1Ms() - sqr(sqr(sin(beta))) *
      displayPhys().t1OV1Ms1loop -
      displayPhys().t2OV2Ms1loop * sqr(sin(beta)) * sqr(cos(beta));

    newMh2sq = sqr(cos(beta)) * (sqr(mAcond) + piaa + sqr(mzRun) - dMA)
      - (sqr(displaySusyMu()) + 0.5 * sqr(mzRun)) -
      displayPhys().t1OV1Ms1loop * sqr(sin(beta)) * sqr(cos(beta)) +
      displayTadpole2Ms() - sqr(sqr(cos(beta))) *
      displayPhys().t2OV2Ms1loop;

    setMh1Squared(newMh1sq);
    setMh2Squared(newMh2sq);

    double m3sqnew;
    if (rewsbM3sq(displaySusyMu(), m3sqnew) == 0) flagB(false);
    else flagB(true);
    setM3Squared(m3sqnew);

    if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
       displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
    flagHiggsufb(true);
  else
    flagHiggsufb(false);

  }
};


/// This boundary condition does *not* set the parameters mH1^2 and mH2^2,
/// since they are set in order to be consistent with a certain mu and MA.
void extendedSugraBcs2(MssmSoftsusy & m,
                       const DoubleVector & inputParameters);
