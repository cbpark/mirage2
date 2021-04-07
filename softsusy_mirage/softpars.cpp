
/** \file softpars.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: softpars.cpp,v $
   Revision 1.5  2006/04/11 13:57:40  allanach
   Better comments in main programs and cleaned up bug in QEWSB non-usage

   Revision 1.5  2006/02/15 18:00:01  allanach
   Bug-fixed 2-loop mER^2 running: thanks to R Ruiz

   Revision 1.4  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.3  2005/05/30 14:22:14  allanach
   Fixed stau mixing in ISAWIG interface to 7.64

   Revision 1.2  2004/12/23 15:29:20  allanach
   Promoted INCLUDE_2_LOOP_SCALAR_CORRECTIONS to a global variable (in
   preparation for its control in the SUSY Les Houches Accord)

   Revision 1.40  2004/05/23 19:17:05  allanach
   Added some quicker coding in the 2-loop beta functions (3rd fam approx)

   Revision 1.39  2004/05/21 16:23:23  allanach
   Fixed bug that prevented compilation of the 2-loop scalar mass/trilinear
   running.

   Revision 1.38  2004/05/19 11:47:33  allanach
   Debugged 2-loop hd beta function

   Revision 1.37  2004/04/21 09:29:25  allanach
   Fixed bug in mQL^2 in GMSB BCs.

   Revision 1.36  2004/01/28 14:05:15  allanach
   Added special cases of x=0,1 to GMSB boundary conditions in order to avoid
   infinities.

   Revision 1.35  2004/01/15 13:54:55  allanach
   New heaer style implemented

   Revision 1.34  2004/01/09 18:47:26  allanach
   Made faster and more efficient

   Revision 1.33  2003/10/27 19:20:32  allanach
   Bug-fix in 2-loop gaugino mass beta function

   Revision 1.32  2003/10/24 16:09:04  allanach
   Implemented running Higgs DRbar vev

   Revision 1.30  2003/08/28 14:40:46  allanach
   Small bug-fix in Higgs mass squared running (two-loop term).

   Revision 1.29  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.28  2003/07/16 11:54:16  allanach
   Added 2-loop RGE terms to higgs mass squared parameters as the default
   (in the dominant third family approximation)

   Revision 1.27  2003/06/25 17:58:40  allanach
   Put in 2-loop higgs mass corrections (in third family approximation) as a
   default in softpars

   Revision 1.25  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.24  2003/05/27 15:05:52  allanach
   Purely efficiency corrections: used variable rather than display() methods
   whenever possible

   Revision 1.23  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.22  2003/02/21 17:59:36  allanach
   Added drbar parameter class and calculation, starting to move to DRbar
   parameters in the 1-loop corrections

   Revision 1.21  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.20  2002/12/19 18:35:20  allanach
   Put missing "double" statements in for two-loop option

   Revision 1.18  2002/11/07 20:10:23  allanach
   Bug-fixed GMSB BCs

   Revision 1.17  2002/09/09 10:42:54  allanach
   TOLERANCE replaces EPS as being more user-friendly

   Revision 1.16  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, TOLERANCE out of def.h to make it quicker to
   compile once they are changed.

   Revision 1.15  2002/08/30 12:31:55  allanach
   Test RPV version that works!

   Revision 1.14  2002/07/30 12:57:31  allanach
   SOFTSUSY1.5

   Revision 1.13  2002/07/03 14:45:41  allanach
   2-loop scalar beta-functions completed

   Revision 1.12  2002/07/01 16:09:10  allanach
   Added 2-loop hd, hu beta functions

   Revision 1.11  2002/06/14 16:26:30  allanach
   Switches included for 2-loop running of scalar masses, and calulating mt
   at mt.

   Revision 1.9  2002/04/18 14:32:05  allanach
   Changed RGEs and anomalous dimensions to be compatible with new notation;
   started implementation of rewsb in R-parity violation

   Revision 1.5  2002/04/14 13:50:41  allanach
   Now use V=m3^2 H1 H2 instead of V=mu B H1 H2. It's more natural!

   Revision 1.4  2002/04/12 16:51:27  allanach
   Added display/set functions to work automatically

   Revision 1.3  2002/04/12 06:24:50  allanach
   Code maintenance - returning a subobject made simpler

   Revision 1.5  2001/10/22 09:19:28  allanach
   Fixed missing fermion contribution in piZZT. Turned some integers to
   decimals for performance on more platforms.

   Revision 1.4  2001/10/04 19:26:34  allanach
   New version deals with AMSB correctly

   Revision 1.3  2001/09/28 14:49:02  allanach
   GMSB boundary conditions added

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#include "softpars.h"
#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

const SoftParsMssm & SoftParsMssm::operator=(const SoftParsMssm & s) {
  if (this == &s) return *this;
  mGaugino = s.mGaugino;
  ua = s.ua;
  da = s.da;
  ea = s.ea;
  m32 = s.m32;
  mQLsq = s.mQLsq;
  mURsq = s.mURsq;
  mDRsq = s.mDRsq;
  mLLsq = s.mLLsq;
  mSEsq = s.mSEsq;
  m3sq = s.m3sq;
  mH1sq = s.mH1sq;
  mH2sq = s.mH2sq;
  setSusy(s.displaySusy());
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  return *this;
}

DoubleMatrix SoftParsMssm::displayTrilinear(trilinears k) const {
  switch(k) {
  case UA: return ua; break;
  case DA: return da; break;
  case EA: return ea; break;
  default:
    ostringstream ii;
    ii << "In SoftParsMssm::displayTrilinear, called with illegal argument";
    ii << " " << k << endl;
    throw ii.str();
    break;
  }
}

double SoftParsMssm::displayTrilinear(trilinears k, int i, int j)
  const {
  switch(k) {
  case UA: return ua.display(i, j); break;
  case DA: return da.display(i, j); break;
  case EA: return ea.display(i, j); break;
  default:
    ostringstream ii;
    ii << "In SoftParsMssm::displayTrilinear, called with illegal argument";
    ii << " " << k << endl;
    throw ii.str();
    break;
  }
}

// Will display zero if that's what the A term is regardless of what the
// Yukawa coupling is (even if it's zero).
double SoftParsMssm::displaySoftA(trilinears k, int i, int j) const {
  double am = displayTrilinear(k, i, j);
  if (fabs(am) < TOLERANCE) return 0.0;

  if (fabs(displayYukawaElement(yukawa(k), i, j)) < 1.0e-10) {
    ostringstream ii;
    ii << "WARNING: asking for SoftParsMssm::displaySoftA(" << int(k) << ","
         << i << "," << j << "), where Yukawa coupling is " <<
      fabs(displayYukawaElement(yukawa(k), i, j)) << endl;
    throw ii.str();
  }
  double temp = displayTrilinear(k, i, j)
    / displayYukawaElement(yukawa(k), i, j);
  return temp;
}

DoubleMatrix SoftParsMssm::displaySoftMassSquared(softMasses k) const {
  switch(k) {
  case mQl: return mQLsq; break;
  case mUr: return mURsq; break;
  case mDr: return mDRsq; break;
  case mLl: return mLLsq; break;
  case mEr: return mSEsq; break;
  default:
    ostringstream ii;
    ii << "SoftParsMssm::displaySoftMassSquared with illegal argument ";
    ii << k << endl;
    throw ii.str();
    break;
  }
}

double SoftParsMssm::displaySoftMassSquared(softMasses k, int i, int j)
  const {
  switch(k) {
  case mQl: return mQLsq.display(i, j); break;
  case mUr: return mURsq.display(i, j); break;
  case mDr: return mDRsq.display(i, j); break;
  case mLl: return mLLsq.display(i, j); break;
  case mEr: return mSEsq.display(i, j); break;
  default:
    ostringstream ii;
    ii << "SoftParsMssm::displaySoftMassSquared with illegal argument ";
    ii << k << endl;
    throw ii.str();
    break;
  }
}

DoubleVector SoftParsMssm::display() const {
  DoubleVector y(MssmSusy::display());
  y.setEnd(numSoftParsMssm);
  int i, j, k=numSusyPars;
  for (i=1; i<=3; i++) {
    k++;
    y(k) = displayGaugino(i);
  }

  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      k++;
      y(k) = displayTrilinear(UA, i, j);
      y(k+9) = displayTrilinear(DA, i, j);
      y(k+18) = displayTrilinear(EA, i, j);
      y(k+27) = displaySoftMassSquared(mQl, i, j);
      y(k+36) = displaySoftMassSquared(mUr, i, j);
      y(k+45) = displaySoftMassSquared(mDr, i, j);
      y(k+54) = displaySoftMassSquared(mLl, i, j);
      y(k+63) = displaySoftMassSquared(mEr, i, j);
    }
  y(k+64) = m3sq;
  y(k+65) = mH1sq;
  y(k+66) = mH2sq;
  return y;
}

void SoftParsMssm::setSoftMassElement(softMasses k, int i, int j,
                                          double f) {
  switch(k) {
  case mQl: mQLsq(i, j) = f; break;
  case mUr: mURsq(i, j) = f; break;
  case mDr: mDRsq(i, j) = f; break;
  case mLl: mLLsq(i, j) = f; break;
  case mEr: mSEsq(i, j) = f; break;
  }
}

void SoftParsMssm::setSoftMassMatrix(softMasses k, const DoubleMatrix & m) {
  switch(k) {
  case mQl: mQLsq = m; break;
  case mUr: mURsq = m; break;
  case mDr: mDRsq = m; break;
  case mLl: mLLsq = m; break;
  case mEr: mSEsq = m; break;
  }
}

void SoftParsMssm::setTrilinearMatrix(trilinears k, const DoubleMatrix & m) {
  switch(k) {
  case UA: ua = m; break;
  case DA: da = m; break;
  case EA: ea = m; break;
  }
}

void SoftParsMssm::setTrilinearElement(trilinears k, int i, int j,
                                          double m) {
  switch(k) {
  case UA: ua(i, j) = m; break;
  case DA: da(i, j) = m; break;
  case EA: ea(i, j) = m; break;
  }
}

void SoftParsMssm::setAllGauginos(const DoubleVector & v) {
  if (v.displayStart() != 1 || v.displayEnd() !=3) {
    ostringstream ii;
    ii << "Initialising SoftParsMssm::setAllGauginos with vector"
         << v;
    throw ii.str();
  }
  mGaugino = v;
}


void SoftParsMssm::set(const DoubleVector & y) {
  MssmSusy::set(y);
  int i, j, k=numSusyPars;
  for (i=1; i<=3; i++) {
    k++;
    mGaugino(i) = y.display(k);
  }

  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      k++;
      ua(i, j) = y.display(k);
      da(i, j) = y.display(k+9);
      ea(i, j) = y.display(k+18);
      mQLsq(i, j) = y.display(k+27);
      mURsq(i, j) = y.display(k+36);
      mDRsq(i, j) = y.display(k+45);
      mLLsq(i, j) = y.display(k+54);
      mSEsq(i, j) = y.display(k+63);
    }
  m3sq = y.display(k+64);
  mH1sq = y.display(k+65);
  mH2sq = y.display(k+66);
}


// Outputs derivatives (DRbar scheme) in the form of dsoft
// thresholds = 0 and NOTHING is decoupled.
// thresholds = 2 and SUSY params/gauge couplings are decoupled at sparticle
// thresholds.
// CHECKED: 24/05/02
SoftParsMssm SoftParsMssm::beta2() const {
  //  susypars beta functions are defined in the
  //  Allanach, Dedes, Dreiner hep-ph/9906209 paper notation
  //  W=  LL Y^E H1 ER + QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1,
  //  (implicit epsilonAb phi^a phi^b a,b=Su(2) indices present where
  //  epsilon12 = +1, epsilon21 = -1, other components zero
  //  The following equations are for our conventions
  //  VTri = QL ua ur H2 + QL da dr H1 + LL ea Er H1 + hc
  //  VBi  = mH1sq |H1|^2 + mH2sq |H2|^2 + Q^* mQLsq Q + L^* mLLsq L +
  //          UR MuR UR^* + DR MDR DR^* + ER MER ER^* + m3sq H2 H1
  // Note in particular that ua(3,3) is NOT AT. AT = ua(3, 3) / ht.

  // Constants for gauge running
  static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3);
  static DoubleMatrix babBeta(3, 3);

  if (bBeta(1) < 1.0e-5) // Constants not set yet
    setBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);

  // For calculational brevity:
  static sBrevity a;
  // convert to beta functions
  static MssmSusy dsb;

  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = MssmSusy::beta(a);

  // To keep this a const function: TIME SAVINGS
  DoubleMatrix &u1=a.u1, &d1=a.d1, &e1=a.e1;
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2,
    &u2t=a.u2t, &d2t=a.d2t, &e2t=a.e2t, &dt=a.dt, &ut=a.ut, &et=a.et;
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq;

  // Best to make these const DoubleMatrix & type
  const DoubleMatrix &hu = displayTrilinear(UA), &hd = displayTrilinear(DA),
    &he = displayTrilinear(EA);

  static DoubleMatrix hut(3, 3), hdt(3, 3), het(3, 3);
  static DoubleMatrix hu2(3, 3), hd2(3, 3), he2(3, 3);
  static DoubleMatrix hu2t(3, 3), hd2t(3, 3), he2t(3, 3);
  const DoubleMatrix &mq = displaySoftMassSquared(mQl),
    &mu = displaySoftMassSquared(mUr), &md = displaySoftMassSquared(mDr),
    &ml = displaySoftMassSquared(mLl), &me = displaySoftMassSquared(mEr);
  static DoubleVector mG(displayGaugino());
  static DoubleVector msq(1, 3), gsqM(1, 3), gMsq(1, 3);

  hut = hu.transpose(); hdt = hd.transpose(); het = he.transpose();
  hu2 = hu * hut; hd2 = hd * hdt; he2 = he * het;
  double huT = hu2.trace(), hdT = hd2.trace(), heT = he2.trace();
  hu2t = hut * hu; hd2t = hdt * hd; he2t = het * he;
  double mqT = mq.trace(), muT = mu.trace(), mdT = md.trace(), meT =
    me.trace(), mlT = ml.trace();

  double huuT = (hu * ut).trace(), hddT = (hd * dt).trace(),
    heeT = (he * et).trace();
  mG = mGaugino; msq = mG * mG; gsqM = gsq * mG; gMsq = gsq * msq;

  double m3sq = displayM3Squared();

  // derivatives of soft parameters
  static DoubleVector dmG(1, 3);
  static double dmH1sq, dmH2sq, dm3sq;
  static DoubleMatrix dmq(3, 3), dmu(3, 3), dmd(3, 3), dme(3, 3), dml(3, 3);
  static DoubleMatrix dhu(3, 3), dhd(3, 3), dhe(3, 3);

  const double ONEO16Pisq = 1.0 / (16.0 * sqr(PI));
  if (displayLoops() > 0) {
    const static double sixteenO3 = 16.0 / 3.0, oneO15 = 1.0 /
      15.0;

    dmG = gsqM * 2.0 * bBeta;

    double curlyS = mH2sq - mH1sq + mqT - mlT - 2.0 * muT + mdT + meT;

    dm3sq = 2.0 * displaySusyMu() * (0.6 * gsqM(1) + 3.0 * gsqM(2) +
                        3.0 * huuT  +  3.0 * hddT + heeT)
      + m3sq * (3.0 * uuT + 3.0 * ddT + eeT - 0.6 * gsq(1) - 3.0 * gsq(2));

    dmH1sq = 2.0 *
      (-0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) +
       3.0 * (mH1sq * ddT + (d2 * mq + d2t * md).trace() + hdT) +
       (mH1sq * eeT + (e2 * ml + e2t * me).trace() + heT));

    dmH2sq = 2.0 *
      (0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) +
       3.0 * (mH2sq * uuT + (u2 * mq + u2t * mu).trace() + huT));

    dmq = 2.0 *
      (0.1 * gsq(1) * curlyS - oneO15 * gMsq(1) - 3.0 * gMsq(2) -
       sixteenO3 * gMsq(3) +
       0.5 * (u2 * mq + mq * u2 + 2.0 * (u1 * mu * ut + mH2sq * u2 + hu2))
       + 0.5 * (d2 * mq + mq * d2 + 2.0 * (d1 * md * dt + mH1sq * d2 +
                                           hd2)));

    dml = 2.0 *
      ( -0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) +
        0.5 * (e2 * ml + ml * e2 + 2.0 * (e1 * me * et + mH1sq * e2 + he2)));

    dmu = 2.0 *
      (-0.4 * gsq(1) * curlyS - 16.0 * oneO15 * gMsq(1) - sixteenO3 *
       gMsq(3) +
       (u2t * mu + mu * u2t + 2.0 * (ut * mq * u1 + mH2sq * u2t + hu2t)));

    dmd = 2.0 * (0.2 * gsq(1) * curlyS
                 - 4.0 * oneO15 * gMsq(1) - sixteenO3 * gMsq(3) +
                 (d2t * md + md * d2t + 2.0 * (dt * mq * d1 + mH1sq * d2t +
                                               hd2t)));

    dme = 2.0 *
      (0.6 * gsq(1) * curlyS - 36.0 * oneO15 * gMsq(1) +
       e2t * me + me * e2t + 2.0 * (et * ml * e1 + mH1sq * e2t + he2t));

    dhu = (3.0 * uuT
           - sixteenO3 * gsq(3) - 3.0 * gsq(2) - 13.0 * oneO15 * gsq(1))
      * hu
      + 2.0 * (3.0 * huuT + 13.0 * oneO15 * gsqM(1) + 3 * gsqM(2) + sixteenO3
             * gsqM(3)) * u1 +
      4.0 * hu * u2t + 5.0 * u2 * hu + 2.0 * hd * dt * u1 + d2 * hu;

    dhd = (eeT + 3.0 * ddT -sixteenO3 * gsq(3) - 3.0 * gsq(2) -
           7.0 * oneO15 * gsq(1)) * hd
      + 2.0 * (7.0 * oneO15 * gsqM(1) + 3.0 * gsqM(2) + sixteenO3
             * gsqM(3) + heeT + 3.0 * hddT) * d1 +
      4.0 * hd * d2t + 5.0 * d2 * hd +
      2.0 * hu * ut * d1 + u2 * hd;

    dhe = (eeT + 3.0 * ddT - 3.0 * gsq(2) - 1.8 * gsq(1)) * he
      + (3.6 * gsqM(1) + 6.0 * gsqM(2) + 2.0 * heeT + 6.0 * hddT) * e1 +
      4.0 * he * e2t + 5.0 * e2 * he;

    // convert to proper derivatives:
    dmG    = ONEO16Pisq * dmG;
    dm3sq  *= ONEO16Pisq;
    dmH1sq *= ONEO16Pisq;
    dmH2sq *= ONEO16Pisq;
    dmq    *= ONEO16Pisq;
    dml    *= ONEO16Pisq;
    dmu    *= ONEO16Pisq;
    dmd    *= ONEO16Pisq;
    dme    *= ONEO16Pisq;
    dhu    *= ONEO16Pisq;
    dhd    *= ONEO16Pisq;
    dhe    *= ONEO16Pisq;
  }

  // two-loop contributions. I got these from hep-ph/9311340. WIth respect to
  // their notation, the Yukawas and h's are TRANSPOSED. Gaugino masses are
  // identical, as are soft masses. B(mV) = mu(BBO) B(BBO)
  if (displayLoops() > 1) {
    static DoubleVector dmG2(1, 3);
    static DoubleMatrix dmq2(3, 3), dmu2(3, 3), dmd2(3, 3), dme2(3, 3),
      dml2(3, 3);
    static DoubleMatrix dhu2(3, 3), dhd2(3, 3), dhe2(3, 3);
    static DoubleVector sigma(3);

    const double oneO16Pif = sqr(ONEO16Pisq);

    DoubleVector &g4=a.g4;

    double mq3 = displaySoftMassSquared(mQl, 3, 3);
    double mu3 = displaySoftMassSquared(mUr, 3, 3);
    double md3 = displaySoftMassSquared(mDr, 3, 3);
    double ml3 = displaySoftMassSquared(mLl, 3, 3);
    double me3 = displaySoftMassSquared(mEr, 3, 3);
    double ht = displayYukawaElement(YU, 3, 3), ht2 = sqr(ht), ht4 = sqr(ht2);
    double htau = displayYukawaElement(YE, 3, 3), htau2 = sqr(htau),
      htau4 = sqr(htau2);
    double hb = displayYukawaElement(YD, 3, 3), hb2 = sqr(hb), hb4 = sqr(hb2);
    double Ut = displayTrilinear(UA, 3, 3), Ut2 = sqr(Ut);
    double Ub = displayTrilinear(DA, 3, 3), Ub2 = sqr(Ub);
    double Utau = displayTrilinear(EA, 3, 3), Utau2 = sqr(Utau);

    double sP =
      -(3.0 * mH2sq + mq3 - 4.0 * mu3) * ht2 +
       (3.0 * mH1sq - mq3 - 2.0 * md3) * hb2 +
       (mH1sq + ml3 - 2.0 * me3) * htau2 +
      (1.5 * gsq(2) + 0.3 * gsq(1)) * (mH2sq - mH1sq - mlT) +
       (8.0 / 3.0 * gsq(3) + 1.5 * gsq(2) + gsq(1) / 30.0) * mqT -
       (16.0 / 3.0 * gsq(3) + 16.0 / 15.0 * gsq(1)) * muT +
       (8.0 / 3.0 * gsq(3) + 2.0 / 15.0 * gsq(1)) * mdT +
      1.2 * gsq(1) * meT;  //checked

    sigma(1) = 0.2 * gsq(1) *
      (3.0 * (mH1sq + mH2sq) + mqT + 3.0 * mlT + 8.0 * muT + 2.0 * mdT + 6.0 *
       meT);
    sigma(2) = gsq(2) * (mH1sq + mH2sq + (3.0 * mqT + mlT));
    sigma(3) = gsq(3) * (2.0 * mqT + muT + mdT); // these 3 checked

    if (INCLUDE_2_LOOP_SCALAR_CORRECTIONS) {
      // There are plenty of time savings that could be made in this section
      // if one were to use dominant third family solely in the 2-loop
      // corrections.
    /* old full 3-family result (slow)
    double sP =
      (-(3.0 * mH2sq + mq) * u2 + 4.0 * u1 * mu * ut +
       (3.0 * mH1sq - mq) * d2 - 2.0 * d1 * md * dt + (mH1sq + ml) * e2 -
       2.0 * e1 * me * et).trace() +
      (1.5 * gsq(2) + 0.3 * gsq(1)) * (mH2sq - mH1sq - mlT) +
       (8.0 / 3.0 * gsq(3) + 1.5 * gsq(2) + gsq(1) / 30.0) * mqT -
       (16.0 / 3.0 * gsq(3) + 16.0 / 15.0 * gsq(1)) * muT +
       (8.0 / 3.0 * gsq(3) + 2.0 / 15.0 * gsq(1)) * mdT +
      1.2 * gsq(1) * meT;  //checked

    double dmH2sq2 = -6.0 *
      (6.0 * (mH2sq + mq) * u2 * u2 + 6.0 * u1 * mu * u2t * ut +
       (mH1sq + mH2sq + mq) * u2 * d2 + u1 * mu * ut * d2 +
       u2 * mq * d2 + u2 * d1 * md * dt + 6.0 * hu2 * u2 +
       6.0 * hu * u2t * hut + hd2 * u2 + d2 * hu2 +
       hd * dt * u1 * hut + d1 * hdt * hu * ut).trace() +
      (32.0 * gsq(3) + 1.6 * gsq(1)) *
      ((mH2sq + mq) * u2 + u1 * mu * ut + hu2).trace() +
      32.0 * gsq(3) * (2.0 * msq(3) * uuT - 2.0 * mG(3) * huuT) +
      1.6 * gsq(1) * (2.0 * msq(1) * uuT - 2.0 * mG(1) * huuT) +
      1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) +
      18.0 / 5.0 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
      621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) +
      0.6 * gsq(1) * sigma(1); // checked

    double dmH1sq2 = -6.0 *
      (6.0 * (mH1sq + mq) * d2 * d2 + 6.0 * dt * md * d2 * d1 +
       (mH1sq + mH2sq + mq) * u2 * d2 + u1 * mu * ut * d2 +
       u2 * mq * d2 + u2 * d1 * md * dt + 2.0 * (mH1sq + ml) * e2 * e2 +
       2.0 * e1 * me * e2t * et + 6.0 * hd2 * d2 + 6.0 * hd * d2t * hdt +
       hu2 * d2 + u2 * hd2 + hu * ut * d1 * hdt + u1 * hut * hd * dt +
       2.0 * he2 * e2 + 2.0 * he * e2t * het).trace() +
      (32.0 * gsq(3) - 0.8 * gsq(1)) * ((mH1sq + mq) * d2 + d1 * md * dt +
                                        hd2).trace() +
      32.0 * gsq(3) * (2.0 * msq(3) * ddT - 2.0 * mG(3) * hddT) -
      0.8 * gsq(1) * (2.0 * msq(1) * ddT - 2.0 * mG(1) * hddT) +
      2.4 * gsq(1) * (((mH1sq + ml) * e2 + e1 * me * et + he2).trace() +
      2.0 * msq(1) * eeT - 2.0 * mG(1) * heeT) - 1.2 * gsq(1) * sP +
      33.0 * g4(2) * msq(2) +
      3.6 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
      621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) *
      sigma(1); // checked+corrected 27/9/05

    double dm3sq2 = m3sq *
      (-3.0 * (3.0 * u2t * u2t + 3.0 * d2t * d2t + 2.0 * ut * d2
               * u1 + e2t * e2t).trace() +
       (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
       (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT +
       1.2 * gsq(1) * eeT + 7.5 * g4(2) + 1.8 * gsq(1) * gsq(2) +
       207. / 50. * g4(1)) +
      displaySusyMu() *
      (-12.0 * (3.0 * hut * u2 * u1 + 3.0 * hdt * d2 * d1 + hut * d2 * u1 +
                hdt * u2 * d1 + het * e2 * e1).trace() +
       (32.0 * gsq(3) + 1.9 * gsq(1)) * huuT +
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT -
       (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
       2.4 * gsq(1) * mG(1) * eeT - 30.0 * g4(2) * mG(2) - 18.0 / 5.0 * gsq(1)
       * gsq(2) * (mG(1) + mG(2)) - 414. / 25.0 * g4(1) * mG(1)); // checked
    */

    /* Old full 3-family version
    dmq2 = -(2.0 * mq + 8.0 * mH2sq) * u2 * u2 - 4.0 * u1 * mu * u2t * ut -
      4.0 * u2 * mq * u2 - 4.0 * u2 * u1 * mu * ut - 2.0 * u2 * u2 * mq -
      (2.0 * mq + 8.0 * mH1sq) * d2 * d2 - 4.0 * d1 * md * d2t * dt -
      4.0 * d2 * mq * d2 - 4.0 * d2 * d1 * md * dt - 2.0 * d2 * d2 * mq -
      ((mq + 4.0 * mH2sq) * u2 + 2.0 * u1 * mu * ut + u2 * mq) * uuT * 3.0 -
      ((mq + 4.0 * mH1sq) * d2 + 2.0 * d1 * md * dt + d2 * mq) * (3.0 * ddT +
                                                                  eeT) -
      6.0 * u2 * (mq * u2 + u1 * mu * ut).trace() -
      d2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 + 2.0 * e1 * me
            * et).trace() -
      4.0 * (u2 * hu2 + hu2 * u2 + u1 * hu2t * ut + hu * u2t * hut) -
      4.0 * (d2 * hd2 + hd2 * d2 + d1 * hd2t * dt + hd * d2t * hdt) -
      hu2 * uuT * 6.0 - u2 * 6.0 * (hu * hut).trace() - hu * ut * 6.0 * huuT -
      u1 * hut * 6.0 * huuT - hd2 * (6.0 * ddT + 2.0 * eeT) -
      d2 * (6.0 * hd2 + 2.0 * he2).trace() -
      hd * dt * (6.0 * hddT + 2.0 * heeT) -
      d1 * hdt * (6.0 * hddT + 2.0 * heeT) +
      0.4 * gsq(1) * ((2.0 * mq +  4.0 * mH2sq) * u2 + 4.0 * u1 * mu * ut +
                      2.0 * u2 * mq + 4.0 * hu2 - 4.0 * mG(1) * hu * ut - 4.0
                      * mG(1) * u1 * hut + 8.0 * msq(1) * u2 +
                      (mq + 2.0 * mH1sq) * d2 + 2.0 * d1 * md * dt +
                      d2 * mq + 2.0 * hd2 - 2.0 * mG(1) * hd * dt -
                      2.0 * mG(1) * d1 * hdt + 4.0 * msq(1) * d2) +
      0.4 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 32.0 * gsq(3) *
      gsq(2) * (msq(3) + msq(2) + mG(2) * mG(3)) +
      32.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + mG(1) * mG(3)) + 33.0
      * g4(2) * msq(2) + 0.4 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) *
                                                  mG(2)) +
      199.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
      3.0 * gsq(2) * sigma(2) + gsq(1) * sigma(1) / 15.0; // checked
    */

    /* new dominant 3-family version */
    dmq2 = -(2.0 * mq + 8.0 * mH2sq) * u2 * u2 - 4.0 * u1 * mu * u2t * ut -
      4.0 * u2 * mq * u2 - 4.0 * u2 * u1 * mu * ut - 2.0 * u2 * u2 * mq -
      (2.0 * mq + 8.0 * mH1sq) * d2 * d2 - 4.0 * d1 * md * d2t * dt -
      4.0 * d2 * mq * d2 - 4.0 * d2 * d1 * md * dt - 2.0 * d2 * d2 * mq -
      ((mq + 4.0 * mH2sq) * u2 + 2.0 * u1 * mu * ut + u2 * mq) * uuT * 3.0 -
      ((mq + 4.0 * mH1sq) * d2 + 2.0 * d1 * md * dt + d2 * mq) * (3.0 * ddT +
                                                                  eeT) -
      6.0 * u2 * (mq * u2 + u1 * mu * ut).trace() -
      d2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 + 2.0 * e1 * me
            * et).trace() -
      4.0 * (u2 * hu2 + hu2 * u2 + u1 * hu2t * ut + hu * u2t * hut) -
      4.0 * (d2 * hd2 + hd2 * d2 + d1 * hd2t * dt + hd * d2t * hdt) -
      hu2 * uuT * 6.0 - u2 * 6.0 * Ut2 - hu * ut * 6.0 * huuT -
      u1 * hut * 6.0 * huuT - hd2 * (6.0 * ddT + 2.0 * eeT) -
      d2 * (6.0 * Ub2 + 2.0 * Utau2) -
      hd * dt * (6.0 * hddT + 2.0 * heeT) -
      d1 * hdt * (6.0 * hddT + 2.0 * heeT) +
      0.4 * gsq(1) * ((2.0 * mq +  4.0 * mH2sq) * u2 + 4.0 * u1 * mu * ut +
                      2.0 * u2 * mq + 4.0 * hu2 - 4.0 * mG(1) * hu * ut - 4.0
                      * mG(1) * u1 * hut + 8.0 * msq(1) * u2 +
                      (mq + 2.0 * mH1sq) * d2 + 2.0 * d1 * md * dt +
                      d2 * mq + 2.0 * hd2 - 2.0 * mG(1) * hd * dt -
                      2.0 * mG(1) * d1 * hdt + 4.0 * msq(1) * d2) +
      0.4 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 32.0 * gsq(3) *
      gsq(2) * (msq(3) + msq(2) + mG(2) * mG(3)) +
      32.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + mG(1) * mG(3)) + 33.0
      * g4(2) * msq(2) + 0.4 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) *
                                                  mG(2)) +
      199.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
      3.0 * gsq(2) * sigma(2) + gsq(1) * sigma(1) / 15.0; // checked

    /*
    dml2 = -(2.0 * ml + 8.0 * mH1sq) * e2 * e2 - 4.0 * e1 * me * e2t * et -
      4.0 * e2 * ml * e2 - 4.0 * e2 * e1 * me * et - 2.0 * e2 * e2 * ml -
      ((ml + 4.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml) *
      (3.0 * ddT + eeT) -
      e2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 +
            2.0 * e1 * me * et).trace() -
      4.0 * (e2 * he2 + he2 * e2 + e1 * he2t * et + he * e2t * het) -
      he2 * (6.0 * ddT + 2.0 * eeT) - e2 * (6.0 * hd2 + 2.0 * he2).trace() -
      he * et * (6.0 * hddT + 2.0 * heeT) -
      e1 * het * (6.0 * hddT + 2.0 * heeT) +
      1.2 * gsq(2) * ((ml + 2.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml +
                      2.0 * he2 - 2.0 * mG(1) * he * et - 2.0 * mG(1) * e1 *
                      het + 4.0 * msq(1) * e2) -
      1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 3.6 * gsq(2) * gsq(1) *
      (msq(2) + msq(1) + mG(1) * mG(2)) + 621.0 / 25.0 * g4(1) * msq(1) +
      3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) * sigma(1); //checked
    */

    dml2 = -(2.0 * ml + 8.0 * mH1sq) * e2 * e2 - 4.0 * e1 * me * e2t * et -
      4.0 * e2 * ml * e2 - 4.0 * e2 * e1 * me * et - 2.0 * e2 * e2 * ml -
      ((ml + 4.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml) *
      (3.0 * ddT + eeT) -
      e2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 +
            2.0 * e1 * me * et).trace() -
      4.0 * (e2 * he2 + he2 * e2 + e1 * he2t * et + he * e2t * het) -
      he2 * (6.0 * ddT + 2.0 * eeT) - e2 * (6.0 * Ub2 + 2.0 * Utau2) -
      he * et * (6.0 * hddT + 2.0 * heeT) -
      e1 * het * (6.0 * hddT + 2.0 * heeT) +
      1.2 * gsq(2) * ((ml + 2.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml +
                      2.0 * he2 - 2.0 * mG(1) * he * et - 2.0 * mG(1) * e1 *
                      het + 4.0 * msq(1) * e2) -
      1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 3.6 * gsq(2) * gsq(1) *
      (msq(2) + msq(1) + mG(1) * mG(2)) + 621.0 / 25.0 * g4(1) * msq(1) +
      3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) * sigma(1); //checked

    /*
    dmu2 = -(2.0 * mu + 8.0 * mH2sq) * u2t * u2t - 3.0 * ut * mq * u2 * u1 -
      4.0 * u2t * mu * u2t - 4.0 * u2t * ut * mq * u1 -
      2.0 * u2t * u2t * mu - (2.0 * mu + 4.0 * mH2sq + 4.0 * mH1sq) * ut * d2
      * u1 - 4.0 * ut * mq * d2 * u1 - 4.0 * ut * d1 * md * dt * u1 -
      4.0 * ut * d2 * mq * u1 - 2.0 * ut * d2 * u1 * mu -
      ((mu + 4.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu) * 6.0 * uuT -
      12.0 * u2t * (mq * u2 + u1 * mu * ut).trace() -
      4.0 * (hu2t * u2t + u2t * hu2t + hut * u2 * hu + ut * hu2 * u1) -
      4.0 * (hut * hd * dt * u1 + ut * d1 * hdt * hu + hut * d2 * hu +
             ut * hd2 * u1) -
      12.0 * (hu2t * uuT + u2t * (hu * hut).trace() + hut * u1 * huuT +
              ut * hu * huuT) +
      (6.0 * gsq(2) - 0.4 * gsq(1)) *
      ((mu + 2.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu + 2.0 * hu2t) +
      12.0 * gsq(2) * (2.0 * msq(2) * u2t - mG(2) * hut * u1 -
                       mG(2) * ut * hu) -
      0.8 * gsq(1) * (2.0 * msq(1) * u2t - mG(1) * hut * u1 - mG(1) * ut * hu)-
      1.6 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) +
      512.0 / 45.0 * gsq(1) * gsq(3) * (msq(3) + msq(1) + mG(1) * mG(3)) +
      3424.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
      16.0 / 15.0 * gsq(1) * sigma(1); // checked
    */

    dmu2 = -(2.0 * mu + 8.0 * mH2sq) * u2t * u2t - 3.0 * ut * mq * u2 * u1 -
      4.0 * u2t * mu * u2t - 4.0 * u2t * ut * mq * u1 -
      2.0 * u2t * u2t * mu - (2.0 * mu + 4.0 * mH2sq + 4.0 * mH1sq) * ut * d2
      * u1 - 4.0 * ut * mq * d2 * u1 - 4.0 * ut * d1 * md * dt * u1 -
      4.0 * ut * d2 * mq * u1 - 2.0 * ut * d2 * u1 * mu -
      ((mu + 4.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu) * 6.0 * uuT -
      12.0 * u2t * (mq * u2 + u1 * mu * ut).trace() -
      4.0 * (hu2t * u2t + u2t * hu2t + hut * u2 * hu + ut * hu2 * u1) -
      4.0 * (hut * hd * dt * u1 + ut * d1 * hdt * hu + hut * d2 * hu +
             ut * hd2 * u1) -
      12.0 * (hu2t * uuT + u2t * Ut2 + hut * u1 * huuT +
              ut * hu * huuT) +
      (6.0 * gsq(2) - 0.4 * gsq(1)) *
      ((mu + 2.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu + 2.0 * hu2t) +
      12.0 * gsq(2) * (2.0 * msq(2) * u2t - mG(2) * hut * u1 -
                       mG(2) * ut * hu) -
      0.8 * gsq(1) * (2.0 * msq(1) * u2t - mG(1) * hut * u1 - mG(1) * ut * hu)-
      1.6 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) +
      512.0 / 45.0 * gsq(1) * gsq(3) * (msq(3) + msq(1) + mG(1) * mG(3)) +
      3424.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
      16.0 / 15.0 * gsq(1) * sigma(1); // checked

    dmd2 = -(2.0 * md + 8.0 * mH1sq) * d2t * d2t - 4.0 * dt * mq * d1 * d2t -
      4.0 * d2t * md * d2t - 4.0 * d2t * dt * mq * d1 - 2.0 * d2t * d2t * md -
      (2.0 * md + 4.0 * mH2sq + 4.0 * mH1sq) * dt * u2 * d1 -
      4.0 * dt * mq * u2 * d1 - 4.0 * dt * u1 * mu * ut * d1 -
      4.0 * dt * u2 * mq * d1 - 2.0 * dt * u2 * d1 * md -
      ((md + 4.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md) *
      (6.0 * ddT + 2.0 * eeT) -
      4.0 * d2t * (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2
                   + e1 * me * et) -
      4.0 * (hd2t * d2t + d2t * hd2t + hdt * d2 * hd + dt * hd2 * d1) -
      4.0 * (hdt * hu * ut * d1 + dt * u1 * hut * hd + hdt * u2 * hd +
             dt * hu2 * d1) -
      4.0 * hd2t * (3.0 * ddT + eeT) - 4.0 * d2t * (3.0 * Ub2 + Utau2) -
      4.0 * hdt * d1 * (3.0 * hddT + heeT) -
      4.0 * dt * hd * (3.0 * hddT + heeT) +
      (6.0 * gsq(2) + 0.4 * gsq(1)) *
      ((md + 2.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md + 2.0 * hd2t) +
      12.0 * gsq(2) * (msq(2) * d2t - mG(2) * hdt * d1 - mG(2) * dt * hd) +
      0.8 * gsq(1) * (2.0 * msq(1) * d2t - mG(1) * hdt * d1 - mG(1) * dt * hd)+
      0.8 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(2) +
      128.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + mG(1) * mG(3)) +
      808.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
      4.0 / 15.0 * gsq(1) * sigma(1); // checked

    dme2 = -(2.0 * me + 8.0 * mH1sq) * e2t * e2t - 4.0 * et * ml * e2 * e1 -
      4.0 * e2t * me * e2t - 4.0 * e2t * et * ml * e1 - 2.0 * e2t * e2t * me -
      ((me + 4.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me) *
      (6.0 * ddT + 2.0 * eeT) - 4.0 * e2t *
      (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 + e1 * me * et).trace() -
      4.0 * (he2t * e2t + e2t * he2t + het * e2 * he + et * he2 * e1) -
      4.0 * he2t * (3.0 * ddT + eeT) - 4.0 * e2t * (3.0 * Ub2 + Utau2) -
      4.0 * het * e1 * (3.0 * hddT + heeT) - 4.0 * et * he *
      (3.0 * hddT + heeT) + (6.0 * gsq(2) - 1.2 * gsq(1)) *
      ((me + 2.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me + 2.0 * he2t) +
      12.0 * gsq(2) * (2.0 * msq(2) * e2t - mG(2) * het * e1
                       - mG(2) * et * he) -
      2.4 * gsq(1) * (2.0 * msq(1) * e2t - mG(1) * het * e1 - mG(1) * et * he)
      + 2.4 * gsq(1) * sP + 2808. / 25.0 * g4(1) * msq(1) +
      2.4 * gsq(1) * sigma(1); // checked

    // old full 3-family version
    /*
      dhu2 =
      (-3.0 * (3.0 * u2 * u2 + ut * d2 * u1).trace() - d2 * (3.0 * ddT + eeT)
       -15.0 * u2 * uuT - 6.0 * u2 * u2 - 2.0 * d2 * d2 -
       4.0 * u2 * d2 + (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
       12.0 * gsq(2) * u2 + 0.4 * gsq(1) * d2 - 16.0 / 9.0 * g4(3) +
       8.0 * gsq(3) * gsq(2) + 136.0 / 45.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
       gsq(2) * gsq(1) + 2743.0 / 450.0 * g4(1)) * hu +
      (-6.0 * (6.0 * hut * u2 * u1 + hut * d2 * u1 + hdt * u2 * d1).trace() -
       18.0 * u2 * huuT - d2 * (6.0 * hddT + 2.0 * heeT) -
       12.0 * hu * ut * uuT - hd * dt * (6.0 * ddT + 2.0 * eeT) -
       6.0 * hu * u2t * ut - 8.0 * u2 * hu * ut - 4.0 * hd * d2t * dt -
       4.0 * d2 * hd * dt - 2.0 * hu * ut * d2 - 4.0 * u2 * hd * dt +
       (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT +
       (6.0 * gsq(2) + 1.2 * gsq(1)) * hu * ut + 0.8 * gsq(1) * hd * dt -
       (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
       (12.0 * gsq(2) * mG(2) + 0.8 * gsq(1) * mG(1)) * u2 -
       0.8 * gsq(1) * mG(1) * d2 + 64.0 / 9.0 * g4(3) * mG(3) -
       16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) -
       272.0 / 45.0 * gsq(3) * gsq(1) * (mG(1) + mG(3)) -
       30.0 * g4(2) * mG(2) - 2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) -
       5486.0 / 225.0 * g4(1) * mG(1)) * u1; // checked
    */

    /* Dominant 3rd family version */
    dhu2 =
      (-3.0 * (3.0 * ht4 + ht2 * hb2) - d2 * (3.0 * ddT + eeT)
       -15.0 * u2 * uuT - 6.0 * u2 * u2 - 2.0 * d2 * d2 -
       4.0 * u2 * d2 + (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
       12.0 * gsq(2) * u2 + 0.4 * gsq(1) * d2 - 16.0 / 9.0 * g4(3) +
       8.0 * gsq(3) * gsq(2) + 136.0 / 45.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
       gsq(2) * gsq(1) + 2743.0 / 450.0 * g4(1)) * hu +
      (-6.0 * (6.0 * Ut * ht2 * ht + Ut * hb2 * ht + Ub * ht2 * hb) -
       18.0 * u2 * huuT - d2 * (6.0 * hddT + 2.0 * heeT) -
       12.0 * hu * ut * uuT - hd * dt * (6.0 * ddT + 2.0 * eeT) -
       6.0 * hu * u2t * ut - 8.0 * u2 * hu * ut - 4.0 * hd * d2t * dt -
       4.0 * d2 * hd * dt - 2.0 * hu * ut * d2 - 4.0 * u2 * hd * dt +
       (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT +
       (6.0 * gsq(2) + 1.2 * gsq(1)) * hu * ut + 0.8 * gsq(1) * hd * dt -
       (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
       (12.0 * gsq(2) * mG(2) + 0.8 * gsq(1) * mG(1)) * u2 -
       0.8 * gsq(1) * mG(1) * d2 + 64.0 / 9.0 * g4(3) * mG(3) -
       16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) -
       272.0 / 45.0 * gsq(3) * gsq(1) * (mG(1) + mG(3)) -
       30.0 * g4(2) * mG(2) - 2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) -
       5486.0 / 225.0 * g4(1) * mG(1)) * u1; // checked

    /* Old full 3-family version
    dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
      (-3.0 * (3.0 * d2t * d2t + ut * d2 * u1 + e2t * e2t).trace() -
       3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d2 * d2 -
       2.0 * u2 * u2 - 4.0 * d2 * u2 + (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT +
       1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 +
       (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) +
       8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
       gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
      (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
               2.0 * het * e2 * e1).trace() - 6.0 * u2 * huuT -
       6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
       4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
       8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
       4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 +
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
       1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
       2.4 * gsq(1) * mG(1) * eeT -
       (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
       1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) -
       16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) -
       16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) - 30.0 * g4(2) * mG(2) -
       2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 574.0 / 45.0 * g4(1) * mG(1))
      * d1;
    */

    // New dominant 3rd family version
    dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
      (-3.0 * (3.0 * hb4 + ht2 * hb2 + htau4) -
       3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d2 * d2 -
       2.0 * u2 * u2 - 4.0 * d2 * u2 + (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT +
       1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 +
       (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) +
       8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
       gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
      (-6.0 * (6.0 * Ub * hb2 * hb + Ut * hb2 * ht + Ub * ht2 * hb +
               2.0 * Utau * htau2 * htau) - 6.0 * u2 * huuT -
       6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
       4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
       8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
       4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 +
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
       1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
       2.4 * gsq(1) * mG(1) * eeT -
       (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
       1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) -
       16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) -
       16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) - 30.0 * g4(2) * mG(2) -
       2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 574.0 / 45.0 * g4(1) * mG(1))
      * d1;

    /* old version */
    dhe2 =
      (-3.0 * (3.0 * d2t * d2t + ut * d2 * u1 + e2t * e2t).trace() -
       5.0 * e2 * (3.0 * ddT + eeT) - 6.0 * e2 * e2 +
       (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT +
       (12.0 * gsq(2) - 1.2 * gsq(1)) * eeT + 7.5 * g4(2) +
       1.8 * gsq(2) * gsq(1) + 13.5 * g4(1)) * he +
      (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
               2.0 * het * e2 * e1).trace() -
       4.0 * he * et * (3.0 * ddT + eeT) - 6.0 * e2 * (3.0 * hddT + heeT) -
       6.0 * he * e2t * et - 8.0 * e2 * he * et +
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT - 2.4 * gsq(1) * heeT +
       (6.0 * gsq(2) + 1.2 * gsq(1)) * he * et -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
       2.4 * gsq(1) * mG(1) * eeT - 12.0 * gsq(2) * mG(2) * e2 -
       30.0 * g4(2) * mG(2) - 3.6 * gsq(2) * gsq(1) * (mG(1) + mG(2)) -
       54.0 * g4(1) * mG(1)) * e1; // checked


    /* 3rd family version */
    dhe2 =
      (-3.0 * (3.0 * hb4 + ht2 * hb2 + htau4) -
       5.0 * e2 * (3.0 * ddT + eeT) - 6.0 * e2 * e2 +
       (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT +
       (12.0 * gsq(2) - 1.2 * gsq(1)) * eeT + 7.5 * g4(2) +
       1.8 * gsq(2) * gsq(1) + 13.5 * g4(1)) * he +
      (-6.0 * (6.0 * Ub * hb2 * hb + Ut * hb2 * ht + Ub * ht2 * hb +
               2.0 * Utau * htau2 * htau) -
       4.0 * he * et * (3.0 * ddT + eeT) - 6.0 * e2 * (3.0 * hddT + heeT) -
       6.0 * he * e2t * et - 8.0 * e2 * he * et +
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT - 2.4 * gsq(1) * heeT +
       (6.0 * gsq(2) + 1.2 * gsq(1)) * he * et -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
       2.4 * gsq(1) * mG(1) * eeT - 12.0 * gsq(2) * mG(2) * e2 -
       30.0 * g4(2) * mG(2) - 3.6 * gsq(2) * gsq(1) * (mG(1) + mG(2)) -
       54.0 * g4(1) * mG(1)) * e1; // checked

    // add onto one-loop beta functions
    dmq2 *= oneO16Pif; dmq += dmq2;
    dml2 *= oneO16Pif; dml += dml2;
    dmu2 *= oneO16Pif; dmu += dmu2;
    dmd2 *= oneO16Pif; dmd += dmd2;
    dme2 *= oneO16Pif; dme += dme2;
    dhu2 *= oneO16Pif; dhu += dhu2;
    dhd2 *= oneO16Pif; dhd += dhd2;
    dhe2 *= oneO16Pif; dhe += dhe2;
    }

    // Default is to include these 2-loop corrections anyhow because they can
    // be so important: gauginos + higgs
    dmG2 = 2.0 * gsq *
      (babBeta * gsqM + mG * (babBeta * gsq) +
       cuBeta * huuT - cuBeta * mG * uuT +
       cdBeta * hddT - cdBeta * mG * ddT +
       ceBeta * heeT - ceBeta * mG * eeT); // checked

    // The following are valid in the third-family approximation
    double dm3sq2 = m3sq *
      (-3.0 * (3.0 * ht4 + 3.0 * hb4 + 2.0 * hb2 * ht2 + htau4) +
       (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
       (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT +
       1.2 * gsq(1) * eeT + 7.5 * g4(2) + 1.8 * gsq(1) * gsq(2) +
       207. / 50. * g4(1)) +
      displaySusyMu() *
      (-12.0 * (3.0 * Ut * ht * ht2 + 3.0 * Ub * hb * hb2 + Ut * hb2 * ht +
                Ub * ht2 * hb + Utau * htau2 * htau) +
       (32.0 * gsq(3) + 1.9 * gsq(1)) * huuT +
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT -
       (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
       2.4 * gsq(1) * mG(1) * eeT - 30.0 * g4(2) * mG(2) - 18.0 / 5.0 * gsq(1)
       * gsq(2) * (mG(1) + mG(2)) - 414. / 25.0 * g4(1) * mG(1)); // checked

    double dmH2sq2 = -6.0 *
      (6.0 * (mH2sq + mq3 + mu3) * ht4 +
      hb2 * ht2 * (mH1sq + mH2sq + 2.0 * mq3 + mu3 + md3) +
       ht2 * 12.0 * Ut2 + ht2 * Ub2 + hb2 * Ut2 + 2.0 * hb * ht * Ut * Ub) +
      (32.0 * gsq(3) + 1.6 * gsq(1)) *
      ((mH2sq + mq3 + mu3) * ht2 + Ut2) +
      32.0 * gsq(3) * (2.0 * msq(3) * ht2 - 2.0 * mG(3) * ht * Ut) +
      1.6 * gsq(1) * (2.0 * msq(1) * ht2 - 2.0 * mG(1) * ht * Ut) +
      1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) +
      18.0 / 5.0 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
      621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) +
      0.6 * gsq(1) * sigma(1); // checked

    double dmH1sq2 = -6.0 *
      (6.0 * (mH1sq + mq3 + md3) * hb4 +
       (mH1sq + mH2sq + 2.0 * mq3 + mu3 + md3) * ht2 * hb2 +
       2.0 * (mH1sq + ml3 + me3) * htau4 +
       12.0 * Ub2 * hb2 + hb2 * Ut2 + ht2 * Ub2 + 2.0 * Ut * ht * Ub * hb +
       4.0 * htau2 * Utau2) +
      (32.0 * gsq(3) - 0.8 * gsq(1)) * ((mH1sq + mq3 + md3) * hb2 + Ub2) +
      32.0 * gsq(3) * (2.0 * msq(3) * ddT - 2.0 * mG(3) * hddT) -
      0.8 * gsq(1) * (2.0 * msq(1) * ddT - 2.0 * mG(1) * hddT) +
      2.4 * gsq(1) * (((mH1sq + ml3 + me3) * htau2 + Utau2) +
                      2.0 * msq(1) * eeT - 2.0 * mG(1) * heeT)
      - 1.2 * gsq(1) * sP +
      33.0 * g4(2) * msq(2) +
      3.6 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
      621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) *
      sigma(1); // checked

    dmG = dmG + dmG2 * oneO16Pif;
    dm3sq = dm3sq + dm3sq2 * oneO16Pif;
    dmH1sq = dmH1sq + dmH1sq2 * oneO16Pif;
    dmH2sq = dmH2sq + dmH2sq2 * oneO16Pif;
  }

  // Check for non-perturbativeness
  /*  const double nonPertLimit = 5.0;
  if (displayYukawaElement(YU, 3, 3) > nonPertLimit ||
      displayYukawaElement(YD, 3, 3)> nonPertLimit) {
    return SoftParsMssm();
    }*/

  SoftParsMssm dsoft(dsb, dmG, dhu, dhd, dhe, dmq, dmu,
                     dmd, dml, dme, dm3sq, dmH1sq, dmH2sq,
                     displayGravitino(), displayMu(),
                     displayLoops(), displayThresholds());

  return dsoft;
}

// Outputs derivatives vector y[109] for SUSY parameters: interfaces to
// integration routines
DoubleVector SoftParsMssm::beta() const {
  // calculate the derivatives
  static SoftParsMssm dsoft; dsoft = beta2();

  return dsoft.display(); // convert to a long vector
}

// Outputs derivatives of anomalous dimensions, from which the running can be
// derived.
void SoftParsMssm::anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
                                  DoubleMatrix & gQQ, DoubleMatrix & gUU,
                                  DoubleMatrix & gDD,
                                  double & gH1H1, double & gH2H2)  const {
  // Constants for gauge running
  static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3);
  static DoubleMatrix babBeta(3, 3);

  if (bBeta(1) < 1.0e-5) // Constants not set yet
    setBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);

  // For calculational brevity:
  static sBrevity a;
  // convert to beta functions
  static MssmSusy dsb;

  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = MssmSusy::beta(a);

  static DoubleVector g1(3);
  g1 = displayGauge();

  // To keep this a const function: TIME SAVINGS
  DoubleMatrix &dt=a.dt, &ut=a.ut, &et=a.et;
  DoubleVector &gsq=a.gsq; // &g3=a.g3, &g4=a.g4;

  static DoubleMatrix hu(3, 3), hd(3, 3), he(3, 3);
  static DoubleMatrix hut(3, 3), hdt(3, 3), het(3, 3);
  static DoubleVector gsqM(1, 3);

  hu = displayTrilinear(UA); hd = displayTrilinear(DA);
  he = displayTrilinear(EA);
  hut = hu.transpose(); hdt = hd.transpose(); het = he.transpose();

  double huuT = (hu * ut).trace(), hddT = (hd * dt).trace(),
    heeT = (he * et).trace();
  gsqM = gsq * mGaugino;

  if (displayLoops() > 0) { // CHECKED: agrees with our conventions
    const static double eightO3 = 8.0 / 3.0, oneO30 = 1.0 / 30.0;
    gLL = - (he * et + 0.3 * gsqM(1) + 1.5 * gsqM(2));
    gEE = - (2.0 * (et * he) + 1.2 * gsqM(1));
    gQQ = - (hd * dt + hu * ut + oneO30 * gsqM(1) + 1.5 * gsqM(2) +
             eightO3 * gsqM(3));
    gDD = - (2.0 * dt * hd + 4.0 * oneO30 * gsqM(1) +
             eightO3 * gsqM(3));
    gUU = - (2.0 * ut * hu + 16.0 * oneO30 * gsqM(1) +
             eightO3 * gsqM(3));
    gH1H1 = - (3.0 * hddT + heeT + 0.3 * gsqM(1) + 1.5 * gsqM(2));
    gH2H2 = - (3.0 * huuT + 0.3 * gsqM(1) + 1.5 * gsqM(2));

    const static double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
    gH1H1 = gH1H1 * oneO16Pisq;
    gH2H2 = gH2H2 * oneO16Pisq;
    gEE  *= oneO16Pisq;
    gLL  *= oneO16Pisq;
    gQQ  *= oneO16Pisq;
    gUU  *= oneO16Pisq;
    gDD  *= oneO16Pisq;
  }
}

// Gives the ytilde terms relevant for the soft mass running: CHECKED 23/5/02
void SoftParsMssm::yTildes(DoubleMatrix & yu, DoubleMatrix & yd, DoubleMatrix
                           &ye) const {
  ye = displaySoftMassSquared(mLl) * displayYukawaMatrix(YE) +
    displayYukawaMatrix(YE) * displayMh1Squared() +
    displayYukawaMatrix(YE) * displaySoftMassSquared(mEr);

  yd = displaySoftMassSquared(mQl) * displayYukawaMatrix(YD) +
    displayYukawaMatrix(YD) * displayMh1Squared() +
    displayYukawaMatrix(YD) * displaySoftMassSquared(mDr);

  yu = displaySoftMassSquared(mQl) * displayYukawaMatrix(YU) +
    displayYukawaMatrix(YU) * displayMh2Squared() +
    displayYukawaMatrix(YU) * displaySoftMassSquared(mUr);
}

/*
   Give it a SUSY object and a value of M3/2, and it will return a soft
   object with AMSB soft breaking terms. Note that the sleptons will be
   tachyonic, ie nothing has been done to fix that problem.
   Note that in the following, we are neglecting all Yukawa couplings except
   that of the third family.

   THE CURRENT STATE OF PLAY:
   trilinear and gauge interactions will generalise to two loops OK
   trilinear interactions have full flavour structure
   soft scalar masses are ONE LOOP ONLY, with 3rd family approximation and no
   particular flavour structure. However, should you require full flavour
   structure (for example for FCNC constraints), this should be relatively
   easy to incorporate. Two loop additions are possible, but a pain.
   */
void SoftParsMssm::addAmsb(double m32) {
  MssmSusy run(displaySusy());
  const double ONEO16pisq = 1.0 / (16. * sqr(PI));

  double
    yt   = run.displayYukawaElement(YU, 3, 3),
    yb   = run.displayYukawaElement(YD, 3, 3),
    ytau = run.displayYukawaElement(YE, 3, 3),
    g1   = run.displayGaugeCoupling(1),
    g2   = run.displayGaugeCoupling(2),
    g3   = run.displayGaugeCoupling(3);

  // gauge parts of beta functions
  double Xu = -13.0 / 15.0 * sqr(g1) - 3.0 * sqr(g2) - 16.0 / 3.0 * sqr(g3);
  double Xd = -7.0 / 15.0 * sqr(g1) - 3.0 * sqr(g2) - 16.0 / 3.0 * sqr(g3);
  double Xe = -9 / 5.0 * sqr(g1) - 3.0 * sqr(g2);

  // Yukawa parts
  double betayt   = yt *   (Xu + 6.0 * sqr(yt) + sqr(yb));
  double betayb   = yb *   (Xd + 6.0 * sqr(yb) + sqr(yt) + sqr(ytau));
  double betaytau = ytau * (Xe + 3.0 * sqr(yb) + 4.0 * sqr(ytau));

  const double ONEO16pif = sqr(ONEO16pisq);

  // Higgs
  double mhusq = ONEO16pif * sqr(m32) *
    (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)) + 3.0 * yt * betayt);
  double mhdsq = ONEO16pif * sqr(m32) *
    (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)) + 3.0 * yb * betayb + ytau
     * betaytau);

  // third generation
  double mstoprsq = ONEO16pif * sqr(m32) *
    (-88. / 25. * sqr(sqr(g1)) + 8.0 * sqr(sqr(g3)) + 2.0 * yt * betayt);
  double msbottomrsq = ONEO16pif * sqr(m32) *
    (-22. / 25. * sqr(sqr(g1)) + 8.0 * sqr(sqr(g3)) + 2.0 * yb * betayb);
  double mQl3sq = ONEO16pif * sqr(m32) *
    (-11. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)) + 8.0 * sqr(sqr(g3)) + yb *
     betayb + yt * betayt);

  double mL3sq = ONEO16pif * sqr(m32) *
    (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)) + ytau
     * betaytau);
  double mstaursq = ONEO16pif * sqr(m32) *
    (-198. / 25. * sqr(sqr(g1)) + 2.0 * ytau * betaytau);

  // other generations
  double mursq = ONEO16pif * sqr(m32) *
    (-88. / 25. * sqr(sqr(g1)) + 8.0 * sqr(sqr(g3)));
  double mdrsq = ONEO16pif * sqr(m32) *
    (-22. / 25. * sqr(sqr(g1)) + 8.0 * sqr(sqr(g3)));
  double mQlsq = ONEO16pif * sqr(m32) *
    (-11. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)) + 8.0 * sqr(sqr(g3)));
  double mLsq = ONEO16pif * sqr(m32) *
    (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)));
  double mersq = ONEO16pif * sqr(m32) *
    (-198. / 25. * sqr(sqr(g1)));

  mQLsq(1, 1) = mQLsq(1, 1) + mQlsq;
  mQLsq(2, 2) = mQLsq(2, 2) + mQlsq;
  mQLsq(3, 3) = mQLsq(3, 3) + mQl3sq;
  mURsq(1, 1) = mURsq(1, 1) + mursq;
  mURsq(2, 2) = mURsq(2, 2) + mursq;
  mURsq(3, 3) = mURsq(3, 3) + mstoprsq;
  mDRsq(1, 1) = mDRsq(1, 1) + mdrsq;
  mDRsq(2, 2) = mDRsq(2, 2) + mdrsq;
  mDRsq(3, 3) = mDRsq(3, 3) + msbottomrsq;
  mLLsq(1, 1) = mLLsq(1, 1) + mLsq;
  mLLsq(2, 2) = mLLsq(2, 2) + mLsq;
  mLLsq(3, 3) = mLLsq(3, 3) + mL3sq;
  mSEsq(1, 1) = mSEsq(1, 1) + mersq;
  mSEsq(2, 2) = mSEsq(2, 2) + mersq;
  mSEsq(3, 3) = mSEsq(3, 3) + mstaursq;

  mH1sq = mH1sq + mhdsq;
  mH2sq = mH2sq + mhusq;

  // For calculational brevity:
  static sBrevity a;
  static MssmSusy dsb;

  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = MssmSusy::beta(a);

  // AMSB spectrum
  DoubleVector amsbGaugino(3);
  DoubleMatrix temp(3, 3);

  int i; for (i=1; i<=3; i++)
    amsbGaugino(i) = dsb.displayGaugeCoupling(i) * m32 /
      run.displayGaugeCoupling(i);
  mGaugino = mGaugino + amsbGaugino;

  ua = ua - dsb.displayYukawaMatrix(YU) * m32;
  da = da - dsb.displayYukawaMatrix(YD) * m32;
  ea = ea - dsb.displayYukawaMatrix(YE) * m32;

  return;
}

// Reads in universal boundary conditions at the current scale:
// m0, M1/2, A0, B and sign of mu
void SoftParsMssm::universal(double m0,  double m12,  double a0,  double mu,
                              double m3sq) {
  standardSugra(m0, m12, a0);
  setSusyMu(mu);
  setM3Squared(m3sq);
}

void SoftParsMssm::universalScalars(double m0) {
  // scalar masses
  DoubleMatrix ID(3, 3), mm0(3, 3);
  int i; for (i=1; i<=3; i++) ID(i, i) = 1.0;
  mm0 = ID * sqr(m0);
  setSoftMassMatrix(mQl, mm0); setSoftMassMatrix(mUr, mm0);
  setSoftMassMatrix(mDr, mm0); setSoftMassMatrix(mLl, mm0);
  setSoftMassMatrix(mEr, mm0);
  setMh1Squared(sqr(m0)); setMh2Squared(sqr(m0));
}

void SoftParsMssm::universalGauginos(double m12) {
  // gaugino masses
  int i; for (i=1; i<=3; i++) setGauginoMass(i, m12);
}

void SoftParsMssm::universalTrilinears(double a0)  {
  // trilinears
  setTrilinearMatrix(UA, a0 * displayYukawaMatrix(YU));
  setTrilinearMatrix(DA, a0 * displayYukawaMatrix(YD));
  setTrilinearMatrix(EA, a0 * displayYukawaMatrix(YE));
}

// Input m0, NOT m0 squared.
void SoftParsMssm::standardSugra(double m0,  double m12, double a0) {
  if (m0 < 0.0) {
    ostringstream ii;
    ii << "m0=" << m0 << " passed to universal boundary" <<
      "conditions illegally negative.";
    throw ii.str();
  }
  universalScalars(m0);
  universalGauginos(m12);
  universalTrilinears(a0);
}

#define HR "---------------------------------------------------------------\n"

ostream & operator <<(ostream &left, const SoftParsMssm &s) {
  left << "SUSY breaking MSSM parameters at Q: " << s.displayMu() << endl;
  left << " UA" << s.displayTrilinear(UA)
       << " UD" << s.displayTrilinear(DA)
       << " UE" << s.displayTrilinear(EA);
  left << " mQLsq" << s.displaySoftMassSquared(mQl)
       << " mURsq" << s.displaySoftMassSquared(mUr)
       << " mDRsq" << s.displaySoftMassSquared(mDr)
       << " mLLsq" << s.displaySoftMassSquared(mLl)
       << " mSEsq" << s.displaySoftMassSquared(mEr);
  left << "m3sq: " << s.displayM3Squared() << " mH1sq: " <<
    s.displayMh1Squared() << " mH2sq: " << s.displayMh2Squared() << '\n';
  left << "Gaugino masses" << s.displayGaugino();
  left << s.displaySusy();
  return left;
}

#undef HR

void SoftParsMssm::inputSoftParsOnly() {
  char c[70];

  cin >> c >> c >> c >> c >> c;
  cin >> c >> ua
       >> c >> da
       >> c >> ea;
  cin >> c >> mQLsq
       >> c >> mURsq
       >> c >> mDRsq
       >> c >> mLLsq
       >> c >> mSEsq;
  cin >> c >> m3sq >> c >> mH1sq >> c >> mH2sq;
  cin >> c >> mGaugino;
}

istream & operator >>(istream &left, SoftParsMssm &s) {
  char c[70];

  left >> c >> c >> c >> c >> c >> c >> c;
  DoubleMatrix ua(3, 3), da(3, 3), ea(3, 3);
  left >> c >> ua
       >> c >> da
       >> c >> ea;
  s.setTrilinearMatrix(UA, ua);
  s.setTrilinearMatrix(DA, da);
  s.setTrilinearMatrix(EA, ea);
  DoubleMatrix mqlsq(3, 3), mursq(3, 3), mdrsq(3, 3), mllsq(3, 3), mersq(3, 3);
  left >> c >> mqlsq
       >> c >> mursq
       >> c >> mdrsq
       >> c >> mllsq
       >> c >> mersq;
  s.setSoftMassMatrix(mQl, mqlsq);
  s.setSoftMassMatrix(mUr, mursq);
  s.setSoftMassMatrix(mDr, mdrsq);
  s.setSoftMassMatrix(mLl, mllsq);
  s.setSoftMassMatrix(mEr, mersq);
  double m3sq, mh1sq, mh2sq;
  left >> c >> m3sq >> c >> mh1sq >> c >> mh2sq;
  s.setM3Squared(m3sq); s.setMh1Squared(mh1sq); s.setMh2Squared(mh2sq);
  DoubleVector mg(3);
  left >> c >> mg;
  s.setAllGauginos(mg);
  MssmSusy ss;
  left >> ss;   s.setSusy(ss);
  return left;
}

// Boundary conditions to be applied at messenger scale for Gauge mediated
// SUSY breaking (see hep-ph/9703211 for example)
void SoftParsMssm::minimalGmsb(int n5, double lambda, double mMess,
                               double cgrav) {

// Modified thresholds by JEL 1-26-04 to accomodate numerical infinities

  const double epstol = 1.0e-4;
  double x = lambda / mMess;

  double f, g;

  if(fabs(x) < epstol) {
    g = 1.0 + x*x/6.0 + sqr(x*x)/15.0;
    f = 1.0 + x*x/36.0 + 11.0*sqr(x*x)/450.0;
  }
  else if(fabs(x-1.0) < 0.0001) {
    g  =  log(4.0);
    f  = -sqr(PI)/6.0 + log(4.0) + 0.5*sqr(log(4.0));
    g -=  0.0008132638905771205626;
    f -= -0.0049563838821509165200;
  }
  else {
    g = 1.0 / sqr(x) *
      ((1.0 + x) * log(1.0 + x) + (1.0 - x) * log(1.0 - x));
    f = (1.0 + x) / sqr(x) *
    (log(1.0 + x) - 2.0 * dilog(x / (1.0 + x)) + 0.5 *
     dilog(2.0 * x / (1.0 + x))) +
     (1.0 - x) / sqr(x) * (log(1.0 - x) - 2.0 * dilog(-x / (1.0 - x)) +
                         0.5 * dilog(-2.0 * x / (1.0 - x)));
  }

  double n5d = double(n5);

  double m1, m2, m3;
  m1 = n5d * sqr(displayGaugeCoupling(1)) / (16.0 * sqr(PI)) * lambda * g;
  m2 = n5d * sqr(displayGaugeCoupling(2)) / (16.0 * sqr(PI)) * lambda * g;
  m3 = n5d * sqr(displayGaugeCoupling(3)) / (16.0 * sqr(PI)) * lambda * g;
  setGauginoMass(1, m1);   setGauginoMass(2, m2);   setGauginoMass(3, m3);

  setM32(2.37e-19 * lambda * mMess * cgrav);

  double g1f = sqr(sqr(displayGaugeCoupling(1)));
  double g2f = sqr(sqr(displayGaugeCoupling(2)));
  double g3f = sqr(sqr(displayGaugeCoupling(3)));

  double mursq, mdrsq, mersq, mqlsq, mllsq;
  mursq = 2.0 * f * sqr(lambda) * n5d *
    (4.0 / 3.0 * g3f + 0.6 * 4.0 / 9.0 * g1f)
    / sqr(16.0 * sqr(PI));
  mdrsq = 2.0 * f * sqr(lambda) * n5d *
    (4.0 / 3.0 * g3f + 0.6 * 1.0 / 9.0 * g1f)
    / sqr(16.0 * sqr(PI));
  mersq = 2.0 * f * sqr(lambda) * n5d *
    (0.6 * g1f)
    / sqr(16.0 * sqr(PI));
  mqlsq = 2.0 * f * sqr(lambda) * n5d *
    (4.0 / 3.0 * g3f + 0.75 * g2f + 0.6 * g1f / 36.0)
    / sqr(16.0 * sqr(PI));
  mllsq = 2.0 * f * sqr(lambda) * n5d *
    (                  0.75 * g2f + 0.6 * 0.25 * g1f)
    / sqr(16.0 * sqr(PI));

  // You need Higgs masses too!

  DoubleMatrix id(3, 3);
  id(1, 1) = 1.0; id(2, 2) = 1.0; id(3, 3) = 1.0;

  setSoftMassMatrix(mQl, mqlsq * id);
  setSoftMassMatrix(mUr, mursq * id);
  setSoftMassMatrix(mDr, mdrsq * id);
  setSoftMassMatrix(mLl, mllsq * id);
  setMh1Squared(mllsq);
  setMh2Squared(mllsq);
  setSoftMassMatrix(mEr, mersq * id);

  universalTrilinears(0.0);
}

struct M32FuncParams {
    double alpha;
    double m0;
    double mpl;
};

double m32Func(double m32, void * p) {
    struct M32FuncParams * params = (struct M32FuncParams *) p;
    double alpha = (params->alpha);
    double m0 = (params->m0);
    double mpl = (params->mpl);

    return alpha * m0 - m32 / log (mpl / m32);
}

// KKLT
// Boundary conditions to be applied at GUT scale for Mirage Mediation soft terms
// Coding at 07/05/14 ref: hep-ph/0703163

void SoftParsMssm::mirage(double alphac, double M0, double aq, double al, double ahu, double ahd, double cq, double cl, double chu, double chd) {
    std::cout << "alpha: " << alphac << ", M0: " << M0 << '\n'
              << "aq: " << aq << ", al: " << al << '\n'
              << "ahu: " << ahu << ", ahd: " << ahd << '\n'
              << "cq: " << cq << ", cl:" << cl << '\n'
              << "chu: " << chu << ", chd: " << chd << '\n';

//    MssmSusy run(displaySusy());
//  double
//    yt   = run.displayYukawaElement(YU, 3, 3),
//    yb   = run.displayYukawaElement(YD, 3, 3),
//    ytau = run.displayYukawaElement(YE, 3, 3),
//    g1   = run.displayGaugeCoupling(1),
//    g2   = run.displayGaugeCoupling(2),
//    g3   = run.displayGaugeCoupling(3);

/// GET Yukawa couplings

  double yt = displayYukawaElement(YU, 3, 3);
  double yb = displayYukawaElement(YD, 3, 3);
  double ytau = displayYukawaElement(YE, 3, 3);

/// GET Gauge coupling squares

  double g1sq = sqr(displayGaugeCoupling(1));
  double g2sq = sqr(displayGaugeCoupling(2));
  double g3sq = sqr(displayGaugeCoupling(3));
  // cout << "g3(mgut)^2 = " << g3sq << '\n';

// SET Gravitino Mass, Planck Mass
// alpha*M0=m32/ln(Mpl/m32)  alpha*M0=500 GeV, Mpl = 2.4*10^18 GeV
  // setM32(16311.2);
  double Mpl = 2.4 * exp10(18);

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double x_lo = 1.0e+2, x_hi = 1.0e+6;
  gsl_function F;
  struct M32FuncParams params = {alphac, M0, Mpl};

  F.function = &m32Func;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  int status;
  int iter = 0, max_iter = 2000;
  double m32_root = 0;
  do {
      ++iter;
      status = gsl_root_fsolver_iterate(s);
      m32_root = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);

      if (status == GSL_SUCCESS) {
          std::cout << "m32: " << m32_root << '\n';
          break;
      }
  } while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free(s);

  setM32(m32_root);
  double m32 = displayGravitino();

// Define Mirage mediation Soft term coefficients

  double gamqp=8.0/3.0*g3sq + 3./2.*g2sq+1.0/30.0*g1sq,
         gamq3=8.0/3.0*g3sq + 3./2.*g2sq+1.0/30.0*g1sq-(sqr(yt)+sqr(yb)),
         gamup=8.0/3.0*g3sq + 8.0/15.0*g1sq,
         gamu3=8.0/3.0*g3sq + 8.0/15.0*g1sq -2.0*sqr(yt),
         gamdp=8.0/3.0*g3sq + 2.0/15.0*g1sq,
         gamd3=8.0/3.0*g3sq + 2.0/15.0*g1sq -2.0*sqr(yb),
         gamlp=3.0/2.0*g2sq + 3.0/10.0*g1sq,
         gaml3=3.0/2.0*g2sq + 3.0/10.0*g1sq - sqr(ytau),
         gamep=6.0/5.0*g1sq,
         game3=6.0/5.0*g1sq - 2.0*sqr(ytau),
         gamh1=3.0/2.0*g2sq + 3.0/10.0*g1sq - 3.0*sqr(yb) - sqr(ytau),
         gamh2=3.0/2.0*g2sq + 3.0/10.0*g1sq - 3.0*sqr(yt),

         bytkklt=6.0*sqr(yt) + sqr(yb) - 16.0/3.0*g3sq - 3.0*g2sq - 13.0/15.0*g1sq,
         bybkklt=sqr(yt) + 6.0*sqr(yb) + sqr(ytau) - 16.0/3.0*g3sq - 3.0*g2sq - 7.0/15.0*g1sq,
         bytaukklt=3.0*sqr(yb) + 4.0*sqr(ytau) - 3.0*g2sq - 9.0/5.0*g1sq,

         dgamqp=-8.0*sqr(g3sq) + 3.0/2.0*sqr(g2sq) + 11.0/50.0*sqr(g1sq),
         dgamq3=-8.0*sqr(g3sq) + 3.0/2.0*sqr(g2sq) + 11.0/50.0*sqr(g1sq) -(sqr(yt)*bytkklt + sqr(yb)*bybkklt),
         dgamup=-8.0*sqr(g3sq) + 88.0/25.0*sqr(g1sq),
         dgamu3=-8.0*sqr(g3sq) + 88.0/25.0*sqr(g1sq) - 2.0*sqr(yt)*bytkklt,
         dgamdp=-8.0*sqr(g3sq) + 22.0/25.0*sqr(g1sq),
         dgamd3=-8.0*sqr(g3sq) + 22.0/25.0*sqr(g1sq) - 2.0*sqr(yb)*bybkklt,
         dgamlp=3.0/2.0*sqr(g2sq) + 99.0/50.0*sqr(g1sq),
         dgaml3=3.0/2.0*sqr(g2sq) + 99.0/50.0*sqr(g1sq) - sqr(ytau)*bytaukklt,
         dgamep=198.0/25.0*sqr(g1sq),
         dgame3=198.0/25.0*sqr(g1sq) - 2.0*sqr(ytau)*bytaukklt,
         dgamh1=3.0/2.0*sqr(g2sq) + 99.0/50.0*sqr(g1sq) -3.0*sqr(yb)*bybkklt - sqr(ytau)*bytaukklt,
         dgamh2=3.0/2.0*sqr(g2sq) + 99.0/50.0*sqr(g1sq) -3.0*sqr(yt)*bytkklt,

         ttgamqp= 16.0/3.0*g3sq + 3.0*g2sq + 1.0/15.0*g1sq,
         ttgamq3= 16.0/3.0*g3sq + 3.0*g2sq + 1.0/15.0*g1sq - 2.0*(sqr(yt)*(ahu+2*aq)+sqr(yb)*(ahd+2*aq)),
         ttgamup= 16.0/3.0*g3sq + 16.0/15.0*g1sq,
         ttgamu3= 16.0/3.0*g3sq + 16.0/15.0*g1sq - 4.0*sqr(yt)*(ahu + 2*aq),
         ttgamdp= 16.0/3.0*g3sq + 4.0/15.0*g1sq,
         ttgamd3= 16.0/3.0*g3sq + 4.0/15.0*g1sq - 4.0*sqr(yb)*(ahd + 2*aq),
         ttgamlp= 3.0*g2sq + 3.0/5.0*g1sq,
         ttgaml3= 3.0*g2sq + 3.0/5.0*g1sq - 2*sqr(ytau)*(ahd + 2*al),
         ttgamep= 12.0/5.0*g1sq,
         ttgame3= 12.0/5.0*g1sq - 4.0*sqr(ytau)*(ahd + 2*al),
         ttgamh1= 3.0*g2sq + 3.0/5.0*g1sq - 6.0*sqr(yb)*(ahd + 2*aq) - 2.0*sqr(ytau)*(ahd + 2*al),
         ttgamh2= 3.0*g2sq + 3.0/5.0*g1sq - 6.0*sqr(yt)*(ahu + 2*aq);


// SET Gaugino Mass
// CAUTION !! Gaugino mass convention for softsusy = (-) * OURS

  double m1, m2, m3;

  m1 =  M0*( 1.0 + log(Mpl/m32)/(16.0*sqr(PI))*33.0/5.0*g1sq*alphac);
  m2 =  M0*( 1.0 + log(Mpl/m32)/(16.0*sqr(PI))*g2sq*alphac );
  m3 =  M0*( 1.0 - log(Mpl/m32)/(16.0*sqr(PI))*3*g3sq*alphac );

  setGauginoMass(1, m1);   setGauginoMass(2, m2);   setGauginoMass(3, m3);

// SET Soft scalar masses
  double mqlsq, mursq, mdrsq, mllsq, mersq;
  double mQlsq, mtrsq, mbrsq, mLlsq, mtaursq;
  double mh1sq, mh2sq;

  mqlsq =sqr(M0)*( cq - log(Mpl/m32)/(16.0*sqr(PI))*ttgamqp*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamqp*sqr(alphac));
  mursq =sqr(M0)*( cq - log(Mpl/m32)/(16.0*sqr(PI))*ttgamup*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamup*sqr(alphac));
  mdrsq =sqr(M0)*( cq - log(Mpl/m32)/(16.0*sqr(PI))*ttgamdp*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamdp*sqr(alphac));
  mllsq =sqr(M0)*( cl - log(Mpl/m32)/(16.0*sqr(PI))*ttgamlp*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamlp*sqr(alphac));
  mersq =sqr(M0)*( cl - log(Mpl/m32)/(16.0*sqr(PI))*ttgamep*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamep*sqr(alphac));

  mQlsq =sqr(M0)*( cq - log(Mpl/m32)/(16.0*sqr(PI))*ttgamq3*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamq3*sqr(alphac));
  mtrsq =sqr(M0)*( cq - log(Mpl/m32)/(16.0*sqr(PI))*ttgamu3*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamu3*sqr(alphac));
  mbrsq =sqr(M0)*( cq - log(Mpl/m32)/(16.0*sqr(PI))*ttgamd3*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamd3*sqr(alphac));
  mLlsq =sqr(M0)*( cl - log(Mpl/m32)/(16.0*sqr(PI))*ttgaml3*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgaml3*sqr(alphac));
  mtaursq =sqr(M0)*( cl - log(Mpl/m32)/(16.0*sqr(PI))*ttgame3*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgame3*sqr(alphac));

  mh1sq =sqr(M0)*( chd - log(Mpl/m32)/(16.0*sqr(PI))*ttgamh1*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamh1*sqr(alphac));
  mh2sq =sqr(M0)*( chu - log(Mpl/m32)/(16.0*sqr(PI))*ttgamh2*alphac
                  - sqr(log(Mpl/m32)/(16.0*sqr(PI)))*dgamh2*sqr(alphac));


  DoubleMatrix ID12(3,3),ID3(3,3),Mmqlsq(3,3),Mmursq(3,3),
                   Mmdrsq(3,3),Mmllsq(3,3),Mmersq(3,3);

  int i; for (i=1; i<=3; i++) {
             if (i<3) { ID12(i,i) = 1.0;  }
             else if (i==3) { ID3(i,i) = 1.0; }
  }

  Mmqlsq = ID12 * mqlsq + ID3 * mQlsq;
  Mmursq = ID12 * mursq + ID3 * mtrsq;
  Mmdrsq = ID12 * mdrsq + ID3 * mbrsq;
  Mmllsq = ID12 * mllsq + ID3 * mLlsq;
  Mmersq = ID12 * mersq + ID3 * mtaursq;

  setSoftMassMatrix(mQl, Mmqlsq);
  setSoftMassMatrix(mUr, Mmursq);
  setSoftMassMatrix(mDr, Mmdrsq);
  setSoftMassMatrix(mLl, Mmllsq);
  setSoftMassMatrix(mEr, Mmersq);
  setMh1Squared(mh1sq);
  setMh2Squared(mh2sq);

//**********************************
// SET Trilinear couplings
//**********************************

   double Atau,Ab,At,Ae,Ad,Au;

   Atau = M0*( (ahd + 2*al) - log(Mpl/m32)/(16.0*sqr(PI))*(gamh1+gaml3+game3)*alphac);
   Ab   = M0*( (ahd + 2*aq) - log(Mpl/m32)/(16.0*sqr(PI))*(gamh1+gamq3+gamd3)*alphac);
   At   = M0*( (ahu + 2*aq) - log(Mpl/m32)/(16.0*sqr(PI))*(gamh2+gamq3+gamu3)*alphac);

   Ae   = M0*( (ahd + 2*al) - log(Mpl/m32)/(16.0*sqr(PI))*(gamh1+gamlp+gamep)*alphac);
   Ad   = M0*( (ahd + 2*aq) - log(Mpl/m32)/(16.0*sqr(PI))*(gamh1+gamqp+gamdp)*alphac);
   Au   = M0*( (ahu + 2*aq) - log(Mpl/m32)/(16.0*sqr(PI))*(gamh2+gamqp+gamup)*alphac);

   DoubleMatrix MUA(3,3),MDA(3,3),MEA(3,3);

   MUA = Au * ID12 * displayYukawaMatrix(YU) + At * ID3 * displayYukawaMatrix(YU);
   MDA = Ad * ID12 * displayYukawaMatrix(YD) + Ab * ID3 * displayYukawaMatrix(YD);
   MEA = Ae * ID12 * displayYukawaMatrix(YE) + Atau * ID3 * displayYukawaMatrix(YE);

  setTrilinearMatrix(UA, MUA);
  setTrilinearMatrix(DA, MDA);
  setTrilinearMatrix(EA, MEA);


}
