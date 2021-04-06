
/** \file lowe.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: lowe.cpp,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.2  2005/05/30 14:22:14  allanach
   Fixed stau mixing in ISAWIG interface to 7.64

   Revision 1.24  2004/05/23 17:50:02  allanach
   Made it so calling toMz automatically calls the pole Mb calculation.
   This is so we can dispense with the annoying input file "massIn".
   Inputs will be specified in the SUSY Les Houches output file in the case of
   a single parameter point run, or in the main program if a scan is preformed.

   Revision 1.23  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.22  2003/12/18 16:30:42  allanach
   Fixed a bug in 3-loop QCD contribution to mass running

   Revision 1.21  2003/11/25 15:52:19  allanach
   If massIn doesn't exist, it initialises SM parameters with defaults

   Revision 1.20  2003/07/28 12:11:37  allanach
   More error trapping, and rearranging rpvsoftsusy to use correct Higgs VEV
   (which is sometimes called at MZ)

   Revision 1.19  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.18  2003/07/21 14:00:17  allanach
   MZ fully implemented as an input now. Kept MZ as the central PDG 2002 value,
   for defaults etc

   Revision 1.17  2003/07/21 09:34:55  allanach
   Changed input from alpha(MZ) to alpha(MZ)^{-1}

   Revision 1.16  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.15  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.14  2002/11/26 13:27:06  allanach
   Cleaned up QCD running so it works with a varying number of loops

   Revision 1.13  2002/09/24 16:15:38  allanach
   Added iteration routines to deal with charged leptons

   Revision 1.12  2002/09/23 15:12:15  allanach
   Quark mixing systematised.

   Revision 1.11  2002/09/20 15:37:10  allanach
   Adding quark-mixing routines

   Revision 1.10  2002/09/19 16:22:03  allanach
   Latest PDG data implemented

   Revision 1.9  2002/09/09 10:42:54  allanach
   TOLERANCE replaces EPS as being more user-friendly

   Revision 1.8  2002/09/03 14:29:38  allanach
   EPS, PRINTOUT, MIXING now all set in massIn

   Revision 1.7  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, EPS out of def.h to make it quicker to compile once
   they are changed.

   Revision 1.6  2002/08/30 12:31:55  allanach
   Test RPV version that works!

   Revision 1.5  2002/08/13 13:52:38  allanach
   Bug-fix in operator >>

   Revision 1.4  2002/06/14 16:26:30  allanach
   Switches included for 2-loop running of scalar masses, and calulating mt at
   mt.

   Revision 1.5  2002/02/22 15:58:05  allanach
   Now number of loops between MSbar and pole mb conversion depends upon the
   number of loops used in the QedQcd object

   Revision 1.4  2002/02/20 17:20:02  allanach
   Treatment of bottom quark mass input changed. Now can input mb(mb) and/or
   pole(mb). If only one is input, with the other being a question mark,
   it will be calculated to three loop QCD from the input.

   Revision 1.3  2002/02/18 19:18:25  allanach
   Input is now bottom POLE mass. Running MSbar mass is now calculated!

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#include "lowe.h"
#include <cstring>

///  external object temp used to get objects into external routines, however:
///  don't use it!
static QedQcd *tempLe;

QedQcd::QedQcd()
  : a(2), mf(9), mtPole(PMTOP), mbPole(PMBOTTOM) {
  setPars(11);
  // Default object: 1998 PDB defined in 'def.h'
  mf(1) = MUP; mf(2) = MCHARM;
  mf(4) = MDOWN; mf(5) = MSTRANGE; mf(6) = MBOTTOM;
  mf(7) = MELECTRON; mf(8) = MMUON; mf(9) = MTAU;
  a(1) = 7.8180608e-3;  a(2) = ALPHASMZ;
  mf(3) = getRunMtFromMz(PMTOP, ALPHASMZ);
  setMu(MZ);
  setLoops(3);
  setThresholds(1);
}

const QedQcd & QedQcd::operator=(const QedQcd & m) {
  if (this == &m) return *this;
  mtPole = m.mtPole;
  mbPole = m.mbPole;
  a = m.a;
  mf = m.mf;
  setLoops(m.displayLoops());
  setThresholds(m.displayThresholds());
  return *this;
}

//For communication with outside routines: sets all data by one vector y=1..11.
void QedQcd::set(const DoubleVector & y) {
  a(ALPHA) = y.display(1);
  a(ALPHAS) = y.display(2);
  int i; for (i=3; i<=11; i++)
    mf(i-2) = y.display(i);
}

DoubleVector QedQcd::display() const {
  DoubleVector y(11);
  y(1) = a.display(ALPHA);
  y(2) = a.display(ALPHAS);
  int i; for (i=3; i<=11; i++)
    y(i) = mf.display(i-2);
  return y;
}

//  Active flavours at energy mu
int QedQcd::flavours(double mu) const {
  int k = 0;
  if (mu > mf.display(mTop)) k++;
  if (mu > mf.display(mCharm)) k++;
  if (mu > mf.display(mUp)) k++;
  if (mu > mf.display(mDown)) k++;
  if (mu > mf.display(mBottom)) k++;
  if (mu > mf.display(mStrange)) k++;
  return k;
}

ostream & operator <<(ostream &left, const QedQcd &m) {
  left << "mU: " << m.displayMass(mUp)
       << "  mC: " << m.displayMass(mCharm)
       << "  mt: " << m.displayMass(mTop)
       << "  mt^pole: " << m.displayPoleMt()
       << endl;
  left << "mD: " << m.displayMass(mDown)
       << "  mS: " << m.displayMass(mStrange)
       << "  mB: " << m.displayMass(mBottom)
       << "  mb^pole: " << m.displayPoleMb()
       << endl;
  left << "mE: " << m.displayMass(mElectron)
       << "  mM: " << m.displayMass(mMuon)
       <<  "  mT: " << m.displayMass(mTau) << endl;
  left << "aE: " << 1.0 / m.displayAlpha(ALPHA) <<
    "  aS: " << m.displayAlpha(ALPHAS)
       << "  scale: " << m.displayMu() << endl;
  left << "loops: " << m.displayLoops()
       << "        thresholds: " << m.displayThresholds() << endl;

  return left;
}

istream & operator >>(istream &left, QedQcd &m) {

  char c[12], cmbmb[12], cmbpole[12];
  double mu, mc, mtpole, md, ms, me, mmu, mtau, invalph,
    alphas, scale;
  int t, l;
  left >> c >> mu >> c >> mc >> c >> c >> c >> mtpole;
  left >> c >> md >> c >> ms >> c >> cmbmb >> c >> cmbpole;
  left >> c >> me >> c >> mmu >> c >> mtau;
  left >> c >> invalph >> c >> alphas >> c >> scale;
  left  >> c >> l >> c >> t;
  m.setMass(mUp, mu);
  m.setMass(mCharm, mc);
  m.setMass(mDown, md);
  m.setMass(mStrange, ms);
  m.setMass(mElectron, me);
  m.setMass(mMuon, mmu);
  m.setMass(mTau, mtau);
  m.setAlpha(ALPHA, 1.0 / invalph);
  m.setAlpha(ALPHAS, alphas);
  m.setMu(scale);
  // y[3] is pole mass
  m.setPoleMt(mtpole);

  // default 3-loop qcd calculation
  m.setLoops(l);
  m.setThresholds(t);

  m.setMass(mTop, getRunMtFromMz(mtpole, alphas));

  if (!strcmp(cmbmb, "?") && !strcmp(cmbpole, "?")) {
    ostringstream ii;
    ii << "Error reading in low energy QCDQED object: must specify ";
    ii << "running AND/OR pole bottom mass" << endl;
    throw ii.str();
  }

  // If you set one of the bottom mass parameters to be "?", it will calculate
  // it from the other one
  if (strcmp(cmbmb, "?") != 0) m.setMass(mBottom, atof(cmbmb));
  if (strcmp(cmbpole, "?") != 0) m.setPoleMb(atof(cmbpole));

  if (strcmp(cmbmb, "?") == 0) m.calcRunningMb();
  if (strcmp(cmbpole, "?") == 0) m.calcPoleMb();

  return left;
}

//  returns qed beta function at energy mu < mtop
double QedQcd::qedBeta() const {
  double x;
  x = 24.0 / 9.0;
  if (displayMu() > mf.display(mCharm)) x += 8.0 / 9.0;
  if (displayMu() > mf.display(mTop)) x += 8.0 / 9.0;
  if (displayMu() > mf.display(mBottom)) x += 2.0 / 9.0;
  if (displayMu() > mf.display(mTau)) x += 2.0 / 3.0;
  if (displayMu() > MW) x += 1.0 / 6.0;
  if (displayMu() > (mtPole + TOLERANCE))  {
    ostringstream ii;
      ii << "qed beta function called at " << displayMu() <<
        "above mt=" << displayPoleMt() <<
        ", outside range of validity";
      ii << " in QedQcd::qedbeta\n";
      throw ii.str();
    }

  return (x * sqr(a.display(ALPHA)) / PI);
}

//  next routine calculates beta function to 3 loops in qcd for The Standard
//  Model. Note that if quark masses are running, the number of active quarks
//  will take this into account. Returns beta
double QedQcd::qcdBeta() const {
  static const double INVPI = 1.0 / PI;
  int quarkFlavours = flavours(this->displayMu());
  double qb0, qb1, qb2;
  qb0 = (11.0e0 - (2.0e0 / 3.0e0 * quarkFlavours)) / 4.0;
  qb1 = (102.0e0 - (38.0e0 * quarkFlavours) / 3.0e0) / 16.0;
  qb2 = (2.857e3 * 0.5 - (5.033e3 * quarkFlavours) / 18.0  +
         (3.25e2 * sqr(quarkFlavours) ) / 5.4e1) / 64;
  double qa0, qa1, qa2;

  if (displayLoops() < 0 || displayLoops() > 3) {
    ostringstream ii;
      ii << "Wrong loops parameter :" << displayLoops() << " in " << *this;
      throw ii.str();
    }
  if (displayLoops() > 0) qa0 = qb0 * INVPI; else qa0 = 0.0;
  if (displayLoops() > 1) qa1 = qb1 * sqr(INVPI); else qa1 = 0.0;
  if (displayLoops() > 2) qa2 = qb2 * sqr(INVPI) * INVPI; else qa2 = 0.0;

  // add contributions of the one, two and three loop constributions resp.
  double beta;
  beta = -2.0 * sqr(displayAlpha(ALPHAS)) *
    (qa0 + qa1 * displayAlpha(ALPHAS) + qa2 *
     sqr(displayAlpha(ALPHAS)));
  return beta;
}

//(See comments for above function). returns a vector x(1..9) of fermion mass
//beta functions -- been checked!
void QedQcd::massBeta(DoubleVector & x) const {
  static const double INVPI = 1.0 / PI, ZETA3 = 1.202056903159594;

  // qcd bits: 1,2,3 loop resp.
  double qg1, qg2, qg3;
  int quarkFlavours = flavours(displayMu());
  if (displayLoops() > 0) qg1 = INVPI; else qg1 = 0.0;
  if (displayLoops() > 1)
    qg2 = (202.0 / 3.0 - (20.0e0 * quarkFlavours) / 9.0) * sqr(INVPI) / 16.0;
  else
    qg2 = 0.0;
  if (displayLoops() > 2)
    qg3 = (1.249e3 - ((2.216e3 * quarkFlavours) / 27.0e0 +
                      1.6e2 * ZETA3 * quarkFlavours / 3.0e0) -
           140.0e0 * quarkFlavours * quarkFlavours / 81.0e0) * sqr(INVPI) *
      INVPI / 64.0;
  else
    qg3 = 0.0;

  double qcd = -2.0 * a.display(ALPHAS) * (qg1  + qg2 * a.display(ALPHAS) +
                                           qg3 * sqr(a.display(ALPHAS)));
  double qed = -a.display(ALPHA) * INVPI / 2;

  int i;
  for (i=1;i<=3;i++)   // up quarks
    x(i) = (qcd + 4.0 * qed / 3.0) * mf.display(i);
  for (i=4;i<=6;i++)   // down quarks
    x(i) = (qcd + qed / 3.0) * mf.display(i);
  for (i=7;i<=9;i++)   // leptons
    x(i) = 3.0 * qed * mf.display(i);
  // switch off relevant beta functions
  if (displayThresholds() > 0)
    for(i=1;i<=9;i++) if (displayMu() < displayMass()(i)) x(i) = 0.0;
  // nowadays, u,d,s masses defined at 2 GeV: don't run them below that
  if (displayMu() < 2.0) x(1) = x(4) = x(5) = 0.0;
}

DoubleVector QedQcd::beta() const {
  DoubleVector dydx(11);
  dydx(1) = qedBeta();
  dydx(2) = qcdBeta();
  DoubleVector y(9);
  massBeta(y);
  int i; for (i=3; i<=11; i++)
    dydx(i) = y(i-2);
  return dydx;
}

void QedQcd::runGauge(double x1, double x2) {
  const double tol = 1.0e-5;

  DoubleVector y(2);
  tempLe = this;
  y(1) = tempLe->displayAlpha(ALPHA);
  y(2) = tempLe->displayAlpha(ALPHAS);

  callRK(x1, x2, y, gaugeDerivs, tol);

  setAlpha(ALPHA, y(1));
  setAlpha(ALPHAS, y(2));
}

// Done at pole mb: extracts running mb(polemb)
double QedQcd::extractRunningMb(double alphasMb) {
  double mbPole = displayPoleMb();

  if (displayMu() != mbPole) {
    ostringstream ii;
    ii << "ERROR: QedQcd::extractRunningMb called at scale "
         << displayMu() << " instead of mbpole\n";
    throw ii.str();
  }

  // Following is the MSbar correction from QCD, hep-ph/9912391 and ZPC48 673
  // (1990)
  double delta = 0.;
  if (displayLoops() > 0) delta = delta + 4.0 / 3.0 * alphasMb / PI;
  if (displayLoops() > 1)
    delta = delta + sqr(alphasMb / PI) *
      (10.1667 + (displayMass(mUp) + displayMass(mDown) +
                  displayMass(mCharm) + displayMass(mStrange)) / mbPole);
  if (displayLoops() > 2)
    delta = delta + 101.45424 * alphasMb / PI * sqr(alphasMb / PI);

  double mbmb = mbPole * (1.0 - delta);

  return mbmb;
}

// Supposed to be done at mb(mb) -- MSbar, calculates pole mass
double QedQcd::extractPoleMb(double alphasMb) {

  if (displayMu() != displayMass(mBottom)) {
    ostringstream ii;
    ii << "ERROR: QedQcd::extractPoleMb called at scale " << displayMu() <<
      " instead of mb(mb)\n";
    throw ii.str();
  }

  // Following is the MSbar correction from QCD, hep-ph/9912391
  double delta = 0.0;
  if (displayLoops() > 0) delta = delta + 4.0 / 3.0 * alphasMb / PI;
  if (displayLoops() > 1) delta = delta + sqr(alphasMb / PI) *
    (9.2778 + (displayMass(mUp) + displayMass(mDown) + displayMass(mCharm) +
               displayMass(mStrange)) / mbPole);
  if (displayLoops() > 2)
    delta = delta + 94.4182 * alphasMb / PI * sqr(alphasMb / PI);

  double mbPole = displayMass(mBottom) * (1.0 + delta);

  return mbPole;
}

// Calculates the running mass from the pole mass:
void QedQcd::calcRunningMb() {

  const double tol = 1.0e-5;

  // Save initial object
  DoubleVector saving(display());
  double saveMu = displayMu();

  // Set arbitrarily low bottom mass to make sure it's included in the RGEs
  setMass(mBottom, 0.);
  runto(displayPoleMb(), tol);
  double mbAtPoleMb = extractRunningMb(displayAlpha(ALPHAS));
  setMass(mBottom, mbAtPoleMb);
  // Now, by running down to 1 GeV, you'll be left with mb(mb) since it will
  // decouple at this scale.
  runto(1.0, tol);
  double mbmb = displayMass(mBottom);

  // restore initial object
  set(saving);
  setMu(saveMu);
  setMass(mBottom, mbmb);
}

// Calculates the pole mass from the running mass, which should be defined at
// mb
void QedQcd::calcPoleMb() {

  double alphasMZ = displayAlpha(ALPHAS);
  double alphaMZ = displayAlpha(ALPHA);
  double saveMu = displayMu();

  runGauge(displayMu(), displayMass(mBottom));
  double poleMb = extractPoleMb(displayAlpha(ALPHAS));
  setPoleMb(poleMb);

  // Reset to erase numerical integration errors.
  setAlpha(ALPHAS, alphasMZ);
  setAlpha(ALPHA, alphaMZ);
  setMu(saveMu);
}

// Takes QedQcd object created at MZ and spits it out at mt
void QedQcd::toMt() {

  const double tol = 1.0e-5;

  double alphasMZ = displayAlpha(ALPHAS);
  double alphaMZ = displayAlpha(ALPHA);

  double mz = displayMu();

  runGauge(mz, 1.0);
  //Run whole lot up to pole top mass
  double mt = this->displayPoleMt();

  run(1.0, mz, tol);

  // Reset alphas to erase numerical integration errors.
  setAlpha(ALPHAS, alphasMZ);
  setAlpha(ALPHA, alphaMZ);
  run(mz, mt, tol);
}

// Takes QedQcd object created at MZ and spits it out at MZ
void QedQcd::toMz() {
  calcPoleMb();

  const double tol = 1.0e-5;

  double alphasMZ = displayAlpha(ALPHAS);
  double alphaMZ = displayAlpha(ALPHA);
  double mz = displayMu();
  runGauge(mz, 1.0);
  run(1.0, mz, tol);
  // Reset alphas to erase numerical integration errors.
  setAlpha(ALPHAS, alphasMZ);
  setAlpha(ALPHA, alphaMZ);
}


// This will calculate the three gauge couplings of the Standard Model at the
// scale m2.
// It's a simple one-loop calculation only and no
// thresholds are assumed. Range of validity is electroweak to top scale.
// alpha1 is in the GUT normalisation. sinth = sin^2 thetaW(Q) in MSbar
// scheme
DoubleVector QedQcd::getGaugeMu(const double m2, const double sinth) const {
  static const double INVPI = 1.0 / PI;
  /*
  if (displayMu() != MZ) {
    ostringstream ii;
    ii << "Error in lowe.cpp:QedQcd::getGaugeMu called with mu=" <<
      displayMu() << " not MZ=" << MZ << endl;
    throw ii.str();
    }
  */
  DoubleVector temp(1, 3);

  double a1, a2, aem = displayAlpha(ALPHA), m1 = displayMu();
  // Set alpha1,2 at scale m1 from data:
  a1 = 5.0 * aem / (3.0 * (1.0 - sinth));
  a2 = aem / sinth;

  // Renormalise a1,a2 to scale m2 assuming topless SM with one light Higgs
  // doublet
  a1 = 1.0 / ( 1.0 / a1 + 4.0 * INVPI * 1.07e2 * log(m1 / m2) / 2.4e2 );
  a2 = 1.0 / ( 1.0 / a2 - 4.0 * INVPI * 2.50e1 * log(m1 / m2) / 4.8e1 );

  temp.set(1, a1);
  temp.set(2, a2);
  // calculate alphas(m2)
  QedQcd oneset(*this); oneset.runto(oneset.displayPoleMt());
  // Set alphas(m) to be what's already calculated.
  temp.set(3, oneset.displayAlpha(ALPHAS));

  return temp;
}

int accessedReadIn; // Should be initialised to zero at start of prog
/*
--------------- read in a qcd-type object ------------------
Call with fname "" if you want it to come from standard input

"massIn" is an example of a data initialisation file:
*/
void readIn(QedQcd &mset, char fname[80]) {
   static QedQcd prevReadIn; // Data will be stored in here for rest of the
                                // run

  // Read in data if it's not been set
  if (accessedReadIn == 0) {
    char c[14];
    if (!strcmp(fname,"")) cin >> prevReadIn >> c >> MIXING >> c >> TOLERANCE
                               >> c >> PRINTOUT; // from standard input
    else {
      // read from filename fname
          fstream fin(fname, ios::in);
          if(!fin) {
            mset = QedQcd();
            return;
            ostringstream ii;
            ii << "Can't find input file " << fname << endl;
            throw ii.str();
          }
          fin >> prevReadIn >> c >> MIXING >> c >> TOLERANCE >> c >> PRINTOUT;
          fin.close();
    }

    if (PRINTOUT) cout << prevReadIn;
    accessedReadIn = 1; // Flag the fact we've read in the data once
  }

  mset = prevReadIn;

}

DoubleVector gaugeDerivs(double x, const DoubleVector & y) {
  tempLe->setMu(exp(x));
  tempLe->setAlpha(ALPHA, y.display(1));
  tempLe->setAlpha(ALPHAS, y.display(2));
  DoubleVector dydx(2);
  dydx(1) = tempLe->qedBeta();
  dydx(2) = tempLe->qcdBeta();

  return dydx;
}

// Given pole mass and alphaS(MZ), returns running top mass -- one loop qcd
double getRunMtFromMz(double poleMt, double asMZ) {
  return getRunMt(poleMt, getAsmt(poleMt, asMZ));
}

// Input pole mass of top and alphaS(mt), outputs running mass mt(mt)
// including one-loop standard model correction only
double getRunMt(double poleMt, double asmt) {
  return poleMt / (1.0 + (4.0 / (3.0 * PI)) * asmt);
}

// Given a value of mt, and alphas(MZ), find alphas(mt) to 1 loops in qcd:
// it's a very good approximation at these scales, better than 10^-3 accuracy
double getAsmt(double mtop, double alphasMz) {
  return alphasMz /
      (1.0 - 23.0 * alphasMz / (6.0 * PI) * log(MZ / mtop));
}

// input diagonal matrices and it'll give you back mixed ones
void doQuarkMixing(DoubleMatrix & vCkm, DoubleMatrix & mDon,
                 DoubleMatrix & mUpq) {
    switch(MIXING) {
    case 1:
      mUpq = vCkm.transpose() * mUpq * vCkm;
      break;
    case 2:
      mDon = vCkm * mDon * vCkm.transpose();
      break;
    default:
      return;
  }
}


// We must first define a down-quark mass matrix: 3 x 3. QedQcd should be at MZ
void massFermions(const QedQcd & r, DoubleMatrix & mDon,
                           DoubleMatrix & mUpq, DoubleMatrix & mEle) {
  /*  if (r.displayMu() != MZ) {
    ostringstream ii;
    ii << "In mQuarks. Called with wrong scale=" << r.displayMu();
    throw ii.str();
  }
  */
  if (MIXING < -1 || MIXING > 2) {
    ostringstream ii;
    ii << "In RgeLowe::massFermion.";
    ii << " MIXING=" << MIXING << " is out of range (-1 -> 2)\n";
    throw ii.str();
  }

  mDon(3, 3) = r.displayMass(mBottom);
  mUpq(3, 3) = r.displayMass(mTop);
  mEle(3, 3) = r.displayMass(mTau);

  if (MIXING > -1) {
    mDon(1, 1) = r.displayMass(mDown);
    mDon(2, 2) = r.displayMass(mStrange);
    mUpq(1, 1) = r.displayMass(mUp);
    mUpq(2, 2) = r.displayMass(mCharm);
    mEle(1, 1) = r.displayMass(mElectron);
    mEle(2, 2) = r.displayMass(mMuon);
  }
}
