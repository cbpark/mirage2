
/** \file lowe.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: QedQcd object contains Standard Model quark and lepton 
   masses. It integrates them using 3 loop qcd x 1 loop qed effective theory.

   $Log: lowe.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.4  2005/08/16 17:22:08  allanach
   Corrected electroweak sbottom corrections

   Revision 1.3  2005/08/16 13:55:16  allanach
   New central mt

   Revision 1.2  2004/12/17 12:06:39  allanach
   PDG 2004 data

   Revision 1.9  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.8  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.7  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.6  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.5  2002/09/23 18:13:50  allanach
   Eigenvalue order with non symmetric diagonalisation setx

   Revision 1.4  2002/09/20 15:37:10  allanach
   Adding quark-mixing routines

   Revision 1.3  2002/07/30 12:57:31  allanach
   SOFTSUSY1.5

   Revision 1.5  2002/02/20 17:20:02  allanach
   Treatment of bottom quark mass input changed. Now can input mb(mb) and/or
   pole(mb). If only one is input, with the other being a question mark, it
   will be calculated to three loop QCD from the input.

   Revision 1.4  2002/02/18 19:18:25  allanach
   Input is now bottom POLE mass. Running MSbar mass is now calculated!

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef LOWE_H
#define LOWE_H

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::ostream;
using std::istream;
#include <fstream>
using std::fstream;
using std::ios;
#include <sstream>
using std::ostringstream;
#include <iostream>
#include "def.h"
#include "utils.h"
#include "linalg.h"
#include "rge.h"

const double MUP = 3.0e-3; ///< default running quark mass from PDG
const double MDOWN = 6.75e-3; ///< default running quark mass from PDG
const double MSTRANGE = 0.1175; ///< default running quark mass from PDG
const double MCHARM = 1.2; ///< default running quark mass from PDG
const double MBOTTOM = 4.25; ///< default running quark mass from PDG
const double MTOP = 165.0; ///< default running quark mass from PDG
/// default pole lepton mass from PDG
const double MELECTRON = 5.10998902e-4; 
const double MMUON = 1.05658357e-1; ///< default pole lepton mass from PDG
const double MTAU = 1.77699; ///< default pole lepton mass from PDG
const double ALPHASMZ = 0.1187; ///< default running mass from PDG
const double ALPHAMZ = 1.0 / 127.918; ///< default running alpha(MZ) from PDG

const double PMTOP = 173.4; ///< default pole mass from PDG
const double PMBOTTOM = 4.9; ///< default pole mass from PDG

/// used to give order of quark masses stored
typedef enum {mUp=1, mCharm, mTop, mDown, mStrange, mBottom, mElectron,
	      mMuon, mTau} mass;
/// order of gauge couplings stored in QedQcd
typedef enum {ALPHA=1, ALPHAS} leGauge;

/// Returns beta functions of alpha, alpha_s only
DoubleVector gaugeDerivs(double, const DoubleVector &);

class QedQcd;
/// Contains data on quark and lepton masses, as well as gauge couplings in an
/// effective QEDxQCD theory.
class QedQcd: public RGE
{
private:
  DoubleVector a;   ///< gauge couplings
  DoubleVector mf;  ///< fermion running masses
  double mtPole, mbPole; ///< pole masses of third family quarks
  
public:
  QedQcd(); ///< Initialises with default values defined in lowe.h
  QedQcd(const QedQcd &); ///< Initialises object with another
  const QedQcd& operator=(const QedQcd & m); ///< Sets two objects equal
  virtual ~QedQcd() {};
  
  void setPoleMt(double mt) { mtPole = mt; }; ///< set pole top mass
  void setPoleMb(double mb) { mbPole = mb; }; ///< set pole bottom mass
  /// sets a running quark mass
  void setMass(mass mno, double m) { mf(mno) = m; }; 
  /// sets QED or QCD structure constant
  void setAlpha(leGauge ai, double ap) { a(ai) = ap; }; 
  /// For exporting beta functions to Runge-Kutta
  void set(const DoubleVector &); 
  
  /// Display pole top mass
  double displayPoleMt() const { return mtPole; };
  /// Returns bottom "pole" mass
  double displayPoleMb() const { return mbPole; };
  /// Returns a vector of running fermion masses
  DoubleVector displayMass() const { return mf; };
  /// Returns a single running mass
  double displayMass(mass mno) const { return mf.display(mno); };
  /// Returns a single gauge structure constant
  double displayAlpha(leGauge ai) const { return a.display(ai); };
  /// Obgligatory: returns vector of all running parameters
  DoubleVector display() const;
  
  int flavours(double) const;  /// returns number of active flavours
  
  double qedBeta() const;   ///< QED beta function
  double qcdBeta() const;   ///< QCD beta function
  void massBeta(DoubleVector &) const; ///< beta functions of masses
  /// Beta functions of both beta-functions and all MSbar masses
  DoubleVector beta() const; 
  
  /// Does not run the masses, just gauge couplings from start to end
  void runGauge(double start, double end);
  /// calculates pole bottom mass given alpha_s(Mb)^{MSbar} from running b mass
  double extractPoleMb(double asMb);
  /// Done at pole mb: extracts running mb(polemb)
  double extractRunningMb(double asMb);
  /// calculates running bottom mass given alpha_s(Mb)^{MSbar} from pole m_b
  void calcRunningMb();
  /// Calculates the pole mass from the running mass, which should be defined
  /// at mb
  void calcPoleMb();

  /// Evolves object to running top mass
  void toMt();
  /// Evolves object to MZ
  void toMz();
  /// This will calculate the three gauge couplings of the Standard Model at
  /// the scale m2.
  /// It's a simple one-loop calculation only and no
  /// thresholds are assumed. Range of validity is electroweak to top scale.
  // alpha1 is in the GUT normalisation. sinth = sin^2 thetaW(Q) in MSbar
  // scheme
  DoubleVector  getGaugeMu(const double m2, const
		     double sinth) const;
};

/// Input numbers into the object: by file stream
ostream & operator <<(ostream &, const QedQcd &);
/// Formatted output from QedQcd object
istream & operator >>(istream &left, QedQcd &m);

/// Reads in a QedQed-type object and returns it in oneset.
/// Call with fname "" if you want it to come from standard input
/// "massIn" is an example of a data initialisation file: 
void readIn(QedQcd & oneset, char fname[80]); 
/// Input pole mass of top and alphaS(mt), outputs running mass mt(mt)
/// including one-loop standard model correction only
double getRunMt(double poleMt, double asmt);
/// Given a value of mt, and alphas(MZ), find alphas(mt) to 1 loops in qcd:
/// it's a very good approximation at these scales, better than 10^-3 accuracy
double getAsmt(double mtop, double alphasMz);
/// Given pole mass and alphaS(MZ), returns running top mass -- one loop qcd
double getRunMtFromMz(double poleMt, double asMZ);

inline QedQcd::QedQcd(const QedQcd &m)
  : a(m.a), mf(m.mf), mtPole(m.mtPole), mbPole(m.mbPole) { 
  setPars(11); 
  setMu(m.displayMu());
  setLoops(m.displayLoops());
  setThresholds(m.displayThresholds());
}

/// Returns diagonal fermion mass matrices given input object r
void massFermions(const QedQcd & r, DoubleMatrix & mDon, 
		  DoubleMatrix & mUpq, DoubleMatrix & mEle);
/// Input diagonal mass matrices and it'll give you back mixed ones, based on
/// the CKM quark mixing matrix you supplied in vCkm
void doQuarkMixing(DoubleMatrix & vCkm, DoubleMatrix & mDon, 
		 DoubleMatrix & mUpq);
#endif

