
/** \file softpars.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: Soft SUSY breaking parameters

   $Log: softpars.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.17  2004/01/15 13:54:55  allanach
   New heaer style implemented

   Revision 1.16  2003/10/24 16:09:04  allanach
   Implemented running Higgs DRbar vev

   Revision 1.14  2003/07/21 14:00:18  allanach
   MZ fully implemented as an input now. Kept MZ as the central PDG 2002 value,
   for defaults etc

   Revision 1.12  2003/05/27 15:05:52  allanach
   Purely efficiency corrections: used variable rather than display() methods
   whenever possible

   Revision 1.11  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.8  2002/07/30 12:57:32  allanach
   SOFTSUSY1.5

   Revision 1.7  2002/06/14 16:26:30  allanach
   Switches included for 2-loop running of scalar masses, and calulating
   mt at mt.

   Revision 1.6  2002/04/26 15:14:44  allanach
   Deleted all translation routines and defined boundary conditions within
   softsusy.h and softsusy.cpp

   Revision 1.5  2002/04/14 13:50:41  allanach
   Now use V=m3^2 H1 H2 instead of V=mu B H1 H2. It's more natural!

   Revision 1.4  2002/04/12 16:51:27  allanach
   Added display/set functions to work automatically

   Revision 1.3  2002/04/12 06:24:50  allanach
   Code maintenance - returning a subobject made simpler

   Revision 1.3  2001/09/28 14:49:02  allanach
   GMSB boundary conditions added

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef SOFTPARS_H
#define SOFTPARS_H

#include <cmath>
#include "susy.h"
#include "def.h"
#include "linalg.h"
#include "utils.h"
#include "numerics.h"

/// SUSY breaking soft mass squared parameters
typedef enum {mQl=1, mUr, mDr, mLl, mEr} softMasses;
/// SUSY breaking trilinear parameters
typedef enum {UA=1, DA, EA} trilinears;

/// Number of parameters contained in RGEs
const static int numSoftParsMssm = 78 + numSusyPars;

/// Soft SUSY breaking parameters and beta functions.
class SoftParsMssm: public MssmSusy
{
private:
  DoubleVector mGaugino; ///< Gaugino masses, see ::beta for definitions
  DoubleMatrix ua, da, ea; ///< Trilinear soft terms..
  /// soft mass squared matrices of \f$ m_Q^2, m_U^2, m_D^2, m_L^2, m_E^2 \f$
  /// respectively.
  DoubleMatrix mQLsq, mURsq, mDRsq, mLLsq, mSEsq;
  /// Bilinear Higgs soft parameters: \f$ m_3^2, m_{H_1}^2, m_{H_2}^2 \f$
  /// respectively.
  double m3sq, mH1sq, mH2sq;
  double m32;         ///< Gravitino mass
public:
  /// Default constructor fills object with zeroes
  SoftParsMssm();
  /// Constructor fills SUSY conserving parts with another object, all
  /// SUSY breaking parameters set to zero
  SoftParsMssm(const MssmSusy &);
  /// Constructor sets all parameters equal to those in another object
  SoftParsMssm(const SoftParsMssm &);
  /// Sets all parameters equal to those in another object
  const SoftParsMssm & operator=(const SoftParsMssm & s);
  /// Constructor sets RPC SUSY parameters to s, gaugino masses to mG,
  /// trilinears to aU, aD, aE for au, ad, ae
  /// trilnears respectively,  \f$m_Q^2\f$=mQl, \f$m_U^2\f$=mUr,
  /// \f$m_D^2\f$=mDr, \f$m_L^2\f$=mLl, \f$m_E^2\f$=mEr, \f$ m_3^2\f$=m3sq,
  /// \f$m_{H_1}^2\f$=mH1sq, \f$m_{H_2}^2\f$=mH2sq, mu parameter, number of
  /// loops=l, and threshold parameter=t
  SoftParsMssm(const MssmSusy & s, const DoubleVector & mG, const
               DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix
               & aE, const DoubleMatrix & mQl, const DoubleMatrix & mUr, const
               DoubleMatrix & mDr, const DoubleMatrix & mLl, const
               DoubleMatrix & mEr, double m3sq, double mH1sq, double mH2sq,
               double mGravitino, double mu, int l, int t);

  /// Returns whole object as a const
  inline SoftParsMssm displaySoftPars() const;

  /// Return a trilinear coupling matrix
  DoubleMatrix displayTrilinear(trilinears) const;
  /// Return a trilinear element
  double displayTrilinear(trilinears, int i, int j) const;
  /// Return a trilinear element in "SUGRA style"
  double displaySoftA(trilinears, int, int) const;
  /// Return a soft mass squared matrix
  DoubleMatrix displaySoftMassSquared(softMasses) const;
  /// Return a soft mass squared element
  double displaySoftMassSquared(softMasses, int i, int j) const;

  double displayGravitino() const; ///< Returns the gravitino mass
  inline double displayM3Squared() const;     ///< Return \f$ m_3^2\f$
  inline double displayMh1Squared() const;    ///< Return \f$m_{H_1}^2\f$
  inline double displayMh2Squared() const;    ///< Return \f$m_{H_2}^2\f$=mH2sq
  inline DoubleVector displayGaugino() const; ///< Return \f$M_{G_i}\f$
  inline double displayGaugino(int i) const;  ///< Return \f$M_{G_i}\f$
  /// Return contents of object in a vector: for RG evolution
  virtual DoubleVector display() const;

  /// Sets gravitino mass
  void setM32(double);
  /// Sets whole thing equal to another object
  void setSoftPars(SoftParsMssm const &);
  /// Set one element of a soft mass squared matrix
  void setSoftMassElement(softMasses, int, int, double);
  /// Set whole of a soft mass squared matrix
  void setSoftMassMatrix(softMasses, const DoubleMatrix &);
  /// Set whole of a trilinear SUSY breaking parameter matrix
  void setTrilinearMatrix(trilinears, const DoubleMatrix &);
  /// Set one element of a trilinear SUSY breaking parameter matrix
  void setTrilinearElement(trilinears k, int i, int j, double a);
  /// Set one gaugino mass
  void setGauginoMass(int, double);
  /// Set all gaugino masses
  void setAllGauginos(const DoubleVector &);
  void setM3Squared(double);  ///< Sets \f$ m_3^2\f$
  void setMh1Squared(double); ///< Sets \f$ m_{H_1}^2\f$
  void setMh2Squared(double); ///< Sets \f$ m_{H_2}^2\f$
  /// Sets total set of RGE parameters equal to elements of a vector
  void set(const DoubleVector &);
  //  void setSusy(const MssmSusy &);

  /// Returns double vector containing numerical beta functions of parameters
  DoubleVector beta() const;
  /// Returns numerical beta functions of parameters
  SoftParsMssm beta2() const;
  /// Returns derivatives of anomalous dimensions of fields with respect to
  /// renormalisation scale in MSSM for: RH leptons, LH leptons, LH quarks, RH
  /// up quarks, RH down quarks, H1 and H2 respectively
  void anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
                      DoubleMatrix & gQQ, DoubleMatrix & gUU,
                      DoubleMatrix & gDD,
                      double & gH1H1, double & gH2H2) const;
  /// Ytilde quantities are for calculational brevity in beta functions.
  void yTildes(DoubleMatrix & yu, DoubleMatrix & yd, DoubleMatrix &ye) const;

  /// Reads in universal boundary conditions at the current scale:
  /// m0, M1/2, A0, B-parameter and mu
  void universal(double m0,  double m12,  double a0,  double mu,
                 double m3sq);
  /// Give it a SUSY object and a value of M3/2, and it will return a soft
  /// object with AMSB soft breaking terms. Note that the sleptons will be
  /// tachyonic, ie nothing has been done to fix that problem.
  /// Note that in the following, we are neglecting all Yukawa couplings except
  /// that of the third family.
  void addAmsb(double m32);
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A0
  void standardSugra(double m0,  double m12, double a0);
  /// Sets all flavour-diagonal SUSY breaking scalar masses to m0
  void universalScalars(double m0);
  /// Sets all flavour-diagonal SUSY breaking gaugino masses to m12
  void universalGauginos(double m12);
  /// Sets all SUSY breaking trilinear couplings to a0
  void universalTrilinears(double a0);
  /// Boundary conditions to be applied at messenger scale for Gauge mediated
  /// SUSY breaking (see hep-ph/9703211 for example), n5 is the number of
  /// 5-plets, mMess is the messenger scale and lambda is the GMSB scale
  void minimalGmsb(int n5, double lambda, double mMess, double cgrav);

  /// *** MIRAGE MEDIATION SOFT TERMS ***
  void mirage(double alphac, double M0, double aq, double al, double ahu, double ahd, double cq, double cl, double chu, double chd);

  /// Reads in soft SUSY breaking parameters from a file
  void inputSoftParsOnly();
};

/// Formatted ouput of whole object
ostream & operator <<(ostream &left, const SoftParsMssm &s);
/// Formatted input of whole object
istream & operator >>(istream &left, SoftParsMssm &s);

inline SoftParsMssm::SoftParsMssm()
  : MssmSusy(), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3),
  mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3), m3sq(0.0),
  mH1sq(0.0), mH2sq(0.0), m32(0.0) {

  setPars(numSoftParsMssm);
  setMu(0.0);
  setLoops(0);
  setThresholds(0);
}

inline SoftParsMssm::SoftParsMssm(const SoftParsMssm & s)
  : MssmSusy(s.displaySusy()),
  mGaugino(s.displayGaugino()), ua(s.displayTrilinear(UA)),
  da(s.displayTrilinear(DA)), ea(s.displayTrilinear(EA)),
  mQLsq(s.displaySoftMassSquared(mQl)),
  mURsq(s.displaySoftMassSquared(mUr)),
  mDRsq(s.displaySoftMassSquared(mDr)),
  mLLsq(s.displaySoftMassSquared(mLl)),
  mSEsq(s.displaySoftMassSquared(mEr)),
  m3sq(s.displayM3Squared()), mH1sq(s.displayMh1Squared()),
  mH2sq(s.displayMh2Squared()), m32(s.displayGravitino()) {

  setPars(numSoftParsMssm);
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}

inline SoftParsMssm::SoftParsMssm(const MssmSusy &s)
  : MssmSusy(s), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3),
    mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3), m3sq(0.0),
    mH1sq(0.0),  mH2sq(0.0), m32(0.0) {
      setPars(numSoftParsMssm);
      setMu(s.displayMu());
      setLoops(s.displayLoops());
      setThresholds(s.displayThresholds());
}

inline SoftParsMssm::SoftParsMssm
(const MssmSusy & s, const DoubleVector & mG, const
 DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE, const
 DoubleMatrix & mQl, const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
 DoubleMatrix & mLl, const DoubleMatrix & mEr, double m3sqn, double mH1sq,
 double mH2sq, double mg, double mu, int l, int t)
  : MssmSusy(s), mGaugino(mG), ua(aU), da(aD), ea(aE),
    mQLsq(mQl), mURsq(mUr), mDRsq(mDr), mLLsq(mLl), mSEsq(mEr), m3sq(m3sqn),
    mH1sq(mH1sq), mH2sq(mH2sq), m32(mg) {
      setPars(numSoftParsMssm);
      setMu(mu);
      setLoops(l);
      setThresholds(t);
}

inline SoftParsMssm SoftParsMssm::displaySoftPars() const { return *this; }

inline double SoftParsMssm::displayM3Squared() const { return m3sq; }

inline double SoftParsMssm::displayMh1Squared() const { return mH1sq; }

inline double SoftParsMssm::displayMh2Squared() const { return mH2sq; }

inline DoubleVector SoftParsMssm::displayGaugino() const { return mGaugino; }

inline double SoftParsMssm::displayGaugino(int i) const {
  return mGaugino.display(i);
}

inline double SoftParsMssm::displayGravitino() const { return m32; }

inline void SoftParsMssm::setGauginoMass(int i, double f) {
  mGaugino(i) = f;
}

inline void SoftParsMssm::setM3Squared(double f) { m3sq = f; }
inline void SoftParsMssm::setMh1Squared(double f) { mH1sq = f; }
inline void SoftParsMssm::setMh2Squared(double f) { mH2sq = f; }
inline void SoftParsMssm::setSoftPars(SoftParsMssm const & s) { *this = s; }
inline void SoftParsMssm::setM32(double a) { m32 = a; }
#endif
