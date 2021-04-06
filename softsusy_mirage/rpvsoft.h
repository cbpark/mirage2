
/** \file rpvsoft.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: Header file for RP violating MSSM object including all (real)
                soft SUSY breaking parameters and (real) SUSY couplings.

*/

#ifndef RPVSOFT_H
#define RPVSOFT_H

#include <iostream>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <def.h>
#include <softsusy.h>
#include <rpvsusypars.h>

/// Number of independent parameters for RGE
static const int numRpvSoftPars = 99 + numSoftParsMssm; 

/// Main class for R-parity violating MSSM. Contains all relevant data,
/// masses, and routines for calculating REWSB and masses etc.
/// Note that \f$ \cos^2 \beta (v_u^2 + v_d^2) = v_u^2 \f$, and 
/// \f$ \sin^2 \beta (v_u^2 + v_d^2) = v_d^2 \f$.
class RpvSoftsusy: public MssmSoftsusy, public RpvSoftPars, public RpvSusyPars
{
private:
  /// Vector of 3 sneutrino vacuum expectation values
  DoubleVector snuVevs;
public:
  /// Default constructor fills object with zeroes
  RpvSoftsusy();
  /// Constructor initialises object equal to another one
  RpvSoftsusy(const RpvSoftsusy &);
  /// All data in object set equal to another one
  const RpvSoftsusy & operator = (const RpvSoftsusy &);

  /// Displays all RGE parameters in a double vector
  DoubleVector display() const;
  /// Sets all RGE parameters from elements of v
  void set(const DoubleVector & v);
  /// Beta functions of RPV MSSM
  DoubleVector beta() const;
  /// Beta functions of RPV MSSM
  RpvSoftsusy beta2() const;

  /// Returns the vacuum expectation values of sneutrinos
  DoubleVector displaySneutrinoVevs() const { return snuVevs; }

  /// Set vacuum expectation values of sneutrinos
  void setSneutrinoVevs(DoubleVector & v) { snuVevs = v; };
  /// Checks to what extent REWSB conditions are satisfied
  void check(const DoubleVector & sneutrinoVevs) const;
  /// Anomalous dimensions of fields in RPV MSSM
  void rpvAnomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
			  DoubleMatrix & gQQ, DoubleMatrix & gUU,
			  DoubleMatrix & gDD, 
			  double & gH1H1, double & gH2H2,
			  DoubleVector & gH1L)  const;

  /// Derivative of anomalous dimensions of fields with respect to
  /// renormalisation scale of RPV parts of RPV MSSM.
  void rpvAnomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
			 DoubleMatrix & gQQ, DoubleMatrix & gUU,
			 DoubleMatrix & gDD, 
			 double & gH1H1, double & gH2H2,
			 DoubleVector & gH1L) const;

  /// Ytilde quantities are for calculational brevity in beta functions. They
  /// are all outputs here. This function returns the RPV parts (only) of the
  /// ytildes
  void rpvyTildes(DoubleMatrix & ye, DoubleMatrix & yd, Tensor & letilde,
		  Tensor & ldtilde, Tensor & lutilde) const;

  /// Performs radiative electroweak symmetry breaking.
  /// IO parameters: sgnMu is +/-1, the sign of mu
  /// mt is the RUNNING DRbar top mass
  void rewsb(int sgnMu, double mt);
  /// Iterative solution to electroweak symmetry breaking boundary conditions.
  /// IO parameters: mu, m3sq are input and output Higgs potential parameters,
  /// iterated along with the sneutrino Vevs. Input sgnMu is the sign of mu,
  /// maxTries is the maximum number of iterations, tol is the desired
  /// fractional accuracy and mt is the DRbar running top mass
  void iterateRewsb(double & mu, double & m3sq, DoubleVector & sneutrinoVevs,
		     int sgnMu, int & numTries, int maxTries, double tol, 
		    double mt);
  /// Returns value of mu consistent with sneutrino vevs given
  /// (which should be a length-3 vector)
  /// IO parameters: sgnMu=+/-1, sign of Higgs potential mu parameter, v1 and
  /// v2 are the two Higgs doublet VEVs and must be input
  double calculateMu(const DoubleVector & sneutrinoVevs, int sgnMu,
		     double v1, double v2);  
  /// Calculates pole Higgs masses and mixings: old feynHiggsfast calculation
  /// IO parameters: piwwt is the W self-energy at M_SUSY, accuracy is number
  /// of loops (0 or 1) to use and pizzt is the Z self-energy at M_SUSY
  virtual void higgs(int accuracy, double piwwtMS, double pizztMS);  
  /// Returns value of m3sq consistent with sneutrino vevs given
  /// (which should be a length-3 vector)
  /// IO parameters: sgnMu=+/-1, sign of Higgs potential mu parameter, v1 and
  /// v2 are the two Higgs doublet VEVs and must be input and snuSq is the sum
  /// of squares of sneutrino VEVs
  double calculateM3sq(const DoubleVector & sneutrinoVevs,
		       double snuSq, double v1, double v2);
  DoubleVector calculateSneutrinoVevs(const DoubleVector & sneutrinoVevs,
				      double tol,
				      double snuSq, double v1, double v2);
  /// Input a set of values for sneutrino VEVs and it returns a more accurate
  /// set - the next step in the iteration. 
  void rotateAwayVevs(DoubleVector & snVevs);

  // Returns some functions of VEVs, gives 0 if there's a problem with the
  // sneutrino VEVs (if they are incompatible with the W and Z masses)
  /// IO parameters: vSM = total VEVs of Higgs+sneutrinos (added in
  // quadrature), v1 and
  /// v2 are the two Higgs doublet VEVs and snuSq is the sum
  /// of squares of sneutrino VEVs: all arguments are outputs
  int usefulVevs(double vSM, const DoubleVector & sneutrinoVevs, 
		 double & snuSq, double & 
		  v1, double & v2) const; 

  /// Theoretical boundary condition upon SUSY breaking terms. SUSY RPV
  /// parameters must be set before this is applied
  /// IO parameters: m0=scalar mass, m12=gaugino mass, a0=trilinear coupling
  void standardSugra(double m0,  double m12, double a0);
  
  /// Returns the 7 by 7 general RPV neutralino mass matrix
  DoubleMatrix neutralinoMassMatrix() const;

  /// Calculates neutralino masses
  void neutralinos() const;
  /// This is used to set both GUT-scale RPV parameters v(i=4,...) and the 
  /// SUSY breaking ones: v(1)=m0, v(2)=m12, v(3)=a0 where 
  /// IO parameters: m0=scalar mass, m12=gaugino mass, a0=trilinear coupling
  void methodBoundaryCondition(const DoubleVector & v);
  /// Calculates DRbar values of gauge and Yukawa couplings from data,
  /// depending upon the spectrum and input value of tan beta=tb
  virtual void sparticleThresholdCorrections(double tb);

  /// Calculates the charged lepton mass matrix once the leptonic Yukawas have
  /// been set. IO parameters: vev is um of slepton+Higgs Vevs in quadrature
  DoubleMatrix chargedLeptons(double vev);
  /// Iterates the charged lepton mass matrix in order to get pole lepton
  /// masses and sneutrino VEVs etc correct.
  /// IO: yeOld is input and a better approximation to it is output.
  /// vev is sum of slepton+Higgs Vevs in quadrature
  /// tol=desired fractional accuracy, maxTries=maximum number of iterations
  /// after which it bails out with zeroes in yeOld and an error flag in
  /// err!=0. mtau=the running DRbar tau mass
  void iterateChargedLeptons(double vev, DoubleMatrix & yeOld, double tol, 
			     int maxTries, int & err, double mtau);

  /// Writes a file for input to ISAWIG, which will write HERWIG information
  /// on a file called herwigInputFile, ISAJET output in isajetOutputFile.
  /// The input file for ISAWIG is called softOutputFile.
  void isawigInterface764(char herwigInputFile [80], 
			  char isajetOutputFile [80],
			  char softOutputFile [80]) const; 
};

/// Formatted input
ostream & operator <<(ostream &left, const RpvSoftsusy & r);

/// Sums up neutrino masses valid for cosmological bound
double neutrinoSum(const RpvSoftsusy & r);

inline RpvSoftsusy::RpvSoftsusy()
  : MssmSoftsusy(), RpvSoftPars(), RpvSusyPars(), snuVevs(3) {
      setPars(numRpvSoftPars);
      setMu(0.0);
      setLoops(0);
      setThresholds(0);
}

inline RpvSoftsusy::RpvSoftsusy(const RpvSoftsusy & s)
  : MssmSoftsusy(s.displayMssmSoft()), RpvSoftPars(s.displayRpvSoft()), 
  RpvSusyPars(s.displayRpvSusy()), snuVevs(s.displaySneutrinoVevs()) {
  // set new parameters here

    setPars(numRpvSoftPars);
    setMu(s.displayMu()); 
    setLoops(s.displayLoops());
    setThresholds(s.displayThresholds());
}

#endif



