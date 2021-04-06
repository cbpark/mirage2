
/** \file softsusy.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: Header file for RP conserving MSSM object including all 
              (real) soft SUSY breaking parameters and (real) SUSY couplings.

   $Log: softsusy.h,v $
   Revision 1.5  2006/04/11 13:57:40  allanach
   Better comments in main programs and cleaned up bug in QEWSB non-usage

   Revision 1.4  2006/01/17 16:40:17  allanach
   Updated to 2.0.4: Higgs mass calculation corrected (problem in the charged
   Higgs contributions)

   Revision 1.7  2005/11/09 14:07:38  allanach
   Added default option to not calculate fine tuning wrt top Yukawa in fineTuning
   calcultion.

   Revision 1.6  2005/09/16 18:16:07  allanach
   Added a function which sets a vector of Higgs masses and couplings.

   Revision 1.5  2005/09/16 16:50:19  allanach
   Added comments to SUSY Les Houches Accord output methods

   Revision 1.4  2005/08/01 11:31:12  allanach
   Added Markus' NLSP routines

   Revision 1.3  2005/07/26 11:00:42  allanach
   Added method nlsp to object MssmSoftsusy by Bernhardt

   Revision 1.2  2005/07/15 15:10:06  allanach
   Added a routine to calculate sin^2 theta_eff

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.65  2004/09/28 10:43:12  allanach
   Name-change to les Houches function - overloaded standard and file output
   functions.

   Revision 1.64  2004/09/28 10:39:07  allanach
   Added Les Houches output function for output to a file

   Revision 1.63  2004/04/20 13:55:08  allanach
   Calculation of tadpoles split up to allow calculation at different scales.
      The situation with mA^2(MZ) < 0 is handled by putting mApole and pretending
      that's the DRbar mass.

   Revision 1.62  2004/03/21 20:43:05  allanach
   Added alternative electroweak symmetry breaking conditions. Added possibility
   of EWSB=SUSY breaking boundary condition scale in SUSY Les Houches Accord.

   Revision 1.61  2004/03/17 13:23:06  allanach
   Solved 1-loop problem in alternative EWSB conditions. 2-loop problem remains.

   Revision 1.60  2004/02/12 14:17:24  allanach
   Added tau^2 corrections and ht hb to higgs masses and Higgs vev calculation

   Revision 1.59  2004/01/28 14:10:38  allanach
   Bug-fixes: for calculation sin theta_w, in the W and Z self-energies, the pole
   top mass is used in the top loop since the BPMZ SM two-loop corrections have
   assume that. Some errors were untrapped before, which is changed and gluino
   mass is allowed to be negative (as it is eg in AMSB).

   Revision 1.58  2004/01/15 13:54:55  allanach
   New heaer style implemented

   Revision 1.57  2003/12/02 19:15:53  allanach
   Added 3-family trilinear information, and added non-universal input

   Revision 1.56  2003/11/06 14:41:47  allanach
   Les Houches interface changed to incorporate running mb

   Revision 1.55  2003/10/27 15:50:39  allanach
   Taken old Higgs routine out of MssmSoftsusy, but added it to RpvSoft (there
   are currently problems using the new routine there -- needs to be fixed)

   Revision 1.52  2003/10/24 10:07:42  allanach
   Any calculation of sin theta_w (DRbar) is now made from the EW gauge
   couplings 

   Revision 1.51  2003/08/19 14:38:59  allanach
   Altered so that mx in arguments to lowOrg is unchanged. It's used as the
   initial guess and lowOrg returns a double number as the calculated mgut.

   Revision 1.50  2003/08/19 14:26:22  allanach
   Changing lowOrg to be more sensible about gauge unification. Should now be
   called with POSITIVE mgut and a flag for gauge unification.

   Revision 1.49  2003/08/12 15:19:54  allanach
   Bug-fixed problem flags, iterateMu. Made iterateMu less precise (was 
   sometimes causing problems)

   Revision 1.48  2003/08/04 16:47:39  allanach
   Added function maxMass to return mass of heaviest superparticle, efficiency
   savings and 2-loop squark/gluino corrections to mt included

   Revision 1.47  2003/07/30 13:39:02  allanach
   Corrected bug in piAA. Started to add piH+H-.

   Revision 1.46  2003/07/28 12:11:37  allanach
   More error trapping, and rearranging rpvsoftsusy to use correct Higgs VEV
   (which is sometimes called at MZ)

   Revision 1.45  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.44  2003/07/24 14:55:28  allanach
   Implemented les Houches input and output properly in the usual command-line
   interface

   Revision 1.43  2003/07/22 10:26:01  allanach
   gdL and gdR used in loop corrections are also 1-loop values now.
   This completes the DR-bar-ness of all 1-loop corrections in SOFTSUSY

   Revision 1.42  2003/07/22 09:20:04  allanach
   Added display functions for running MW/MZ

   Revision 1.41  2003/07/21 16:08:55  allanach
   displayHiggsVevMs now only called at MSUSY, and is checked @ that scale.

   Revision 1.40  2003/07/21 14:00:18  allanach
   MZ fully implemented as an input now. Kept MZ as the central PDG 2002 value,
   for defaults etc

   Revision 1.39  2003/07/18 15:39:14  allanach
   Added prediction of MW to definition of a softsusy object

   Revision 1.38  2003/07/18 14:39:20  allanach
   Implemented MW as a global variable (in preparation for predicting it),
   and also speed corrections in getVev and rhohat: allowing input of
   self-energies to remove their calculation several times

   Revision 1.37  2003/07/16 11:07:06  allanach
   Changed isajet number to 764

   Revision 1.35  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.34  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.22  2002/11/27 18:19:16  allanach
   New higgs mass calculation and bugfix in piZZ (included neutrino self
   energy diagrams) 

   Revision 1.21  2002/11/26 13:25:59  allanach
   New mb calculation

   Revision 1.20  2002/11/19 17:00:00  allanach
   Added FCNC routine

   Revision 1.18  2002/10/22 13:12:10  allanach
   Introduced new problem flag for infra-red quasi fixed points

   Revision 1.17  2002/10/14 14:17:30  allanach
   Added "runto" command in softpoint.x to get micromegas inputs to a
   different scale

   Revision 1.16  2002/10/10 13:37:12  allanach
   Added physical parameters to Micromegas output

   Revision 1.15  2002/10/09 17:18:51  allanach
   Added new softpoint option to specify mu and mA instead of mH1, mH2. Some
   fine-tuning still necessary

   Revision 1.14  2002/09/20 15:37:10  allanach
   Adding quark-mixing routines

   Revision 1.12  2002/09/04 13:59:45  allanach
   Added gauge unification possibility

   Revision 1.11  2002/08/12 16:29:49  allanach
   display and set whole-object routines added

   Revision 1.10  2002/07/30 12:57:32  allanach
   SOFTSUSY1.5

   Revision 1.9  2002/06/14 16:26:30  allanach
   Switches included for 2-loop running of scalar masses, and calulating mt at
   mt. 

   Revision 1.7  2002/04/26 15:14:44  allanach
   Deleted all translation routines and defined boundary conditions within
   softsusy.h and softsusy.cpp

   Revision 1.6  2002/04/14 13:50:41  allanach
   Now use V=m3^2 H1 H2 instead of V=mu B H1 H2. It's more natural!

   Revision 1.5  2002/04/12 16:51:27  allanach
   Added display/set functions to work automatically

   Revision 1.4  2002/04/12 06:24:50  allanach
   Code maintenance - returning a subobject made simpler

   Revision 1.11  2002/02/04 15:07:22  allanach
   Latest version

   Revision 1.10  2001/10/31 09:12:08  allanach
   Altered so that an early no rho convergence can disappear

   Revision 1.9  2001/10/23 12:53:58  allanach
   Easier user interface utilised

   Revision 1.8  2001/10/03 13:34:16  allanach
   Changed name of complex.h to avoid conflict with STD libraries

   Revision 1.7  2001/09/28 13:59:58  allanach
   Split rhohat determination up into sensible pieces - small bugs fixed there.

   Revision 1.6  2001/07/30 14:08:47  allanach
   Added ISAWIG and SSRUN interface

   Revision 1.5  2001/07/26 14:20:34  allanach
   Added isajet interface

   Revision 1.4  2001/07/18 15:54:50  allanach
   Put MIXING switch into def.h

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

/** \mainpage Detailed SOFTSUSY Documentation 

    \section install Installation or downloads
    For installation instructions or a download, please go to the 
    <a href="http://allanach.home.cern.ch/allanach/softsusy.html">
    SOFTSUSY Homepage</a>

    \section manual Official manual
    If you use SOFTSUSY to write a paper, please cite 
    <a href="http://xxx.soton.ac.uk/abs/hep-ph/0104145">
    B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145, 
    </a> which is the SOFTSUSY manual.
    
    \section documentation Documentation
    These web-pages contain the documentation of the latest SOFTSUSY code.
    There are class diagrams and cross-referenced links a la doxygen to help 
    you navigate.

    \section updates Official Updates
    Updates will be posted on the    
    <a href="http://allanach.home.cern.ch/allanach/softsusy.html">
    SOFTSUSY Homepage</a>, and the 
    <a href="http://xxx.soton.ac.uk/abs/hep-ph/0104145">manual</a>
    will also be updated.
 */

#ifndef SOFTSUSY_H
#define SOFTSUSY_H

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::istream;
using std::ostream;
#include <fstream>
using std::fstream;
using std::ios;
#include <sstream>
using std::ostringstream;
#include <string>
using std::string;
#include <cstdlib>
using std::abs;
#include <cmath>
using std::fabs;
#include <fstream>
#include "def.h"
#include "utils.h"
#include "numerics.h"
#include "physpars.h"
#include "lowe.h"
#include "softpars.h"
#include "twoloophiggs.h"

class MssmSoftsusy; 
std::istream & operator >>(std::istream &left, MssmSoftsusy &s);

/// Contains all supersymmetric MSSM parameters, and is basically the SOFTSUSY
/// kernel that incorporates everything
/// - Soft terms
/// - DRbar masses
/// - RGEs
/// - Higgs VEVs, tadpoles
/// - Info on any problems with the parameter space point
/// - main routines for driving the calculation
/// - loop corrections for sparticle masses and Standard Model couplings
/// - fine-tuning and flavour changing neutral currents calculation
class MssmSoftsusy: public SoftParsMssm {
private:
  sPhysical physpars; ///< Contains pole masses and mixings of sparticles
  drBarPars forLoops; ///< Contains DRbar tree-level masses
  sProblem problem;   ///< Contains problem flags 
  double msusy;       ///< Scale at which Higgs potential is minimised
  double minV;        ///< Value of Higgs potential at minimum
  double mw;          ///< Pole W mass prediction
  QedQcd dataSet;     ///< contains low energy data on quark masses etc

public:
  /// Default constructor fills object with zeroes
  MssmSoftsusy();
  /// Constructor sets SUSY parameters only from another object
  MssmSoftsusy(const MssmSusy &);
  /// Constructor copies another object
  MssmSoftsusy(const MssmSoftsusy &);
  /// Sets all parameters from s, sp, mu is the mu superpotential parameter, l
  /// is the number of loops used for RG evolution, t is the thresholds
  /// accuracy parameter, mg is the gravitino mass, hv is the Higgs VEV
  /// parameter. 
  MssmSoftsusy(const SoftParsMssm & s, const sPhysical & sp, double mu, int l, 
	       int t, double hv);
  /// Set all data in the object equal to another
  const MssmSoftsusy & operator=(const MssmSoftsusy & s);
  
  /// Displays whole object as a const
  inline MssmSoftsusy displayMssmSoft() const;
  /// Displays physical parameters only
  inline sPhysical displayPhys() const;
  /// Displays tree-level masses and mixings of sparticles and third
  /// generation fermions
  inline drBarPars displayDrBarPars() const;
  /// Returns any problem flags associated with the object
  sProblem displayProblem() const {return problem; };
  /// Gives the low energy Standard Model data set used for the object
  inline QedQcd displayDataSet() const;
  double displayMinpot() const;    ///< Returns minimum of Higgs potential
  double displayMsusy() const; ///< Returns Higgs minimisation scale
  double displayMw() const; ///< Returns predicted pole MW
  /// Returns DRbar MW, must be calculated by calcDrBarPars first
  double displayMwRun() const; 
  /// Returns DRbar MZ, must be calculated by calcDrBarPars first
  double displayMzRun() const; 
  double displayTadpole1Ms() const; ///< displays t_1/v_1 tadpole
  double displayTadpole2Ms() const; ///< displays t_2/v_2 tadpole
  /// Returns object as a const
  MssmSoftsusy displaySoftsusy() const { return *this; }
  /// Returns value of pole MZ being used
  double displayMz() const { return displayDataSet().displayMu(); }
  
  /// Flags Infra-red quasi fixed point breach problem
  void flagIrqfp(bool a) { problem.irqfp = a; };
  /// Flags non-perturbative RG evolution
  void flagNonperturbative(bool a) { problem.nonperturbative = a; };
  /// Flags a negative-mass squared scalar (really a CCB problem)
  void flagTachyon(bool a) { problem.tachyon = a; };
  /// Flags problem with Higgs potential minimum
  void flagB(bool a) { problem.b = a; };
  /// Flags a really bad convergence: nowhere near a solution
  void flagBadConvergenve(bool a) {problem.badConvergence = a; };
  /// Flags fact that calculation hasn't acheived required accuracy
  void flagNoConvergence(bool a) { problem.noConvergence = a; };
  /// Flags fact that mu sub iteration didn't converge
  void flagNoMuConvergence(bool a) { problem.noMuConvergence = a; };
  /// Flags fact that rho parameter sub iteration didn't converge
  void flagNoRhoConvergence(bool a) { problem.noRhoConvergence = a; };
  /// Flags point inconsistent with electroweak symmetry breaking
  void flagMusqwrongsign(bool a) { problem.muSqWrongSign = a; };
  /// Flags an inconsistent Higgs minimum
  void flagHiggsufb(bool a) { problem.higgsUfb = a; };
  /// Sets all problems equal to either true or false (contained in a)
  void flagAllProblems(bool a) { problem.irqfp = a; 
  problem.tachyon = a; problem.b = a; problem.badConvergence = a;
  problem.noConvergence = a; problem.higgsUfb = a;
  problem.nonperturbative = a; problem.noRhoConvergence = a; 
  problem.noMuConvergence = a; problem.muSqWrongSign = a; }

  /// Sets whole object equal to another  
  void setSoftsusy(const MssmSoftsusy & s) { *this = s; };
  /// Sets low energy Standard Model fermion mass and gauge coupling data
  void setData(const QedQcd & r) { dataSet = r; };
  /// Sets potential value at minimum of Higgs potential
  void setMinpot(double);
  /// Sets scale of Higgs potential minimisation and sparticle mass calculation
  void setMsusy(double);
  /// sets pole MW prediction
  void setMw(double);
  /// Sets all physical parameters
  void setPhys(const sPhysical & s) { physpars = s; };
  /// Sets tree-level DRbar parameters
  void setDrBarPars(const drBarPars & s) { forLoops = s; };
  /// Does the full 2-loop calculation of both tadpoles and sets them
  void doTadpoles(double mt, double sinthDRbar);
  /// Does the calculation of one-loop pieces of \f$ t_1 / v_1 \f$ 
  double doCalcTadpole1oneLoop(double mt, double sinthDRbar);
  /// Does the calculation of one-loop pieces of \f$ t_2 / v_2 \f$ 
  double doCalcTadpole2oneLoop(double mt, double sinthDRbar);
  /// Calculates and sets the one-loop pieces of \f$ t_1 / v_1 \f$ 
  virtual void calcTadpole1Ms1loop(double mt, double sinthDRbar);
  /// Calculates then sets the one-loop pieces of \f$ t_2 / v_2 \f$ 
  virtual void calcTadpole2Ms1loop(double mt, double sinthDRbar);
  /// Adds one-loop corrections to stop mass matrix
  /// IO parameters: p=external momentum, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  void addStopCorrection(double p, DoubleMatrix & mass, double mt);
  /// Adds one-loop corrections to sbottom mass matrix at p=root(mb1 mb2)
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  void addSdownCorrection(DoubleMatrix & mass, int family);
  /// Adds one-loop corrections to sbottom mass matrix at p=root(mb1 mb2)
  /// IO parameters: p=external momentum scale, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  void addSbotCorrection(double p, DoubleMatrix & mass, double mb);
  /// Adds one-loop corrections to sel_fam mass matrix at p=root(msel1 msel2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  void addSlepCorrection(DoubleMatrix & mass, int family);
  /// Adds one-loop corrections to stau mass matrix at p=root(mtau1 mtau2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  void addStauCorrection(DoubleMatrix & mass, double mtau);
  /// Adds one-loop corrections to stau mass matrix at p=root(mtau1 mtau2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  void addSupCorrection(DoubleMatrix & mass, int family);
  /// Adds one-loop corrections to tau sneutrino mass 
  /// IO parameters: p=external momentum, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  void addSnuTauCorrection(double & mass);
  /// Adds one-loop corrections to sneutrino mass of family "family"
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added
  void addSnuCorrection(double & mass, int family);
  /// Adds approximate one-loop corrections to squark mass matrix for first
  /// two families.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added
  void addSquarkCorrection(DoubleMatrix & mass);
  /// Organises calculation of all up squark masses.
  /// IO parameters: mt=DRbar top mass, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w, accuracy=0,1: number of loops
  /// to add to tree-level squark mass matrix
  void doUpSquarks(double mt, double pizztMS, double sinthDRbarMS, int
		      accuracy); 
  /// Organises calculation of all down squark masses.
  /// IO parameters: mb=DRbar bottom mass, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w, accuracy=0,1: number of loops
  /// to add to tree-level squark mass matrix
  void doDownSquarks(double mb, double pizztMS, double sinthDRbarMS, int
		      accuracy, double mt);
  /// Organises calculation of all slepton masses.
  /// IO parameters: mT=DRbar tau mass, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w, accuracy=0,1: number of loops
  /// to add to tree-level squark mass matrix
  void doChargedSleptons(double mT, double pizztMS, double sinthDRbarMS, int
		      accuracy);
  /// Organises calculation of all sneutrino masses, pizztMS=Z self energy at
  /// Q=M_SUSY
  void doSnu(double pizztMS, int accuracy = 0);
  /// Returns tree-level up squark mass matrix in "mass".
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mtrun=DR bar top
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
  void treeUpSquark(DoubleMatrix & mass, double mtrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// Returns tree-level down squark mass matrix in "mass".
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mbrun=DR bar bottom
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
  void treeDownSquark(DoubleMatrix & mass, double mbrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// Returns tree-level down squark mass matrix in "mass".
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mTrun=DR bar tau
  /// mass, family=generation of slepton, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w
  void treeChargedSlepton(DoubleMatrix & mass, double mTrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// Organises calculation of all sneutrino masses, pizztMS=Z self energy at
  /// Q=M_SUSY, mSq=mass of sneutrino, family=generation of sneutrino
  void treeSnu(double & mSq, double pizztMS, int family);

  /// Organises tree-level calculation of all sparticle masses and mixings
  void calcDrBarPars();
  /// For an input tan beta=tb, sets gauge and Yukawa couplings according to
  /// the tree-level spectrum and data set 
  virtual void sparticleThresholdCorrections(double tb);
  /// Does SUSY (and other) threshold corrections to alphaEm - returns alpha in
  /// DRbar scheme at scale Q. From hep-ph/9606211. Input empirical value of
  /// alpha at 0 momentum....
  double qedSusythresh(double alphaEm, double Q) const;
  /// Does SUSY (and other) threshold corrections to alphaS
  /// Input alphas in MSbar and it returns it in DRbar scheme. 
  /// From hep-ph/9606211
  double qcdSusythresh(double alphasMSbar, double Q) const;
  /// Calculates the best scale at which to do symmetry breaking:
  /// \f$ M_{SUSY}=Q_{EWSB} \sqrt{m_{{\tilde t}_1 {\tilde t}_2}} \f$. 
  /// Should only be called after calcDrBarPars.
  double calcMs() const;
  /// Calculates physical sparticle masses to accuracy number of loops. Should
  /// be called at M_{SUSY}.
  virtual void physical(int accuracy);
  /// Applies 1-loop SUSY and 2-loop QCD corrections to pole mt in order to
  /// return the DRbar running value at the current scale
  double calcRunningMt();
  /// Applies approximate 1-loop SUSY corrections to pole mtau in order to
  /// return the DRbar running value at the current scale
  double calcRunningMtau() const;  
  /// Applies approximate 1-loop SUSY corrections to mb(MZ) in order to
  /// return the DRbar running value
  double calcRunningMb() const;
  /*
  /// Calculates top Yukawa coupling, supply Higgs vev parameter at current
  /// scale 
  double calcHt(double vev);
  /// Calculates bottom Yukawa coupling, supply Higgs vev parameter at current
  /// scale 
  double calcHb(double vev) const;
  /// Calculates tau Yukawa coupling, supply Higgs vev parameter at current
  /// scale 
  double calcHtau(double vev) const;
  */
  /// Calculates DRbar sin theta_w at the current scale from gauge couplings 
  double calcSinthdrbar() const;
  /// Calculates Higgs VEV parameter from gauge couplings and MZ
  double getVev() const;
  /// Input for this one (saves time, possibly) is to give the self-energy of
  /// the Z at the current scale
  double getVev(double pizzt) const;
  /// Calculates pole chargino masses and mixing using approximate 1-loop SUSY
  /// corrections. IO parameters: piwwt is the W self-energy at the current,
  /// accuracy is the number of loops required (0 or 1 currently)
  virtual void charginos(int accuracy, double piwwt);
  /// Adds the loop corrections on to an input tree-level chargino mass
  virtual void addCharginoLoop(DoubleMatrix &);
  /// Calculates pole neutralino masses and mixingusing approximate 1-loop SUSY
  /// corrections. IO parameters: piwwt is the W self-energy at M_SUSY,
  /// accuracy is the number of loops required (0 or 1 currently), pizzt is
  /// the Z self-energy at M_SUSY
  virtual void neutralinos(int accuracy, double piwwt, double pizzt);
  /// Adds the loop corrections on to an input tree-level neutralino mass
  virtual void addNeutralinoLoop(DoubleMatrix &);
  void addNeutralinoLoopNew(DoubleMatrix &); // DEBUG
  /// Calculates pole gluino mass to 1-loop SUSY corrections
  virtual void gluino(int accuracy);
  /// Calculates pole Higgs masses and mixings: full 1-loop SUSY corrections
  /// and 2-loop alpha_t (alpha_s + alpha_t) + alpha_s alpha_b effective
  /// potential corrections. 
  /// IO parameters: piwwt is the W self-energy at M_SUSY, accuracy is number
  /// of loops (0 or 1) to use and pizzt is the Z self-energy at M_SUSY
  virtual void higgs(int accuracy, double piwwt, double pizzt);
  /// Calculates pole Higgs masses and mixings: full 1-loop SUSY corrections
  /// and 2-loop alpha_t (alpha_s + alpha_t) + alpha_s alpha_b effective
  /// potential corrections. 
  /// IO parameters: piwwt is the W self-energy at M_SUSY, accuracy is number
  /// of loops (0 or 1) to use and pizzt is the Z self-energy at M_SUSY
  //  virtual void newhiggs(int accuracy, double piwwt, double pizzt);
  /// Tree-level REWSB calculation, returning mu at correct value. sgnMu is the
  /// required sign (+/- 1). Returns 1 if mu^2<0, indicating an inconsistent
  /// minimum 
  virtual int rewsbMu(int sgnMu, double & mu) const;
  /// Tree-level REWSB calculation, returning m3sq at correct value consistent
  /// with mu
  virtual int rewsbM3sq(double, double &) const;
  /// Organises high accuracy rewsb: call it at the low scale M_{SUSY}
  /// IO parameters: sgnMu is +/-1 (desired sign of mu), mt is DRbar top mass
  virtual void rewsb(int sgnMu, double mt);
  /// Organises tree-level rewsb: call it at the low scale M_{SUSY}
  /// IO parameters: sgnMu is +/-1 (desired sign of mu)
  virtual void rewsbTreeLevel(int sgnMu);
  /// Obtains solution of one-loop effective potential minimisation via
  /// iteration technique. Currently includes: all 1-loop SUSY tadpoles, plus
  /// 2-loop alpha_t (alpha_t + alpha_s) + alpha_b alpha_s corrections
  /// IO parameters: 
  /// munew=current value of mu for this iteration, sgnMu=desired sign of mu,
  /// mt=DRbar mtop, maxTries=maximum number of iterations before it bails
  /// out, pizzMS=self-energy of MZ at current scale, sinthDRbar=DRbar value
  /// of sin theta_w, tol=desired fractional accuracy on mu, err=error flag:
  /// err=1 if no iteration reached, 2 if incorrect rewsb
  void iterateMu(double & munew, int sgnMu, double mt, 
		  int maxTries, double pizztMS, double sinthDRbar, double tol,
		 int  & err);
  /// This is a check: predicts tan beta from the values of soft parameters
  /// and mu that we have
  double predTanb() const;
  /// Predicts value of MZ(pole) from values of soft parameters and mu that we
  /// have. tanb=tan beta is also predicted
  double predMzsq(double & tanb) const;
  /// Calculates fine-tuning for soft parameters and mu, m_3^2, top Yukawa. 
  /// IO parameters: bcPars 
  /// should be a vector giving the high-scale SUSY breaking boundary
  /// condition parameters, MX is the high-scale, boundaryCondition is the
  /// user-supplied function that sets the SUSY breaking BCs.
  /// If doTop is true, it also calculates the fine tuning associated with the
  /// top Yukawa coupling. 
  DoubleVector fineTune(void (*boundaryCondition)(MssmSoftsusy &, 
						  const DoubleVector &), 
			const DoubleVector & bcPars, double MX, 
			bool doTop = false);
  /// Give it a GUT scale object consistent with rewsb
  /// and it'll return the fine tuning of one parameter specified by numPar
  /// ht, mu  and m3sq at the high 
  double it1par(int numPar, const DoubleVector & bcPars);
  /// Input mx the scale up to which you search for minima.
  /// Returns minimum value of potential along UFB3 direction.
  /// Does ufbs truly properly but takes a long time.
  double ufb3sl(double);
  /// You should evaluate this at a scale MSusy average of stops.
  /// Returns depth of electroweak minimum
  double realMinMs() const;
  
  /// Calculates transverse part of Z self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piZZT(double p, double Q, bool usePoleMt = false) const;
  /// Calculates transverse part of W self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piWWT(double p, double Q, bool usePoleMt = false) const;
  /// Calculates transverse part of H^+H^- self-energy:  
  /// for p=external momentum, Q=renormalisation scale
  virtual double piHpHm(double p, double Q) const;
  /// Calculates Z gamma self-energy:  
  /// for p=external momentum, Q=renormalisation scale
  double piZGT(double p, double Q) const;
  /// Calculates transverse part of A^0 self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piAA(double p, double Q) const;
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s1(double p, double q) const;
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s2(double p, double q) const;
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis2s2(double p, double q) const;
  /// Calculates sin^2 theta^l_eff
  double sinSqThetaEff();

  /// Iterative determination of rho parameter consistent with muon decay
  /// constant, MZ and alpha_0. 
  /// IO parameters: outrho and outsin are the current DRbar values of sin
  /// theta_w and the rho parameter. alphaMZDRbar=already calculated value of
  /// fine-structure constant, pizztMZ=Z self-energy at Q=MZ, tol=desired
  /// piwwt0=self-energy of the W at p=0, piwwtMW=self-energy of the W at p=MW
  /// accuracy and maxTries is the maximum number of iterations the routines
  /// will allow
  virtual void rhohat(double & outrho, double & outsin, double alphaMZDRbar,
		      double pizztMZ, double piwwt0, double piwwtMW, 
		      double tol, int maxTries);
  /// Calculates delta_v corrections for outrho=DRbar rho parameter,
  /// outsin=DRbar sin theta_w, alphaDRbar=alpha(Q) in the DR bar scheme,
  /// pizztMZ=self-energy of the Z at MZ
  double deltaVb(double outrho, double outsin, double alphaDRbar, 
		double pizztMZ) const;
  /// Calculates delta rho corrections for outrho=DRbar rho parameter,
  /// outsin=DRbar sin theta_w, alphaDRbar=alpha(Q) in the DR bar scheme,
  /// pizztMZ=self-energy of the Z at MZ, piwwtMW=self-energy of the W at p=MW
  double dRho(double outrho, double outsin, double alphaDRbar, 
		 double pizztMZ, double piwwtMW);
  /// Calculates delta r corrections for outrho=DRbar rho parameter,
  /// outsin=DRbar sin theta_w, alphaDRbar=alpha(Q) in the DR bar scheme,
  /// pizztMZ=self-energy of the Z at p=MZ, pizzt0=self-energy of the W at p=0
  double dR(double outrho, double outsin, double alphaDRbar, double pizztMZ,
	    double piwwt0);
  /// Returns the mass of the heaviest SUSY particle, excluding gravitino
  double maxMass() const;
  /// Returns lsp mass in mass and function return labels which particle is 
  /// lsp=0 means LSP is neutralino, 1=up squark, 2=down squark, 3=sleptons,
  /// 4=charginos, 5=sneutrinos, 6=gluino. posi, posj give the 
  /// "handedness" (ie 1 or 2) for scalars and family respectively. 
  int lsp(double & mass, int & posi, int & posj) const;
  /// Returns nlsp mass in mass and function return labels which particle is 
  /// nlsp=0 means NLSP is neutralino, 1=up squark, 2=down squark, 3=sleptons,
  /// 4=charginos, 5=sneutrinos, 6=gluino. posi, posj give the 
  /// "handedness" (ie 1 or 2) for scalars and family respectively.
  int nlsp(double & mass, int & posi, int & posj) const;
  
  /// Prints a list of important sparticle/Higgs masses to standard output
  void printShort() const;
  /// Prints a list of all sparticle/Higgs masses to standard output
  void printLong();
  /// Prints whols object to standard output
  virtual void printObj() { cout << *this; };
  
  /// log(max(a^2, b^2, c^2) / Q^2)
  double thet(double a, double b, double c);
 
  /// Driver calculation to determine all sparticle masses and parameters.
  /// Returns low energy softsusy object consistent with BC's m0 etc at MGUT.
  /// oneset should be at MZ and contains the SM data to fit the model to.
  /// If the running comes into difficulty, eg if a Landau pole is reached, it
  /// returns a ZERO object: no result is possible!
  /// Boundary condition is the theoretical condition on parameters at the high
  /// energy scale mx: the parameters themselves are contained within the
  /// vector. IO parameters:
  /// sgnMu is the desired sign of mu: + or - 1. If mu is 0, mu is set
  /// initially as a boundary condition. tanb = desired value of DR bar tan
  /// beta(MZ). boundaryCondition is the function which sets to SUSY BCs at
  /// the high scale, mxGuess is the GUT scale
  /// gaugeUnification=true if the user requests true gauge unification, in
  /// which case mxGuess is the first guess. Returns actual GUT scale.
  /// ewsbBCscale = false (or omitted) means that the boundary condition on
  /// SUSY breaking is set in the usual way. If it is true, the boundary
  /// condition is set to \f$\sqrt{m_{\tilde t}_1 m_{\tilde t}_2} \f$, ie like
  /// in the "pheno MSSM".
  double lowOrg(void (*boundaryCondition)
		(MssmSoftsusy &, const DoubleVector &),
		double mxGuess, 
		const DoubleVector & pars, int sgnMu, double tanb,
		const QedQcd & oneset, bool gaugeUnification, 
		bool ewsbBCscale =  false); 
  /// Main iteration routine: 
  /// Boundary condition is the theoretical condition on parameters at the high
  /// energy scale mx: the parameters themselves are contained within the
  /// vector. IO parameters:  
  /// maxTries is the maximum number of iterations allowed, mx is the GUT
  /// scale (negative if you require gauge unification),
  /// sgnMu is the desired sign of mu: + or - 1. If mu is 0, mu is set
  /// initially as a boundary condition. tanb = desired value of DR bar tan
  /// beta(MZ).
  void itLowsoft(int maxTries, double & mx, int sgnMu, double tol, 
		 double tanb, void (*boundaryCondition)(MssmSoftsusy &, 
							const DoubleVector &), 
		 const DoubleVector & pars, bool gaugeUnification, 
		 bool ewsbBCscale);
  
  /// Dummy function to allow users to re-define it in user supplied
  /// objects. Pars will contain necessary parameters to describe high-scale
  /// boundayr conditions on SUSY breaking terms
  virtual void methodBoundaryCondition(const DoubleVector & pars);

  /// Works out how best to fit the isajet numbers to the spectrum.
  /// There are problems with the Higgs and sbottoms because ISAJET assumes
  /// certain tree-level relations between masses that are broken by SOFTSUSY's
  /// higher accuracy. The differences get large for high tan beta around 50,
  /// at around 10 they're typically only a percent.
  /// Output parameters: mtopPole is the pole mass of the top quark, mGPole is
  /// the gluino pole mass, smu is the superpotential Higgs parameter, mA is
  /// the pseudoscalar Higgs mass, tanb is tan beta, mq1l is the mass of the
  /// LH first family squark, mdr is the mass of the RH first family squark,
  /// meL is the LH selectron mass, meR is the RH selectron mass, mql3 is the
  /// LH stop mass, mdr3 is the RH sbottom mass, mur3 is the RH stop mass,
  /// mtauL is the LH stau mass, mtauR is the RH stau mass, at is the stop
  /// trilinear term, ab is the sbottom trilinear term, atau is the stau
  /// trilinear term, mq2l, msr, mcr are 2nd family squark masses, mmuL, mmuR
  /// are the smu masses, m1 and m2 are the gaugino mass parameters.
  void isajetNumbers764
  (double & mtopPole, double & mGPole, double & smu, double & mA, 
   double & tanb, double & mq1l, double & mdr, double & mur, double & meL, 
   double & meR, double & mql3, double & mdr3, double & mur3, double &  mtauL, 
   double & mtauR, double & at, double & ab, double & atau, double & mq2l, 
   double & msr, double & mcr, double & mmuL, double & mmuR, double & m1, 
   double & m2) const; 
  /// prints a file into fname which acts as an input to isajet
  void isajetInterface764(char fname[80]) const;
  /// prints a file into fname which acts as an input to isassrun: fstream
  /// should be opened before calling
  void ssrunInterface764Inside(char fname [80], fstream & ) const;
  /// prints a file into fname which acts as an input to isajet
  /// First name input is the name of an OUTPUT file from ssrun, the second
  /// name is the name of the interface file for INPUT to ssrun
  void ssrunInterface764(char fname [80], char softfname [80]) const;
  /// Prints a file into fnamesoft that can be input into isawig. fnamein
  /// gives the ISASSRUN output file name, fnameout is the filename for the
  /// ISAWIG input file
  void isawigInterface764(char fnamein [80], char fnameout [80], 
			  char fnamesoft[80]) const;
  /// Outputs to softoutput a micromegas file (the name of which is input)
  void microMegasInterface(char fname[80]) const;
  /// Outputs with Les Houches accord conventions to standard output. 
  /// Inputs:
  /// model contains what form of model is used for the SUSY breaking terms
  /// (eg sugra, gmsb, amsb, nonUniversal). qMax is only relevant if you want
  /// a gridded output of running parameters up to some scale qMax. Put
  /// numPoints = 1 if you don't want to use this option - then qMaz is
  /// immaterial. mb is mb(mb) in the MSbar scheme used to produce the output,
  /// whereas mtau is the pole mass used (eg 1.777). mgut is the GUT scale
  /// that has been determined, and altEwsb is true if you specified mu and mA
  /// as input parameters (not tan beta and mH1, mH2).
  void lesHouchesAccordOutput(char model[], const DoubleVector & pars, 
			      int sgnMu, double tanb, double qMax, 
			      int numPoints, double mb, 
			      double mtau, double mgut, 
			      bool altEwsb = false);
  /// Outputs with Les Houches accord conventions to a file called
  /// fileName. Inputs: 
  /// model contains what form of model is used for the SUSY breaking terms
  /// (eg sugra, gmsb, amsb, nonUniversal). qMax is only relevant if you want
  /// a gridded output of running parameters up to some scale qMax. Put
  /// numPoints = 1 if you don't want to use this option - then qMaz is
  /// immaterial. mb is mb(mb) in the MSbar scheme used to produce the output,
  /// whereas mtau is the pole mass used (eg 1.777). mgut is the GUT scale
  /// that has been determined, and altEwsb is true if you specified mu and mA
  /// as input parameters (not tan beta and mH1, mH2).
  void lesHouchesAccordOutput(char model[], const DoubleVector & pars, 
			      int sgnMu, double tanb, double qMax, 
			      char fileName[], int numPoints, double mb, 
			      double mtau, double mgut, 
			      bool altEwsb = false);
  /// Inputs a micromegas filename into a MssmSoftsusy object
  void microMegasIn(char fname[80]);

  /// all matrices in inputs should be 3x3: they're output as the quark
  /// FCNC matrices.
  void sCkmRotation
  (DoubleMatrix & deltaULL, DoubleMatrix & deltaURR, DoubleMatrix & deltaULR,  
   DoubleMatrix & deltaDLL, DoubleMatrix & deltaDRR, DoubleMatrix & deltaDLR) 
    const;


  /// Prints out quark FCNC delta parameters (defined by Masiero and co)
  void outputFcncs() const;

  /// Sets the minimum of potential to be the difference between the UFB-3
  /// direction minimum and the standard EW breaking minimum. mgut is
  /// obviously the high-scale at which boundary conditions are employed
  void doUfb3(double mgut);

  /// Utility function: sets Higgs masses of neutral Higgs' (higgsm) and
  /// charged (higgsc). They should be of dimension 4 and 2 respectively.
  /// Also sets couplings dnu(4), dnd(4) and cn(4). beta is from tan beta.
  void assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc, 
		   DoubleVector & dnu, DoubleVector & dnd, 
		   DoubleVector & cn, double beta) const;

  /// Utility function: sets Higgs masses of neutral Higgs' (higgsm) and
  /// charged (higgsc). They should be of dimension 4 and 2 respectively.
  void assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc) const;
};

inline MssmSoftsusy::MssmSoftsusy()
  : SoftParsMssm(), physpars(), forLoops(), 
    problem(), msusy(0.0), minV(6.66e66), 
    mw(0.0), dataSet() { 
      setPars(110);
      setMu(0.0);
      setLoops(0);
      setThresholds(0);
}


inline MssmSoftsusy::MssmSoftsusy(const MssmSoftsusy & s)
  : SoftParsMssm(s.displaySoftPars()), physpars(s.displayPhys()), 
    forLoops(s.displayDrBarPars()), 
    problem(s.problem), msusy(s.msusy), minV(s.minV), 
    mw(s.mw), dataSet(s.displayDataSet()) {
    setPars(110);
    setMu(s.displayMu()); 
    setLoops(s.displayLoops());
    setThresholds(s.displayThresholds());
}

inline MssmSoftsusy::MssmSoftsusy(const MssmSusy &s)
  : SoftParsMssm(s), physpars(), forLoops(), problem(), 
    msusy(0.0), minV(6.66e66), mw(0.0), dataSet() { 
      setPars(110);
      setMu(s.displayMu()); 
      setLoops(s.displayLoops());
      setThresholds(s.displayThresholds());
}

inline MssmSoftsusy::MssmSoftsusy
(const SoftParsMssm & s, const sPhysical & sp, double mu, int l, int t, 
 double hv) 
  : SoftParsMssm(s), physpars(sp), forLoops(), problem(), msusy(0.0),
    minV(6.66e66), mw(0.0), dataSet() {
      setPars(110);
      setMu(mu);
      setLoops(l);
      setThresholds(t);
}

inline MssmSoftsusy MssmSoftsusy::displayMssmSoft() const { return *this; }

inline QedQcd MssmSoftsusy::displayDataSet() const { return dataSet; }

inline sPhysical MssmSoftsusy::displayPhys() const { return physpars; }

inline drBarPars MssmSoftsusy::displayDrBarPars() const { return forLoops; }


inline double MssmSoftsusy::displayMinpot() const { return minV; } 
inline double MssmSoftsusy::displayMsusy() const { return msusy; } 
inline double MssmSoftsusy::displayMw() const { return mw; } 

inline double MssmSoftsusy::displayTadpole1Ms() const {
  return physpars.t1OV1Ms; 
}

inline double MssmSoftsusy::displayTadpole2Ms() const {
  return physpars.t2OV2Ms; 
}

inline void MssmSoftsusy::setMinpot(double f) { minV = f; }
inline void MssmSoftsusy::setMsusy(double f) { msusy = f; }
inline void MssmSoftsusy::setMw(double f) { mw = f; }
inline void MssmSoftsusy::doUfb3(double mgut) { setMinpot(ufb3sl(mgut) -
							  realMinMs()); } 
/// Allows user to specify a boundary condition where ALL SUSY breaking
/// parameters are specified in inputParameters
void generalBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Prints out header line for print-short output
void printShortInitialise();
/// Formatted output
ostream & operator <<(ostream &, const MssmSoftsusy &); 
/// Calculates fractional difference in Drbar masses between in and out
double sumTol(const MssmSoftsusy & in, const MssmSoftsusy & out);
/// returns the square root of the absolute value of the argument
// returns sqrt(f) for f>0 or -sqrt(|f|) for f<0
inline double ccbSqrt(double f){ return sqrt(fabs(f)); }
/// Prints out the identity of the LSP to standard output. 
/// temp, j are defined from lsp function 
void recogLsp(int temp, int j);
/// Two-loop Standard Model corrections to rho parameter
double rho2(double r);
/// Provides the first guess at a SUSY object at mt, inputting tanb and oneset
/// (should be at MZ) - it's very crude, doesn't take radiative corrections
/// into account etc. oneset provides low energy data and tanb=tan beta
MssmSusy guessAtSusyMt(double tanb, const QedQcd & oneset);

/// For a given trial value of the log of field H2, gives the value of the
/// potential at the minimum. The following global variables must be set before
/// it is called: unificationScale=high BC scale, minTol=fractional accuracy
/// with which minimum is found
extern double minimufb3(double);
/// Given mu paramer, tau Yukawa htau, family number examined, finds height of
/// potential for temp at |H_2|=h2.
double ufb3fn(double mu, double htau, double h2, int family, const MssmSoftsusy
	      & temp);
/// For UFB-3direction, returns scale at which one-loop corrections are
/// smallest. IO parameters: inminTol is fractional accuracy with which
/// minimum is found, eR is value of RH selectron field, h2 is value of H2
/// field, Lisq=|L_i|^2 slepton VEV value, mx=high BC-scale
double getQhat(double inminTol,double eR, double h2, double Lisq, double mx,
	       MssmSoftsusy & temp);

/// non-universal mSUGRA boundary conditions
void extendedSugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// universal mSUGRA boundary conditions
void sugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Adds 2-loop AMSB boundary conditions onto m
void amsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// One-loop GMSB boundary conditions
void gmsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
void userDefinedBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Sets all soft parameters in m except for mh1sq or mh2sq: it is intended
/// for the case where mu and M_A^0(pole) is specified
void generalBcs2(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// initialises factors (gfL,R) for loop calculations
void initialise();

/// function used for calculating sin theta_eff
double gEff(double x);
/// function used for calculating sin theta_eff
double fEff(double x);
#endif



