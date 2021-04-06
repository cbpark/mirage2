
/** \file softpoint.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: main calling program: command line interface

   $Log: softpoint.cpp,v $
   Revision 1.7  2006/04/11 13:57:40  allanach
   Better comments in main programs and cleaned up bug in QEWSB non-usage

   Revision 1.6  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.10  2005/09/16 13:03:07  allanach
   Modified to use new SLHA convention for phenomenological MSSM: supply
   mu(MSUSY) and mA(pole). Other bug-fixes related to phenoMSSM.

   Revision 1.9  2005/08/16 17:22:08  allanach
   Corrected electroweak sbottom corrections

   Revision 1.8  2005/08/16 14:31:00  allanach
   Removed dependence upon massIn and made full 2-loop RGE running the default

   Revision 1.7  2005/08/16 13:54:02  allanach
   Got rid of testing code

   Revision 1.6  2005/07/15 15:10:06  allanach
   Added a routine to calculate sin^2 theta_eff

   Revision 1.5  2005/06/09 10:14:04  allanach
   Deleted general boundary conditions -- they were confusing.

   Revision 1.4  2005/05/30 14:22:14  allanach
   Fixed stau mixing in ISAWIG interface to 7.64

   Revision 1.3  2004/12/23 15:38:14  allanach
   Implementation of 2-loop scalar terms in SUSY Les Houches Accord

   Revision 1.2  2004/12/23 15:29:20  allanach
   Promoted INCLUDE_2_LOOP_SCALAR_CORRECTIONS to a global variable (in
   preparation for its control in the SUSY Les Houches Accord)

   Revision 1.62  2004/09/28 10:43:12  allanach
   Name-change to les Houches function - overloaded standard and file output
   functions.

   Revision 1.61  2004/07/30 13:32:47  allanach
   xxx renamed QEWSB and added to Les Houches accord in SOFTSUSY block

   Revision 1.60  2004/05/23 19:18:05  allanach
   Initial propaganda printed in error output channel to avoid spoiling les
   Houches output format files

   Revision 1.59  2004/05/23 18:33:28  allanach
   Added a SOFTSUSY Block to SUSY Les Houches accord: you can now specify (within
   this block):
     1   TOLERANCE
     2   MIXING
     3   PRINTOUT

   Revision 1.58  2004/05/21 16:11:19  allanach
   Bug-fixed SLHA with respect to GMSB. Wasn't working previously

   Revision 1.57  2004/03/26 17:33:38  allanach
   Added additional option to make pheno MSSM easier: same mass for scalars
   except Higgs' and gauginos/A0 parameter.

   Revision 1.56  2004/03/23 00:16:24  allanach
   Added option to set MGUT=mEWSB

   Revision 1.55  2004/03/21 20:43:05  allanach
   Added alternative electroweak symmetry breaking conditions. Added possibility
   of EWSB=SUSY breaking boundary condition scale in SUSY Les Houches Accord.

   Revision 1.54  2004/03/17 13:23:06  allanach
   Solved 1-loop problem in alternative EWSB conditions. 2-loop problem remains.

   Revision 1.53  2004/02/26 16:33:35  allanach
   Removed random interfaces

   Revision 1.52  2004/02/09 16:54:55  allanach
   mgut now set to gauge unification

   Revision 1.51  2004/01/15 13:54:55  allanach
   New header style implemented

   Revision 1.50  2003/12/02 19:15:53  allanach
   Added 3-family trilinear information, and added non-universal input

   Revision 1.49  2003/11/25 15:52:49  allanach
   massIn does not have to exist - in which case it will use defaults

   Revision 1.48  2003/11/19 17:07:44  allanach
   Running mb mass implemented in SUSY Les Houches Accord

   Revision 1.47  2003/10/24 16:30:26  allanach
   Deleted old higgsVevMs variable, using running variable instead

   Revision 1.46  2003/08/19 14:38:59  allanach
   Altered so that mx in arguments to lowOrg is unchanged. It's used as the
   initial guess and lowOrg returns a double number as the calculated mgut.

   Revision 1.45  2003/08/19 14:26:22  allanach
   Changing lowOrg to be more sensible about gauge unification. Should now be
   called with POSITIVE mgut and a flag for gauge unification.

   Revision 1.44  2003/08/19 14:02:03  allanach
   Can now use "unified" to flag gauge unification

   Revision 1.43  2003/08/07 12:27:49  allanach
   Version relevant for producing numbers for Pietro

   Revision 1.42  2003/08/04 16:45:47  allanach
   Conforms with latest SUSY Les Houches draft

   Revision 1.41  2003/07/28 17:12:42  allanach
   More error trapping

   Revision 1.40  2003/07/24 15:30:42  allanach
   More changes to Les Houches interfacing: LH accord output only happens if LH
   accord stuff is input

   Revision 1.39  2003/07/24 14:55:28  allanach
   Implemented les Houches input and output properly in the usual command-line
   interface

   Revision 1.38  2003/07/22 09:16:04  allanach
   Added running MW, MZ to definition of drbar parameters

   Revision 1.37  2003/07/21 14:00:18  allanach
   MZ fully implemented as an input now. Kept MZ as the central PDG 2002 value,
   for defaults etc

   Revision 1.36  2003/07/21 13:10:39  allanach
   Removed MZ as constant and added it as an input

   Revision 1.35  2003/07/21 08:47:28  allanach
   MW now part of MssmSoftsusy class, and prediction is implemented fully

   Revision 1.34  2003/07/18 14:39:20  allanach
   Implemented MW as a global variable (in preparation for predicting it),
   and also speed corrections in getVev and rhohat: allowing input of
   self-energies to remove their calculation several times

   Revision 1.33  2003/07/16 11:07:06  allanach
   Changed isajet number to 764

   Revision 1.31  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.30  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.29  2003/03/28 15:59:35  allanach
   Two-loop stuff relegated to seperate file

   Revision 1.28  2003/03/09 19:04:08  allanach
   General specification of mA now works to 2-loop order

   Revision 1.27  2003/02/24 18:17:55  allanach
   More work on DRbar loop corrections. Delete physical calculation each
   iteration, check fine-tuning and also pole fermion masses used in loops.

   Revision 1.26  2003/02/24 14:26:39  allanach
   Implementing DRbar parameters in loop corrections. Half way though:
   drbarpars now changes to MPZ notation, but need to get rid of pole masses
   in loop corrections (and re-calculate tadpoles).

   Revision 1.25  2003/02/21 17:59:36  allanach
   Added drbar parameter class and calculation, starting to move to DRbar
   parameters in the 1-loop corrections

   Revision 1.24  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.20  2002/10/17 13:13:43  allanach
   Error output changed to be more helpful

   Revision 1.19  2002/10/15 13:52:32  allanach
   mu is set beforehand in the case of "general" boundary conditions,
   in order to get a decent solution to REWSB.

   Revision 1.18  2002/10/14 14:17:30  allanach
   Added "runto" command in softpoint.x to get micromegas inputs to a
   different scale

   Revision 1.17  2002/10/14 12:41:34  allanach
   Added tachyon information when PRINTOUT > 2. Higgs tachyons now trapped.

   Revision 1.16  2002/10/11 15:26:26  allanach
   Micromegas interface changed

   Revision 1.15  2002/10/10 13:37:12  allanach
   Added physical parameters to Micromegas output

   Revision 1.14  2002/10/09 17:18:51  allanach
   Added new softpoint option to specify mu and mA instead of mH1, mH2. Some
   fine-tuning still necessary

   Revision 1.13  2002/09/11 13:36:52  allanach
   Working version: RPV works with one-loop RPC corrections

   Revision 1.12  2002/09/11 09:52:34  allanach
   TOLERANCE etc now used in preference.

   Revision 1.11  2002/09/09 10:42:54  allanach
   TOLERANCE replaces EPS as being more user-friendly

   Revision 1.10  2002/09/09 10:21:12  allanach
   For gauge unification, enter "?" for mgut

   Revision 1.9  2002/09/03 14:29:38  allanach
   EPS, PRINTOUT, MIXING now all set in massIn

   Revision 1.8  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, EPS out of def.h to make it quicker to compile once
   they are changed.

   Revision 1.6  2002/07/19 15:54:07  allanach
   SOFTSUSY1.5 version

   Revision 1.4  2002/04/26 15:14:44  allanach
   Deleted all translation routines and defined boundary conditions within
   softsusy.h and softsusy.cpp

   Revision 1.3  2002/04/25 12:19:11  allanach
   Cleared bug that happened if there was no argument to the function

   Revision 1.1  2002/03/01 17:33:58  allanach
   Added softpoint: allows you to call softsusy just from the command line.
   Also, interfaces to micromegas added. All interfaces write to input files:
   SoftsusyMicromegasInterface, SoftsusyIsawigInterface,
   SoftsusyIsajetInterface, SoftsusySsrunInterface respectively.
*/

#include <softpoint.h>

void extendedSugraBcs2(MssmSoftsusy & m,
                       const DoubleVector & inputParameters) {
  int i;
  for (i=1; i<=3; i++) m.setGauginoMass(i, inputParameters.display(i));
  m.setTrilinearElement(UA, 1, 1, m.displayYukawaElement(YU, 1, 1) *
                        inputParameters.display(11));
  m.setTrilinearElement(UA, 2, 2, m.displayYukawaElement(YU, 2, 2) *
                        inputParameters.display(11));
  m.setTrilinearElement(UA, 3, 3, m.displayYukawaElement(YU, 3, 3) *
                        inputParameters.display(11));
  m.setTrilinearElement(DA, 1, 1, m.displayYukawaElement(YD, 1, 1) *
                        inputParameters.display(12));
  m.setTrilinearElement(DA, 2, 2, m.displayYukawaElement(YD, 2, 2) *
                        inputParameters.display(12));
  m.setTrilinearElement(DA, 3, 3, m.displayYukawaElement(YD, 3, 3) *
                        inputParameters.display(12));
  m.setTrilinearElement(EA, 1, 1, m.displayYukawaElement(YE, 1, 1) *
                        inputParameters.display(13));
  m.setTrilinearElement(EA, 2, 2, m.displayYukawaElement(YE, 2, 2) *
                        inputParameters.display(13));
  m.setTrilinearElement(EA, 3, 3, m.displayYukawaElement(YE, 3, 3) *
                        inputParameters.display(13));
  m.setSoftMassElement(mLl, 1, 1, sqr(inputParameters.display(31)));
  m.setSoftMassElement(mLl, 2, 2, sqr(inputParameters.display(32)));
  m.setSoftMassElement(mLl, 3, 3, sqr(inputParameters.display(33)));
  m.setSoftMassElement(mEr, 1, 1, sqr(inputParameters.display(34)));
  m.setSoftMassElement(mEr, 2, 2, sqr(inputParameters.display(35)));
  m.setSoftMassElement(mEr, 3, 3, sqr(inputParameters.display(36)));
  m.setSoftMassElement(mQl, 1, 1, sqr(inputParameters.display(41)));
  m.setSoftMassElement(mQl, 2, 2, sqr(inputParameters.display(42)));
  m.setSoftMassElement(mQl, 3, 3, sqr(inputParameters.display(43)));
  m.setSoftMassElement(mUr, 1, 1, sqr(inputParameters.display(44)));
  m.setSoftMassElement(mUr, 2, 2, sqr(inputParameters.display(45)));
  m.setSoftMassElement(mUr, 3, 3, sqr(inputParameters.display(46)));
  m.setSoftMassElement(mDr, 1, 1, sqr(inputParameters.display(47)));
  m.setSoftMassElement(mDr, 2, 2, sqr(inputParameters.display(48)));
  m.setSoftMassElement(mDr, 3, 3, sqr(inputParameters.display(49)));
}

// Returns a string with all characters in upper case: very handy
string ToUpper(const string & s) {
        string result;
        unsigned int index;
        for (index = 0; index < s.length(); index++) {
          char a = s[index];
          a = toupper(a);
          result = result + a;
        }

        return result;
    }

void errorCall() {
  ostringstream ii;
  ii << "SOFTSUSY called with incorrect arguments. Need to put either:\n";
  ii << "softpoint.x sugra <m0> <m12> <a0> <tanb> <mgut> <sgnMu>\n";
  ii << "softpoint.x amsb <m0> <m32> <tanb> <mgut> <sgnMu>\n";
  ii << "softpoint.x gmsb <n5> <mMess> <lambda> <cgrav> <tanb> <sgnMu> \n";
  ii << "softpoint.x mirage <alpha> <M0> <tanb> <sgnMu> <am> <ah> <cm> <chu> <chd> <mgut>\n";
  ii << "softpoint.x leshouches < lesHouchesInput \n";
  ii << "where bracketed entries are numerical values.\n";
  ii << "<mgut> denotes the scale at which the SUSY breaking ";
  ii << "parameters are specified. \n";
  ii << "Enter <mgut>=unified to define MGUT by g1(MGUT)=g2(MGUT), otherwise";
  ii << "\nit will be fixed at the numerical value input.\n";
  ii << "For SUSY breaking terms set at MSUSY, enter <mgut>=msusy.\n";
  ii << "lesHouchesInput contains the SUSY Les Houches accord";
  ii << " input.\n";
  throw ii.str();
}

/// start of definition of global variables
/// no quark mixing (dominant third family approx), and no verbose output
int MIXING = -1, PRINTOUT = 0;
/// Do we include 2-loop RGEs of *all* scalar masses and A-terms, or only the
/// scalar mass Higgs parameters? (Other quantities all 2-loop anyway): the
/// default in SOFTSUSY 2.x is to include all 2-loop terms
bool INCLUDE_2_LOOP_SCALAR_CORRECTIONS = true;
/// fractional accuracy required
double TOLERANCE = 1.0e-3;
/// there are two possible conventions: if QEWSB > MZ, its value is assumed
/// in GeV and used as a constant MSUSY. Otherwise, it MULTIPLIES the usual
/// MSUSY value, of root(mstop1 mstop2)
double QEWSB = 1.0;
/// number of loops used to calculate Higgs mass and tadpoles. They should be
/// identical for a consistent calculation
int numRewsbLoops = 2, numHiggsMassLoops = 2;
/// decay constant of muon
double GMU = 1.16637e-5;
/// end of global variable declaration

int main(int argc, char *argv[]) {

  int numPoints = 1;

  double qMax = 0.;

  // Sets format of output: 4 decimal places
  outputCharacteristics(6);

  void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &)=sugraBcs;

  MssmSoftsusy m; MssmSoftsusy2 l;
  MssmSoftsusy * r = &m;

  try {
    QedQcd oneset;

  if (argc !=1 && strcmp(argv[1],"leshouches") != 0) {
    cerr << "SOFTSUSY" << VERSION << endl;
    cerr << "B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331,";
    cerr << " hep-ph/0104145\n";
    cerr << "Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE="
         << TOLERANCE << endl;
    cerr << "G_F=" << GMU << " GeV^2" << endl;
  }

  double mgutGuess = 2.0e16, tanb = 10., mbmb = MBOTTOM, mtau = MTAU;
  int sgnMu = 1;
  bool gaugeUnification = true, ewsbBCscale = false, altEwsb = false;

    // If there are no arguments, give error message,
    // or if none of the options are called, then go to error message
    if (argc == 1 || (strcmp(argv[1], "sugra") && strcmp(argv[1], "amsb") &&
                      strcmp(argv[1], "gmsb") &&
                      strcmp(argv[1], "mirage") &&
                      strcmp(argv[1], "runto") &&
                      strcmp(argv[1], "leshouches")))
      errorCall();

    DoubleVector pars(3);

    char * modelIdent = "";

    if (!strcmp(argv[1], "runto")) {
      if (argc == 3) {
        MssmSoftsusy r;
        r.microMegasIn("SoftSusyMicromegasInterface");
        double mu = atof(argv[2]);
        r.runto(mu);
        r.microMegasInterface("SoftSusyMicromegasInterface");
      }
      else
        errorCall();
    }

    if (!strcmp(argv[1], "sugra")) {
      cout << "SOFTSUSY SUGRA calculation" << endl;
      boundaryCondition = &sugraBcs;
      if (argc == 8) {
        double m0 = atof(argv[2]);
        double m12 = atof(argv[3]);
        double a0 = atof(argv[4]);
        tanb = atof(argv[5]);
        mgutGuess = mgutCheck(argv[6], gaugeUnification, ewsbBCscale);
        sgnMu =  atoi(argv[7]);
        pars(1) = m0; pars(2) = m12; pars(3) = a0;
        r = &m;
      } else if (argc == 9) {
        double m0 = atof(argv[2]);
        double m12 = atof(argv[3]);
        double a0 = atof(argv[4]);
        tanb = atof(argv[5]);
        mgutGuess = mgutCheck(argv[6], gaugeUnification, ewsbBCscale);
        sgnMu =  atoi(argv[7]);
        pars(1) = m0; pars(2) = m12; pars(3) = a0;
        QEWSB = atof(argv[8]);
        r = &m;
      }
      else
        errorCall();
      // end of SUGRA option
    }

    if (!strcmp(argv[1], "amsb")) {
      cout << "SOFTSUSY mAMSB calculation" << endl;
      boundaryCondition = &amsbBcs;
      if (argc == 7) {
        double m0 = atof(argv[2]);
        double m32 = atof(argv[3]);
        tanb = atof(argv[4]);
        mgutGuess = mgutCheck(argv[5], gaugeUnification, ewsbBCscale);
        sgnMu =  atoi(argv[6]);
        pars(1) = m32; pars(2) = m0;
        r = &m;
      }
      else
        errorCall();
    }

    if (!strcmp(argv[1], "gmsb")) {
      cout << "SOFTSUSY mGMSB calculation" << endl;

      boundaryCondition = &gmsbBcs;

      if (argc == 8) {
        double n5 = atof(argv[2]);
        double mMess = atof(argv[3]);
        double lambda = atof(argv[4]);
        double cgrav = atof(argv[5]);
        tanb = atof(argv[6]);
        sgnMu =  atoi(argv[7]);
        mgutGuess = mMess;
        gaugeUnification = false;
        pars.setEnd(4);
        pars(1) = n5; pars(2) = mMess; pars(3) = lambda; pars(4) = cgrav;
        r = &m;
        if (lambda > mMess) {
          ostringstream ii;
          ii << "Input lambda=" << lambda << " should be greater than mMess="
             << mMess << endl;
          throw ii.str();
        }
        if (cgrav > 1.0) {
          ostringstream ii;
          ii << "Input cgrav=" << cgrav << " a real number bigger than or "
             << " equal to 1 (you can use 1 as a default value).\n";
          throw ii.str();
        }
      }
      else
        errorCall();
    }

    if (!strcmp(argv[1], "mirage")) {
      cout<< "SOFTSUSY Mirage calculation" <<endl;

      boundaryCondition = &userDefinedBcs;

      if (argc == 11) {
        double alphac = atof(argv[2]);
        double M0 = atof(argv[3]);
        tanb = atof(argv[4]);
        sgnMu= atoi(argv[5]);
        double am = atof(argv[6]);
        double ah = atof(argv[7]);
        double cm = atof(argv[8]);
        double ch = atof(argv[9]);
        mgutGuess = mgutCheck(argv[10], gaugeUnification, ewsbBCscale);
        pars.setEnd(6);
        pars(1) = alphac; pars(2) = M0; pars(3) = am; pars(4) = ah;
        pars(5) = cm; pars(6) = ch;
        r = &m;
      }
      else
        //cout<<"softpoint.x mirage alpha M0 tanb sgnMu am ah cm ch mgut"<<endl;
        errorCall();
    }

    bool flag = false;
    if (!strcmp(argv[1], "leshouches")) {
      if (argc == 2) {
        string line, block;
        int model;

        while (getline(cin,line)) {
          //      mgutGuess = mgutCheck("unified", gaugeUnification);

          //    cout << line << endl;
          istringstream input(line);
          string word1, word2;
          input >> word1;

          if (word1.find("#") == string::npos) {
            // read in another word if there's no comment
            input >> word2;

            if (ToUpper(word1) == "BLOCK")  {
              block = ToUpper(word2);

            } else { // ought to be data
              istringstream kk(line);
              if (block == "MODSEL") {
                int i; kk >> i;

                switch(i) {
                case 1: kk >> model;
                  switch(model) {
                  case 0: boundaryCondition = &extendedSugraBcs;
                    modelIdent = "nonUniversal"; r=&m;
                    break;
                  case 1: boundaryCondition = &sugraBcs; pars.setEnd(3);
                    modelIdent = "sugra"; r=&m;
                    break;
                  case 2: boundaryCondition = &gmsbBcs; pars.setEnd(4); r=&m;
                    modelIdent = "gmsb"; r=&m;
                    break;
                  case 3: boundaryCondition = &amsbBcs; pars.setEnd(2); r=&m;
                    modelIdent = "amsb"; r=&m;
                    break;
                  case 4: boundaryCondition = &userDefinedBcs; pars.setEnd(6); r=&m;
                    modelIdent = "mirage"; r=&m;
                    break;
                  default:
                    ostringstream ii;
                    ii << "SOFTSUSY" << VERSION << " cannot yet do model "
                       << model << ": terminal error\n";
                    throw ii.str();
                  }
                  break;
                case 11: kk >> numPoints;
                  if (numPoints < 1) {
                    ostringstream ii;
                    ii << "MODEL 10 selecting silly number of points"
                       << "(" << numPoints << ") to output" << endl;
                    throw ii.str();
                  }
                  break;
                case 12: double d; kk >> d;
                  if (d < MZ) {
                    ostringstream ii;
                    ii << "MODSEL 12 selecting silly scale Qmax"
                       << "(" << d << ") < MZ to output" << endl;
                    throw ii.str();
                  }
                  qMax = d; break;
                default:
                  ostringstream ii;
                  ii << "# WARNING: don't understand first integer " << word1
                     << " " << word2 << " in block " << block
                     << ": ignoring it\n";
                  throw ii.str();
                  break;
                }
              }
              else if (block == "MINPAR") {
                int i; double d; kk >> i >> d;
                switch (i) {
                case 3: tanb = d; break;
                case 4: sgnMu = int(d); break;
                default:
                  switch(model) {
                  case 0:
                    // SUGRA inputs to fill out the pheno MSSM case
                    switch(i) {
                    case 1: pars(1) = d; break;
                    case 2: pars(2) = d; break;
                    case 5: pars(3) = d; break;
                    default:
                      ostringstream ii;
                      ii << "Didn't understand pheno MSSM input " << i << endl;
                      break;
                    } break;
                  case 1: // SUGRA inputs
                    switch(i) {
                    case 1: pars(1) = d; break;
                    case 2: pars(2) = d; break;
                    case 5: pars(3) = d; break;
                    default:
                      ostringstream ii;
                      ii << "Didn't understand SUGRA input " << i << endl;
                      break;
                    } break;
                  case 2: // GMSB inputs
                    switch(i) {
                    case 1: pars(3) = d; break;
                    case 2: pars(2) = d; mgutGuess = d;
                      gaugeUnification = false; break;
                    case 5: pars(1) = d; break;
                    case 6: pars(4) = d; break;
                    default:
                      ostringstream ii;
                        ii << "Didn't understand GMSB input " << i << endl;
                        break;
                    } break;
                  case 3: // AMSB inputs
                    switch(i) {
                    case 1: pars(2) = d; break;
                    case 2: pars(1) = d; break;
                    default:
                      ostringstream ii;
                      ii << "Didn't understand AMSB input " << i << endl;
                      break;
                    } break;
                  case 4: //mirage inputs
                     switch(i) {
                     case 1: pars(1) = d; break;
                     case 2: pars(2) = d; break;
                     case 5: pars(3) = d; break;
                     case 6: pars(4) = d; break;
                     case 7: pars(5) = d; break;
                     case 8: pars(6) = d; break;
                     case 9: pars(7) = d; break;
                     default:
                       ostringstream ii;
                       ii <<"Didn't understand mirage model input"<<i<<endl;
                       break;
                     } break;
                  default:
                    ostringstream ii;
                    ii << "Didn't understand model input " << i << endl;
                    break;
                  }
                  break;
                }
              }
              // Adding non-minimal options. However, mA and mu option is not
              // yet supported. Also, we assume the initial model was mSUGRA
              // (for now).
              else if (block == "EXTPAR") {
                if (modelIdent == "sugra" || modelIdent == "nonUniversal") {
                  if (pars.displayEnd() != 49) { // initialise vector
                    modelIdent = "nonUniversal"; r=&m;
                    boundaryCondition = &extendedSugraBcs;
                    double m0 = pars(1), m12 = pars(2), a0 = pars(3);
                    pars.setEnd(49);
                    int i; for (i=1; i<=3; i++) pars(i) = m12;
                    for (i=11; i<=13; i++) pars(i) = a0;
                    pars(21) = m0 * m0; pars(22) = m0 * m0;
                    for (i=31; i<=36; i++) pars(i) = m0;
                    for (i=41; i<=49; i++) pars(i) = m0;
                  }
                  int i; double d; kk >> i >> d;
                  if ((i > 0 && i <=  3) || (i >= 11 && i <= 13) ||
                      (i >= 21 && i <= 23) || (i == 26)
                      || (i >= 31 && i <= 36) ||
                      (i >= 41 && i <= 49)) pars(i) = d;
                  if (i == 0) {
                    mgutGuess = d;
                    gaugeUnification = false;
                    // setting Minput=-1 should yield MSSM BCs at MSUSY
                    if ((d - 1.0) < EPSTOL) {
                      mgutGuess = 1.0e3;
                      ewsbBCscale = true;
                      QEWSB = 1.0;
                      if (gaugeUnification)
                        cout << "# Gauge unification ignored since pheno MSSM"
                             << " assumes BC set at QEWSB\n";
                      gaugeUnification = false;
                    }
                  }
                  if (i == 23 || i == 26) altEwsb = true;
                  if (!((i >= -1 && i <=  3) || (i >= 11 && i <= 13) ||
                      (i >= 21 && i <= 23) || (i == 26)
                        || (i >= 31 && i <= 36) ||
                        (i >= 41 && i <= 49)))
                    cout << "# Didn't understand extra parameter " << i
                            << " - ignoring it" << endl;
                }
                else cout << "# Can't yet handle extra parameters in "
                          << modelIdent << ": ignoring them" << endl;
              }
              else if (block == "SMINPUTS") {
                int i; double d; kk >> i >> d;
                switch (i) {
                case 1: oneset.setAlpha(ALPHA, 1.0 / d); break;
                case 2: GMU = d; break;
                case 3: oneset.setAlpha(ALPHAS, d); break;
                case 4: oneset.setMu(d); break;
                case 5: oneset.setMass(mBottom, d); mbmb = d; flag = true;
                  break;
                case 6: oneset.setPoleMt(d); break;
                case 7: oneset.setMass(mTau, d); mtau = d; break;
                default:
                  ostringstream ii;
                  ii << "# WARNING: Don't understand data input " << i
                       << " " << d << " in block "
                       << block << ": ignoring it\n"; break;
                  throw ii.str();
                }
              }
              else if (block == "SOFTSUSY") {
                int i; double d; kk >> i >> d;
                switch(i) {
                case 1: TOLERANCE = d; break;
                case 2: MIXING = int(d); break;
                case 3: PRINTOUT = int(d); break;
                case 4: QEWSB = d; break;
                case 5: INCLUDE_2_LOOP_SCALAR_CORRECTIONS = bool(int(d));
                  break;
                default:
                  ostringstream ii;
                  ii << "# WARNING: Don't understand data input " << i
                       << " " << d << " in block "
                       << block << ": ignoring it\n"; break;
                  throw ii.str();
                }
              }
              //              else if (block == "SPCFILE")
              //        kk >> fileName;
              else {
                ostringstream ii;
                ii << "# WARNING: don't recognise block " << block
                   << ": ignoring all data in it" << endl;
                throw ii.str();
              }
              // end if blocks

            } // end of data
          } // end of no-comment

        } // end of file

      }
      else errorCall();
    }

    if (altEwsb) {
      boundaryCondition = &extendedSugraBcs2;
      l.muCond = pars(23);
      l.mAcond = pars(26);
      l.setSusyMu(l.muCond);
      sgnMu = 0; // Flags different BCs
      r = &l;
    }

    // intput error checking
    if (sgnMu != 1 && sgnMu != -1 && sgnMu != 0) {
      ostringstream ii;
      ii << "Incorrect input for sign(mu)=" << sgnMu <<endl;
      throw ii.str();
    }
    if (tanb < 1.0 || tanb > 1.0e2)  {
      ostringstream ii;
      ii << "Incorrect input for tan beta=" << tanb <<endl;
      throw ii.str();
    }

    if (flag) oneset.calcPoleMb();
    oneset.toMz();

    double mgut =  r->lowOrg(boundaryCondition, mgutGuess, pars, sgnMu,
                             tanb, oneset, gaugeUnification, ewsbBCscale);

    if (strcmp(argv[1],"leshouches") != 0) {
      cout << "mgut = " << mgut << endl << *r;
    }
    else {
      r->lesHouchesAccordOutput(modelIdent, pars, sgnMu, tanb, qMax,
                                numPoints, mbmb, mtau, mgut, altEwsb);
    }
  }
  catch(const string & a) {
    cout << a;
  }

  if (r->displayProblem().test()) {
    cout << "# SOFTSUSY problem with point: " << r->displayProblem();
  }

  exit(0);
}
