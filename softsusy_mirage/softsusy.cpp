/** \file softsusy.cpp
   Project:     SOFTSUSY
   Author:      Ben Allanach
   Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: softsusy.cpp,v $
   Revision 1.7  2006/04/11 13:57:40  allanach
   Better comments in main programs and cleaned up bug in QEWSB non-usage

   Revision 1.6  2006/01/17 16:40:17  allanach
   Updated to 2.0.4: Higgs mass calculation corrected (problem in the charged
   Higgs contributions)

   Revision 1.40  2005/11/10 18:00:48  allanach
   Fixed potential ordering problem with stop/sbottom masses

   Revision 1.39  2005/11/09 14:07:38  allanach
   Added default option to not calculate fine tuning wrt top Yukawa in fineTuning
   calcultion.

   Revision 1.38  2005/10/27 11:49:24  allanach
   Allowed sbottom mass calculation at p=msbot1, msbot2. Not sure if this will
   cause problems yet. Corrected q-scanning behaviour in SLHA routines.

   Revision 1.37  2005/10/04 11:28:59  allanach
   Cleaned up some of loop corrections, made simpler etc

   Revision 1.36  2005/10/03 15:35:19  allanach
   Debugged corrections to RH selectron mass.

   Revision 1.35  2005/10/03 14:31:48  allanach
   Corrected some sparticle radiative corrections: charged higgs
   contributions. Only makes a small numerical difference.

   Revision 1.34  2005/09/16 18:16:07  allanach
   Added a function which sets a vector of Higgs masses and couplings.

   Revision 1.33  2005/09/16 13:03:37  allanach
   Changes to SLHA implemented with respect to mA(pole).

   Revision 1.32  2005/09/16 12:58:18  allanach
   Modified to account for phenomenological MSSM. Used displayMsusy rather than
   re-calculating it in order to use MSUSY(MSUSY) rather than MSUSY(MZ), which
   was a bug. Changes Les Houches Accord output to be more compliant.

   Revision 1.31  2005/08/25 15:31:20  allanach
   Checked and corrected neutralino corrections

   Revision 1.30  2005/08/25 11:58:24  allanach
   Added neutralino full loop corrections - they need checking

   Revision 1.29  2005/08/23 17:24:19  allanach
   Checked charginos

   Revision 1.28  2005/08/23 14:14:46  allanach
   Added chargino corrections

   Revision 1.27  2005/08/19 16:08:21  allanach
   Corrected sdown corrections

   Revision 1.26  2005/08/19 15:28:55  allanach
   Completed all lighter generations of sfermions

   Revision 1.25  2005/08/19 13:50:02  allanach
   Added up squark loops

   Revision 1.24  2005/08/18 17:32:54  allanach
   Added loop corrections to electron/muon sneutrinos

   Revision 1.23  2005/08/18 16:50:56  allanach
   Checked snu_tau corrections + finalised

   Revision 1.22  2005/08/18 16:39:55  allanach
   Checked some of snu_tau corrections

   Revision 1.21  2005/08/18 14:57:05  allanach
   Added stau sneutrino loop corrections

   Revision 1.20  2005/08/18 11:16:13  allanach
   Corrected and checked stau corrections

   Revision 1.19  2005/08/18 10:13:40  allanach
   Correcting stau loop corrections

   Revision 1.18  2005/08/17 16:09:50  allanach
   Implemented 1-loop corrections to mstau

   Revision 1.17  2005/08/17 13:50:56  allanach
   Added sneutrino contributions to sbottom mass loops

   Revision 1.16  2005/08/16 17:22:08  allanach
   Corrected electroweak sbottom corrections

   Revision 1.15  2005/08/16 16:53:20  allanach
   Added electroweak loops to sbottom loop corrections

   Revision 1.14  2005/08/01 11:31:12  allanach
   Added Markus' NLSP routines

   Revision 1.13  2005/07/26 11:00:42  allanach
   Added method nlsp to object MssmSoftsusy by Bernhardt

   Revision 1.12  2005/07/26 10:51:22  allanach
   new method nmin added to DoubleVector and DoubleMatrix by Bernhardt

   Revision 1.11  2005/07/15 15:10:06  allanach
   Added a routine to calculate sin^2 theta_eff

   Revision 1.10  2005/06/09 10:14:04  allanach
   Deleted general boundary conditions -- they were confusing. Also,

   Revision 1.9  2005/05/30 16:46:37  allanach
   Corrected neutralinos somewhat in the isajet interface

   Revision 1.8  2005/05/30 15:15:19  allanach
   mgluino now calculated at MSUSY for a fairer calculation of the theoretical
   error associated with it

   Revision 1.7  2005/05/30 14:22:14  allanach
   Fixed stau mixing in ISAWIG interface to 7.64

   Revision 1.6  2005/05/24 09:24:36  allanach
   Improved the speed of lsp function slightly

   Revision 1.5  2005/05/13 15:53:35  allanach
   Fixed small bug in SLHA output

   Revision 1.4  2005/04/11 14:10:47  allanach
   Fixed small bug in SLHA of mA^2 running.

   Revision 1.3  2005/03/02 16:43:47  allanach
   Fixed minor bug printing out mA^2(MSUSY) in Les Houches output.

   Revision 1.2  2004/12/23 15:29:20  allanach
   Promoted INCLUDE_2_LOOP_SCALAR_CORRECTIONS to a global variable (in
   preparation for its control in the SUSY Les Houches Accord)

   Revision 1.159  2004/09/28 10:43:12  allanach
   Name-change to les Houches function - overloaded standard and file output
   functions.

   Revision 1.158  2004/09/28 10:39:07  allanach
   Added Les Houches output function for output to a file

   Revision 1.157  2004/09/17 10:13:27  allanach
   If there's no convergence, it still calculates the physical spectrum
   properly

   Revision 1.156  2004/07/13 10:15:06  allanach
   Formally added global variable QEWSB as the possible fixed EWSB scale.
   Bug-fixed tau Yukawa coupling threshold corrections.

   Revision 1.155  2004/05/28 21:10:36  allanach
   Implemented full MSSM corrections to mtau at MZ

   Revision 1.154  2004/05/28 19:18:14  allanach
   Promoted hb Yukawa to full 1-loop corrections

   Revision 1.153  2004/05/20 19:10:25  allanach
   Removed bogus "resummation" of gluino mass 1-loop FX

   Revision 1.152  2004/05/17 15:55:24  allanach
   Added extra throw/catches that could otherwise have screwed up the
   iterative routines when the calling program uses a loop.

   Revision 1.151  2004/05/06 10:16:44  allanach
   Input of mH1 and mH2 @ MSUSY in SLHA format now updated (squared values)

   Revision 1.150  2004/04/20 13:53:21  allanach
   Fixed the spurious error messages from the mA calculation

   Revision 1.149  2004/04/15 18:12:22  allanach
   Calculation of tadpoles split up to allow calculation at different scales.
   The situation with mA^2(MZ) < 0 is handled by putting mApole and pretending
   that's the DRbar mass.

   Revision 1.148  2004/04/07 16:19:09  allanach
   Added proper NaN tests

   Revision 1.147  2004/04/07 15:11:40  allanach
   Improved error checking, specifically wrt Pietro's routines

   Revision 1.146  2004/03/23 00:16:24  allanach
   Added option to set MGUT=mEWSB

   Revision 1.145  2004/03/21 20:43:05  allanach
   Added alternative electroweak symmetry breaking conditions. Added
   possibility of EWSB=SUSY breaking boundary condition scale in SUSY Les
   Houches Accord.

   Revision 1.144 2004/03/19 Mistake

   Revision 1.143  2004/03/19 21:27:25  allanach
   Dispensed with global variable storing 1-loop value of
   tadpoles. Intermediate step is just to include them in physpars structure.

   Revision 1.142  2004/03/19 20:42:02  allanach
   Bug-fixed neutralino and chargino mass calcation.

   Revision 1.141  2004/03/19 16:26:19  allanach
   Added 2-loop tau contribution to mA and printed out mH1 an mH2 running
   parameters in the SUSY Les Houches Accord

   Revision 1.140  2004/03/17 13:23:06  allanach
   Solved 1-loop problem in alternative EWSB conditions. 2-loop problem
   remains.

   Revision 1.139  2004/03/11 18:53:50  allanach
   Bug-fixed conversion of mb(MZ)MSbar -> mb(MZ)DRbar

   Revision 1.138  2004/03/03 11:23:03  allanach
   Fixed bug in alpha -- added finite DRbar-MSbar conversion

   Revision 1.137  2004/02/26 16:33:00  allanach
   Corrected MW computation: re-using the DRbar value of mt in piWW for that
   computation.

   Revision 1.136  2004/02/26 13:54:00  allanach
   Pietro's routines changed so they take negative gluino mass sometimes

   Revision 1.135  2004/02/17 10:36:28  allanach
   Fine-tuning Les Houches accord

   Revision 1.134  2004/02/16 15:51:08  allanach
   Updated les Houches accord to latest standard

   Revision 1.133  2004/02/12 14:17:24  allanach
   Added tau^2 corrections and ht hb to higgs masses and Higgs vev calculation

   Revision 1.132  2004/02/09 16:55:59  allanach
   Bug-fixed, stopped running to zero renormalisation scale.

   Revision 1.131  2004/02/09 09:43:21  allanach
   Bail-out added for fatal problems which can have NaN's

   Revision 1.130  2004/01/28 14:08:53  allanach
   Bug-fixes: for calculation sin theta_w, in the W and Z self-energies, the
   pole top mass is used in the top loop since the BPMZ SM two-loop
   corrections have
   assume that. Some errors were untrapped before, which is changed and gluino
   mass is allowed to be negative (as it is eg in AMSB).

   Revision 1.129  2004/01/20 13:04:03  allanach
   setting Higgs VEV has moved to provide proper convergence

   Revision 1.128  2004/01/09 18:46:52  allanach
   Resummed log in piZZ at one point and correction to alpha(MZ)
   interpretations

   Revision 1.127  2003/12/04 20:10:55  allanach
   Added 3-fam information setting in non-minimal models

   Revision 1.126  2003/12/04 18:35:23  allanach
   Added non-minimal SUSY les Houches options to output

   Revision 1.125  2003/12/02 19:15:53  allanach
   Added 3-family trilinear information, and added non-universal input

   Revision 1.124  2003/11/25 15:53:37  allanach
   Les Houches accord output corrected

   Revision 1.123  2003/11/21 12:08:25  allanach
   Bug-fix in calcDrBarPars

   Revision 1.122  2003/11/06 14:41:12  allanach
   Small compiler error and les Houches output error fixed

   Revision 1.121  2003/10/27 15:50:39  allanach
   Taken old Higgs routine out of MssmSoftsusy, but added it to RpvSoft (there
   are currently problems using the new routine there -- needs to be fixed)

   Revision 1.118  2003/10/24 16:30:25  allanach
   Deleted old higgsVevMs variable, using running variable instead

   Revision 1.117  2003/10/24 10:07:42  allanach
   Any calculation of sin theta_w (DRbar) is now made from the EW gauge
   couplings

   Revision 1.116  2003/09/23 16:12:47  allanach
   Fixed bug in pseudoscalar Higgs couplings in the stop and sbottom loop
   corrections.

   Revision 1.115  2003/09/09 11:48:37  allanach
   Corrected resummed calculation of mb (NB if you put in any higher-order
   corrections, you'll have to change it). See hep-ph/9912251.

   Revision 1.114  2003/09/03 12:39:38  allanach
   Small bug-fixes in 2-loop Higgs calls, and pisisi

   Revision 1.113  2003/08/19 14:38:59  allanach
   Altered so that mx in arguments to lowOrg is unchanged. It's used as the
   initial guess and lowOrg returns a double number as the calculated mgut.

   Revision 1.112  2003/08/19 14:26:22  allanach
   Changing lowOrg to be more sensible about gauge unification. Should now be
   called with POSITIVE mgut and a flag for gauge unification.

   Revision 1.111  2003/08/12 15:19:54  allanach
   Bug-fixed problem flags, iterateMu. Made iterateMu less precise (was
   sometimes causing problems)

   Revision 1.110  2003/08/08 15:36:33  allanach
   Bug-fixed calcRunningMt: static initialisation problem

   Revision 1.109  2003/08/07 12:27:16  allanach
   Some more work on Les Houches accord

   Revision 1.108  2003/08/04 16:47:39  allanach
   Added function maxMass to return mass of heaviest superparticle, efficiency
   savings and 2-loop squark/gluino corrections to mt included

   Revision 1.107  2003/07/31 12:50:25  allanach
   Added rest of piH+H- (unchecked as yet) and Higgs calculation doesn't flag a
   tachyon if only the tree-level mA^2 value is negative

   Revision 1.106  2003/07/30 13:39:01  allanach
   Corrected bug in piAA. Started to add piH+H-.

   Revision 1.105  2003/07/30 09:01:36  allanach
   Cleaning up piAA

   Revision 1.104  2003/07/30 08:51:02  allanach
   Checking piAA routine

   Revision 1.103  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.102  2003/07/24 14:55:28  allanach
   Implemented les Houches input and output properly in the usual command-line
   interface

   Revision 1.101  2003/07/22 10:26:01  allanach
   gdL and gdR used in loop corrections are also 1-loop values now.
   This completes the DR-bar-ness of all 1-loop corrections in SOFTSUSY

   Revision 1.100  2003/07/22 09:47:41  allanach
   Implemented running MW and MZ in all loop calculations, although for
   external momenta, the pole value is used.

   Revision 1.99  2003/07/22 09:20:04  allanach
   Added display functions for running MW/MZ

   Revision 1.98  2003/07/22 09:16:04  allanach
   Added running MW, MZ to definition of drbar parameters

   Revision 1.97  i2003/07/21 16:08:54  allanach
   displayHiggsVevMs now only called at MSUSY, and is checked @ that scale.

   Revision 1.96  2003/07/21 14:00:18  allanach
   MZ fully implemented as an input now. Kept MZ as the central PDG 2002 value,
   for defaults etc

   Revision 1.95  2003/07/21 13:08:37  allanach
   Changed so that MZ comes from the input file "massIn"

   Revision 1.94  2003/07/21 08:47:28  allanach
   MW now part of MssmSoftsusy class, and prediction is implemented fully

   Revision 1.93  2003/07/18 15:39:14  allanach
   Added prediction of MW to definition of a softsusy object

   Revision 1.92  2003/07/18 14:57:46  allanach
   MW prediction implemented

   Revision 1.91  2003/07/18 14:39:20  allanach
   Implemented MW as a global variable (in preparation for predicting it),
   and also speed corrections in getVev and rhohat: allowing input of
   self-energies to remove their calculation several times

   Revision 1.90  2003/07/16 11:07:06  allanach
   Changed isajet number to 764

   Revision 1.89  2003/06/25 17:58:40  allanach
   Put in 2-loop higgs mass corrections (in third family approximation) as a
   default in softpars

   Revision 1.88  2003/06/24 08:29:06  allanach
   Setting MSUSY=root(mt1(mt1) mt2(mt2)) rather than root(mt1(MZ) mt2(MZ)),
   also have fixed the implementation of mh 2-loop corrections and
   a small sneutrino correction in pis1s1 (negligible)

   Revision 1.86  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.85  2003/05/27 15:05:52  allanach
   Purely efficiency corrections: used variable rather than display() methods
   whenever possible

   Revision 1.84  2003/05/22 12:45:41  allanach
   Fortran version of 2-loop Higgs corrections included with corrections:
   alpha_s (alpha_t + alpha_b) + alpha_t^2 + alpha_b^2 + alpha_t alpha_b

   Revision 1.83  2003/05/14 15:01:20  allanach
   Made sure that initial object overwritten with zeroes in order to get
   identical results if the parameter is scanned over or not

   Revision 1.82  2003/03/31 12:08:15  allanach
   Bug-fixed alpha_s alpha_b corrections to the tadpoles

   Revision 1.81  2003/03/28 15:56:07  allanach
   Seperated 2-loop Higgs/EWSB symmetry breaking bits off into a different
   file, also 2-loop EWSB all controlled by doTadpoles now

   Revision 1.80  2003/03/25 16:27:24  allanach
   Added routines to protect against nan's during physical calculation

   Revision 1.77  2003/03/24 09:53:09  allanach
   Work progressing on isawig interface

   Revision 1.76  2003/03/09 18:53:13  allanach
   Fixed slight bug in 2-loop corrections

   Revision 1.75  2003/03/07 11:32:48  allanach
   Added tested 2-loop alpha_t^2, alpha_s alpha_t, alpha_s alpha_b corrections
   to REWSB and bug-fixed the Higgs mass calculation

   Revision 1.74  2003/03/05 10:26:26  allanach
   Added alpha_s alpha_b corrections to Higgs masses

   Revision 1.73  2003/03/04 12:54:39  allanach
   Removed Slavich's routine again.

   Revision 1.72: 2003/02/28 18:11:57 Added routines to convert to Slavich's
   stop conventions in order to calculate the 2-loop Higgs mass

   Revision 1.71: 2003/02/25 14:06:11 Final working DRbar corrected version,
   but must do: 2-loop higgs and newhiggs calculation

   Revision 1.70: 2003/02/25 13:21:10 Working version with DRbar parameters,
   there was a bug in 1.69 of calcDrBarPars

   Revision 1.69: 2003/02/24 18:17:55 More work on DRbar loop
   corrections. Delete physical calculation each iteration, check fine-tuning
   and also pole fermion masses used in loops.

   Revision 1.68 2003/02/24 14:26:39
   Implementing DRbar parameters in loop corrections.
   Half way though: drbarpars now changes to MPZ notation, but need to get rid
   of pole masses in loop corrections (and re-calculate tadpoles).

   Revision 1.67 2003/02/21 17:59:36
   Added drbar parameter class and calculation,
   starting to move to DRbar parameters in the 1-loop corrections

   Revision 1.66 2003/02/19 12:50:02
   Changed headings to new conventions

   Revision 1.65  2003/02/19 12:50:02  allanach
   New stop routine with varying p's implemented

   Revision 1.64  2003/02/18 14:53:09  allanach
   Deleted increased accuracy (forever) of running in the event of
   non-convergence, and made quark mixing more precise

   Revision 1.63  2003/01/31 14:56:45  allanach
   Debugged 2-loop QCD contribution to mt

   Revision 1.62  2002/12/10 17:39:11  allanach
   Added 2-loop higgs calls

   Revision 1.61  2002/12/09 17:09:44  allanach
   Added electroweak corrections to stop masses (unchecked)

   Revision 1.60  2002/12/06 17:03:52  allanach
   Added various safety checks vs nan in physical mass determinations

   Revision 1.58  2002/12/05 17:26:50  allanach
   Fixed PRINTOUT dependence and performed stop calculation at proper values
   of p

   Revision 1.57  2002/11/27 18:19:16  allanach
   New higgs mass calculation and bugfix in piZZ (included neutrino self
   energy diagrams)

   Revision 1.56  2002/11/26 13:28:21  allanach
   New mb calculation

   Revision 1.55  2002/11/19 17:00:00  allanach
   Added FCNC routine

   Revision 1.53  2002/11/11 16:58:50  allanach
   checking with Abdel...

   Revision 1.52  2002/11/07 13:15:26  allanach
   Bug-fixed (but not implemented) piAA

   Revision 1.50  2002/11/04 14:35:18  allanach
   Have reverted to old Higgs calculation while new one is debugged

   Revision 1.48  2002/10/28 14:42:46  allanach
   Bugchecked the pisisj functions

   Revision 1.47  2002/10/24 16:27:12  allanach
   Added pis1s2 and pis2s2 (they're not checked yet!)

   Revision 1.46  2002/10/24 14:19:36  allanach
   Wrote pi_s1s1 as a start to a better Higgs calculation

   Revision 1.45  2002/10/24 09:43:40  allanach
   Bug-fixed piAA: (s)lepton contribution added, and third family mixing
   incorporated properly

   Revision 1.43  2002/10/22 13:12:10  allanach
   Introduced new problem flag for infra-red quasi fixed points

   Revision 1.42  2002/10/17 13:11:39  allanach
   Been more careful with mgut fix: bounded by 10^15,10^17 GeV now.
   (you'll need to change if you want intermediate unification)

   Revision 1.41  2002/10/16 16:54:09  allanach
   Altered calculation of next mgut so that tachyons don't give a ridiculous
   one

   Revision 1.40  2002/10/16 10:17:53  allanach
   Added a nan trap in squark mass loop corrections in the case that gluino is
   degenerate with one of the squarks.

   Revision 1.39  2002/10/15 15:31:59  allanach
   Squark mass calculation now undertaken near msquark, to avoid problems with
   large logs

   Revision 1.38  2002/10/15 13:52:57  allanach
   m3sq set to begin with to the same as mu in the hope of better REWSB

   Revision 1.37  2002/10/14 14:17:30  allanach
   Added "runto" command in softpoint.x to get micromegas inputs to a
   different scale

   Revision 1.36  2002/10/14 12:41:34  allanach
   Added tachyon information when PRINTOUT > 2. Higgs tachyons now trapped.

   Revision 1.35  2002/10/11 15:26:26  allanach
   Micromegas interface changed

   Revision 1.34  2002/10/10 13:37:12  allanach
   Added physical parameters to Micromegas output

   Revision 1.33  2002/10/09 17:18:51  allanach
   Added new softpoint option to specify mu and mA instead of mH1, mH2. Some
   fine-tuning still necessary

   Revision 1.31  2002/10/01 11:52:16  allanach
   Bug-fixed bound-finding routines. MGUT always determined.

   Revision 1.30  2002/09/23 18:13:51  allanach
   Eigenvalue order with non symmetric diagonalisation setx

   Revision 1.29  2002/09/23 15:12:15  allanach
   Quark mixing systematised.

   Revision 1.27  2002/09/19 16:22:03  allanach
   Latest PDG data implemented

   Revision 1.25  2002/09/11 13:36:52  allanach
   Working version: RPV works with one-loop RPC corrections

   Revision 1.24  2002/09/11 09:52:34  allanach
   TOLERANCE etc now used in preference.

   Revision 1.23  2002/09/09 10:42:54  allanach
   TOLERANCE replaces EPS as being more user-friendly

   Revision 1.22  2002/09/04 13:54:54  allanach
   Gauge unification condition now imposed if mx is input as negative

   Revision 1.21  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, TOLERANCE out of def.h to make it quicker to
   compile once they are changed.

   Revision 1.20  2002/09/03 12:57:03  allanach
   Stopped possible nan's in chargino mass calculation

   Revision 1.19  2002/09/03 06:05:13  allanach
   Added trap for non-perturbativity: higgs mass no longer nan

   Revision 1.18  2002/08/30 12:31:56  allanach
   Test RPV version that works!

   Revision 1.17  2002/08/12 16:29:49  allanach
   display and set whole-object routines added

   Revision 1.16  2002/08/08 17:22:59  allanach
   Bug-fixed

   Revision 1.15  2002/07/30 12:57:32  allanach
   SOFTSUSY1.5

   Revision 1.14  2002/07/19 15:54:07  allanach
   SOFTSUSY1.5 version

   Revision 1.13  2002/07/03 14:45:59  allanach
   Negative mA^2 problem fixed

   Revision 1.12  2002/06/14 16:26:30  allanach
   Switches included for 2-loop running of scalar masses, and calulating mt
   at mt.

   Revision 1.10  2002/04/26 15:14:44  allanach
   Deleted all translation routines and defined boundary conditions within
   softsusy.h and softsusy.cpp

   Revision 1.8  2002/04/14 13:50:41  allanach
   Now use V=m3^2 H1 H2 instead of V=mu B H1 H2. It's more natural!

   Revision 1.5  2002/03/27 20:37:30  allanach
   pole-ized stop corrections

   Revision 1.4  2002/03/27 20:17:13  allanach
   pole-ized running top mass

   Revision 1.3  2002/03/27 18:01:09  allanach
   Pole-ized gluino corrections

   Revision 1.2  2002/03/27 17:53:30  allanach
   tadpole contributions pole-ized

   Revision 1.21  2002/03/18 18:48:15  allanach
   mb resummation made more accurate

   Revision 1.20  2002/03/01 17:33:58  allanach
   Added softpoint: allows you to call softsusy just from the command line.
   Also, interfaces to micromegas added. All interfaces write to input files:
   SoftsusyMicromegasInterface, SoftsusyIsawigInterface,
   SoftsusyIsajetInterface, SoftsusySsrunInterface respectively.

   Revision 1.19  2002/02/25 14:33:30  allanach
   Added conversion from MSbar to DRbar of bottom and tau Yukawa couplings.

   Revision 1.18  2002/02/14 17:16:17  allanach
   Resummed bottom mass corrections (important for high tan beta)

   Revision 1.16  2001/10/31 10:30:39  allanach
   Error trapping: checking for squark/slepton tachyons improved

   Revision 1.15  2001/10/25 11:51:30  allanach
   Taking sqrt(abs(slepton mass squareds)) in addCharginoLoops
   and addNeutralinoLoops because sometimes in AMSB, the slepton mass squareds
   are negative.

   Revision 1.14  2001/10/23 12:53:58  allanach
   Easier user interface utilised

   Revision 1.13  2001/10/22 09:19:28  allanach
   Fixed missing fermion contribution in piZZT. Turned some integers to
   decimals for performance on more platforms.

   Revision 1.12  2001/10/04 19:26:34  allanach
   New version deals with AMSB correctly

   Revision 1.10  2001/09/29 08:38:17  allanach
   Re-ordered REWSB constraints to be after calculation of physical masses
   fixed high tan beta bug

   Revision 1.9  2001/09/28 13:59:58  allanach
   Split rhohat determination up into sensible pieces - small bugs fixed there.

   Revision 1.8  2001/08/01 13:27:46  allanach
   SOFTSUSY1.2 revision

   Revision 1.7  2001/07/31 10:40:11  allanach
   Bug fix: radiative corrections to tachyons are now trapped

   Revision 1.6  2001/07/30 14:08:47  allanach
   Added ISAWIG and SSRUN interface

   Revision 1.5  2001/07/26 14:20:22  allanach
   Added isajet interface

   Revision 1.4  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#include "softsusy.h"
#include <cstring>

static double sw2 = 1.0 - sqr(MW / MZ),
  gnuL = 0.5,
  guL = 0.5 - 2.0 * sw2 / 3.0,
  gdL = -0.5 + sw2 / 3.0,
  geL = -0.5 + sw2,
  guR = 2.0 * sw2 / 3.0,
  gdR = -sw2 / 3.0,
  geR = -sw2, gsnuL = 0.5,
  yuL = 1.0 / 3.0,
  yuR = -4.0 / 3.0,
  ydL = 1.0 / 3.0,
  ydR = 2.0 / 3.0,
  yeL = -1.0,
  yeR = 2.0,
  ynuL = -1.0;

void initialise() {
  sw2 = 1.0 - sqr(MW / MZ);
  guL = 0.5 - 2.0 * sw2 / 3.0;
  gdL = -0.5 + sw2 / 3.0;
  geL = -0.5 + sw2;
  guR = 2.0 * sw2 / 3.0;
  gdR = -sw2 / 3.0;
  geR = -sw2;
}

const MssmSoftsusy & MssmSoftsusy::operator=(const MssmSoftsusy & s) {
  if (this == &s) return *this;
  setSoftPars(s.displaySoftPars());
  physpars = s.physpars;
  problem = s.problem;
  mw = s.mw;
  minV = s.minV;
  forLoops = s.forLoops;
  msusy = s.msusy;
  dataSet = s.dataSet;
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  return *this;
}

// Returns mu from rewsb requirement (TREE LEVEL is default, without deltaT2)
// returns 1 if there's a problem. Call at MSusy
int MssmSoftsusy::rewsbMu(int sgnMu, double & mu) const {
  int flag = 0;
   if (abs(sgnMu) != 1) {
    ostringstream ii;
    ii << "Error: sign mu = " << sgnMu << "\n";
    throw ii.str();
  }
  double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(), tanb =
    displayTanb(), tanb2 =  sqr(tanb);

  // Tree-level relation
  double musq = (mH1sq - mH2sq * tanb2)
    / (tanb2 - 1.0) - 0.5 * sqr(displayMz());

  if (musq < 0.0) {
    mu = sgnMu * sqrt(fabs(musq)); flag = 1; // mu has incorrect sign
  }
  else
    mu = sgnMu * sqrt(musq);

  return flag;
}

// returns 1 if mu < 1.0e-9
int MssmSoftsusy::rewsbM3sq(double mu, double & m3sq) const {
  int flag = 0;

  if (fabs(mu) < 1.0e-9)
    { flag = 1; m3sq = 0.0; }
  else
    m3sq =  0.5 *
         (displayMh1Squared() + displayMh2Squared() - displayTadpole2Ms() -
          displayTadpole1Ms() + 2.0 * sqr(mu)) *
      sin(2 * atan(displayTanb()));

  // Following means no good rewsb
  if (m3sq < 0.0) flag = 1;

  return flag;
}


// Predicts tan beta once mu and soft terms are predicted at low energy
// Useful for fine-tuning calculation. Call at MSusy only.
double MssmSoftsusy::predTanb() const  {
  double sin2t = 2.0 * displayM3Squared() /
    (displayMh1Squared() - displayTadpole1Ms() +
     displayMh2Squared() - displayTadpole2Ms() + 2.0 *
     sqr(displaySusyMu()));

  // Note: we want to take inverse sine so that fundamental domain is greater
  // than pi/4. sin(pi - 2 beta)=sin 2 beta should achieve this.
  // we also use tan (pi/2 - theta) = 1/tan(theta)
  double theta;
  if (fabs(sin2t) < 1.0) theta = asin(sin2t) * 0.5;
  else return 0.0;

  return 1.0 / tan(theta);
}

void MssmSoftsusy::doTadpoles(double mt, double sinthDRbar) {

    calcTadpole1Ms1loop(mt, sinthDRbar);
    calcTadpole2Ms1loop(mt, sinthDRbar);

    physpars.t1OV1Ms1loop = displayTadpole1Ms();
    physpars.t2OV2Ms1loop = displayTadpole2Ms();

    if (numRewsbLoops > 1) {
      // add the two-loop terms, prepare inputs
      double s1s = 0., s2s = 0., s1t = 0., s2t = 0.,
        gs = displayGaugeCoupling(3),
        rmtsq = sqr(forLoops.mt), scalesq = sqr(displayMu()),
        vev2 = sqr(displayHvev()), tanb = displayTanb(),
        amu = -displaySusyMu(), mg = displayGaugino(3),
        mAsq = sqr(forLoops.mhiggs(2));

      double sxt = sin(forLoops.thetat), cxt = cos(forLoops.thetat);
      double mst1sq = sqr(forLoops.mu(1, 3)),
        mst2sq = sqr(forLoops.mu(2, 3));
      // two-loop Higgs corrections: alpha_s alpha_b
      double sxb = sin(forLoops.thetab),
        cxb = cos(forLoops.thetab);
      double sintau = sin(forLoops.thetatau),
        costau = cos(forLoops.thetatau);
      double msnusq = sqr(forLoops.msnu(3));
      double msb1sq = sqr(forLoops.md(1, 3)),
        msb2sq = sqr(forLoops.md(2, 3));
      double mstau1sq = sqr(forLoops.me(1, 3)),
        mstau2sq = sqr(forLoops.me(2, 3));
      double cotbeta = 1.0 / tanb;
      double rmbsq = sqr(forLoops.mb);
      double rmtausq = sqr(forLoops.mtau);
      double s1b = 0.0, s2b = 0.0, s1tau = 0.0, s2tau = 0.0;

      ewsb2loop_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq,
                 &amu, &tanb, &vev2, &gs, &s1s, &s2s);
      ddstad_(&rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s1t,
              &s2t);
      ewsb2loop_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq,
                 &amu, &cotbeta, &vev2, &gs, &s2b, &s1b);
      tausqtad_(&rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
                &costau, &scalesq, &amu, &tanb, &vev2, &s1tau, &s2tau);

      if (!testNan(s1s * s1t * s1b * s1tau * s2s * s2t * s2b * s2tau)) {
        physpars.t1OV1Ms = physpars.t1OV1Ms1loop - s1s - s1t - s1b - s1tau;
        physpars.t2OV2Ms = physpars.t2OV2Ms1loop - s2s - s2t - s2b - s2tau;
        // end of 2-loop bit
      }
      else  {
        flagNoMuConvergence(true);
        if (PRINTOUT > 1) cout << "2-loop tadpoles are nans\n";
      }
    }
}



// From hep-ph/9606211's appendix. It should be done at MSusy to minimize the
// 1-loop contributions. Only call if you've calculated drbarpars.
// inputs are running top/bottom masses: call at MSusy only
double MssmSoftsusy::doCalcTadpole1oneLoop(double mt, double sinthDRbar) {
 if (forLoops.mu(1, 3) == 0.0 || forLoops.mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
         << " the DRbar masses.\n";
   return 0.0;
  }

  double g = displayGaugeCoupling(2),
    costhDRbar = cos(asin(sinthDRbar)),
    tanb = displayTanb(), cosb = cos(atan(tanb)),
    ht = displayYukawaElement(YU, 3, 3),
    mu = -displaySusyMu(), q = displayMu(),
    hb = displayYukawaElement(YD, 3, 3),
    hbsq = sqr(hb), mz = displayMzRun();
  double beta = atan(displayTanb());
  double v1 = displayHvev() * cos(beta);

  // sneutrino coupling
  double lSnu = g * mz / costhDRbar * 0.5 * cosb;

  // stop couplings
  DoubleMatrix lTS1Lr(2, 2), lTS112(2, 2), rotate(2, 2);
  lTS1Lr(1, 1) = g * mz / costhDRbar * guL * cosb;
  lTS1Lr(1, 2) = ht * mu / root2;
  lTS1Lr(2, 1) = lTS1Lr(1, 2);
  lTS1Lr(2, 2) = g * mz / costhDRbar * guR * cosb;
  rotate = rot2d(forLoops.thetat);

  lTS112 = rotate * lTS1Lr * rotate.transpose();

  // sbottom couplings
  DoubleMatrix lBS1Lr(2, 2), lBS112(2, 2);
  double mA = forLoops.mhiggs(2), mH0 = forLoops.mhiggs(3), mHp =
    forLoops.mhiggs(4);
  double mb = forLoops.mb;
  double mtau = forLoops.mtau;

  lBS1Lr(1, 1) = g * mz / costhDRbar * gdL * cosb + hbsq * v1;
  lBS1Lr(1, 2) = displayTrilinear(DA, 3, 3) / root2;
  lBS1Lr(2, 1) = lBS1Lr(1, 2);
  lBS1Lr(2, 2) = g * mz / costhDRbar * gdR * cosb + hbsq * v1;

  rotate = rot2d(forLoops.thetab);

  lBS112 = rotate * lBS1Lr * rotate.transpose();

  // stau couplings
  DoubleMatrix lTauS1Lr(2, 2), lTauS112(2, 2);
  double htau = displayYukawaElement(YE, 3, 3), htausq = sqr(htau);

  lTauS1Lr(1, 1) = g * mz / costhDRbar * geL * cosb + htausq * v1;
  lTauS1Lr(1, 2) = displayTrilinear(EA, 3, 3) / root2;
  lTauS1Lr(2, 1) = lTauS1Lr(1, 2);
  lTauS1Lr(2, 2) = g * mz / costhDRbar * geR * cosb + htausq * v1;

  rotate = rot2d(forLoops.thetatau);

  lTauS112 = rotate * lTauS1Lr * rotate.transpose();

  // bottom quark and tau, ignore others - factor (10^-2)^3 down
  // I have included the bottom pole mass in the propagators and the Yukawa
  // for the coupling, hence BPMZ's hb is written mb * root2 / v1
  double fermions = - 6.0 * hb * mb * root2 / v1 * a0(mb, q)
    - 2.0 * htau * mtau * root2 / v1 * a0(mtau, q);

  double gO2mwcosb = g / (2.0 * displayMwRun() * cosb);
  double stops = 0., sbots = 0.;
  // third generation squarks
  stops = stops + 3.0 * gO2mwcosb * lTS112(1, 1) *
    a0(forLoops.mu(1, 3), q);
  stops = stops + 3.0 * gO2mwcosb * lTS112(2, 2) *
    a0(forLoops.mu(2, 3), q);
  sbots = sbots + 3.0 * gO2mwcosb * lBS112(1, 1) *
    a0(forLoops.md(1, 3), q);
  sbots = sbots + 3.0 * gO2mwcosb * lBS112(2, 2) *
    a0(forLoops.md(2, 3), q);

  // third generation sleptons
  double staus = gO2mwcosb * lTauS112(1, 1) *
    a0(forLoops.me(1, 3), q);
  staus = staus + gO2mwcosb * lTauS112(2, 2) *
    a0(forLoops.me(2, 3), q);

  double sneuts = 0.0;
  sneuts = sneuts + gO2mwcosb * lSnu * a0(forLoops.msnu(3), q);

  int family;
  double gsqMZO2mwcw = sqr(g) * mz * 0.5 / (displayMwRun() * costhDRbar);
  double sups = 0.0, sdowns = 0.0;
  // first two families of squarks
  for (family = 1; family <=2; family++) {
    sups = sups + 3.0 * gsqMZO2mwcw *
      ( guL * a0(forLoops.mu(1, family), q) +
        guR * a0(forLoops.mu(2, family), q));
    sdowns = sdowns + 3.0 * gsqMZO2mwcw *
      ( gdL * a0(forLoops.md(1, family), q) +
        gdR * a0(forLoops.md(2, family), q));
  }

  double sleps = 0.;
  // sleptons
  for (family = 1; family <=2; family++) {
    sleps = sleps + gsqMZO2mwcw *
      (geL * a0(forLoops.me(1, family), q) +
       geR * a0(forLoops.me(2, family), q));
    sneuts = sneuts + gO2mwcosb * lSnu * a0(forLoops.msnu(family), q);
  }

  // Higgs
  double alpha = forLoops.thetaH, sina2 = sqr(sin(alpha)),
    cosa2 = 1.0 - sina2, cos2b = cos(2.0 * atan(tanb)),
    costhDRbar2 = sqr(costhDRbar),
    mh = forLoops.mhiggs(1), sin2a = sin(2.0 * alpha);
  double higgs = 0.0;
  higgs = higgs - sqr(g) * cos2b / (8.0 * costhDRbar2) *
    (a0(mA, q) + 2.0 * a0(mHp, q)) +
    sqr(g) * a0(mHp, q) * 0.5 +
    sqr(g) / (8.0 * costhDRbar2) * a0(mh, q) *
    (3.0 * sina2 - cosa2 + sin2a * tanb) +
    sqr(g) / (8.0 * costhDRbar2) * a0(mH0, q) *
    (3.0 * cosa2 - sina2 - sin2a * tanb);

  // Neutralinos
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  double neutralinos = 0.0;
  double tanthDRbar = tan(acos(costhDRbar));
  for (family = 1; family <= 4; family++)
    neutralinos = neutralinos -
      sqr(g) * mneut(family) / (displayMwRun() * cosb) *
      (n(family, 3) * (n(family, 2) - n(family, 1) * tanthDRbar)).real() *
      a0(mneut(family), q);

  // Charginos
  double charginos = 0.0;
  for (family=1; family<=2; family++)
    charginos = charginos - root2 * sqr(g) / (displayMwRun() * cosb)
      * mch(family) * (v(family, 1) * u(family, 2)).real()
         * a0(mch(family), q);

  // Weak bosons
  double gaugeBosons = 0.0;
  gaugeBosons = gaugeBosons + 3.0 * sqr(g) / 4.0 *
    (2.0 * a0(displayMwRun(), q) + a0(mz, q) / costhDRbar2)
    + sqr(g) * cos2b / (8.0 * costhDRbar2) * (2.0 * a0(displayMwRun(), q) +
                                              a0(mz, q));

  double sfermions = stops + sbots + staus + sneuts + sleps + sups + sdowns;

  double delta = fermions + sfermions + higgs + charginos + neutralinos +
    gaugeBosons;

  return delta / (16.0 * sqr(PI));
}

void MssmSoftsusy::calcTadpole1Ms1loop(double mt, double sinthDRbar) { // checked

  physpars.t1OV1Ms1loop = doCalcTadpole1oneLoop(mt, sinthDRbar);
  if (testNan(physpars.t1OV1Ms1loop)) {
    physpars.t1OV1Ms1loop = 0.0;
    flagNoMuConvergence(true);
  }

  physpars.t1OV1Ms = physpars.t1OV1Ms1loop;
}

double MssmSoftsusy::displayMwRun() const {
  double costhDRbar = cos(asin(calcSinthdrbar()));
  return displayMzRun() * costhDRbar;
}

// From hep-ph/9311269's appendix. It should be done at MSusy to minimize the
// 1-loop contributions. Only call if you've calculated physpars
// inputs are running top/bottom masses. Call at MSusy
double MssmSoftsusy::doCalcTadpole2oneLoop(double mt, double sinthDRbar) {
// CHECKED
 if (forLoops.mu(1, 3) == 0.0 || forLoops.mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
         << " the DRbar masses.\n";
   return 0.0;
  }

  double g = displayGaugeCoupling(2),
    costhDRbar = cos(asin(sinthDRbar)),
    tanb = displayTanb(),
    sinb = sin(atan(tanb)),
    ht = displayYukawaElement(YU, 3, 3),
    htsq = sqr(ht),
    mu = -displaySusyMu(), q = displayMu(),
    hb = displayYukawaElement(YD, 3, 3), mz = displayMzRun();
  double beta = atan(displayTanb());
  double v2 = displayHvev() * sin(beta);
  double mtop = forLoops.mt;

  // Stop couplings
  DoubleMatrix lTS2Lr(2, 2), lTS212(2, 2), rotate(2, 2);

  lTS2Lr(1, 1) = - g * mz / costhDRbar * guL * sinb + htsq * v2;
  lTS2Lr(1, 2) = displayTrilinear(UA, 3, 3) / root2;
  lTS2Lr(2, 1) = lTS2Lr(1, 2);
  lTS2Lr(2, 2) = - g * mz / costhDRbar * guR * sinb + htsq * v2;

  rotate = rot2d(forLoops.thetat);

  lTS212 = rotate * lTS2Lr * rotate.transpose();

  // sneutrino coupling
  double lSnu = - g * mz / costhDRbar * 0.5 * sinb;

  // sbottom couplings
  DoubleMatrix lBS2Lr(2, 2), lBS212(2, 2);
  double mA = forLoops.mhiggs(2), mH0 = forLoops.mhiggs(3), mHp =
    forLoops.mhiggs(4);

  lBS2Lr(1, 1) = - g * mz / costhDRbar * gdL * sinb;
  lBS2Lr(1, 2) = hb * mu / root2;
  lBS2Lr(2, 1) = lBS2Lr(1, 2);
  lBS2Lr(2, 2) = - g * mz / costhDRbar * gdR * sinb;

  rotate = rot2d(forLoops.thetab);

  lBS212 = rotate * lBS2Lr * rotate.transpose();

  // stau couplings
  DoubleMatrix lTauS2Lr(2, 2), lTauS212(2, 2);
  double htau = displayYukawaElement(YE, 3, 3);

  lTauS2Lr(1, 1) = - g * mz / costhDRbar * geL * sinb;
  lTauS2Lr(1, 2) = htau * mu / root2;
  lTauS2Lr(2, 1) = lTauS2Lr(1, 2);
  lTauS2Lr(2, 2) = - g * mz / costhDRbar * geR * sinb;

  rotate = rot2d(forLoops.thetatau);

  lTauS212 = rotate * lTauS2Lr * rotate.transpose();

  double delta = 0.0;

  // top quark, ignore others - factor (10^-2)^3 down
  double fermions = - 6.0 * sqr(ht) * a0(mtop, q);

  // third generation squarks
  double sfermions = 3.0 * g * lTS212(1, 1) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.mu(1, 3), q);
  sfermions = sfermions + 3.0 * g * lTS212(2, 2) /
    (2.0 * displayMwRun() * sinb) *
    a0(forLoops.mu(2, 3), q);
  sfermions = sfermions + 3.0 * g * lBS212(1, 1) /
    (2.0 * displayMwRun() * sinb) *
    a0(forLoops.md(1, 3), q);
  sfermions = sfermions + 3.0 * g * lBS212(2, 2) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.md(2, 3), q);
  // third generation sleptons
  sfermions = sfermions + g * lTauS212(1, 1) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.me(1, 3), q);
  sfermions = sfermions + g * lTauS212(2, 2) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.me(2, 3), q);
  sfermions = sfermions +
    g * lSnu / (2.0 * displayMwRun() * sinb) * a0(forLoops.msnu(3), q);

  int family;
  // first two families of squarks
  for (family = 1; family <=2; family++)
    {
      sfermions = sfermions + 3.0 * g *
        (- g * mz / costhDRbar * guL * sinb * a0(forLoops.mu(1, family), q)
         - g * mz / costhDRbar * guR * sinb * a0(forLoops.mu(2, family), q))
        / (2.0 * displayMwRun() * sinb);
      sfermions = sfermions + 3.0 * g *
        (- g * mz / costhDRbar * gdL * sinb * a0(forLoops.md(1, family), q)
         - g * mz / costhDRbar * gdR * sinb * a0(forLoops.md(2, family), q))
        / (2.0 * displayMwRun() * sinb);
    }

  // sleptons
  for (family = 1; family <=2; family++)
    {
      sfermions = sfermions + g *
        (- g * mz / costhDRbar * geL * sinb * a0(forLoops.me(1, family), q)
         - g * mz / costhDRbar * geR * sinb * a0(forLoops.me(2, family), q))
        / (2.0 * displayMwRun() * sinb);
      sfermions = sfermions + g * lSnu * a0(forLoops.msnu(family), q)
        / (2.0 * displayMwRun() * sinb);
    }

  // Higgs
  double alpha = forLoops.thetaH, sina2 = sqr(sin(alpha)), cosa2 = 1.0 -
    sina2, cos2b = cos(2.0 * atan(tanb)), costhDRbar2 = sqr(costhDRbar),
    mh = forLoops.mhiggs(1), sin2a = sin(2.0 * alpha);
  double higgs = 0.0;
  higgs = sqr(g) * cos2b / (8.0 * costhDRbar2) *
    (a0(mA, q) + 2.0 * a0(mHp, q)) +
    sqr(g) * a0(mHp, q) * 0.5 +
    sqr(g) / (8.0 * costhDRbar2) * a0(mh, q) *
    (3.0 * cosa2 - sina2 + sin2a / tanb) +
    sqr(g) / (8.0 * costhDRbar2) * a0(mH0, q) *
    (3.0 * sina2 - cosa2 - sin2a / tanb);

  // Neutralinos
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  double tanthDRbar = tan(acos(costhDRbar));
  double neutralinos = 0.0;
  for (family = 1; family <= 4; family++)
    neutralinos = neutralinos + sqr(g) * mneut(family) / (displayMwRun() * sinb) *
      (n(family, 4) * (n(family, 2) - n(family, 1) * tanthDRbar)).real() *
      a0(mneut(family), q);
  double charginos = 0.0;
  for (family = 1; family <= 2; family++)
    charginos = charginos - root2 * sqr(g) * mch(family) /
      (displayMwRun() * sinb)
      * (v(family, 2) * u(family, 1)).real() * a0(mch(family), q);

  // Weak bosons
  double gaugeBosons = 3.0 * sqr(g) / 4.0 *
    (2.0 * a0(displayMwRun(), q) + a0(mz, q) / costhDRbar2)
    - sqr(g) * cos2b / (8.0 * costhDRbar2) *
    (2.0 * a0(displayMwRun(), q) + a0(mz, q));

  delta = fermions + sfermions + higgs + charginos + neutralinos +
    gaugeBosons;

  return delta / (16.0 * sqr(PI));
}

void MssmSoftsusy::calcTadpole2Ms1loop(double mt, double sinthDRbar) {// CHECKED
  physpars.t2OV2Ms1loop = doCalcTadpole2oneLoop(mt, sinthDRbar);
  if (testNan(physpars.t2OV2Ms1loop)) {
    flagNoMuConvergence(true);
    physpars.t2OV2Ms1loop = 0.0;
  }
  physpars.t2OV2Ms = physpars.t2OV2Ms1loop;
}

// Apply at scale MSusy: checked 19.12.2000
// Displays PHYSICAL MZ, ie MZ(q) - piZz^T(q)
double MssmSoftsusy::predMzsq(double & tanb) const {
  tanb = predTanb();

  double MZsq = 2.0 *
    ((displayMh1Squared() - displayTadpole1Ms() -
      (displayMh2Squared() - displayTadpole2Ms()) *
      sqr(tanb)) /
     (sqr(tanb) - 1.0) - sqr(displaySusyMu())) -
    piZZT(displayMz(), displayMu());

  return MZsq;
}

// Used to get useful information into ftCalc
static MssmSoftsusy *tempSoft1;
static int ftFunctionality;
static DoubleVector ftPars(3);
void (*ftBoundaryCondition)(MssmSoftsusy &, const DoubleVector &);

double ftCalc(double x) {
  // Stores running parameters in a vector
  DoubleVector storeObject(tempSoft1->display());
  double initialMu = tempSoft1->displayMu();
  drBarPars saveDrBar(tempSoft1->displayDrBarPars());

  if (PRINTOUT > 1) cout << "#";

  if (ftFunctionality <= ftPars.displayEnd()) {
    ftPars(ftFunctionality) = x;
    ftBoundaryCondition(*tempSoft1, ftPars);
    if (PRINTOUT > 1) cout << "p" << ftFunctionality << "=";
  }
  else if (ftFunctionality == ftPars.displayEnd() + 1) {
    tempSoft1->setSusyMu(x); if (PRINTOUT > 1) cout << "mu= ";
  }
  else if (ftFunctionality == ftPars.displayEnd() + 2) {
    tempSoft1->setM3Squared(x); if (PRINTOUT > 1) cout << "m3sq= ";
  }
  else if (ftFunctionality == ftPars.displayEnd() + 3) {
    tempSoft1->setYukawaElement(YU, 3, 3, x); if (PRINTOUT > 1) cout << "ht= ";
  }
  else {
    ostringstream ii;
    ii << "MssmSoftsusy:ftCalc called with incorrect functionality=" <<
      ftFunctionality << endl;
    throw ii.str();
  }

  double referenceMzsq, predTanb;

  tempSoft1->runto(tempSoft1->calcMs());
  tempSoft1->calcDrBarPars();
  tempSoft1->runto(tempSoft1->calcMs());

  tempSoft1->setHvev(tempSoft1->getVev());
  double mt = tempSoft1->displayDrBarPars().mt;
  double sinthDRbar = tempSoft1->calcSinthdrbar();
  // We miss two-loop terms in our calculation of fine-tuning...
  tempSoft1->calcTadpole2Ms1loop(mt, sinthDRbar);
  tempSoft1->calcTadpole1Ms1loop(mt, sinthDRbar);

  referenceMzsq = tempSoft1->predMzsq(predTanb);

  if (PRINTOUT > 1) cout << x << " MZ=" << sqrt(fabs(referenceMzsq))
                         << " tanb=" << predTanb << "\n";

  // Restore initial parameters at correct scale
  tempSoft1->setMu(initialMu);
  tempSoft1->set(storeObject);
  tempSoft1->setDrBarPars(saveDrBar);

  return referenceMzsq;
}

// Give it a GUT scale object consistent with rewsb
// and it'll return the fine tuning by varying each parameter in turn then mu
// and m3sq at the high scale
double MssmSoftsusy::it1par(int numPar, const DoubleVector & bcPars) {

  double ftParameter = 0.0, err, h = 0.01;

  tempSoft1 = this;

  // Stores running parameters in a vector
  DoubleVector storeObject(display());
  double initialMu = displayMu();
  sPhysical savePhys(displayPhys());

  ftFunctionality = numPar;
  double x;

  // Defines starting value to calculate derivative from
  if (numPar > 0 && numPar <= bcPars.displayEnd()) {
    x = bcPars.display(numPar); h = 0.01 * x;
  }
  else if (numPar == bcPars.displayEnd() + 1) {//mu
    x = displaySusyMu(); h = 0.01 * x;
  }
  else if (numPar == bcPars.displayEnd() + 2) {//m3sq
    x = displayM3Squared(); h = 0.01 * x;
  }
  else if (numPar == bcPars.displayEnd() + 3) {//ht
    x = displayYukawaElement(YU, 3, 3); h = x * 0.0005;
  }
  else {
    ostringstream ii;
    ii << "it1par called with functionality " << ftFunctionality <<
      " out of range.\n";
    throw ii.str();
  }

  ftPars.setEnd(bcPars.displayEnd()); ftPars = bcPars;

  // Fine tuning is a / MZ^2 d MZ^2 / d a for any parameter a
  double derivative;
  if (fabs(x) < 1.0e-10) {
    ftParameter = 0.0; derivative = 0.0; err = 0.0;
  }
  else {
    derivative = calcDerivative(ftCalc, x, h, &err);
    ftParameter = x * derivative / sqr(displayMz());
  }

  if (PRINTOUT > 1)
    cout << "derivative=" << derivative << " error=" << err << "\n";

  // High error: if can't find a derivative, error comes back with 1.0e30
  if (ftParameter > TOLERANCE && fabs(err / derivative) > 1.0)
    return 6.66e66;

  // Restore initial parameters at correct scale
  setMu(initialMu);
  set(storeObject);
  setPhys(savePhys);

  return fabs(ftParameter);
}

// Pass it an object and it'll return the fine tuning parameters
DoubleVector MssmSoftsusy::fineTune
(void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &),
 const DoubleVector  & bcPars, double mx, bool doTop) {

  // Stores running parameters in a vector
  DoubleVector storeObject(display());
  double initialMu = displayMu();
  sPhysical savePhys(displayPhys());

  runto(mx);

  ftBoundaryCondition = boundaryCondition;

  int numPars = bcPars.displayEnd() + 2;
  if (doTop) numPars = numPars + 1;

  DoubleVector tempFineTuning(numPars);

  int i; for (i = 1; i<= numPars; i++) {
    tempFineTuning(i) = it1par(i, bcPars);
    // flag problem FT calculation with NaN...too inaccurate
    if (tempFineTuning(i) > 1.0e66) tempFineTuning(i) = asin(2.0);
  }

  // Restore initial parameters at correct scale
  setMu(initialMu);
  set(storeObject);
  setPhys(savePhys);

  return tempFineTuning;
}

// Obtains solution of one-loop effective potential minimisation via iteration
// technique
// err is 1 if no iteration reached
// 2 if incorrect rewsb
void MssmSoftsusy::iterateMu(double & muold, int sgnMu,
                             double mt, int maxTries, double pizzMS,
                             double sinthDRbar, double tol, int & err) {
  static int numTries = 0;
  static double munew = 0.;

  if (numTries - 1 > maxTries) {
    if (PRINTOUT) cout << "iterateMu reached maxtries\n";
    numTries = 0; munew = 0.0;
    err = 1; return;
  }
  // How close to convergence are we?
  double c = 1.0 - minimum(fabs(muold), fabs(munew)) /
    maximum(fabs(muold), fabs(munew));

  if (PRINTOUT > 2) cout << " diff=" << c;

  if (c < tol) {
    muold = munew; //err = 0;
    numTries = 0; munew = 0.0;
    if (PRINTOUT > 2) cout << " mu converged\n";
    return;
  }

  numTries = numTries + 1;

  muold = munew;

  double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(),
    tanb = displayTanb(), tanb2 = sqr(tanb);

  double treeLevelMusq = (mH1sq - mH2sq * tanb2) / (tanb2 - 1.0) - 0.5 *
    sqr(displayMz());

  try {
  calcDrBarPars();
  double oneLoopMusq = treeLevelMusq;
  // calculate the new one-loop tadpoles with old value of mu
  if (numRewsbLoops > 0) {
    doTadpoles(mt, sinthDRbar);
    oneLoopMusq = treeLevelMusq - 0.5 * pizzMS +
      (displayTadpole2Ms() * sqr(tanb) - displayTadpole1Ms()) /
      (sqr(tanb) - 1.0);
  }
  // Error in rewsb
  if (oneLoopMusq < 0.0) {
    err = 2;
    if (PRINTOUT > 1) cout << " mu^2<0 ";
  }

  munew = sgnMu * sqrt(fabs(oneLoopMusq));
  setSusyMu(munew); double m3sqnew;
  if (rewsbM3sq(munew, m3sqnew) == 0) flagB(false);
  else flagB(true);
  setM3Squared(m3sqnew);
  }
  catch(const char *a) {
    numTries = 0;
    throw a;
  }
  catch(const string &a) {
    numTries = 0;
    throw a;
  }

  if (PRINTOUT > 2) cout << " mu=" << munew;

  iterateMu(muold, sgnMu, mt, maxTries, pizzMS, sinthDRbar, tol, err);
}

void MssmSoftsusy::rewsbTreeLevel(int sgnMu) {
  double mu, m3sq;
  if (rewsbMu(sgnMu, mu)) flagMusqwrongsign(true);
  else flagMusqwrongsign(false);
  setSusyMu(mu);
  if (rewsbM3sq(mu, m3sq)) flagB(true);
  else flagB(false);

  setM3Squared(m3sq);

  if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
       displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
    flagHiggsufb(true);
  else
    flagHiggsufb(false);

  return;
}

// Organises rewsb: call it at the low scale MS^2 = sqrt(0.5 * (mT1^2 +
// mT2^2)) is best, or below if it's decoupled from there.
// Call with zero, or no mt if you want tree level
void MssmSoftsusy::rewsb(int sgnMu, double mt) {
  double munew, m3sqnew;

  double sinthDRbarMS = calcSinthdrbar();

  calcTadpole1Ms1loop(mt, sinthDRbarMS);
  calcTadpole2Ms1loop(mt, sinthDRbarMS);

  munew = displaySusyMu();

  // Iterate to get a self-consistent mu solution
  int maxTries = 20, err = 0;
  double tol = TOLERANCE * 1.0e-4;

  double pizztMS = sqr(displayMzRun()) - sqr(displayMz()); // resums logs
  //double pizztMS = piZZT(displayMz(), displayMu());

  iterateMu(munew, sgnMu, mt, maxTries, pizztMS, sinthDRbarMS,
            tol, err);

  if (err == 2) flagMusqwrongsign(true);
  else flagMusqwrongsign(false);
  if (err == 1) flagNoMuConvergence(true);
  else
    setSusyMu(munew);

  if (rewsbM3sq(munew, m3sqnew) == 0) flagB(false);
  else flagB(true);

  setM3Squared(m3sqnew);

  if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
       displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
    flagHiggsufb(true);
  else
    flagHiggsufb(false);
}

#define HR "----------------------------------------------------------"

ostream & operator <<(ostream &left, const MssmSoftsusy &s) {
  left << HR << endl;
  left << "Gravitino mass M3/2: " << s.displayGravitino() << endl;
  left << "Msusy: " << s.displayMsusy() << " MW: " << s.displayMw() << endl;
  left << "Data set:\n" << s.displayDataSet();
  left << HR << endl;
  left << s.displaySoftPars();
  left << HR << "\nPhysical MSSM parameters:\n";
  left << s.displayPhys();
  double mass; int posi, posj, id;
  id = s.lsp(mass, posi, posj);
  left << "lsp is ";
  recogLsp(id, posj);
  left << " of mass " << mass << " GeV\n";
  if (s.displayProblem().test()) left << "***** PROBLEM *****" <<
                                   s.displayProblem() << " *****" << endl;
  left << HR << endl;
  return left;
}

#undef HR

// Gives a summary of important properties of a SUSY object: mu, m3sq, mH0,
// mChi0, lightest stau, lightest mGluino, lightest stop, lightest chargino,
// lightest tau sneutrino, lightest sbottom, down squark, up squark and
// selectron masses, minimum of potential (if calculated) and fine-tuning
// parameter (if passed)
void printShortInitialise() {
  cout <<
    "     mu     " << "   m3sq     " <<
    "   mstau1   " << "   msbott   " << "  mstop1    " <<
    " mSquark   " << "   mSlep   " << "  msnu1     " <<
    "    mhl0    " << "   mchi10   " << " mchargino1 " <<
    "  mGluino  " << flush;
  }

// Prints mu, B and important spectral information
void MssmSoftsusy::printShort() const {

  const double problemFlag = -1.0;

  double mu, m3sq;
  mu = displaySusyMu(); m3sq = displayM3Squared();
  int pos;
  sPhysical s(this->displayPhys());

  if (displayProblem().muSqWrongSign || displayProblem().b ||
      displayProblem().higgsUfb || displayProblem().noRhoConvergence)
    mu = -1.0 *  mu;
  if (displayProblem().tachyon) m3sq = -1.0 * m3sq;

  if (displayProblem().nonperturbative || displayProblem().noMuConvergence) {
    cout << " ";
    int i; for (i=1; i<=12; i++) cout << problemFlag << " ";
  }
    else {
      cout << " " << mu << " " << m3sq << " "
           << minimum(s.me(1, 3), s.me(2, 3)) << " "
           << minimum(s.md(1, 3), s.md(2, 3)) << " "
           << minimum(s.mu(1, 3), s.mu(2, 3)) << " "
         << minimum(minimum(s.md(1, 1), s.md(2, 1)),
                    minimum(s.mu(1, 1), s.mu(2, 1)))
           << " " << minimum (s.me(1, 1), s.me(2, 1)) << " "
           << s.msnu.min(pos) << " "
           << s.mhiggs(1) << " "
           << (s.mneut.apply(fabs)).min(pos) << " "
           << (s.mch.apply(fabs)).min(pos) << " "
           << s.mGluino << " ";
    }

  cout << flush;
}

void MssmSoftsusy::printLong() {
  // output:
  //  1  2     3      4   5   6   7    8  9  10 11   12    13
  // mu  m3sq mH1sq mH2sq g1 g2 mt(mt) mh mA mH mH+ alphaH msnu3
  //  14     15     16     17    18     19    20
  // msnu1 mstopL mstopR msupL msupR msbotL msbotR
  //  21   22   23     24     25     26     27    28     29
  // msdL msdR mstauL mstauR mselL mselR thetat thetab thetatau
  // 30   31   32   33   34      35      36     37     38
  // mgl mch1 mch2 thetaL thetaR mneut1 mneut2 mneut3 mneut4
  //   39    40     41
  // sinthW t1ov1 t2ov2
  double mu = displaySusyMu();
  if (displayProblem().muSqWrongSign || displayProblem().b ||
      displayProblem().higgsUfb)
    mu = -1.0 *  mu;
  //  if (displayProblem().noConvergence)
  // mu = asin(-2.0);

  if (displayProblem().tachyon) {
    drBarPars s(displayDrBarPars());
    s.me = -1.0 * s.me;
    s.md = -1.0 * s.md;
    s.mu = -1.0 * s.mu;
    s.msnu = -1.0 * s.msnu;
    s.mhiggs = -1.0 * s.mhiggs;
    setDrBarPars(s);
  }

  cout <<" " << mu << " "
       << displayM3Squared() << " "
       << displayMh1Squared() << " " <<
    displayMh2Squared() << " " <<
    displayGaugeCoupling(1) << " " <<
    displayGaugeCoupling(2) << " " <<
    calcRunningMt() << " " <<
    displayPhys().mhiggs(1) << " " <<
    displayPhys().mhiggs(2) << " " <<
    displayPhys().mhiggs(3) << " " <<
    displayPhys().mhiggs(4) << " " <<
    displayPhys().thetaH << " " <<
    displayPhys().msnu(3) << " " <<
    displayPhys().msnu(1) << " " <<
    displayPhys().mu(1, 3) << " " <<
    displayPhys().mu(2, 3) << " " <<
    displayPhys().mu(1, 1) << " " <<
    displayPhys().mu(2, 1) << " " <<
    displayPhys().md(1, 3) << " " <<
    displayPhys().md(2, 3) << " " <<
    displayPhys().md(1, 1) << " " <<
    displayPhys().md(2, 1) << " " <<
    displayPhys().me(1, 3) << " " <<
    displayPhys().me(2, 3) << " " <<
    displayPhys().me(1, 1) << " " <<
    displayPhys().me(2, 1) << " " <<
    displayPhys().thetat << " " <<
    displayPhys().thetab << " " <<
    displayPhys().thetatau << " " <<
    displayPhys().mGluino << " " <<
    displayPhys().mch(1) << " " <<
    displayPhys().mch(2) << " " <<
    displayPhys().thetaL << " " <<
    displayPhys().thetaR << " " <<
    displayPhys().mneut(1) << " " <<
    displayPhys().mneut(2) << " " <<
    displayPhys().mneut(3) << " " <<
    displayPhys().mneut(4) << " " <<
    calcSinthdrbar() << " " <<
    displayTadpole1Ms() << " " <<
    displayTadpole2Ms() << " ";

  if (displayProblem().test()) cout << "%" << displayProblem();
  cout << endl << flush;
}



void MssmSoftsusy::higgs(int accuracy, double piwwtMS, double pizztMS) {

  double tanb = displayTanb();
  double beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  double sinb2 = sqr(sinb), cosb2 = sqr(cosb), mz = displayMz(),
    mz2 = sqr(mz) + pizztMS;
  double mApole = physpars.mhiggs(2); // physical value
  //  double mApole2 = sqr(mApole);

  // tree level: there'll be trouble if B has the opp sign to mu
  double mAsq =  (displayMh2Squared() - displayTadpole2Ms() -
                  displayMh1Squared() + displayTadpole1Ms()) /
    cos(2.0 * beta) - mz2;

  DoubleMatrix mHtree(2, 2);
  mHtree(1, 1) = mAsq * sinb2 + mz2 * cosb2;
  mHtree(1, 2) = - sinb * cosb * (mAsq + mz2);
  mHtree(2, 2) = mAsq * cosb2 + mz2 * sinb2;
  mHtree(2, 1) = mHtree(1 ,2);

  DoubleMatrix mhAtmh(mHtree), mhAtmH(mHtree);
  DoubleMatrix sigmaMh(2, 2), sigmaMH(2, 2);

  double q = displayMu(), p;
  double p2s = 0.0, p2w = 0.0, p2b = 0.0, p2tau = 0.0, dMA = 0.;
  // radiative corrections:
  if (accuracy > 0) {
    // one-loop radiative corrections included in sigma
    p = physpars.mhiggs(1);
    sigmaMh(1, 1) = pis1s1(p, q);
    sigmaMh(1, 2) = pis1s2(p, q);
    sigmaMh(2, 2) = pis2s2(p, q);

    p = physpars.mhiggs(3);
    sigmaMH(1, 1) = pis1s1(p, q);
    sigmaMH(1, 2) = pis1s2(p, q);
    sigmaMH(2, 2) = pis2s2(p, q);

    if (numHiggsMassLoops > 1) {
      double s11s = 0., s22s = 0., s12s = 0.,
        gstrong = displayGaugeCoupling(3),
        rmtsq = sqr(forLoops.mt), scalesq = sqr(displayMu()),
        vev2 = sqr(displayHvev()), tbeta = displayTanb(),
        amu = -displaySusyMu(), mg = displayGaugino(3);

      double sxt = sin(forLoops.thetat),
        cxt = cos(forLoops.thetat);
      double mst1sq = sqr(forLoops.mu(1, 3)),
        mst2sq = sqr(forLoops.mu(2, 3));
      // two-loop Higgs corrections: alpha_s alpha_b
      double sxb = sin(forLoops.thetab),
      cxb = cos(forLoops.thetab);
      double msb1sq = sqr(forLoops.md(1, 3)),
        msb2sq = sqr(forLoops.md(2, 3));
      double cotbeta = 1.0 / tbeta;
      double rmbsq = sqr(forLoops.mb);

      double s11b = 0.0, s12b = 0.0, s22b = 0.0;
      double s11tau = 0.0, s12tau = 0.0, s22tau = 0.0;
      double s11w = 0.0, s22w = 0.0, s12w = 0.0;
      int kkk = 0; // chooses DR-bar scheme from slavich et al

      double fmasq = fabs(mAsq);
      // two-loop Higgs corrections: alpha_s alpha_t, alpha_s alpha_b and
      // alpha_b^2, alpha_t*2, alpha_b alpha_t
      dszhiggs_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
                &tbeta, &vev2, &gstrong, &kkk, &s11s, &s22s, &s12s);
      dszodd_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
              &tbeta, &vev2, &gstrong, &p2s);
      dszhiggs_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
                &cotbeta, &vev2, &gstrong, &kkk, &s22b, &s11b, &s12b);
      dszodd_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
              &cotbeta, &vev2, &gstrong, &p2b);
      ddshiggs_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s11w,
              &s12w, &s22w);
      ddsodd_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2,
              &p2w);
      /*
      // In hep-ph/0406277 found the lambda_tau^2 and lambda_tau lambda_b
      // corrections to be negligible. Have commented them here for
      // calculational efficiency. Uncomment to see their effect.
      double sintau = sin(forLoops.thetatau),
        costau = cos(forLoops.thetatau);
      double rmtausq = sqr(forLoops.mtau);
      int OS = 0;
      double mstau1sq = sqr(forLoops.me(1, 3)),
        mstau2sq = sqr(forLoops.me(2, 3));
      double msnusq = sqr(forLoops.msnu(3));
      tausqhiggs_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
                  &costau, &scalesq, &amu, &tanb, &vev2, &OS, &s11tau,
                  &s22tau, &s12tau);
      tausqodd_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
                &costau, &scalesq, &amu, &tanb, &vev2, &p2tau);
      */

      sigmaMh(1, 1) = sigmaMh(1, 1) - s11s - s11w - s11b - s11tau;
      sigmaMH(1, 1) = sigmaMH(1, 1) - s11s - s11w - s11b - s11tau;
      sigmaMh(1, 2) = sigmaMh(1, 2) - s12s - s12w - s12b - s12tau;
      sigmaMH(1, 2) = sigmaMH(1, 2) - s12s - s12w - s12b - s12tau;
      sigmaMh(2, 2) = sigmaMh(2, 2) - s22s - s22w - s22b - s22tau;
      sigmaMH(2, 2) = sigmaMH(2, 2) - s22s - s22w - s22b - s22tau;
      }

    sigmaMh(2, 1) = sigmaMh(1, 2);
    sigmaMH(2, 1) = sigmaMH(1, 2);

    /*
      As the calculation stands without the two-loop terms, BPMZ have
      obviously organised PI_Sij (CP-even loop corrections) so that their pole
      masses come out correctly to one-loop order, hence there is no need to
      add any one-loop self energies of the PI^AA terms to DMA. Secondly, the
      one loop terms involving mA in the PI_Sij are not proportional to either
      alpha_s, alpha_t or alpha_b. So, within PI_Sij(1 loop), if I consider
      the 2-loop order difference due to not using mA(pole), it will be (some
      other coupling) x (either alpha_s OR alpha_t OR alpha_b) ie it will not
      be of order alpha_s alpha_t, alpha_b alpha_s or alpha_t^2, thus
      remaining consistent with your two-loop terms. They have also performed
      the calculation where their stuff "adds on" to the 1-loop BPMZ stuff.
      The tadpoles therefore appear only as 1-loop pieces, except where
      minimisation conditions are explicitly used (in the calculation of mA)
     */
    dMA = p2s + p2w + p2b + p2tau;
    mhAtmh(1, 1) = mHtree(1, 1) + physpars.t1OV1Ms1loop + dMA * sqr(sin(beta));
    mhAtmh(1, 2) = mHtree(1, 2) - dMA * sin(beta) * cos(beta);
    mhAtmh(2, 2) = mHtree(2, 2) + physpars.t2OV2Ms1loop + dMA * sqr(cos(beta));
    mhAtmh(2, 1) = mhAtmh(1 ,2);
    mhAtmh = mhAtmh - sigmaMh;

    mhAtmH(1, 1) = mHtree(1, 1) + physpars.t1OV1Ms1loop + dMA * sqr(sin(beta));
    mhAtmH(1, 2) = mHtree(1, 2) - dMA * sin(beta) * cos(beta);
    mhAtmH(2, 2) = mHtree(2, 2) + physpars.t2OV2Ms1loop + dMA * sqr(cos(beta));
    mhAtmH(2, 1) = mhAtmH(1 ,2);
    mhAtmH = mhAtmH - sigmaMH;
  }

  DoubleVector temp(2);
  double theta;

  temp = mhAtmh.sym2by2(theta);

  if (temp(1) < 0.0 || temp(2) < 0.0) {
    flagTachyon(true);
    if (PRINTOUT > 2) cout << " h0/H tachyon ";
  }
  temp = temp.apply(ccbSqrt);

  // Definitions are such that theta should diagonalise the matrix like
  // O = [ cos  sin ]
  //     [-sin  cos ]
  // and
  // O m O^T = [ m2^2      ]
  //           [      m1^2 ]
  // where m1 < m2, therefore if they come out in the wrong order,
  // add pi/2 onto theta.
  if (temp(2) > temp(1)) theta = theta + PI * 0.5;

  physpars.thetaH = theta; // theta defined for p=mh
  int i; double littleMh = temp.apply(fabs).min(i);

  temp = mhAtmH.sym2by2(theta);
  if (temp(1) < 0.0 && temp(2) < 0.0) {
    flagTachyon(true);
    if (PRINTOUT > 2) cout << " h0/H tachyon ";
  }
  temp = temp.apply(ccbSqrt);
  double bigMh = temp.max();

  double piaa = piAA(mApole, displayMu());
  //  double piaa = piAA(displayDrBarPars().mhiggs(2), displayMu());
  double poleMasq = (displayMh2Squared() - displayMh1Squared() )
    / cos(2.0 * beta) - sqr(mz);

  if (accuracy > 0) {
      poleMasq =
        (displayMh2Squared() - displayTadpole2Ms() -
         displayMh1Squared() + displayTadpole1Ms()) /
        cos(2.0 * beta) - mz2 - piaa +
        sqr(sin(beta)) * physpars.t1OV1Ms1loop + sqr(cos(beta)) *
        physpars.t2OV2Ms1loop + dMA;
    }

  if (poleMasq < 0.0) { // This definitely *is* a tachyon.
    flagTachyon(true);
    if (PRINTOUT) cout << " mA(phys)^2=" << poleMasq;
    poleMasq = fabs(poleMasq);
  }

  double pihphm = piHpHm(physpars.mhiggs(4), displayMu());

  double poleMhc = poleMasq + sqr(displayMw()) + piaa + piwwtMS - pihphm;

  physpars.mhiggs(1) = littleMh;
  physpars.mhiggs(2) = ccbSqrt(poleMasq);
  physpars.mhiggs(3) = bigMh;
  physpars.mhiggs(4) = ccbSqrt(poleMhc);
}


 void MssmSoftsusy::addCharginoLoop(DoubleMatrix & mass) {
  double g = displayGaugeCoupling(2),
    gp = displayGaugeCoupling(1) * sqrt(0.6),
    ht = displayYukawaElement(YU, 3, 3),
    htau = displayYukawaElement(YE, 3, 3),
    hb = displayYukawaElement(YD, 3, 3),
    q = displayMu(),
    tanb = displayTanb(),
    mz = displayMzRun(), mw = displayMwRun();
  double beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  double mt = ht * displayHvev() * sin(beta) / root2,
    mb = hb * displayHvev() * cos(beta) / root2,
    mtau = htau * displayHvev() * cos(beta) / root2;
  double e = g * calcSinthdrbar();

  double p = sqrt(forLoops.mchBpmz(1) * forLoops.mchBpmz(2));
  DoubleMatrix sigmaL(2, 2), sigmaR(2, 2), sigmaS(2, 2);

  // basis is (psi1+ psi2+, L R)
  DoubleMatrix aPsicDu(2, 2), aPsicUd(2, 2), aPsicENu(2, 2), aPsicNue(2, 2);
  DoubleMatrix aPsicBt(2, 2), aPsicTb(2, 2), aPsicTauNu(2, 2),
    aPsicNuTau(2, 2);
  DoubleMatrix aPsicBtm(2, 2), aPsicTbm(2, 2), aPsicTauNum(2, 2),
    aPsicNuTaum(2, 2);

  DoubleMatrix bPsicDu(2, 2), bPsicUd(2, 2), bPsicENu(2, 2), bPsicNue(2, 2);
  DoubleMatrix bPsicBt(2, 2), bPsicTb(2, 2), bPsicTauNu(2, 2),
    bPsicNuTau(2, 2);
  DoubleMatrix bPsicBtm(2, 2), bPsicTbm(2, 2), bPsicTauNum(2, 2),
    bPsicNuTaum(2, 2);

  ComplexMatrix aChicDu(2, 2), aChicUd(2, 2), aChicENu(2, 2), aChicNue(2, 2);
  ComplexMatrix aChicBt(2, 2), aChicTb(2, 2), aChicTauNu(2, 2),
    aChicNuTau(2, 2);

  ComplexMatrix bChicDu(2, 2), bChicUd(2, 2), bChicENu(2, 2), bChicNue(2, 2);
  ComplexMatrix bChicBt(2, 2), bChicTb(2, 2), bChicTauNu(2, 2),
    bChicNuTau(2, 2);

  aPsicDu(1, 1) = g;       aPsicENu(1, 1) = g;
  bPsicUd(1, 1) = g;       bPsicNue(1, 1) = g;
  aPsicBt(1, 1) = g;       aPsicTauNu(1, 1) = g;
  bPsicTb(1, 1) = g;       bPsicNuTau(1, 1) = g;
  bPsicBt(2, 1) = - hb;    bPsicTauNu(2, 1) = - htau;
  bPsicTb(2, 2) = - hb;    bPsicNuTau(2, 2) = - htau;
  aPsicTb(2, 1) = - ht;
  aPsicBt(2, 2) = - ht;

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  // mix 3rd generation sfermions
  aPsicBtm = aPsicBt * rot2d(forLoops.thetat).transpose();
  aPsicTbm = aPsicTb * rot2d(forLoops.thetab).transpose();
  aPsicNuTaum = aPsicNuTau * rot2d(forLoops.thetatau).transpose();
  bPsicBtm = bPsicBt * rot2d(forLoops.thetat).transpose();
  bPsicTbm = bPsicTb * rot2d(forLoops.thetab).transpose();
  bPsicNuTaum = bPsicNuTau * rot2d(forLoops.thetatau).transpose();

  // mix up to chargino basis
  aChicDu = v * aPsicDu;
  aChicUd = v * aPsicUd;
  aChicBt = v * aPsicBtm;
  aChicTb = v * aPsicTbm;
  aChicNue = v * aPsicNue;
  aChicNue = v * aPsicENu;
  aChicNuTau = v * aPsicNuTaum;
  bChicDu = u * bPsicDu;
  bChicUd = u * bPsicUd;
  bChicBt = u * bPsicBtm;
  bChicTb = u * bPsicTbm;
  bChicNue = u * bPsicNue;
  bChicNue = u * bPsicENu;
  bChicNuTau = u * bPsicNuTaum;

  DoubleVector msup(2), msdown(2), msel(2), mscharm(2), msstrange(2),
    msmuon(2), msnumu(2), mstop(2), msbot(2), mstau(2), msnutau(2), msnue(2);
  msup(1) = forLoops.mu(1, 1);      msup(2) = forLoops.mu(2, 1);
  mscharm(1) = forLoops.mu(1, 2);   mscharm(2) = forLoops.mu(2, 2);
  mstop(1) = forLoops.mu(1, 3);     mstop(2) = forLoops.mu(2, 3);
  msdown(1) = forLoops.md(1, 1);    msdown(2) = forLoops.md(2, 1);
  msstrange(1) = forLoops.md(1, 2); msstrange(2) = forLoops.md(2, 2);
  msbot(1) = forLoops.md(1, 3);     msbot(2) = forLoops.md(2, 3);
  msel(1) = forLoops.me(1, 1);      msel(2) = forLoops.me(2, 1);
  msmuon(1) = forLoops.me(1, 2);    msmuon(2) = forLoops.me(2, 2);
  mstau(1) = forLoops.me(1, 3);     mstau(2) = forLoops.me(2, 3);
  msnue(1) = forLoops.msnu(1);
  msnumu(1) = forLoops.msnu(2);
  msnutau(1) = forLoops.msnu(3);

  ComplexMatrix aPsi0PsicW(4, 2), bPsi0PsicW(4, 2), aChi0PsicW(4, 2),
    bChi0PsicW(4, 2);
  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;
  bPsi0PsicW(3, 2) = -g / root2;

  aChi0PsicW = n * aPsi0PsicW;
  bChi0PsicW = n.complexConjugate() * bPsi0PsicW;

  ComplexMatrix aPsiPsiZ(2, 2), bPsiPsiZ(2, 2), aPsiChiZ(2, 2), bPsiChiZ(2, 2);
  double sinthW = calcSinthdrbar();
  double thetaWDRbar = asin(sinthW);
  aPsiPsiZ(1, 1) = g * cos(thetaWDRbar);
  aPsiPsiZ(2, 2) = g * cos(2.0 * thetaWDRbar) / (2.0 * cos(thetaWDRbar));
  bPsiPsiZ = aPsiPsiZ;

  aPsiChiZ = aPsiPsiZ * v.transpose();
  bPsiChiZ = bPsiPsiZ * u.hermitianConjugate();

  ComplexMatrix aPsiChiGam(2, 2), bPsiChiGam(2, 2);
  aPsiChiGam = e * v.complexConjugate();
  bPsiChiGam = e * u;

  DoubleMatrix aPsiPsiHc1(4, 2), bPsiPsiHc1(4, 2);
  DoubleMatrix aPsiPsiHc2(4, 2), bPsiPsiHc2(4, 2);
  aPsiPsiHc1(1, 2) = gp / root2;
  bPsiPsiHc2(1, 2) = aPsiPsiHc1(1, 2);
  aPsiPsiHc1(2, 2) = g / root2;
  bPsiPsiHc2(2, 2) = g / root2;
  aPsiPsiHc1(3, 1) = -g;
  bPsiPsiHc2(4, 1) = g;

  ComplexMatrix aChiPsiHc1(4, 2), bChiPsiHc1(4, 2);
  ComplexMatrix aChiPsiHc2(4, 2), bChiPsiHc2(4, 2);
  aChiPsiHc1 = n.complexConjugate() * aPsiPsiHc1;
  bChiPsiHc1 = n * bPsiPsiHc1;
  aChiPsiHc2 = n.complexConjugate() * aPsiPsiHc2;
  bChiPsiHc2 = n * bPsiPsiHc2;

  ComplexMatrix aChiPsiHHp(4, 2), bChiPsiHHp(4, 2);
  ComplexMatrix aChiPsiHGp(4, 2), bChiPsiHGp(4, 2);
  int i,j,k; for (i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      aChiPsiHGp(i, j) = cosb * aChiPsiHc1(i, j) + sinb * aChiPsiHc2(i, j);
      bChiPsiHGp(i, j) = cosb * bChiPsiHc1(i, j) + sinb * bChiPsiHc2(i, j);
      aChiPsiHHp(i, j) =-sinb * aChiPsiHc1(i, j) + cosb * aChiPsiHc2(i, j);
      bChiPsiHHp(i, j) =-sinb * bChiPsiHc1(i, j) + cosb * bChiPsiHc2(i, j);
    }

  ComplexMatrix aPsiPsis1(2, 2), aPsiPsis2(2, 2),
    aPsiPsip1(2, 2), aPsiPsip2(2, 2);
  ComplexMatrix bPsiPsis1(2, 2), bPsiPsis2(2, 2),
    bPsiPsip1(2, 2), bPsiPsip2(2, 2);
  aPsiPsis1(1, 2) = g / root2;
  aPsiPsis2(2, 1) = g / root2;
  aPsiPsip1(1, 2) = g / root2;
  aPsiPsip2(2, 1) =-g / root2;
  bPsiPsis1 = aPsiPsis1;
  bPsiPsis2 = aPsiPsis2;
  bPsiPsip1 = -1.0 * aPsiPsip1.transpose();
  bPsiPsip2 = -1.0 * aPsiPsip2.transpose();

  ComplexMatrix aPsiChis1(2, 2), aPsiChis2(2, 2),
    aPsiChip1(2, 2), aPsiChip2(2, 2);
  ComplexMatrix bPsiChis1(2, 2), bPsiChis2(2, 2),
    bPsiChip1(2, 2), bPsiChip2(2, 2);
  aPsiChis1 = aPsiPsis1 * u.hermitianConjugate();
  aPsiChis2 = aPsiPsis2 * u.hermitianConjugate();
  aPsiChip1 = aPsiPsip1 * u.hermitianConjugate();
  aPsiChip2 = aPsiPsip2 * u.hermitianConjugate();
  bPsiChis1 = bPsiPsis1 * v.transpose();
  bPsiChis2 = bPsiPsis2 * v.transpose();
  bPsiChip1 = bPsiPsip1 * v.transpose();
  bPsiChip2 = bPsiPsip2 * v.transpose();

  ComplexMatrix aPsiChiH(2, 2), aPsiChih(2, 2), aPsiChiG(2, 2), aPsiChiA(2, 2);
  ComplexMatrix bPsiChiH(2, 2), bPsiChih(2, 2), bPsiChiG(2, 2),
    bPsiChiA(2, 2);
  double cosa = cos(forLoops.thetaH), sina = sin(forLoops.thetaH);
  aPsiChiH = cosa * aPsiChis1 + sina * aPsiChis2;
  aPsiChih =-sina * aPsiChis1 + cosa * aPsiChis2;
  aPsiChiG = cosb * aPsiChip1 + sinb * aPsiChip2;
  aPsiChiA =-sinb * aPsiChip1 + cosb * aPsiChip2;
  bPsiChiH = cosa * bPsiChis1 + sina * bPsiChis2;
  bPsiChih =-sina * bPsiChis1 + cosa * bPsiChis2;
  bPsiChiG = cosb * bPsiChip1 + sinb * bPsiChip2;
  bPsiChiA =-sinb * bPsiChip1 + cosb * bPsiChip2;

  // actual contributions start here
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++)
      for (k=1; k<=2; k++) {
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (3.0 * aPsicDu(i, k) * aPsicDu(j, k) * b1(p, 0., msup(k), q) +
           3.0 * aPsicUd(i, k) * aPsicUd(j, k) * b1(p, 0., msdown(k), q) +
           aPsicNue(i, k) * aPsicNue(j, k) * b1(p, 0., msel(k), q) +
           aPsicENu(i, k) * aPsicENu(j, k) * b1(p, 0., msnue(k), q));
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (3.0 * aPsicDu(i, k) * aPsicDu(j, k) * b1(p, 0., mscharm(k), q) +
           3.0 * aPsicUd(i, k) * aPsicUd(j, k) * b1(p, 0., msstrange(k), q) +
           aPsicNue(i, k) * aPsicNue(j, k) * b1(p, 0., msmuon(k), q) +
           aPsicENu(i, k) * aPsicENu(j, k) * b1(p, 0., msnumu(k), q));
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (3.0 * aPsicBt(i, k) * aPsicBt(j, k) * b1(p, mt, mstop(k), q) +
           3.0 * aPsicTb(i, k) * aPsicTb(j, k) * b1(p, mb, msbot(k), q) +
           aPsicNuTau(i, k) * aPsicNuTau(j, k) * b1(p, mtau, mstau(k), q) +
           aPsicTauNu(i, k) * aPsicTauNu(j, k) * b1(p, 0., msnutau(k), q));
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (3.0 * bPsicDu(i, k) * bPsicDu(j, k) * b1(p, 0., msup(k), q) +
           3.0 * bPsicUd(i, k) * bPsicUd(j, k) * b1(p, 0., msdown(k), q) +
           bPsicNue(i, k) * bPsicNue(j, k) * b1(p, 0., msel(k), q) +
           bPsicENu(i, k) * bPsicENu(j, k) * b1(p, 0., msnue(k), q));
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (3.0 * bPsicDu(i, k) * bPsicDu(j, k) * b1(p, 0., mscharm(k), q) +
           3.0 * bPsicUd(i, k) * bPsicUd(j, k) * b1(p, 0., msstrange(k), q) +
           bPsicNue(i, k) * bPsicNue(j, k) * b1(p, 0., msmuon(k), q) +
           bPsicENu(i, k) * bPsicENu(j, k) * b1(p, 0., msnumu(k), q));
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (3.0 * bPsicBt(i, k) * bPsicBt(j, k) * b1(p, mt, mstop(k), q) +
           3.0 * bPsicTb(i, k) * bPsicTb(j, k) * b1(p, mb, msbot(k), q) +
           bPsicNuTau(i, k) * bPsicNuTau(j, k) * b1(p, mtau, mstau(k), q) +
           bPsicTauNu(i, k) * bPsicTauNu(j, k) * b1(p, 0., msnutau(k), q));
        sigmaS(i, j) = sigmaS(i, j) +
          (3.0 * bPsicBt(i, k) * aPsicBt(j, k) * mt * b0(p, mt, mstop(k), q) +
           3.0 * bPsicTb(i, k) * aPsicTb(j, k) * mb * b0(p, mb, msbot(k), q) +
           bPsicNuTau(i, k) * aPsicNuTau(j, k) * mtau *
           b0(p, mtau, mstau(k), q));
      }

  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++)
      for (k=1; k<=4; k++) {
        double b1p = b1(p, mneut(k), mw, q);
        double b0p = b0(p, mneut(k), mw, q);
        sigmaL(i, j) = sigmaL(i, j) +
          (aChi0PsicW(k, i).conj() * aChi0PsicW(k, j) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) +
          (bChi0PsicW(k, i).conj() * bChi0PsicW(k, j) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) - 4.0 * mneut(k) *
          (bChi0PsicW(k, i).conj() * aChi0PsicW(k, j) * b0p).real();

        // G+
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (aChiPsiHGp(k, i).conj() * aChiPsiHGp(k, j) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (bChiPsiHGp(k, i).conj() * bChiPsiHGp(k, j) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + mneut(k) *
          (bChiPsiHGp(k, i).conj() * aChiPsiHGp(k, j) * b0p).real();

        // H+
        b1p = b1(p, mneut(k), forLoops.mhiggs(4), q);
        b0p = b0(p, mneut(k), forLoops.mhiggs(4), q);
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (aChiPsiHHp(k, i).conj() * aChiPsiHHp(k, j) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (bChiPsiHHp(k, i).conj() * bChiPsiHHp(k, j) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + mneut(k) *
          (bChiPsiHHp(k, i).conj() * aChiPsiHHp(k, j) * b0p).real();

        if (k <= 2) {
          b1p = b1(p, mch(k), mz, q);
          b0p = b0(p, mch(k), mz, q);
          sigmaL(i, j) = sigmaL(i, j) +
            (aPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b1p).real();
          sigmaR(i, j) = sigmaR(i, j) +
            (bPsiChiZ(i, k).conj() * bPsiChiZ(j, k) * b1p).real();
          sigmaS(i, j) = sigmaS(i, j) - 4.0 * mch(k) *
            (bPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b0p).real();
          b1p = b1(p, mch(k), 0., q);
          b0p = b0(p, mch(k), 0., q);
          sigmaL(i, j) = sigmaL(i, j) +
            (aPsiChiGam(i, k).conj() * aPsiChiGam(j, k) * b1p).real();
          sigmaR(i, j) = sigmaR(i, j) +
            (bPsiChiGam(i, k).conj() * bPsiChiGam(j, k) * b1p).real();
          sigmaS(i, j) = sigmaS(i, j) - 4.0 * mch(k) *
            (bPsiChiGam(i, k).conj() * aPsiChiGam(j, k) * b0p).real();

          // H
          b1p = b1(p, mneut(k), forLoops.mhiggs(3), q);
          b0p = b0(p, mneut(k), forLoops.mhiggs(3), q);
          sigmaL(i, j) = sigmaL(i, j) + 0.5 *
            (aPsiChiH(i, k).conj() * aPsiChiH(j, k) * b1p).real();
          sigmaR(i, j) = sigmaR(i, j) + 0.5 *
            (bPsiChiH(i, k).conj() * bPsiChiH(j, k) * b1p).real();
          sigmaS(i, j) = sigmaS(i, j) + mch(k) *
            (bPsiChiH(i, k).conj() * aPsiChiH(j, k) * b0p).real();

          // h
          b1p = b1(p, mneut(k), forLoops.mhiggs(1), q);
          b0p = b0(p, mneut(k), forLoops.mhiggs(1), q);
          sigmaL(i, j) = sigmaL(i, j) + 0.5 *
            (aPsiChih(i, k).conj() * aPsiChih(j, k) * b1p).real();
          sigmaR(i, j) = sigmaR(i, j) + 0.5 *
            (bPsiChih(i, k).conj() * bPsiChih(j, k) * b1p).real();
          sigmaS(i, j) = sigmaS(i, j) + mch(k) *
            (bPsiChih(i, k).conj() * aPsiChih(j, k) * b0p).real();

          // G0
          b1p = b1(p, mneut(k), mz, q);
          b0p = b0(p, mneut(k), mz, q);
          sigmaL(i, j) = sigmaL(i, j) + 0.5 *
            (aPsiChiG(i, k).conj() * aPsiChiG(j, k) * b1p).real();
          sigmaR(i, j) = sigmaR(i, j) + 0.5 *
            (bPsiChiG(i, k).conj() * bPsiChiG(j, k) * b1p).real();
          sigmaS(i, j) = sigmaS(i, j) + mch(k) *
            (bPsiChiG(i, k).conj() * aPsiChiG(j, k) * b0p).real();

          // A0
          b1p = b1(p, mneut(k), forLoops.mhiggs(2), q);
          b0p = b0(p, mneut(k), forLoops.mhiggs(2), q);
          sigmaL(i, j) = sigmaL(i, j) + 0.5 *
            (aPsiChiA(i, k).conj() * aPsiChiA(j, k) * b1p).real();
          sigmaR(i, j) = sigmaR(i, j) + 0.5 *
            (bPsiChiA(i, k).conj() * bPsiChiA(j, k) * b1p).real();
          sigmaS(i, j) = sigmaS(i, j) + mch(k) *
            (bPsiChiA(i, k).conj() * aPsiChiA(j, k) * b0p).real();
        }
      }
  mass = mass - 1.0 / (16.0 * sqr(PI)) *
    (sigmaR * mass + mass * sigmaL + sigmaS);
}

/*
void MssmSoftsusy::addCharginoLoop(DoubleMatrix & mass) {
  double g2 = displayGaugeCoupling(2),
    ht = displayYukawaElement(YU, 3, 3),
    hb = displayYukawaElement(YD, 3, 3),
    q = displayMu(), mu = displaySusyMu(), M2 = displayGaugino(2),
    MA = forLoops.mhiggs(2), tanb = displayTanb(), mz = displayMzRun();

  // From hep-ph/9606211, approx formulae only, both checked! Results are a
  // bit inaccurate.... Note that a0,b0,b1 do NOT depend on the sign of the
  // masses/p input, so don't have to change sign of mu in these terms.
  // quark/squark loops
  double quarkSquarkM2 = - sqr(g2) / (16.0 * sqr(PI)) *
    (9.0 * b1(M2, 0.0, ccbSqrt(displaySoftMassSquared(mQl, 1, 1)), q) +
     3.0 * b1(M2, 0.0, ccbSqrt(displaySoftMassSquared(mLl, 1, 1)), q) -
     mu / M2 * sin(2.0 * atan(tanb)) *
     (b0(M2, mu, MA, q) - b0(M2, mu, mz, q)) +
     b1(M2, mu, MA, q) + b1(M2, mu, mz, q));

  // Gauge loops
  double gaugeBosonM2 =
    sqr(g2) / (4.0 * sqr(PI)) *
    (2.0 * b0(M2, M2, displayMwRun(), q) - b1(M2, M2, displayMwRun(), q));

  double deltaM2oM2 = quarkSquarkM2 + gaugeBosonM2;

  double deltaMuoMu = - 3.0 / (32.0 * sqr(PI)) *
    ((sqr(hb) + sqr(ht)) *
     b1(mu, 0.0, ccbSqrt(displaySoftMassSquared(mQl, 3, 3)), q) +
     sqr(ht) *
     b1(mu, 0.0, ccbSqrt(displaySoftMassSquared(mUr, 3, 3)), q) + sqr(hb) *
     b1(mu, 0.0, ccbSqrt(displaySoftMassSquared(mDr, 3, 3)), q))
    - 3.0 * sqr(g2) / (64.0 * sqr(PI)) *
    (b1(mu, M2, MA, q) + b1(mu, M2, mz, q) + 2.0 * b1(mu, mu, mz, q) - 4.0 *
     b0(mu, mu, mz, q));

  if (testNan(deltaM2oM2 * deltaMuoMu)) {
    if (PRINTOUT > 2) cout << "Chargino loop nan\n";
    flagTachyon(true);
    return;
  }

  mass(1, 1) = mass(1, 1) * (1.0 + deltaM2oM2);
  mass(2, 2) = mass(2, 2) * (1.0 + deltaMuoMu);
}
*/

// checked
void MssmSoftsusy::charginos(int accuracy, double piwwtMS) {
  double tanb = displayTanb(), smu = displaySusyMu();
  DoubleMatrix mCh(2, 2);
    double m2 = displayGaugino(2);

  double mwOneLarg = sqr(displayMw()) + piwwtMS;

  // tree level mass matrix
  mCh(1, 1) = m2;
  mCh(2, 1) = root2 * sqrt(fabs(mwOneLarg)) * cos(atan(tanb));
  mCh(1, 2) = mCh(2, 1) * tanb;
  mCh(2, 2) = smu;

  if (accuracy > 0) addCharginoLoop(mCh);

  physpars.mch = mCh.asy2by2(physpars.thetaL, physpars.thetaR);
}

double MssmSoftsusy::thet(double a, double b, double c = 0.0) {
  double yy = maximum(sqr(a), sqr(b));
  yy = maximum(yy, sqr(c));
  return log(yy / sqr(displayMu()));
}

/*void MssmSoftsusy::addNeutralinoLoop(DoubleMatrix & mass) {
  double gp = displayGaugeCoupling(1) * sqrt(0.6),
    g2 = displayGaugeCoupling(2), ht = displayYukawaElement(YU, 3, 3),
    hb = displayYukawaElement(YD, 3, 3),
    q = displayMu(), mu = -displaySusyMu(), M2 = displayGaugino(2), M1 =
    displayGaugino(1), MA = forLoops.mhiggs(2), tanb = displayTanb(),
    mz = displayMzRun();
  // brevity

  double msq1 = ccbSqrt(displaySoftMassSquared(mQl, 1, 1)),
    mslep1 = ccbSqrt(displaySoftMassSquared(mLl, 1, 1));

  // From hep-ph/9606211, approx formulae only, both checked!
  double deltaM1 = -sqr(gp) / (16.0 * sqr(PI)) *
    (11.0 * b1(M1, 0.0, msq1, q) + 9.0 * b1(M1, 0.0, mslep1, q) + mu / M1 * sin
     (2.0 * atan(tanb)) * (b0(M1, mu, MA, q) - b0(M1, mu, mz, q))
     + b1(M1, mu, MA, q) + b1(M1, mu, mz, q));

  double deltaM2 = - sqr(g2) / (16.0 * sqr(PI)) *
    (9.0 * b1(M2, 0.0, ccbSqrt(displaySoftMassSquared(mQl, 1, 1)), q) +
     3.0 * b1(M2, 0.0, ccbSqrt(displaySoftMassSquared(mLl, 1, 1)), q) +
     mu / M2 * sin(2.0 * atan(tanb)) *
     (b0(M2, mu, MA, q) - b0(M2, mu, mz, q)) +
     b1(M2, mu, MA, q) + b1(M2, mu, mz, q)
     + 4.0 * b1(M2, M2, displayMwRun(), q)
     - 8.0 * b0(M2, M2, displayMwRun(), q));

  double deltaMu = - 3.0 / (32.0 * sqr(PI)) *
    ((sqr(hb) + sqr(ht)) *
     b1(mu, 0.0, ccbSqrt(displaySoftMassSquared(mQl, 3, 3)), q) +
     sqr(ht) *
     b1(mu, 0.0, ccbSqrt(displaySoftMassSquared(mUr, 3, 3)), q) + sqr(hb) *
     b1(mu, 0.0, ccbSqrt(displaySoftMassSquared(mDr, 3, 3)), q))
    - 3.0 * sqr(g2) / (64.0 * sqr(PI)) *
    (b1(mu, M2, MA, q) + b1(mu, M2, mz, q) + 2.0 * b1(mu, mu, mz, q) - 4.0 *
     b0(mu, mu, mz, q));

  if (testNan(deltaM1 * deltaM2 * deltaMu)) {
    if (PRINTOUT > 2) cout << "Neutralino loop tachyon";
    flagTachyon(true);
    return;
  }

  mass(1, 1) = mass(1, 1) * (1.0 + deltaM1);
  mass(2, 2) = mass(2, 2) * (1.0 + deltaM2);
  mass(3, 4) = mass(3, 4) * (1.0 + deltaMu);
  }
*/

void MssmSoftsusy::addNeutralinoLoop(DoubleMatrix & mass) {
  double g = displayGaugeCoupling(2),
    gp = displayGaugeCoupling(1) * sqrt(0.6),
    ht = displayYukawaElement(YU, 3, 3),
    htau = displayYukawaElement(YE, 3, 3),
    hb = displayYukawaElement(YD, 3, 3),
    q = displayMu(),
    tanb = displayTanb(),
    mz = displayMzRun(), mw = displayMwRun();
  double beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  double mt = ht * displayHvev() * sin(beta) / root2,
    mb = hb * displayHvev() * cos(beta) / root2,
    mtau = htau * displayHvev() * cos(beta) / root2;
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;

  double p = sqrt(forLoops.mneut(1) * forLoops.mneut(4));
  DoubleMatrix sigmaL(4, 4), sigmaR(4, 4), sigmaS(4, 4);

  DoubleVector msup(2), msdown(2), msel(2), mscharm(2), msstrange(2),
    msmuon(2), msnumu(2), mstop(2), msbot(2), mstau(2), msnutau(2), msnue(2);
  msup(1) = forLoops.mu(1, 1);      msup(2) = forLoops.mu(2, 1);
  mscharm(1) = forLoops.mu(1, 2);   mscharm(2) = forLoops.mu(2, 2);
  mstop(1) = forLoops.mu(1, 3);     mstop(2) = forLoops.mu(2, 3);
  msdown(1) = forLoops.md(1, 1);    msdown(2) = forLoops.md(2, 1);
  msstrange(1) = forLoops.md(1, 2); msstrange(2) = forLoops.md(2, 2);
  msbot(1) = forLoops.md(1, 3);     msbot(2) = forLoops.md(2, 3);
  msel(1) = forLoops.me(1, 1);      msel(2) = forLoops.me(2, 1);
  msmuon(1) = forLoops.me(1, 2);    msmuon(2) = forLoops.me(2, 2);
  mstau(1) = forLoops.me(1, 3);     mstau(2) = forLoops.me(2, 3);
  msnue(1) = forLoops.msnu(1);
  msnumu(1) = forLoops.msnu(2);
  msnutau(1) = forLoops.msnu(3);

  DoubleMatrix aPsiUu(4, 2), aPsiDd(4, 2), aPsiEe(4, 2), aPsiNuNu(4, 2);
  DoubleMatrix bPsiUu(4, 2), bPsiDd(4, 2), bPsiEe(4, 2), bPsiNuNu(4, 2);
  DoubleMatrix aPsiTt(4, 2), aPsiBb(4, 2), aPsiTauTau(4, 2), aPsiNutNut(4, 2);
  DoubleMatrix bPsiTt(4, 2), bPsiBb(4, 2), bPsiTauTau(4, 2), bPsiNutNut(4, 2);
  aPsiUu(1, 2) = gp / root2 * yuR;
  bPsiUu(1, 1) = gp / root2 * yuL;
  bPsiUu(2, 1) = g * root2 * 0.5;
  aPsiTt(1, 2) = gp / root2 * yuR;
  bPsiTt(1, 1) = gp / root2 * yuL;
  bPsiTt(2, 1) = g * root2 * 0.5;
  aPsiTt(4, 1) = ht;
  bPsiTt(4, 2) = ht;

  aPsiDd(1, 2) = gp / root2 * ydR;
  bPsiDd(1, 1) = gp / root2 * ydL;
  bPsiDd(2, 1) = -g * root2 * 0.5;
  aPsiBb(1, 2) = gp / root2 * ydR;
  bPsiBb(1, 1) = gp / root2 * ydL;
  bPsiBb(2, 1) = -g * root2 * 0.5;
  aPsiBb(3, 1) = hb;
  bPsiBb(3, 2) = hb;

  bPsiNuNu(1, 1) = gp / root2 * ynuL;
  bPsiNuNu(2, 1) = g * root2 * 0.5;
  bPsiNutNut(1, 1) = gp / root2 * ynuL;
  bPsiNutNut(2, 1) = g * root2 * 0.5;

  aPsiEe(1, 2) = gp / root2 * yeR;
  bPsiEe(1, 1) = gp / root2 * yeL;
  bPsiEe(2, 1) = -g * root2 * 0.5;
  aPsiTauTau(1, 2) = gp / root2 * yeR;
  bPsiTauTau(1, 1) = gp / root2 * yeL;
  bPsiTauTau(2, 1) = -g * root2 * 0.5;
  aPsiTauTau(3, 1) = htau;
  bPsiTauTau(3, 2) = htau;

  /// mix up the third family sfermions
  DoubleMatrix O(2, 2);
  O = rot2d(thetat);
  DoubleVector t1(2), t2(2), tt(2);
  int i; for (i=1; i<=4; i++) {
    tt(1) = aPsiTt(i, 1); tt(2) = aPsiTt(i, 2);
    t1 = O * tt;
    tt(1) = bPsiTt(i, 1); tt(2) = bPsiTt(i, 2);
    t2 = O * tt;
    aPsiTt(i, 1) = t1(1);     aPsiTt(i, 2) = t1(2);
    bPsiTt(i, 1) = t2(1);     bPsiTt(i, 2) = t2(2);
  }

  O = rot2d(thetab);
  for (i=1; i<=4; i++) {
    tt(1) = aPsiBb(i, 1); tt(2) = aPsiBb(i, 2);
    t1 = O * tt;
    tt(1) = bPsiBb(i, 1); tt(2) = bPsiBb(i, 2);
    t2 = O * tt;
    aPsiBb(i, 1) = t1(1);     aPsiBb(i, 2) = t1(2);
    bPsiBb(i, 1) = t2(1);     bPsiBb(i, 2) = t2(2);
  }

  O = rot2d(thetatau);
  for (i=1; i<=4; i++) {
    tt(1) = aPsiTauTau(i, 1); tt(2) = aPsiTauTau(i, 2);
    t1 = O * tt;
    tt(1) = bPsiTauTau(i, 1); tt(2) = bPsiTauTau(i, 2);
    t2 = O * tt;
    aPsiTauTau(i, 1) = t1(1);     aPsiTauTau(i, 2) = t1(2);
    bPsiTauTau(i, 1) = t2(1);     bPsiTauTau(i, 2) = t2(2);
  }

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  ComplexMatrix aPsi0PsicW(4, 2), bPsi0PsicW(4, 2), aPsi0ChicW(4, 2),
    bPsi0ChicW(4, 2);
  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;
  bPsi0PsicW(3, 2) = -g / root2;

  aPsi0ChicW = aPsi0PsicW * v.transpose();
  bPsi0ChicW = bPsi0PsicW * u.hermitianConjugate();

  ComplexMatrix aPsiPsiZ(4, 4), bPsiPsiZ(4, 4), aPsiChiZ(4, 4), bPsiChiZ(4, 4);
  double sinthW = calcSinthdrbar();
  double thetaWDRbar = asin(sinthW);
  double costh = cos(thetaWDRbar);
  aPsiPsiZ(3, 3) = g / (2.0 * costh);
  aPsiPsiZ(4, 4) =-g / (2.0 * costh);
  bPsiPsiZ = -1.0 * aPsiPsiZ;
  aPsiChiZ = aPsiPsiZ * n.transpose();
  bPsiChiZ = aPsiChiZ * n.hermitianConjugate();

  DoubleMatrix aPsiPsiHc1(4, 2), bPsiPsiHc1(4, 2);
  DoubleMatrix aPsiPsiHc2(4, 2), bPsiPsiHc2(4, 2);
  aPsiPsiHc1(1, 2) = gp / root2;
  bPsiPsiHc2(1, 2) = aPsiPsiHc1(1, 2);
  aPsiPsiHc1(2, 2) = g / root2;
  bPsiPsiHc2(2, 2) = g / root2;
  aPsiPsiHc1(3, 1) = -g;
  bPsiPsiHc2(4, 1) = g;

  ComplexMatrix aPsiChiHc1(4, 2), bPsiChiHc1(4, 2);
  ComplexMatrix aPsiChiHc2(4, 2), bPsiChiHc2(4, 2);
  aPsiChiHc1 = aPsiPsiHc1 * u.hermitianConjugate();
  aPsiChiHc2 = aPsiPsiHc2 * u.hermitianConjugate();
  bPsiChiHc1 = bPsiPsiHc1 * v.transpose();
  bPsiChiHc2 = bPsiPsiHc2 * v.transpose();

  ComplexMatrix aPsiChiHHp(4, 2), bPsiChiHHp(4, 2);
  ComplexMatrix aPsiChiHGp(4, 2), bPsiChiHGp(4, 2);
  int j,k; for (i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      aPsiChiHGp(i, j) = cosb * aPsiChiHc1(i, j) + sinb * aPsiChiHc2(i, j);
      bPsiChiHGp(i, j) = cosb * bPsiChiHc1(i, j) + sinb * bPsiChiHc2(i, j);
      aPsiChiHHp(i, j) =-sinb * aPsiChiHc1(i, j) + cosb * aPsiChiHc2(i, j);
      bPsiChiHHp(i, j) =-sinb * bPsiChiHc1(i, j) + cosb * bPsiChiHc2(i, j);
    }

  ComplexMatrix aPsiPsis1(4, 4), aPsiPsis2(4, 4),
    aPsiPsip1(4, 4), aPsiPsip2(4, 4);
  ComplexMatrix bPsiPsis1(4, 4), bPsiPsis2(4, 4),
    bPsiPsip1(4, 4), bPsiPsip2(4, 4);
  ComplexMatrix aPsiChis1(4, 4), aPsiChis2(4, 4),
    aPsiChip1(4, 4), aPsiChip2(4, 4);
  ComplexMatrix bPsiChis1(4, 4), bPsiChis2(4, 4),
    bPsiChip1(4, 4), bPsiChip2(4, 4);
  aPsiPsis1(1, 3) = -gp * 0.5;
  aPsiPsis1(2, 3) = g * 0.5;
  aPsiPsis2(2, 4) = -g * 0.5;
  aPsiPsis2(1, 4) = gp * 0.5;
  aPsiPsip1(1, 3) = -gp * 0.5;
  aPsiPsip1(2, 3) = g * 0.5;
  aPsiPsip2(2, 4) = g * 0.5;
  aPsiPsip2(1, 4) = -gp * 0.5;
  aPsiPsis1.symmetrise();
  aPsiPsis2.symmetrise();
  aPsiPsip1.symmetrise();
  aPsiPsip2.symmetrise();
  bPsiPsis1 = aPsiPsis1;
  bPsiPsis2 = aPsiPsis2;
  bPsiPsip1 =-1.0 * aPsiPsip1;
  bPsiPsip2 =-1.0 * aPsiPsip2;

  aPsiChis1 = aPsiPsis1 * n.hermitianConjugate();
  aPsiChis2 = aPsiPsis2 * n.hermitianConjugate();
  aPsiChip1 = aPsiPsip1 * n.hermitianConjugate();
  aPsiChip2 = aPsiPsip2 * n.hermitianConjugate();
  bPsiChis1 = bPsiPsis1 * n.transpose();
  bPsiChis2 = bPsiPsis2 * n.transpose();
  bPsiChip1 = bPsiPsip1 * n.transpose();
  bPsiChip2 = bPsiPsip2 * n.transpose();

  ComplexMatrix aPsiChiH(4, 4), aPsiChih(4, 4), aPsiChiG(4, 4), aPsiChiA(4, 4);
  ComplexMatrix bPsiChiH(4, 4), bPsiChih(4, 4), bPsiChiG(4, 4),
    bPsiChiA(4, 4);
  double cosa = cos(forLoops.thetaH), sina = sin(forLoops.thetaH);
  aPsiChiH = cosa * aPsiChis1 + sina * aPsiChis2;
  aPsiChih =-sina * aPsiChis1 + cosa * aPsiChis2;
  aPsiChiG = cosb * aPsiChip1 + sinb * aPsiChip2;
  aPsiChiA =-sinb * aPsiChip1 + cosb * aPsiChip2;
  bPsiChiH = cosa * bPsiChis1 + sina * bPsiChis2;
  bPsiChih =-sina * bPsiChis1 + cosa * bPsiChis2;
  bPsiChiG = cosb * bPsiChip1 + sinb * bPsiChip2;
  bPsiChiA =-sinb * bPsiChip1 + cosb * bPsiChip2;

  // corrections here
  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++)
      for (k=1; k<=2; k++) {
        sigmaL(i, j) = sigmaL(i, j) +
          (3.0 * aPsiUu(i, k) * aPsiUu(j, k) * b1(p, 0., msup(k), q) +
           3.0 * aPsiDd(i, k) * aPsiDd(j, k) * b1(p, 0., msdown(k), q) +
           aPsiEe(i, k) * aPsiEe(j, k) * b1(p, 0., msel(k), q) +
           aPsiNuNu(i, k) * aPsiNuNu(j, k) * b1(p, 0., msnue(k), q));
        sigmaL(i, j) = sigmaL(i, j) +
          (3.0 * aPsiUu(i, k) * aPsiUu(j, k) * b1(p, 0., mscharm(k), q) +
           3.0 * aPsiDd(i, k) * aPsiDd(j, k) * b1(p, 0., msstrange(k), q) +
           aPsiEe(i, k) * aPsiEe(j, k) * b1(p, 0., msmuon(k), q) +
           aPsiNuNu(i, k) * aPsiNuNu(j, k) * b1(p, 0., msnumu(k), q));
        sigmaL(i, j) = sigmaL(i, j) +
          (3.0 * aPsiTt(i, k) * aPsiTt(j, k) * b1(p, mt, mstop(k), q) +
           3.0 * aPsiBb(i, k) * aPsiBb(j, k) * b1(p, mb, msbot(k), q) +
           aPsiTauTau(i, k) * aPsiTauTau(j, k) * b1(p, mtau, mstau(k), q) +
           aPsiNutNut(i, k) * aPsiNutNut(j, k) * b1(p, 0., msnutau(k), q));
        sigmaR(i, j) = sigmaR(i, j) +
          (3.0 * bPsiUu(i, k) * bPsiUu(j, k) * b1(p, 0., msup(k), q) +
           3.0 * bPsiDd(i, k) * bPsiDd(j, k) * b1(p, 0., msdown(k), q) +
           bPsiEe(i, k) * bPsiEe(j, k) * b1(p, 0., msel(k), q) +
           bPsiNuNu(i, k) * bPsiNuNu(j, k) * b1(p, 0., msnue(k), q));
        sigmaR(i, j) = sigmaR(i, j) +
          (3.0 * bPsiUu(i, k) * bPsiUu(j, k) * b1(p, 0., mscharm(k), q) +
           3.0 * bPsiDd(i, k) * bPsiDd(j, k) * b1(p, 0., msstrange(k), q) +
           bPsiEe(i, k) * bPsiEe(j, k) * b1(p, 0., msmuon(k), q) +
           bPsiNuNu(i, k) * bPsiNuNu(j, k) * b1(p, 0., msnumu(k), q));
        sigmaR(i, j) = sigmaR(i, j) +
          (3.0 * bPsiTt(i, k) * bPsiTt(j, k) * b1(p, mt, mstop(k), q) +
           3.0 * bPsiBb(i, k) * bPsiBb(j, k) * b1(p, mb, msbot(k), q) +
           bPsiTauTau(i, k) * bPsiTauTau(j, k) * b1(p, mtau, mstau(k), q) +
           bPsiNutNut(i, k) * bPsiNutNut(j, k) * b1(p, 0., msnutau(k), q));
        sigmaS(i, j) = sigmaS(i, j) + 2.0 *
          (3.0 * bPsiTt(i, k) * aPsiTt(j, k) * mt * b0(p, mt, mstop(k), q) +
           3.0 * bPsiBb(i, k) * aPsiBb(j, k) * mb * b0(p, mb, msbot(k), q) +
           bPsiTauTau(i, k) * aPsiTauTau(j, k) * mtau *
           b0(p, mtau, mstau(k), q));
      }

  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++)
      for (k=1; k<=4; k++) {
        double b1p=0.;
        double b0p=0.;
        if (k<=2) {
          b1p = b1(p, mch(k), mw, q);
          b0p = b0(p, mch(k), mw, q);
        sigmaL(i, j) = sigmaL(i, j) + 2.0 *
          (aPsi0ChicW(i, k).conj() * aPsi0ChicW(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) + 2.0 *
          (bPsi0ChicW(i, k).conj() * bPsi0ChicW(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) - 8.0 * mch(k) *
          (bPsi0ChicW(i, k).conj() * aPsi0ChicW(j, k) * b0p).real();

        // G+
        sigmaL(i, j) = sigmaL(i, j) +
          (aPsiChiHGp(i, k).conj() * aPsiChiHGp(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) +
          (bPsiChiHGp(i, k).conj() * bPsiChiHGp(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + 2.0 * mch(k) *
          (bPsiChiHGp(i, k).conj() * aPsiChiHGp(j, k) * b0p).real();

        // H+
        b1p = b1(p, mch(k), forLoops.mhiggs(4), q);
        b0p = b0(p, mch(k), forLoops.mhiggs(4), q);
        sigmaL(i, j) = sigmaL(i, j) +
          (aPsiChiHHp(i, k).conj() * aPsiChiHHp(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) +
          (bPsiChiHHp(i, k).conj() * bPsiChiHHp(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + 2.0 * mch(k) *
          (bPsiChiHHp(i, k).conj() * aPsiChiHHp(j, k) * b0p).real();
        }

        b1p = b1(p, mneut(k), mz, q);
        b0p = b0(p, mneut(k), mz, q);
        sigmaL(i, j) = sigmaL(i, j) +
          (aPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) +
          (bPsiChiZ(i, k).conj() * bPsiChiZ(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) - 4.0 * mneut(k) *
          (bPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b0p).real();

        // H
        b1p = b1(p, mneut(k), forLoops.mhiggs(3), q);
        b0p = b0(p, mneut(k), forLoops.mhiggs(3), q);
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (aPsiChiH(i, k).conj() * aPsiChiH(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (bPsiChiH(i, k).conj() * bPsiChiH(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + mneut(k) *
          (bPsiChiH(i, k).conj() * aPsiChiH(j, k) * b0p).real();

        // h
        b1p = b1(p, mneut(k), forLoops.mhiggs(1), q);
        b0p = b0(p, mneut(k), forLoops.mhiggs(1), q);
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (aPsiChih(i, k).conj() * aPsiChih(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (bPsiChih(i, k).conj() * bPsiChih(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + mneut(k) *
          (bPsiChih(i, k).conj() * aPsiChih(j, k) * b0p).real();

        // G0
        b1p = b1(p, mneut(k), mz, q);
        b0p = b0(p, mneut(k), mz, q);
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (aPsiChiG(i, k).conj() * aPsiChiG(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (bPsiChiG(i, k).conj() * bPsiChiG(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + mneut(k) *
          (bPsiChiG(i, k).conj() * aPsiChiG(j, k) * b0p).real();

        // A0
        b1p = b1(p, mneut(k), forLoops.mhiggs(2), q);
        b0p = b0(p, mneut(k), forLoops.mhiggs(2), q);
        sigmaL(i, j) = sigmaL(i, j) + 0.5 *
          (aPsiChiA(i, k).conj() * aPsiChiA(j, k) * b1p).real();
        sigmaR(i, j) = sigmaR(i, j) + 0.5 *
          (bPsiChiA(i, k).conj() * bPsiChiA(j, k) * b1p).real();
        sigmaS(i, j) = sigmaS(i, j) + mneut(k) *
          (bPsiChiA(i, k).conj() * aPsiChiA(j, k) * b0p).real();
      }

  DoubleMatrix deltaM(4, 4);
  deltaM = -sigmaR * mass - mass * sigmaL - sigmaS;

  deltaM = (deltaM + deltaM.transpose()) / (32.0 * sqr(PI));

  mass = mass + deltaM;
}



// mixNeut set to diagonal = mixNeut^T mNeutralino mixNeut: checked
void MssmSoftsusy::neutralinos(int accuracy, double piwwtMS, double pizztMS) {
  double tanb = displayTanb();
  double cosb = cos(atan(tanb));
  DoubleMatrix mNeut(4, 4);

  double m1 = displayGaugino(1);
  double m2 = displayGaugino(2);
  double smu = displaySusyMu();

  // tree level
  mNeut(1, 1) = m1;
  mNeut(2, 2) = m2;
  mNeut(1, 3) = - displayMzRun() * cosb * calcSinthdrbar();
  mNeut(1, 4) = - mNeut(1, 3) * tanb;
  mNeut(2, 3) = displayMwRun() * cosb;

  mNeut(2, 4) = - mNeut(2, 3) * tanb;
  mNeut(3, 4) = - smu;

  if (accuracy > 0) addNeutralinoLoop(mNeut);

  mNeut.symmetrise();
  if (mNeut.diagonaliseSym(physpars.mixNeut, physpars.mneut) > TOLERANCE *
      1.0e-3) {
    ostringstream ii;
    ii << "accuracy bad in neutralino diagonalisation"<< flush;
    ii << "diagonalising " << physpars.mneut << " with "
       << physpars.mixNeut << "giving " << physpars.mixNeut.transpose() *
      physpars.mneut * physpars.mixNeut << flush;
    throw ii.str();
    }
}

// One loop corrections to gluino pole mass: hep-ph/9606211
// BUG fixed to use g3 at current scale 8.1.2001
// Changed to resummed version 11.05.2001
void MssmSoftsusy::gluino(int accuracy) {

  if (accuracy == 0) {
    physpars.mGluino = displayGaugino(3);
    return;
  }

  if (forLoops.mGluino == 0.0) forLoops.mGluino = displayGaugino(3);

  double p = fabs(displayGaugino(3)), Q = displayMu();

  // SUSY QCD radiative correction (gluon/gluino)
  double delta = 15.0 + 9.0 * log(sqr(Q / p));

  // Quark/squark correction
  delta = delta -
    (b1(p, dataSet.displayMass(mUp), forLoops.mu(1, 1), Q) +
     b1(p, dataSet.displayMass(mCharm), forLoops.mu(1, 2), Q) +
     b1(p, forLoops.mt, forLoops.mu(1, 3), Q) +
     b1(p, dataSet.displayMass(mUp), forLoops.mu(2, 1), Q) +
     b1(p, dataSet.displayMass(mCharm), forLoops.mu(2, 2), Q) +
     b1(p, forLoops.mt, forLoops.mu(2, 3), Q) +
     b1(p, dataSet.displayMass(mDown), forLoops.md(1, 1), Q) +
     b1(p, dataSet.displayMass(mStrange), forLoops.md(1, 2), Q) +
     b1(p, forLoops.mb, forLoops.md(1, 3), Q) +
     b1(p, dataSet.displayMass(mDown), forLoops.md(2, 1), Q) +
     b1(p, dataSet.displayMass(mStrange), forLoops.md(2, 2), Q) +
     b1(p, forLoops.mb, forLoops.md(2, 3), Q) );

  // Third family mixing contribution:
  delta = delta - forLoops.mt * sin(2.0 * forLoops.thetat) *
    (b0(p, forLoops.mt, forLoops.mu(1, 3), Q) -
     b0(p, forLoops.mt, forLoops.mu(2, 3), Q)) / p;
  delta = delta - forLoops.mb * sin(2.0 * forLoops.thetab) *
    (b0(p, forLoops.mb, forLoops.md(1, 3), Q) -
     b0(p, forLoops.mb, forLoops.md(2, 3), Q)) / p;

  delta = delta * sqr(displayGaugeCoupling(3)) / (16.0 * sqr(PI));

  if (testNan(delta)) {
    if (PRINTOUT > 2) cout << "Gluino loop tachyon\n";
    flagTachyon(true);
    physpars.mGluino = displayGaugino(3);
    return;
  }

  physpars.mGluino = displayGaugino(3) * (1.0 + delta);
}

//  Formulae from hep-ph/9801365: checked but should be checked again!
// Implicitly calculates at the current scale.
double MssmSoftsusy::calcRunningMt() {

  // For brevity
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    mtpole      = dataSet.displayPoleMt();
  double    mt          = forLoops.mt;
  const double  costh   = (displayMw() / displayMz());
  const double    cw2   = sqr(costh) ;
  const double    sw2   = (1.0 - cw2);

  double    ht      = displayYukawaElement(YU, 3, 3);
  double    hb      = displayYukawaElement(YD, 3, 3);
  double    mstop1  = forLoops.mu(1,3);
  double    mstop2  = forLoops.mu(2,3);
  double    mg      = forLoops.mGluino;
  double    thetat  = forLoops.thetat ;
  double    mH      = forLoops.mhiggs(3);
  double    alpha   = forLoops.thetaH ;
  double    g       = displayGaugeCoupling(2);
  double    e       = g * calcSinthdrbar();
  double    mh0     = forLoops.mhiggs(1);
  double    mA      = forLoops.mhiggs(2);
  double    beta    = atan(displayTanb());
  double    mb      = forLoops.mb;
  double    mHc     = forLoops.mhiggs(4);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    thetab  = forLoops.thetab;
  double    mz      = displayMzRun();
  double    q       = displayMu();
  double resigmat = 0.0;

  double qcd = 0.0, stopGluino = 0.0, twoLoopQcd = 0.0 , higgs = 0.0;

  // 1 loop QCD only -- DRbar 10% correction
  qcd = - (5.0 + 6.0 * log(displayMu() / mt)) * 4.0 *
    sqr(displayGaugeCoupling(3)) / 3.0;
  resigmat = resigmat + qcd;

  double p = mtpole;
  // stop/gluino correction 6% correction
  stopGluino = 4.0 * sqr(displayGaugeCoupling(3)) / 3.0 *
    (b1(p, mg, mstop1, displayMu()) +
     b1(p, mg, mstop2, displayMu()) -
     sin(2.0 * thetat) * mg / mtpole *
     (b0(p, mg, mstop1, displayMu()) -
      b0(p, mg, mstop2, displayMu())));

  resigmat = resigmat + stopGluino;

  // 2 loop QCD: hep-ph/0210258 -- debugged 15-6-03
  double l = 2.0 * log(mt / displayMu());
  twoLoopQcd = sqr(sqr(displayGaugeCoupling(3))) *
    (-0.538314 + 0.181534*l - 0.0379954*sqr(l));
  resigmat = resigmat + twoLoopQcd;

  // 2 loop QCD involving MSSM sparticles -- hep-ph/0210258, in the
  // approximation that all squarks and the gluino
  // have mass mSUSY: a few per mille error induced at SPS1a.
  /*
  const static double cf = 4.0 / 3.0, ca = 3.0;
  double m = sqrt(forLoops.mu(1, 3) * forLoops.mu(2, 3));
  double aq = displaySoftA(UA, 3, 3) - displaySusyMu() / displayTanb();
  double logMoQsq = 2.0 * log(m / q);
  double twoLoopMssm = -cf * sqr(sqr(displayGaugeCoupling(3))) /
    (16.0 * sqr(PI)) *
    (47.0 / 3.0 + 20.0 * logMoQsq + 12.0 * logMoQsq * log(m / mt) +
     cf * (23.0 / 24.0 - 13.0 / 6.0 * logMoQsq + 0.5 * sqr(logMoQsq) -
           6.0 * logMoQsq * log(mt / q)) +
     ca * (175.0 / 72.0 + 41.0 / 6.0 * logMoQsq - 0.5 * sqr(logMoQsq) -
           4.0 * logMoQsq * log(mt / q)) +
     aq / m * (-4.0 - 8.0 * logMoQsq) +
     cf * aq / m * (7.0 / 3.0 - 11.0 / 3.0 * logMoQsq + 6.0 * log(mt / q)) +
     ca * aq / m * (-8.0 / 3.0 + 4.0 * logMoQsq));

  resigmat = resigmat + twoLoopMssm;
  */

  // rest are extra bits from Matchev et al: 2% corrections
  double gtL = 0.5 - 2.0 * sw2 / 3.0, gtR = 2.0 * sw2 / 3.0;

  higgs = sqr(ht) / 2.0 *
    (sqr(sin(alpha)) * (b1(p, mt, mH, q) + b0(p, mt, mH, q))
     + sqr(cos(alpha)) * (b1(p, mt, mh0, q) +
                          b0(p, mt, mh0, q))
     + sqr(cos(beta)) * (b1(p, mt, mA, q) -
                         b0(p, mt, mA, q))
     + sqr(sin(beta)) * (b1(p, mt, mz, q) -
                         b0(p, mt, mz, q))) +
    0.5 * ( (sqr(hb) * sqr(sin(beta)) + sqr(ht) * sqr(cos(beta))) *
            b1(p, mb, mHc, q) +
            (sqr(g) + sqr(hb) * sqr(cos(beta)) + sqr(ht) * sqr(sin(beta))) *
            b1(p, mb, displayMwRun(), q)) +
    sqr(hb) * sqr(cos(beta)) * (b0(p, mb, mHc, q)
                                - b0(p, mb, displayMwRun(), q)) -
    sqr(e) * 4.0 / 9.0 * (5.0 + 6.0 * log(q / mt)) +
    sqr(g) / cw2DRbar * ( (sqr(gtL) + sqr(gtR)) * b1(p, mt, mz, q) +
                          4.0 * gtL * gtR * b0(p, mt, mz, q) );
  resigmat = resigmat + higgs;

  DoubleMatrix neutralinoContribution(4, 2);
  // Neutralino contribution
  DoubleVector aPsi0TStopr(4), bPsi0TStopr(4), aPsi0TStopl(4),
    bPsi0TStopl(4);
  aPsi0TStopr(1) = - 4 * gp / (3.0 * root2);
  bPsi0TStopl(1) = gp / (3.0 * root2);
  bPsi0TStopl(2) = g / root2;
  aPsi0TStopl(4) = ht;
  bPsi0TStopr(4) = ht;

  ComplexVector aChi0TStopl(4), bChi0TStopl(4), aChi0TStopr(4),
    bChi0TStopr(4);

  // Neutralinos
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  ComplexMatrix aNeutTStop(4, 2), bNeutTStop(4, 2);
  DoubleMatrix fNeutTStop(4, 2), gNeutTStop(4, 2);

  int i, j;
  DoubleMatrix O(2, 2);
  ComplexVector tt(2), t1(2), t2(2);
  O = rot2d(thetat);
  for (i=1; i<=4; i++) {
    tt(1) = aChi0TStopl(i); tt(2) = aChi0TStopr(i);
    t1 = O * tt;

    tt(1) = bChi0TStopl(i); tt(2) = bChi0TStopr(i);
    t2 = O * tt;
    for (j=1; j<=2; j++) {
      aNeutTStop(i, j) = t1(j);
      bNeutTStop(i, j) = t2(j);
      // functions of couplings needed for loops
      fNeutTStop(i, j) = sqr(aNeutTStop(i, j).mod()) +
        sqr(bNeutTStop(i, j).mod());

      gNeutTStop(i, j) = 2.0 *
        (aNeutTStop(i, j) * bNeutTStop(i, j).conj()).real();

      neutralinoContribution(i, j) = (fNeutTStop(i, j) *
         b1(p, mneut(i), forLoops.mu(j, 3), q) +
         gNeutTStop(i, j) * mneut(i) /  mtpole *
         b0(p, mneut(i), forLoops.mu(j, 3), q)) * 0.5;

      resigmat = resigmat + neutralinoContribution(i, j);
    }
  }

  // Chargino contribution
  DoubleVector bPsicTSbotl(2), bPsicTSbotr(2), aPsicTSbotl(2);

  bPsicTSbotl(1) = g;
  bPsicTSbotr(2) = -hb;
  aPsicTSbotl(2) = -ht;

  DoubleVector aPsicTSbotr(2), aPsicCSbotl(2);
  ComplexVector aChicTSbotr(2), aChicTSbotl(2), bChicTSbotl(2),
      bChicTSbotr(2);
  ComplexMatrix aChTSbot(2, 2), bChTSbot(2, 2);
  DoubleMatrix fChTSbot(2, 2), gChTSbot(2, 2);

  aChicTSbotl = v.complexConjugate() * aPsicTSbotl;
  bChicTSbotl = u * bPsicTSbotl;
  aChicTSbotr = v.complexConjugate() * aPsicTSbotr;
  bChicTSbotr = u * bPsicTSbotr;

  double charginoContribution = 0.0;

  for(i=1; i<=2; i++) {
    O = rot2d(thetab);
    tt(1) = aChicTSbotl(i); tt(2) = aChicTSbotr(i);
    t1 = O * tt;
    tt(1) = bChicTSbotl(i); tt(2) = bChicTSbotr(i);
    t2 = O * tt;
    for (j=1; j<=2; j++) {
        aChTSbot(i, j) = t1(j);
        bChTSbot(i, j) = t2(j);

        fChTSbot(i, j) = sqr(aChTSbot(i, j).mod()) +
          sqr(bChTSbot(i, j).mod());
        gChTSbot(i, j) = 2.0 * (aChTSbot(i, j) *
                                bChTSbot(i, j).conj()).real();

        charginoContribution = charginoContribution +
          (fChTSbot(i, j) *
           b1(p, mch(i), fabs(forLoops.md(j, 3)),
              q) +
           gChTSbot(i, j) * mch(i) / mtpole *
           b0(p, mch(i), forLoops.md(j, 3), q)) * 0.5;
    }
  }

  resigmat = resigmat + charginoContribution;

  resigmat = resigmat * mtpole / (16.0 * sqr(PI));

  return mtpole + resigmat;
}

double MssmSoftsusy::calcRunningMb() const {

  if (displayMu() != displayMz()) {
    ostringstream ii;
    ii << "MssmSoftsusy::calcRunningMb called with mu=" <<
      displayMu() << endl;
    throw ii.str();
  }

  drBarPars forLoops(displayDrBarPars());

  double mbMZ = dataSet.displayMass(mBottom),
    alphasMZ = sqr(displayGaugeCoupling(3)) / (4.0 * PI);

  double    msbot1  = forLoops.md(1,3);
  double    msbot2  = forLoops.md(2,3);
  double    mg      = forLoops.mGluino;
  double    thetab  = forLoops.thetab;
  double    thetat  = forLoops.thetat;
  double    g       = displayGaugeCoupling(2);
  double    g1      = displayGaugeCoupling(1);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    mbMSSM  = forLoops.mb;
  double    hb      = displayYukawaElement(YD, 3, 3);
  double    ht      = displayYukawaElement(YU, 3, 3);
  double    mb      = forLoops.mb;
  double    mt      = forLoops.mt;
  double    mh      = forLoops.mhiggs(1);
  double    mA      = forLoops.mhiggs(2);
  double    mH      = forLoops.mhiggs(3);
  double    mHp     = forLoops.mhiggs(4);
  double    mz = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    ca      = cos(forLoops.thetaH);
  double    sa      = sin(forLoops.thetaH);
  double    cosb    = cos(atan(displayTanb()));
  double    sinb    = sin(atan(displayTanb()));

  double p = mbMZ;
  double q = displayMu();

  // First convert mbMZ into DRbar value from hep-ph/9703293,0207126,9701308
  // (SM gauge boson contributions)
  mbMZ = mbMZ *
    (1.0 - alphasMZ / (3.0 * PI) - 23.0 / 72.0 * sqr(alphasMZ) / sqr(PI) +
  3.0 * sqr(g) / (128.0 * sqr(PI)) +
  13.0 * sqr(g1) / (1152. * sqr(PI)));

  double deltaSquarkGluino = - alphasMZ / (3.0 * PI) *
    (b1(p, mg, msbot1, displayMu()) +
     b1(p, mg, msbot2, displayMu()) -
     sin(2.0 * thetab) * mg / mbMSSM *
     (b0(p, mg, msbot1, displayMu()) -
      b0(p, mg, msbot2, displayMu())));

  // Chargino contribution
  DoubleVector bPsicBstopl(2), bPsicBstopr(2),
    aPsicBstopl(2), aPsicBstopr(2);

  aPsicBstopl(1) = g;
  aPsicBstopr(2) = -displayYukawaElement(YU, 3, 3);
  bPsicBstopl(2) = -displayYukawaElement(YD, 3, 3);

  DoubleVector aPsicCStopl(2);
  ComplexVector aChicBstopr(2), aChicBstopl(2), bChicBstopl(2),
      bChicBstopr(2);
  ComplexMatrix aChBstop(2, 2), bChBstop(2, 2);
  DoubleMatrix fChBstop(2, 2), gChBstop(2, 2);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChicBstopl = v.complexConjugate() * aPsicBstopl;
  bChicBstopl = u * bPsicBstopl;
  aChicBstopr = v.complexConjugate() * aPsicBstopr;
  bChicBstopr = u * bPsicBstopr;

  double charginoContribution = 0.0;

  int i, j; DoubleMatrix O(2, 2); ComplexVector tt(2), t1(2), t2(2);
  for(i=1; i<=2; i++) {
    O = rot2d(thetat);
    tt(1) = aChicBstopl(i); tt(2) = aChicBstopr(i);
    t1 = O * tt;
    tt(1) = bChicBstopl(i); tt(2) = bChicBstopr(i);
    t2 = O * tt;
    for (j=1; j<=2; j++) {
        aChBstop(i, j) = t1(j);
        bChBstop(i, j) = t2(j);

        fChBstop(i, j) = sqr(aChBstop(i, j).mod()) +
          sqr(bChBstop(i, j).mod());
        gChBstop(i, j) = 2.0 * (aChBstop(i, j) *
                                bChBstop(i, j).conj()).real();

        charginoContribution = charginoContribution +
          (fChBstop(i, j) *
           b1(p, mch(i), fabs(forLoops.mu(j, 3)),
              q) +
           gChBstop(i, j) * mch(i) / mbMSSM *
           b0(p, mch(i), forLoops.mu(j, 3), q)) * 0.5;
    }
  }
  double deltaSquarkChargino = -charginoContribution / (16.0 * sqr(PI));

  double deltaHiggs = 0.;
  // new corrections follow and they must be checked! Neutralinos to follow...
  // Higgs
  deltaHiggs = 0.5 * sqr(hb) *
    (sqr(ca) * (b1(p, mb, mH, q) + b0(p, mb, mH, q)) +
     sqr(sa) * (b1(p, mb, mh, q) + b0(p, mb, mh, q)) +
     sqr(sinb) * (b1(p, mb, mA, q) - b0(p, mb, mA, q)) +
     sqr(cosb) * (b1(p, mb, mz, q) - b0(p, mb, mz, q))) +
    0.5 * ((sqr(ht) * sqr(cosb) + sqr(hb) * sqr(sinb)) * b1(p, mt, mHp, q) +
           (sqr(g) + sqr(ht) * sqr(sinb) + sqr(hb) * sqr(cosb)) *
           b1(p, mt, mw, q)) +
    sqr(ht) * sqr(sinb) * (b0(p, mt, mHp, q) - b0(p, mt, mw, q)) +
    sqr(g) / cw2DRbar * ((sqr(gdL) + sqr(gdR)) * b1(p, mb, mz, q) +
                                    4.0 * gdL * gdR * b0(p, mb, mz, q));

  deltaHiggs = - deltaHiggs / (16.0 * sqr(PI));

  // Neutralinos
  DoubleVector aPsi0Bsbotr(4), bPsi0Bsbotr(4), aPsi0Bsbotl(4),
    bPsi0Bsbotl(4);
  aPsi0Bsbotr(1) = gp / (root2 * 3.0) * 2.0;
  bPsi0Bsbotl(1) = gp / (root2 * 3.0);
  bPsi0Bsbotl(2) = -root2 * g * 0.5;
  aPsi0Bsbotl(3) = hb;
  bPsi0Bsbotr(3) = hb;

  ComplexVector aChi0Bsbotl(4), bChi0Bsbotl(4), aChi0Bsbotr(4),
    bChi0Bsbotr(4);

  double deltaNeutralino = 0.;
  aChi0Bsbotl = n.complexConjugate() * aPsi0Bsbotl;
  bChi0Bsbotl = n * bPsi0Bsbotl;
  aChi0Bsbotr = n.complexConjugate() * aPsi0Bsbotr;
  bChi0Bsbotr = n * bPsi0Bsbotr;

  ComplexMatrix aNeutBsbot(4, 2), bNeutBsbot(4, 2);
  DoubleMatrix fNeutBsbot(4, 2), gNeutBsbot(4, 2),
    neutralinoContribution(4, 2);

  O = rot2d(thetab);
  for (i=1; i<=4; i++) {
    tt(1) = aChi0Bsbotl(i); tt(2) = aChi0Bsbotr(i);
    t1 = O * tt;

    tt(1) = bChi0Bsbotl(i); tt(2) = bChi0Bsbotr(i);
    t2 = O * tt;
    for (j=1; j<=2; j++) {
      aNeutBsbot(i, j) = t1(j);
      bNeutBsbot(i, j) = t2(j);
      // functions of couplings needed for loops
      fNeutBsbot(i, j) = sqr(aNeutBsbot(i, j).mod()) +
        sqr(bNeutBsbot(i, j).mod());

      gNeutBsbot(i, j) = 2.0 *
        (aNeutBsbot(i, j) * bNeutBsbot(i, j).conj()).real();

      neutralinoContribution(i, j) = (fNeutBsbot(i, j) *
         b1(p, mneut(i), forLoops.md(j, 3), q) +
         gNeutBsbot(i, j) * mneut(i) /  mbMSSM *
         b0(p, mneut(i), forLoops.md(j, 3), q)) * 0.5;

      deltaNeutralino = deltaNeutralino + neutralinoContribution(i, j);
    }
  }

  deltaNeutralino = -deltaNeutralino / (16.0 * sqr(PI));

  // it's NOT clear if this resummation is reliable in the full 1-loop scheme
  // but it's at least valid to 1 loop. Warning though: if you add higher
  // loops, you'll have to re-arrange.
  return mbMZ / (1.0 + deltaSquarkGluino + deltaSquarkChargino + deltaHiggs
                 + deltaNeutralino);
}


// Full BPMZ expression
double MssmSoftsusy::calcRunningMtau() const {

  drBarPars forLoops(displayDrBarPars());

  // MSbar value
  double mTauSMMZ = displayDataSet().displayMass(mTau);
  double mTauPole = MTAU;

  // conversion to DRbar
  mTauSMMZ = mTauSMMZ *
    (1.0 - 3.0 * (sqr(displayGaugeCoupling(1)) - sqr(displayGaugeCoupling(2)))
     / (128.0 * sqr(PI)));

  double    msnutau = fabs(forLoops.msnu(3));
  double    thetatau= forLoops.thetatau;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    htau    = displayYukawaElement(YE, 3, 3);
  double    mtau    = forLoops.mtau;
  double    mh      = forLoops.mhiggs(1);
  double    mA      = forLoops.mhiggs(2);
  double    mH      = forLoops.mhiggs(3);
  double    mHp     = forLoops.mhiggs(4);
  double    mz = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    ca      = cos(forLoops.thetaH);
  double    sa      = sin(forLoops.thetaH);
  double    cosb    = cos(atan(displayTanb()));
  double    sinb    = sin(atan(displayTanb()));

  double p = mTauPole;
  double q = displayMu();

  // Chargino contribution
  DoubleVector aPsicTauSnul(2), bPsicTauSnul(2);

  aPsicTauSnul(1) = g;
  bPsicTauSnul(2) = -displayYukawaElement(YE, 3, 3);

  ComplexVector aChicTauSnul(2), bChicTauSnul(2);
  DoubleVector fChiTauSnu(2), gChiTauSnu(2);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  // Mass eignebasis of charginos
  aChicTauSnul = v.complexConjugate() * aPsicTauSnul;
  bChicTauSnul = u * bPsicTauSnul;

  DoubleVector charg(2);
  int i;
  for(i=1; i<=2; i++) {
    fChiTauSnu(i) = sqr(aChicTauSnul(i).mod()) +
      sqr(bChicTauSnul(i).mod());
    gChiTauSnu(i) = 2.0 * (aChicTauSnul(i) * bChicTauSnul(i).conj()).real();

    charg(i) =
      (fChiTauSnu(i) * b1(p, mch(i), msnutau, q) +
       gChiTauSnu(i) * mch(i) / mTauSMMZ * b0(p, mch(i), msnutau, q)) * 0.5;
    }

  double sigmaChargino = (charg(1) + charg(2)) / (16.0 * sqr(PI));
  // checked charginos

  // Higgs
  double mnu = 0.;
  double sigmaHiggs = 0.5 * sqr(htau) *
    (sqr(ca) * (b1(p, mtau, mH, q) + b0(p, mtau, mH, q)) +
     sqr(sa) * (b1(p, mtau, mh, q) + b0(p, mtau, mh, q)) +
     sqr(sinb) * (b1(p, mtau, mA, q) - b0(p, mtau, mA, q)) +
     sqr(cosb) * (b1(p, mtau, mz, q) - b0(p, mtau, mz, q))) +
    0.5 * (sqr(htau) * sqr(sinb) * b1(p, mnu, mHp, q) +
           (sqr(g) + sqr(htau) * sqr(cosb)) * b1(p, mnu, mw, q)) +
    sqr(g) / cw2DRbar * ((sqr(geL) + sqr(geR)) * b1(p, mtau, mz, q) +
    4.0 * geL * geR * b0(p, mtau, mz, q));

  sigmaHiggs = sigmaHiggs / (16.0 * sqr(PI));

  // Neutralinos
  DoubleVector aPsi0TauStaur(4), bPsi0TauStaur(4), aPsi0TauStaul(4),
    bPsi0TauStaul(4);
  aPsi0TauStaur(1) = gp / root2 * 2.0;
  bPsi0TauStaul(1) = -gp / root2;
  bPsi0TauStaul(2) = -root2 * g * 0.5;
  aPsi0TauStaul(3) = htau;
  bPsi0TauStaur(3) = htau;

  ComplexVector aChi0TauStaul(4), bChi0TauStaul(4), aChi0TauStaur(4),
    bChi0TauStaur(4);

  double sigmaNeutralino = 0.;
  aChi0TauStaul = n.complexConjugate() * aPsi0TauStaul;
  bChi0TauStaul = n * bPsi0TauStaul;
  aChi0TauStaur = n.complexConjugate() * aPsi0TauStaur;
  bChi0TauStaur = n * bPsi0TauStaur;

  ComplexMatrix aNeutTauStau(4, 2), bNeutTauStau(4, 2);
  DoubleMatrix fNeutTauStau(4, 2), gNeutTauStau(4, 2),
    neutralinoContribution(4, 2);

  DoubleMatrix O(2, 2); O = rot2d(thetatau);
  ComplexVector t1(2), t2(2), tt(2);
  for (i=1; i<=4; i++) {
    tt(1) = aChi0TauStaul(i); tt(2) = aChi0TauStaur(i);
    t1 = O * tt;

    tt(1) = bChi0TauStaul(i); tt(2) = bChi0TauStaur(i);
    t2 = O * tt;
    int j;
    for (j=1; j<=2; j++) {
      aNeutTauStau(i, j) = t1(j);
      bNeutTauStau(i, j) = t2(j);
      // functions of couplings needed for loops
      fNeutTauStau(i, j) = sqr(aNeutTauStau(i, j).mod()) +
        sqr(bNeutTauStau(i, j).mod());

      gNeutTauStau(i, j) = 2.0 *
        (aNeutTauStau(i, j) * bNeutTauStau(i, j).conj()).real();

      neutralinoContribution(i, j) = (fNeutTauStau(i, j) *
         b1(p, mneut(i), forLoops.me(j, 3), q) +
         gNeutTauStau(i, j) * mneut(i) /  mTauSMMZ *
         b0(p, mneut(i), forLoops.me(j, 3), q)) * 0.5;

      sigmaNeutralino = sigmaNeutralino + neutralinoContribution(i, j);
    }
  }

  sigmaNeutralino = sigmaNeutralino / (16.0 * sqr(PI));

  // old calculation of tau mass
  /*  double delta = sqr(displayGaugeCoupling(2)) / (16 * sqr(PI)) *
    (-displaySusyMu()) * displayGaugino(2) * displayTanb() /
    (sqr(displaySusyMu()) - sqr(displayGaugino(2))) *
    (b0(mTauPole, displayGaugino(2),
        forLoops.msnu(3), displayMu()) -
        b0(mTauPole, -displaySusyMu(), forLoops.msnu(3), displayMu()));*/

  // From hep-ph/9912516
  return mTauSMMZ * (1.0 + sigmaNeutralino + sigmaChargino + sigmaHiggs);
}

void MssmSoftsusy::treeUpSquark(DoubleMatrix & mass, double mtrun,
                                double pizztMS, double sinthDRbarMS,
                                int family) {
  const double cu = 2.0 / 3.0;
  double mz2 = sqr(displayMzRun()), mt2 = sqr(mtrun);
  double beta = atan(displayTanb()), mu = displaySusyMu(),
    c2b = cos(2 * beta), tanb = displayTanb(), sinth2 = sqr(sinthDRbarMS);

  mass(1, 1) = displaySoftMassSquared(mQl, family, family) +
    mz2 * (0.5 - cu * sinth2) * c2b;
  mass(2, 2) = displaySoftMassSquared(mUr, family, family) +
    mz2 * cu * sinth2 * c2b;

  if (family != 3) mass(1, 2) = 0.0;
  else {
    mass(1, 2) = mtrun * (displaySoftA(UA, 3, 3) - mu / tanb);
    mass(1, 1) = mass(1, 1) + mt2;
    mass(2, 2) = mass(2, 2) + mt2;
  }

  mass(2, 1) = mass(1, 2);
}

// Now these are calculated at the squark scale
void MssmSoftsusy::addSquarkCorrection(DoubleMatrix & mass) {

// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) {
    if (PRINTOUT > 2) cout << " squark tachyonic ";
    flagTachyon(true);
    return;
  }

  DoubleVector x(2), delta(2);
  int i; for (i=1; i<=2; i++)
    {
      x(i)= sqr(forLoops.mGluino) / mass(i, i);

      if (close(x(i), 1.0, EPSTOL))
        delta(i) = sqr(displayGaugeCoupling(3)) / (6.0 * sqr(PI)) *
          (1.0 + 3.0 * x(i)
           - sqr(x(i)) * log(fabs(x(i))) + 2.0 * x(i) *
           log(sqr(displayMu()) / mass(i, i)));
      else
        delta(i) = sqr(displayGaugeCoupling(3)) / (6.0 * sqr(PI)) *
          (1.0 + 3.0 * x(i) + sqr(x(i) - 1.0) * log(fabs(x(i) - 1.0))
           - sqr(x(i)) * log(fabs(x(i))) + 2.0 * x(i) *
           log(sqr(displayMu()) / mass(i, i)));

      mass(i, i) = mass(i, i) * (1.0 + delta(i));
    }
}

void MssmSoftsusy::addSnuTauCorrection(double & mass) {

  // No point adding radiative corrections to tachyonic particles
  if (mass < 0.0) {
    if (PRINTOUT > 2) cout << " snu_tau tachyonic ";
    flagTachyon(true);
    return;
  }

  double p = sqrt(mass);

  // one-loop correction matrix
  double piSq; // Self-energy matrix

  // brevity
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());
  double mtau = displayYukawaElement(YE, 3, 3) * displayHvev() / root2 *
    cos(beta);

  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double   thetatau = forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double  ctau      = cos(thetatau);
  double  stau      = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  double    smu     = -displaySusyMu();
  double q = displayMu(),
    htau = displayYukawaElement(YE, 3, 3),
    htausq = sqr(htau),
    sinb = sin(beta), cosb = cos(beta),
    v1 = displayHvev() * cos(beta),
    mz = displayMzRun();

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2), dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleVector lsSnuLSnuLR(4);
  // Order (s1 s2 G A, L R)
  lsSnuLSnuLR(1) = g * mz * gnuL * cosb / costhDrbar;
  lsSnuLSnuLR(2) = -g * mz * gnuL * sinb / costhDrbar;

  DoubleVector lHSnuLSnu12(lsSnuLSnuLR);
  DoubleVector temp(2), temp2(2);
  // Mix CP-even Higgses up
  temp(1) = lsSnuLSnuLR(1);
  temp(2) = lsSnuLSnuLR(2);
  temp2 = rot2d(alpha) * temp;
  lHSnuLSnu12(1) = temp2(1);
  lHSnuLSnu12(2) = temp2(2);

  // Charged Higgs Feynman rules
  DoubleMatrix lChHsnuLstauLR(2, 2); // (G+ H+, L R) basis
  lChHsnuLstauLR(2, 1) = (g * displayMwRun() * sin(2.0 * beta) / root2
    - htau * mtau * sinb);
  lChHsnuLstauLR(2, 2) = (smu * htau * cosb -
                           displayTrilinear(EA, 3, 3) * sinb);
  lChHsnuLstauLR(1, 1) = (-g * displayMwRun() * cos(2.0 * beta)
    + htausq * v1 * cosb) / root2;
  lChHsnuLstauLR(1, 2) = htau * smu * sinb + displayTrilinear(EA, 3, 3) * cosb;

  DoubleMatrix lChHsnuLstau12(2, 2);
  temp(1) = lChHsnuLstauLR(1, 1);
  temp(2) = lChHsnuLstauLR(1, 2);
  temp2 = rot2d(thetatau) * temp;
  lChHsnuLstau12(1, 1) = temp2(1);
  lChHsnuLstau12(1, 2) = temp2(2);
  temp(1) = lChHsnuLstauLR(2, 1);
  temp(2) = lChHsnuLstauLR(2, 2);
  temp2 = rot2d(thetatau) * temp;
  lChHsnuLstau12(2, 1) = temp2(1);
  lChHsnuLstau12(2, 2) = temp2(2);

  // Neutralino Feynman rules
  DoubleVector aPsi0TSnul(4), bPsi0TSnul(4);
  bPsi0TSnul(1) = gp * ynuL / root2;
  bPsi0TSnul(2) = g / root2;

  ComplexVector aChi0TSnul(4), bChi0TSnul(4);
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0TSnul = n.complexConjugate() * aPsi0TSnul;
  bChi0TSnul = n * bPsi0TSnul;

  DoubleVector gChi0NuSnuLL(4), fChi0NuSnuLL(4);
  int i; for (i=1; i<=4; i++) {
    fChi0NuSnuLL(i) = (aChi0TSnul(i) * aChi0TSnul(i).conj() +
      bChi0TSnul(i) * bChi0TSnul(i).conj()).real();
    gChi0NuSnuLL(i) = (bChi0TSnul(i).conj() * aChi0TSnul(i) +
      bChi0TSnul(i) * aChi0TSnul(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicBSnul(2), aPsicBSnul(2);
  aPsicBSnul(1) = g;
  bPsicBSnul(2) = -htau;

  DoubleVector aPsicCStaul(2);
  ComplexVector aChicBSnur(2), aChicBSnul(2), bChicBSnul(2),
      bChicBSnur(2);
  ComplexMatrix aChBSnu(2, 2), bChBSnu(2, 2);

  aChicBSnul = v.complexConjugate() * aPsicBSnul;
  bChicBSnul = u * bPsicBSnul;

  DoubleVector fChBSnuLL(2), gChBSnuLL(2) ;
  for (i=1; i<=2; i++) {
    fChBSnuLL(i) = (aChicBSnul(i).conj() * aChicBSnul(i) +
                      bChicBSnul(i).conj() * bChicBSnul(i)).real();
    gChBSnuLL(i) = (bChicBSnul(i).conj() * aChicBSnul(i) +
                      aChicBSnul(i).conj() * bChicBSnul(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  double stop = 0., sbottom = 0., higgs = 0., chargino = 0.,
    neutralino = 0.;

  sbottom =
    htausq * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q));

  for (i=1; i<=4; i++) {
    higgs = higgs +
      0.5 * (- sqr(g) * gnuL * 0.5 / sqr(costhDrbar) * cn(i))
      * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs = higgs +
      (htausq * dnu(i) + sqr(g) *
       (gnuL * 0.5 / costhDrbar2 - 0.5) * cn(i))* a0p;
  }
  int j; for(i=1; i<=4; i++) {
    double b0p = b0(p, higgsm(i), msnu(3), q);
   higgs = higgs + sqr(lHSnuLSnu12(i)) * b0p;
  }

 for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), mstau(j), q);
      higgs = higgs + sqr(lChHsnuLstau12(i, j)) * b0p;
    }

  // EW bosons
  higgs = higgs +
    4.0 * sqr(g) / costhDrbar2 * sqr(gnuL) * a0(mz, q) +
    2.0 * sqr(g) * a0(displayMwRun(), q) +
    sqr(g * gnuL / costhDrbar) * ffn(p, msnu(3), mz, q) +
    sqr(g) * 0.5 * (sqr(ctau) * ffn(p, mstau(1), displayMwRun(), q) +
                    sqr(stau) * ffn(p, mstau(2), displayMwRun(), q)) +
    sqr(g) * 0.25 *
    (a0(msnu(3), q) + 2.0 *
     (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))) +
    sqr(g) * 0.5 *
    (1.5 * a0(forLoops.mu(1, 1), q) +
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) +
            sqr(st) * a0(forLoops.mu(2, 2), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
            sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) +
            a0(forLoops.msnu(2), q) +
            a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(ynuL) * a0(msnu(3), q) +
    sqr(gp) * 0.25 * ynuL *
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) +
                  a0(forLoops.mu(1, 2), q) +
                  sqr(ct) * a0(forLoops.mu(1, 3), q) +
                  sqr(st) * a0(forLoops.mu(2, 2), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) +
                  a0(forLoops.mu(2, 2), q) +
                  sqr(st) * a0(forLoops.mu(1, 3), q) +
                  sqr(ct) * a0(forLoops.mu(2, 2), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) +
                  a0(forLoops.md(1, 2), q) +
                  sqr(cb) * a0(forLoops.md(1, 3), q) +
                  sqr(sb) * a0(forLoops.md(2, 2), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) +
                  a0(forLoops.md(2, 2), q) +
                  sqr(sb) * a0(forLoops.md(1, 3), q) +
                  sqr(cb) * a0(forLoops.md(2, 2), q)) +
     yeL * (a0(forLoops.me(1, 1), q) +
            a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 2), q)) +
     yeR * (a0(forLoops.me(2, 1), q) +
            a0(forLoops.me(2, 2), q) +
            sqr(stau) * a0(forLoops.me(1, 3), q) +
            sqr(ctau) * a0(forLoops.me(2, 2), q)) +
     ynuL * (a0(forLoops.msnu(1), q) +
             a0(forLoops.msnu(2), q) +
             a0(forLoops.msnu(3), q)));

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), mtau, q);
    double two = mch(i) * mtau * b0(p, mch(i), mtau, q) * 2.0;
    chargino = chargino + fChBSnuLL(i) * one - gChBSnuLL(i) * two;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p, mneut(i), 0., q);
    neutralino = neutralino + fChi0NuSnuLL(i) * one;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (stop + sbottom + higgs + chargino + neutralino);

  mass = mass - piSq;
}

/// Found+fixed bug 7/09/06. Thanks to J Kersten.
void MssmSoftsusy::addSnuCorrection(double & mass, int family) {

  // No point adding radiative corrections to tachyonic particles
  if (mass < 0.0) {
    if (PRINTOUT > 2) cout << " snu_tau tachyonic ";
    flagTachyon(true);
    return;
  }

  double p = sqrt(mass);

  // one-loop correction matrix
  double piSq; // Self-energy

  // brevity
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());

  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double   thetatau = forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double  ctau      = cos(thetatau);
  double  stau      = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  double q = displayMu(),
    sinb = sin(beta), cosb = cos(beta),
    mz = displayMzRun();

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2), dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleVector lsSnuLSnuLR(4);
  // Order (s1 s2 G A, L R)
  lsSnuLSnuLR(1) = g * mz * gnuL * cosb / costhDrbar;
  lsSnuLSnuLR(2) = -g * mz * gnuL * sinb / costhDrbar;

  DoubleVector lHSnuLSnu12(lsSnuLSnuLR);
  DoubleVector temp(2), temp2(2);
  // Mix CP-even Higgses up
  temp(1) = lsSnuLSnuLR(1);
  temp(2) = lsSnuLSnuLR(2);
  temp2 = rot2d(alpha) * temp;
  lHSnuLSnu12(1) = temp2(1);
  lHSnuLSnu12(2) = temp2(2);

  // Charged Higgs Feynman rules
  DoubleMatrix lChHsnuLstauLR(2, 2); // (G+ H+, L R) basis
  lChHsnuLstauLR(2, 1) = (g * displayMwRun() * sin(2.0 * beta) / root2);
  lChHsnuLstauLR(1, 1) = -g * displayMwRun() * cos(2.0 * beta);

  DoubleMatrix lChHsnuLstau12(2, 2);
  temp(1) = lChHsnuLstauLR(1, 1);
  temp(2) = lChHsnuLstauLR(1, 2);
  temp2 = temp;
  lChHsnuLstau12(1, 1) = temp2(1);
  lChHsnuLstau12(1, 2) = temp2(2);
  temp(1) = lChHsnuLstauLR(2, 1);
  temp(2) = lChHsnuLstauLR(2, 2);
  temp2 = temp;
  lChHsnuLstau12(2, 1) = temp2(1);
  lChHsnuLstau12(2, 2) = temp2(2);

  // Neutralino Feynman rules
  DoubleVector aPsi0TSnul(4), bPsi0TSnul(4);
  bPsi0TSnul(1) = gp * ynuL / root2;
  bPsi0TSnul(2) = g / root2;

  ComplexVector aChi0TSnul(4), bChi0TSnul(4);
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0TSnul = n.complexConjugate() * aPsi0TSnul;
  bChi0TSnul = n * bPsi0TSnul;

  DoubleVector gChi0NuSnuLL(4), fChi0NuSnuLL(4);
  int i; for (i=1; i<=4; i++) {
    fChi0NuSnuLL(i) = (aChi0TSnul(i) * aChi0TSnul(i).conj() +
      bChi0TSnul(i) * bChi0TSnul(i).conj()).real();
    gChi0NuSnuLL(i) = (bChi0TSnul(i).conj() * aChi0TSnul(i) +
      bChi0TSnul(i) * aChi0TSnul(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicBSnul(2), aPsicBSnul(2);
  aPsicBSnul(1) = g;

  DoubleVector aPsicCStaul(2);
  ComplexVector aChicBSnur(2), aChicBSnul(2), bChicBSnul(2),
      bChicBSnur(2);
  ComplexMatrix aChBSnu(2, 2), bChBSnu(2, 2);

  aChicBSnul = v.complexConjugate() * aPsicBSnul;
  bChicBSnul = u * bPsicBSnul;

  DoubleVector fChBSnuLL(2), gChBSnuLL(2) ;
  for (i=1; i<=2; i++) {
    fChBSnuLL(i) = (aChicBSnul(i).conj() * aChicBSnul(i) +
                      bChicBSnul(i).conj() * bChicBSnul(i)).real();
    gChBSnuLL(i) = (bChicBSnul(i).conj() * aChicBSnul(i) +
                      aChicBSnul(i).conj() * bChicBSnul(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  double higgs = 0., chargino = 0., neutralino = 0., sfermions = 0.;

  for (i=1; i<=4; i++) {
    higgs = higgs +
      0.5 * (- sqr(g) * gnuL * 0.5 / sqr(costhDrbar) * cn(i))
      * a0(higgsm(i), q);
  }

  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs = higgs +
      (gnuL * 0.5 / costhDrbar2 - 0.5) * cn(i)* a0p * sqr(g);
  }

  int j; for(i=1; i<=4; i++) {
    double b0p = b0(p, higgsm(i), msnu(family), q);
   higgs = higgs + sqr(lHSnuLSnu12(i)) * b0p;
  }

  double meL = forLoops.me(1, family);
  DoubleVector msel(2);
  msel(1) = forLoops.me(1, family);
  msel(2) = forLoops.me(2, family);
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), msel(j), q);
      higgs = higgs + sqr(lChHsnuLstau12(i, j)) * b0p;
  }

  // EW bosons
  higgs = higgs +
    4.0 * sqr(g) / costhDrbar2 * sqr(gnuL) * a0(mz, q) +
    2.0 * sqr(g) * a0(displayMwRun(), q) +
    sqr(g * gnuL / costhDrbar) * ffn(p, msnu(family), mz, q) +
    sqr(g) * 0.5 * ffn(p, meL, displayMwRun(), q);

  sfermions = sfermions +
    sqr(g) * 0.25 *
    (a0(msnu(family), q) + 2.0 * a0(meL, q)) + sqr(g) * 0.5 *
    (1.5 * a0(forLoops.mu(1, 1), q) +
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) +
            sqr(st) * a0(forLoops.mu(2, 2), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
            sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) +
            a0(forLoops.msnu(2), q) +
            a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(ynuL) * a0(msnu(family), q) +
    sqr(gp) * 0.25 * ynuL *
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) +
                  a0(forLoops.mu(1, 2), q) +
                  sqr(ct) * a0(forLoops.mu(1, 3), q) +
                  sqr(st) * a0(forLoops.mu(2, 2), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) +
                  a0(forLoops.mu(2, 2), q) +
                  sqr(st) * a0(forLoops.mu(1, 3), q) +
                  sqr(ct) * a0(forLoops.mu(2, 2), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) +
                  a0(forLoops.md(1, 2), q) +
                  sqr(cb) * a0(forLoops.md(1, 3), q) +
                  sqr(sb) * a0(forLoops.md(2, 2), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) +
                  a0(forLoops.md(2, 2), q) +
                  sqr(sb) * a0(forLoops.md(1, 3), q) +
                  sqr(cb) * a0(forLoops.md(2, 2), q)) +
     yeL * (a0(forLoops.me(1, 1), q) +
            a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 2), q)) +
     yeR * (a0(forLoops.me(2, 1), q) +
            a0(forLoops.me(2, 2), q) +
            sqr(stau) * a0(forLoops.me(1, 3), q) +
            sqr(ctau) * a0(forLoops.me(2, 2), q)) +
     ynuL * (a0(forLoops.msnu(1), q) +
             a0(forLoops.msnu(2), q) +
             a0(forLoops.msnu(3), q)));

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), 0., q);
    chargino = chargino + fChBSnuLL(i) * one;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p, mneut(i), 0., q);
    neutralino = neutralino + fChi0NuSnuLL(i) * one;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (higgs + chargino + neutralino + sfermions);

  mass = mass - piSq;
}

// As in BPMZ appendix, INCLUDING weak boson loops.
void MssmSoftsusy::addStopCorrection(double p, DoubleMatrix & mass,
                                     double mt) {

// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) {
    if (PRINTOUT > 2) cout << " stop tachyonic ";
    flagTachyon(true);
    return;
  }

  // one-loop correction matrix
  DoubleMatrix piSq(2, 2); // Self-energy matrix

  // brevity
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());
  double    e       = g * sinthDrbar; // DRbar value of e

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double   thetatau = forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double  ctau      = cos(thetatau);
  double  stau      = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  double msbot1 = forLoops.md(1, 3), msbot2 = forLoops.md(2, 3);
  double    mg      = forLoops.mGluino;
  double    smu     = -displaySusyMu();
  double q = displayMu(), g3sq = sqr(displayGaugeCoupling(3)),
    ht = displayYukawaElement(YU, 3, 3),
    hb = displayYukawaElement(YD, 3, 3),
    htsq = sqr(displayYukawaElement(YU, 3, 3)),
    sinb = sin(beta), cosb = cos(beta),
    hbsq = sqr(displayYukawaElement(YD, 3, 3)),
    v1 = displayHvev() * cos(beta),
    v2 = displayHvev() * sin(beta),
    mb = forLoops.mb, mz = displayMzRun();

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);
  /*
  cn(1)  = - cos(2.0 * forLoops.thetaH);
  cn(2)  = - cn(1);
  cn(3)  = - cos(2.0 * beta);
  cn(4)  = - cn(3);
  dnu(1) = sqr(sin(forLoops.thetaH));
  dnu(2) = sqr(cos(forLoops.thetaH));
  dnu(3) = sqr(sinb);
  dnu(4) = sqr(cosb);
  dnd(1) = dnu(2);
  dnd(2) = dnu(1);
  dnd(3) = dnu(4);
  dnd(4) = dnu(3);
  */
  DoubleMatrix lsStopLStopLR(4, 2), lsStopLStop12(4, 2);
  DoubleMatrix lsStopRStopLR(4, 2), lsStopRStop12(4, 2);
  // Order (s1 s2 G A, L R)
  lsStopLStopLR(1, 1) = g * mz * guL * cosb / costhDrbar;
  lsStopLStopLR(1, 2) = ht * smu / root2;
  lsStopLStopLR(3, 2) = 1.0 / root2 *
    (smu * cosb * ht + displayTrilinear(UA, 3, 3) * sinb);
  lsStopLStopLR(4, 2) = -1.0 / root2 *
    (smu * sinb * ht - displayTrilinear(UA, 3, 3) * cosb);
  lsStopLStopLR(2, 1) = -g * mz * guL * sinb / costhDrbar
    + htsq * v2;
  lsStopLStopLR(2, 2) = displayTrilinear(UA, 3, 3) / root2;

  lsStopRStopLR(1, 1) = lsStopLStopLR(1, 2);
  lsStopRStopLR(1, 2) = g * mz * guR * cosb / costhDrbar;
  lsStopRStopLR(2, 1) = lsStopLStopLR(2, 2);
  lsStopRStopLR(2, 2) = -g * mz * guR * sinb / costhDrbar
    + htsq * v2;
  lsStopRStopLR(3, 1) = -lsStopLStopLR(3, 2);
  lsStopRStopLR(4, 1) = -lsStopLStopLR(4, 2);

  // Mix stops up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsStopLStopLR(i, 1);
    temp(2) = lsStopLStopLR(i, 2);
    temp2 = rot2d(thetat) * temp;
    lsStopLStop12(i, 1) = temp2(1);
    lsStopLStop12(i, 2) = temp2(2);
    temp(1) = lsStopRStopLR(i, 1);
    temp(2) = lsStopRStopLR(i, 2);
    temp2 = rot2d(thetat) * temp;
    lsStopRStop12(i, 1) = temp2(1);
    lsStopRStop12(i, 2) = temp2(2);
 }

  DoubleMatrix lHStopLStop12(lsStopLStop12), lHStopRStop12(lsStopRStop12);
  // Mix CP-even Higgses up
  for (i=1; i<=2; i++) { // i is the s1/s2 label
    temp(1) = lsStopLStop12(1, i);
    temp(2) = lsStopLStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopLStop12(1, i) = temp2(1);
    lHStopLStop12(2, i) = temp2(2);
    temp(1) = lsStopRStop12(1, i);
    temp(2) = lsStopRStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopRStop12(1, i) = temp2(1);
    lHStopRStop12(2, i) = temp2(2);
  }

  // Charged Higgs Feynman rules
  DoubleMatrix lChHstopLsbotLR(2, 2); // (G+ H+, L R) basis
  lChHstopLsbotLR(2, 1) = (g * displayMwRun() * sin(2.0 * beta)
    - htsq * v2 * cosb - hbsq * v1 * sinb) / root2;
  lChHstopLsbotLR(2, 2) = (smu * hb * cosb -
                           displayTrilinear(DA, 3, 3) * sinb);
  lChHstopLsbotLR(1, 1) = (-g * displayMwRun() * cos(2.0 * beta)
    - htsq * v2 * sinb + hbsq * v1 * cosb) / root2;
  lChHstopLsbotLR(1, 2) = hb * smu * sinb + displayTrilinear(DA, 3, 3) * cosb;

  DoubleMatrix lChHstopLsbot12(2, 2);
  temp(1) = lChHstopLsbotLR(1, 1);
  temp(2) = lChHstopLsbotLR(1, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopLsbot12(1, 1) = temp2(1);
  lChHstopLsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopLsbotLR(2, 1);
  temp(2) = lChHstopLsbotLR(2, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopLsbot12(2, 1) = temp2(1);
  lChHstopLsbot12(2, 2) = temp2(2);

  DoubleMatrix lChHstopRsbotLR(2, 2); // (G+ H+, L R) basis
  lChHstopRsbotLR(2, 1) = ht * smu * sinb - displayTrilinear(UA, 3, 3) * cosb;
  lChHstopRsbotLR(2, 2) = ht * hb * (- v1 * cosb - v2 * sinb) / root2;
  lChHstopRsbotLR(1, 1) = -ht * smu * cosb - displayTrilinear(UA, 3, 3) * sinb;
  DoubleMatrix lChHstopRsbot12(2, 2);
  temp(1) = lChHstopRsbotLR(1, 1);
  temp(2) = lChHstopRsbotLR(1, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopRsbot12(1, 1) = temp2(1);
  lChHstopRsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopRsbotLR(2, 1);
  temp(2) = lChHstopRsbotLR(2, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopRsbot12(2, 1) = temp2(1);
  lChHstopRsbot12(2, 2) = temp2(2);

  // Neutralino Feynman rules
  DoubleVector aPsi0TStopr(4), bPsi0TStopr(4), aPsi0TStopl(4),
    bPsi0TStopl(4);
  aPsi0TStopr(1) = - 4.0 * gp / (3.0 * root2);
  bPsi0TStopl(1) = gp / (3.0 * root2);
  bPsi0TStopl(2) = g / root2;
  aPsi0TStopl(4) = ht;
  bPsi0TStopr(4) = ht;

  ComplexVector aChi0TStopl(4), bChi0TStopl(4), aChi0TStopr(4),
    bChi0TStopr(4);
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  DoubleVector gChi0TopStopLL(4), fChi0TopStopLL(4);
  DoubleVector gChi0TopStopLR(4), fChi0TopStopLR(4);
  DoubleVector gChi0TopStopRR(4), fChi0TopStopRR(4);
  for (i=1; i<=4; i++) {
    fChi0TopStopLL(i) = (aChi0TStopl(i) * aChi0TStopl(i).conj() +
      bChi0TStopl(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLL(i) = (bChi0TStopl(i).conj() * aChi0TStopl(i) +
      bChi0TStopl(i) * aChi0TStopl(i).conj()).real();
    fChi0TopStopRR(i) = (aChi0TStopr(i) * aChi0TStopr(i).conj() +
      bChi0TStopr(i) * bChi0TStopr(i).conj()).real();
    gChi0TopStopRR(i) = (bChi0TStopr(i).conj() * aChi0TStopr(i) +
      bChi0TStopr(i) * aChi0TStopr(i).conj()).real();
    fChi0TopStopLR(i) = (aChi0TStopr(i) * aChi0TStopl(i).conj() +
      bChi0TStopr(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLR(i) = (bChi0TStopl(i).conj() * aChi0TStopr(i) +
      bChi0TStopr(i) * aChi0TStopl(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicBStopl(2), bPsicBStopr(2), aPsicBStopl(2),
    aPsicBStopr(2);

  aPsicBStopl(1) = g;
  aPsicBStopr(2) = -ht;
  bPsicBStopl(2) = -hb;

  DoubleVector aPsicCSbotl(2);
  ComplexVector aChicBStopr(2), aChicBStopl(2), bChicBStopl(2),
      bChicBStopr(2);
  ComplexMatrix aChBStop(2, 2), bChBStop(2, 2);

  aChicBStopl = v.complexConjugate() * aPsicBStopl;
  bChicBStopl = u * bPsicBStopl;
  aChicBStopr = v.complexConjugate() * aPsicBStopr;
  bChicBStopr = u * bPsicBStopr;

  DoubleVector fChBStopLL(2), gChBStopLL(2) ;
  DoubleVector fChBStopLR(2), gChBStopLR(2);
  DoubleVector fChBStopRR(2), gChBStopRR(2);
  for (i=1; i<=2; i++) {
    fChBStopLL(i) = (aChicBStopl(i).conj() * aChicBStopl(i) +
                      bChicBStopl(i).conj() * bChicBStopl(i)).real();
    gChBStopLL(i) = (bChicBStopl(i).conj() * aChicBStopl(i) +
                      aChicBStopl(i).conj() * bChicBStopl(i)).real();
    fChBStopLR(i) = (aChicBStopl(i).conj() * aChicBStopr(i) +
                      bChicBStopl(i).conj() * bChicBStopr(i)).real();
    gChBStopLR(i) = (bChicBStopl(i).conj() * aChicBStopr(i) +
                      aChicBStopl(i).conj() * bChicBStopr(i)).real();
    fChBStopRR(i) = (aChicBStopr(i).conj() * aChicBStopr(i) +
                      bChicBStopr(i).conj() * bChicBStopr(i)).real();
    gChBStopRR(i) = (bChicBStopr(i).conj() * aChicBStopr(i) +
                      aChicBStopr(i).conj() * bChicBStopr(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2),
    higgs(2, 2), chargino(2, 2), neutralino(2, 2);
  double a0t1 = a0(mstop(1), q), a0t2 = a0(mstop(2), q);
  double ft1 = ffn(p, mstop(1), 0.0, q), ft2 = ffn(p, mstop(2), 0.0, q);
  double ggt = gfn(p, mg, mt, q);
  strong(1, 1) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(ct) * (ft1 + a0t1) + sqr(st) * (ft2 + a0t2));
  strong(2, 2) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(st) * (ft1 + a0t1) + sqr(ct) * (ft2 + a0t2));
  strong(1, 2) = 4.0 * g3sq / 3.0 *
    (4.0 * mg * mt * b0(p, mg, mt, q) +
     st * ct * (ft1 - a0t1 - ft2 + a0t2));

  stop(1, 1) = htsq * (sqr(st) * a0t1 + sqr(ct) * a0t2);
  stop(2, 2) = htsq * (sqr(ct) * a0t1 + sqr(st) * a0t2);
  stop(1, 2) = htsq * ct * st * 3.0 * (a0t1 - a0t2);

  sbottom(1, 1) =
    hbsq * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q));
  sbottom(2, 2) =
    htsq * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q));

  for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) +
      0.5 * (htsq * dnu(i) - sqr(g) * guL * 0.5 / sqr(costhDrbar) * cn(i))
      * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) +
      0.5 * (htsq * dnu(i) - sqr(g) * guR * 0.5 / sqr(costhDrbar) * cn(i))
      * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) +
      (hbsq * dnu(i) + sqr(g) * (guL * 0.5 / costhDrbar2 - 0.5) * cn(i))* a0p;
    higgs(2, 2) = higgs(2, 2) +
      (htsq * dnd(i) + sqr(g) * guR * 0.5 / costhDrbar2 * cn(i))* a0p;
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsm(i), mstop(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHStopLStop12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) +
        lHStopLStop12(i, j) * lHStopRStop12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lHStopRStop12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), msbot(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lChHstopLsbot12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) +
        lChHstopLsbot12(i, j) * lChHstopRsbot12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lChHstopRsbot12(i, j)) * b0p;
    }

  // EW bosons
  higgs(1, 1) = higgs(1, 1) +
    4.0 * sqr(g) / costhDrbar2 * sqr(guL) * a0(mz, q) +
    2.0 * sqr(g) * a0(displayMwRun(), q) + sqr(2.0 / 3.0 * e) *
    (sqr(ct) * ffn(p, mstop(1), 0.0, q) + sqr(st) * ffn(p, mstop(2), 0.0, q))
    + sqr(g * guL / costhDrbar) *
    (sqr(ct) * ffn(p, mstop(1), mz, q) + sqr(st) * ffn(p, mstop(2), mz, q)) +
    sqr(g) * 0.5 * (sqr(cb) * ffn(p, msbot1, displayMwRun(), q) + sqr(sb) *
                    ffn(p, msbot2, displayMwRun(), q)) +
    sqr(g) * 0.25 *
    (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q) + 2.0 *
     (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q))) +
    sqr(g) * 0.5 *
    (1.5 * a0(forLoops.mu(1, 1), q) +
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) +
            sqr(st) * a0(forLoops.mu(2, 2), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
            sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) +
            a0(forLoops.msnu(2), q) +
            a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(yuL) *
    (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q)) +
    sqr(gp) * 0.25 * yuL *
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) +
                  a0(forLoops.mu(1, 2), q) +
                  sqr(ct) * a0(forLoops.mu(1, 3), q) +
                  sqr(st) * a0(forLoops.mu(2, 2), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) +
                  a0(forLoops.mu(2, 2), q) +
                  sqr(st) * a0(forLoops.mu(1, 3), q) +
                  sqr(ct) * a0(forLoops.mu(2, 2), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) +
                  a0(forLoops.md(1, 2), q) +
                  sqr(cb) * a0(forLoops.md(1, 3), q) +
                  sqr(sb) * a0(forLoops.md(2, 2), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) +
                  a0(forLoops.md(2, 2), q) +
                  sqr(sb) * a0(forLoops.md(1, 3), q) +
                  sqr(cb) * a0(forLoops.md(2, 2), q)) +
     yeL * (a0(forLoops.me(1, 1), q) +
            a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 2), q)) +
     yeR * (a0(forLoops.me(2, 1), q) +
            a0(forLoops.me(2, 2), q) +
            sqr(stau) * a0(forLoops.me(1, 3), q) +
            sqr(ctau) * a0(forLoops.me(2, 2), q)) +
     ynuL * (a0(forLoops.msnu(1), q) +
             a0(forLoops.msnu(2), q) +
             a0(forLoops.msnu(3), q)));

  higgs(2, 2) = higgs(2, 2) +
    4.0 * sqr(g) / costhDrbar2 * sqr(guR) * a0(mz, q) +
    sqr(2.0 / 3.0 * e) *
    (sqr(st) * ffn(p, mstop(1), 0.0, q) + sqr(ct) * ffn(p, mstop(2), 0.0, q))
    + sqr(g * guR / costhDrbar) *
    (sqr(st) * ffn(p, mstop(1), mz, q) + sqr(ct) * ffn(p, mstop(2), mz, q)) +
    sqr(gp) * 0.25 * sqr(yuR) *
    (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q)) +
    sqr(gp) * 0.25 * yuR *
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) +
                  a0(forLoops.mu(1, 2), q) +
                  sqr(ct) * a0(forLoops.mu(1, 3), q) +
                  sqr(st) * a0(forLoops.mu(2, 2), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) +
                  a0(forLoops.mu(2, 2), q) +
                  sqr(st) * a0(forLoops.mu(1, 3), q) +
                  sqr(ct) * a0(forLoops.mu(2, 2), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) +
                  a0(forLoops.md(1, 2), q) +
                  sqr(cb) * a0(forLoops.md(1, 3), q) +
                  sqr(sb) * a0(forLoops.md(2, 2), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) +
                  a0(forLoops.md(2, 2), q) +
                  sqr(sb) * a0(forLoops.md(1, 3), q) +
                  sqr(cb) * a0(forLoops.md(2, 2), q)) +
     yeL * (a0(forLoops.me(1, 1), q) +
            a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 2), q)) +
     yeR * (a0(forLoops.me(2, 1), q) +
            a0(forLoops.me(2, 2), q) +
            sqr(stau) * a0(forLoops.me(1, 3), q) +
            sqr(ctau) * a0(forLoops.me(2, 2), q)) +
     ynuL * (a0(forLoops.msnu(1), q) +
             a0(forLoops.msnu(2), q) +
             a0(forLoops.msnu(3), q)));

  higgs(1, 2) = higgs(1, 2) +
    sqr(gp) * yuL * yuR * st * ct * (a0(mstop(1), q) - a0(mstop(2), q)) +
    sqr(2.0 / 3.0 * e) * st * ct *
    (ffn(p, mstop(1), 0.0, q) - ffn(p, mstop(2), 0.0, q)) -
    sqr(g) * 0.25 * guL * guR * st * ct *
    (ffn(p, mstop(1), mz, q) - ffn(p, mstop(2), mz, q));

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), mb, q);
    double two = mch(i) * mb * b0(p, mch(i), mb, q) * 2.0;
    chargino(1, 1) = chargino(1, 1) + fChBStopLL(i) * one -
      gChBStopLL(i) * two;
    chargino(1, 2) = chargino(1, 2) + fChBStopLR(i) * one -
      gChBStopLR(i) * two;
    chargino(2, 2) = chargino(2, 2) + fChBStopRR(i) * one -
      gChBStopRR(i) * two;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p, mneut(i), mt, q);
    double two = 2.0 * mneut(i) * mt * b0(p, mneut(i), mt, q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0TopStopLL(i) * one - gChi0TopStopLL(i) * two;
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0TopStopRR(i) * one - gChi0TopStopRR(i) * two;
    neutralino(1, 2) = neutralino(1, 2) +
      fChi0TopStopLR(i) * one - gChi0TopStopLR(i) * two;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (strong + stop + sbottom + higgs + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);

  mass = mass - piSq;
}

void MssmSoftsusy::assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc)
  const {
  drBarPars f(displayDrBarPars());

  higgsm(1) = f.mhiggs(3);
  higgsm(2) = f.mhiggs(1);
  higgsm(3) = displayMzRun();
  higgsm(4) = f.mhiggs(2);
  higgsc(1) = displayMwRun();
  higgsc(2) = f.mhiggs(4);
}

void MssmSoftsusy::assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc,
                               DoubleVector & dnu, DoubleVector & dnd,
                               DoubleVector & cn, double beta) const {
  double sinb = sin(beta), cosb = cos(beta);

  assignHiggs(higgsm, higgsc);

  cn(1)  = - cos(2.0 * forLoops.thetaH);
  cn(2)  = - cn(1);
  cn(3)  = - cos(2.0 * beta);
  cn(4)  = - cn(3);
  dnu(1) = sqr(sin(displayDrBarPars().thetaH));
  dnu(2) = sqr(cos(displayDrBarPars().thetaH));
  dnu(3) = sqr(sinb);
  dnu(4) = sqr(cosb);
  dnd(1) = dnu(2);
  dnd(2) = dnu(1);
  dnd(3) = dnu(4);
  dnd(4) = dnu(3);
}

// 16.09.05 checked.
void MssmSoftsusy::addSlepCorrection(DoubleMatrix & mass, int family) {

// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) {
    if (PRINTOUT > 2) cout << " slepton tachyonic ";
    flagTachyon(true);
    return;
  }

  // one-loop correction matrix
  DoubleMatrix piSq(2, 2); // Self-energy matrix

  // brevity
  double    mz      = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    e       = g * sinthDrbar;
  double    beta    = atan(displayTanb());

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double    ctau    = cos(thetatau);
  double    stau    = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msel(2);
  msel(1)           = forLoops.me(1, family);
  msel(2)           = forLoops.me(2, family);
  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);

  double q = displayMu(), sinb = sin(beta), cosb = cos(beta);
  double    p       = sqrt(msel(1) * msel(2));

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2), dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  // stau now stands for generation 1 or 2
  DoubleMatrix lsStauLStauLR(4, 2), lsStauLStau12(4, 2);
  DoubleMatrix lsStauRStauLR(4, 2), lsStauRStau12(4, 2);
  // Order (s1 s2 G A, L R)
  lsStauLStauLR(1, 1) = g * mz * geL * cosb / costhDrbar;
  lsStauLStauLR(2, 1) = -g * mz * geL * sinb / costhDrbar;

  lsStauRStauLR(1, 2) = g * mz * geR * cosb / costhDrbar;
  lsStauRStauLR(2, 2) = -g * mz * geR * sinb / costhDrbar;
  lsStauRStauLR(3, 1) = -lsStauLStauLR(3, 2);
  lsStauRStauLR(4, 1) = -lsStauLStauLR(4, 2);

  DoubleMatrix lHStauLStau12(lsStauLStauLR), lHStauRStau12(lsStauRStauLR);
  DoubleVector temp(2), temp2(2);
  // Mix CP-even Higgses up
  int i; for (i=1; i<=2; i++) { // i is the s1/s2 label
    temp(1) = lsStauLStauLR(1, i);
    temp(2) = lsStauLStauLR(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauLStau12(1, i) = temp2(1);
    lHStauLStau12(2, i) = temp2(2);
    temp(1) = lsStauRStauLR(1, i);
    temp(2) = lsStauRStauLR(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauRStau12(1, i) = temp2(1);
    lHStauRStau12(2, i) = temp2(2);
  }

  // Charged Higgs Feynman rules
  DoubleMatrix lChHstauLsnu12(2, 2); // (G+ H+, L R) basis
  lChHstauLsnu12(1, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHstauLsnu12(2, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;

  // Neutralino Feynman rules
  DoubleVector aPsi0TauStaur(4), bPsi0TauStaur(4), aPsi0TauStaul(4),
    bPsi0TauStaul(4);
  aPsi0TauStaur(1) = gp * root2;
  bPsi0TauStaul(1) = -gp / root2;
  bPsi0TauStaul(2) = -g / root2;

  ComplexVector aChi0TauStaul(4), bChi0TauStaul(4), aChi0TauStaur(4),
    bChi0TauStaur(4);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0TauStaul = n.complexConjugate() * aPsi0TauStaul;
  bChi0TauStaul = n * bPsi0TauStaul;
  aChi0TauStaur = n.complexConjugate() * aPsi0TauStaur;
  bChi0TauStaur = n * bPsi0TauStaur;

  DoubleVector gChi0TauStauLL(4), fChi0TauStauLL(4);
  DoubleVector gChi0TauStauRR(4), fChi0TauStauRR(4);
  for (i=1; i<=4; i++) {
    fChi0TauStauLL(i) = (aChi0TauStaul(i) * aChi0TauStaul(i).conj() +
      bChi0TauStaul(i) * bChi0TauStaul(i).conj()).real();
    gChi0TauStauLL(i) = (bChi0TauStaul(i).conj() * aChi0TauStaul(i) +
      bChi0TauStaul(i) * aChi0TauStaul(i).conj()).real();
    fChi0TauStauRR(i) = (aChi0TauStaur(i) * aChi0TauStaur(i).conj() +
      bChi0TauStaur(i) * bChi0TauStaur(i).conj()).real();
    gChi0TauStauRR(i) = (bChi0TauStaur(i).conj() * aChi0TauStaur(i) +
      bChi0TauStaur(i) * aChi0TauStaur(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicNuStaul(2);//, bPsicNuStaur(2);

  bPsicNuStaul(1) = g;

  ComplexVector bChicNuStaul(2);
  ComplexMatrix bChNuStau(2, 2);

  bChicNuStaul = u * bPsicNuStaul;

  DoubleVector fChNuStauLL(2) ;
  DoubleVector fChNuStauRR(2);
  for (i=1; i<=2; i++) {
    fChNuStauLL(i) = (bChicNuStaul(i).conj() * bChicNuStaul(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  DoubleMatrix higgs(2, 2), chargino(2, 2), neutralino(2, 2);
  for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) +
      0.5 * ( - sqr(g) * geL * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) +
      0.5 * (- sqr(g) * geR * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }

  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) + a0p *
          (sqr(g) * (geL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) = higgs(2, 2) + a0p *
      (sqr(g) * geR / (2.0 * sqr(costhDrbar)) * cn(i));
          (sqr(g) * (geL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
  }

  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsm(i), msel(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHStauLStau12(i, j)) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lHStauRStau12(i, j)) * b0p;
    }


  for(i=1; i<=2; i++) {
      double b0p = b0(p, msnu(family), higgsc(i), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lChHstauLsnu12(i, 1)) * b0p;
    }

  DoubleMatrix electroweak(2, 2);
  // line by line....
  higgs(1, 1) = higgs(1, 1) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e) *
    ffn(p, msel(1), 0., q);
  higgs(1, 1) = higgs(1, 1) +
    sqr(g) / sqr(costhDrbar) * sqr(geL) *
    ffn(p, msel(family), mz, q) +
    sqr(g) * 0.5 * ffn(p, msnu(family), mw, q);

  electroweak(1, 1) = electroweak(1, 1) +
    sqr(g) * 0.25 * (a0(msel(1), q) + 2.0 * a0(msnu(family), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) *
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(yeL) * a0(msel(1), q);
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * yeL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));


  electroweak(2, 2) = electroweak(2, 2) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geR) * a0(mz, q) +
    sqr(e) * ffn(p, msel(2), 0., q);
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(g) / sqr(costhDrbar) * sqr(geR) * ffn(p, msel(2), mz, q);
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(gp) * 0.25 * sqr(yeR) * a0(msel(2), q);
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * yeR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), 0., q);
    chargino(1, 1) = chargino(1, 1) + fChNuStauLL(i) * one;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p, mneut(i), 0., q);
    neutralino(1, 1) = neutralino(1, 1) + fChi0TauStauLL(i) * one;
    neutralino(2, 2) = neutralino(2, 2) + fChi0TauStauRR(i) * one;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (higgs + chargino + neutralino + electroweak);

  mass = mass - piSq;
}

void MssmSoftsusy::addStauCorrection(DoubleMatrix & mass, double mtau) {

// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) {
    if (PRINTOUT > 2) cout << " sbottom tachyonic ";
    flagTachyon(true);
    return;
  }

  // one-loop correction matrix
  DoubleMatrix piSq(2, 2); // Self-energy matrix

  // brevity
  double    mz      = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    e       = g * sinthDrbar;
  double    beta    = atan(displayTanb());

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double    ctau    = cos(thetatau);
  double    stau    = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);

  double    smu     = -displaySusyMu();
  double q = displayMu(),
    htau = displayYukawaElement(YE, 3, 3),
    sinb = sin(beta), cosb = cos(beta),
    htausq = sqr(htau);
  double    p       = sqrt(mstau(1) * mstau(2));

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleMatrix lsStauLStauLR(4, 2), lsStauLStau12(4, 2);
  DoubleMatrix lsStauRStauLR(4, 2), lsStauRStau12(4, 2);
  // Order (s1 s2 G A, L R)
  lsStauLStauLR(1, 1) = g * mz * geL * cosb / costhDrbar +
    root2 * htau * mtau;
  lsStauLStauLR(1, 2) = displayTrilinear(EA, 3, 3) / root2;
  lsStauLStauLR(3, 2) = -1.0 / root2 *
    (smu * sinb * htau + displayTrilinear(EA, 3, 3) * cosb);
  lsStauLStauLR(4, 2) = -1.0 / root2 *
    (smu * cosb * htau - displayTrilinear(EA, 3, 3) * sinb);
  lsStauLStauLR(2, 1) = -g * mz * geL * sinb / costhDrbar;
  lsStauLStauLR(2, 2) = htau * smu / root2;

  lsStauRStauLR(1, 1) = lsStauLStauLR(1, 2);
  lsStauRStauLR(1, 2) = g * mz * geR * cosb / costhDrbar +
    root2 * htau * mtau;
  lsStauRStauLR(2, 1) = lsStauLStauLR(2, 2);
  lsStauRStauLR(2, 2) = -g * mz * geR * sinb / costhDrbar;
  lsStauRStauLR(3, 1) = -lsStauLStauLR(3, 2);
  lsStauRStauLR(4, 1) = -lsStauLStauLR(4, 2);

  // Mix staus up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsStauLStauLR(i, 1);
    temp(2) = lsStauLStauLR(i, 2);
    temp2 = rot2d(thetatau) * temp;
    lsStauLStau12(i, 1) = temp2(1);
    lsStauLStau12(i, 2) = temp2(2);
    temp(1) = lsStauRStauLR(i, 1);
    temp(2) = lsStauRStauLR(i, 2);
    temp2 = rot2d(thetatau) * temp;
    lsStauRStau12(i, 1) = temp2(1);
    lsStauRStau12(i, 2) = temp2(2);
 }

  DoubleMatrix lHStauLStau12(lsStauLStau12), lHStauRStau12(lsStauRStau12);
  // Mix CP-even Higgses up
  for (i=1; i<=2; i++) { // i is the s1/s2 label
    temp(1) = lsStauLStau12(1, i);
    temp(2) = lsStauLStau12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauLStau12(1, i) = temp2(1);
    lHStauLStau12(2, i) = temp2(2);
    temp(1) = lsStauRStau12(1, i);
    temp(2) = lsStauRStau12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauRStau12(1, i) = temp2(1);
    lHStauRStau12(2, i) = temp2(2);
  }

  // Charged Higgs Feynman rules
  DoubleMatrix lChHstauLsnu12(2, 2); // (G+ H+, L R) basis
  lChHstauLsnu12(1, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2
    + htau * mtau * cosb;
  lChHstauLsnu12(2, 1) = g * displayMwRun() * sin(2.0 * beta) / root2
    - htau * mtau * sinb;

  DoubleMatrix lChHstauRsnu12(2, 2); // (G+ H+, L R) basis
  lChHstauRsnu12(1, 1) = htau * smu * sinb + displayTrilinear(EA, 3, 3) * cosb;
  lChHstauRsnu12(2, 1) = htau * smu * cosb - displayTrilinear(EA, 3, 3) * sinb;

  // Neutralino Feynman rules
  DoubleVector aPsi0TauStaur(4), bPsi0TauStaur(4), aPsi0TauStaul(4),
    bPsi0TauStaul(4);
  aPsi0TauStaur(1) = gp * root2;
  bPsi0TauStaul(1) = -gp / root2;
  bPsi0TauStaul(2) = -g / root2;
  aPsi0TauStaul(3) = htau;
  bPsi0TauStaur(3) = htau;

  ComplexVector aChi0TauStaul(4), bChi0TauStaul(4), aChi0TauStaur(4),
    bChi0TauStaur(4);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0TauStaul = n.complexConjugate() * aPsi0TauStaul;
  bChi0TauStaul = n * bPsi0TauStaul;
  aChi0TauStaur = n.complexConjugate() * aPsi0TauStaur;
  bChi0TauStaur = n * bPsi0TauStaur;

  DoubleVector gChi0TauStauLL(4), fChi0TauStauLL(4);
  DoubleVector gChi0TauStauLR(4), fChi0TauStauLR(4);
  DoubleVector gChi0TauStauRR(4), fChi0TauStauRR(4);
  for (i=1; i<=4; i++) {
    fChi0TauStauLL(i) = (aChi0TauStaul(i) * aChi0TauStaul(i).conj() +
      bChi0TauStaul(i) * bChi0TauStaul(i).conj()).real();
    gChi0TauStauLL(i) = (bChi0TauStaul(i).conj() * aChi0TauStaul(i) +
      bChi0TauStaul(i) * aChi0TauStaul(i).conj()).real();
    fChi0TauStauRR(i) = (aChi0TauStaur(i) * aChi0TauStaur(i).conj() +
      bChi0TauStaur(i) * bChi0TauStaur(i).conj()).real();
    gChi0TauStauRR(i) = (bChi0TauStaur(i).conj() * aChi0TauStaur(i) +
      bChi0TauStaur(i) * aChi0TauStaur(i).conj()).real();
    fChi0TauStauLR(i) = (aChi0TauStaur(i) * aChi0TauStaul(i).conj() +
      bChi0TauStaur(i) * bChi0TauStaul(i).conj()).real();
    gChi0TauStauLR(i) = (bChi0TauStaul(i).conj() * aChi0TauStaur(i) +
      bChi0TauStaur(i) * aChi0TauStaul(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicNuStaul(2), bPsicNuStaur(2);

  bPsicNuStaul(1) = g;
  bPsicNuStaur(2) = -htau;

  ComplexVector bChicNuStaul(2), bChicNuStaur(2);
  ComplexMatrix bChNuStau(2, 2);

  bChicNuStaul = u * bPsicNuStaul;
  bChicNuStaur = u * bPsicNuStaur;

  DoubleVector fChNuStauLL(2) ;
  DoubleVector fChNuStauLR(2);
  DoubleVector fChNuStauRR(2);
  for (i=1; i<=2; i++) {
    fChNuStauLL(i) = (bChicNuStaul(i).conj() * bChicNuStaul(i)).real();
    fChNuStauLR(i) = (bChicNuStaul(i).conj() * bChicNuStaur(i)).real();
    fChNuStauRR(i) = (bChicNuStaur(i).conj() * bChicNuStaur(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  DoubleMatrix stop(2, 2), sbottom(2, 2),
    higgs(2, 2), chargino(2, 2), neutralino(2, 2);
  double a0t1 = a0(mstau(1), q), a0t2 = a0(mstau(2), q);

  sbottom(2, 2) = htausq * a0(msnu(3), q);

  // start
  stop(1, 1) = htausq * (sqr(stau) * a0t1 + sqr(ctau) * a0t2);
  stop(2, 2) = htausq * (sqr(ctau) * a0t1 + sqr(stau) * a0t2);
  stop(1, 2) = htausq * ctau * stau * (a0t1 - a0t2);

  for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) +
      0.5 * (htausq * dnd(i) - sqr(g) * geL * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) +
      0.5 * (htausq * dnd(i) - sqr(g) * geR * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) + a0p *
          (sqr(g) * (geL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) = higgs(2, 2) + a0p *
      (htausq * dnu(i) + sqr(g) * geR / (2.0 * sqr(costhDrbar)) * cn(i));
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsm(i), mstau(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHStauLStau12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) +
        lHStauLStau12(i, j) * lHStauRStau12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lHStauRStau12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++) {
    double b0p = b0(p, higgsc(i), msnu(3), q);
    higgs(1, 1) = higgs(1, 1) + sqr(lChHstauLsnu12(i, 1)) * b0p;
    higgs(1, 2) = higgs(1, 2) +
      lChHstauLsnu12(i, 1) * lChHstauRsnu12(i, 1) * b0p;
    higgs(2, 2) = higgs(2, 2) + sqr(lChHstauRsnu12(i, 1)) * b0p;
  }

  DoubleMatrix electroweak(2, 2);
  // line by line....
  electroweak(1, 1) = electroweak(1, 1) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e) *
    (sqr(ctau) * ffn(p, mstau(1), 0., q) +
     sqr(stau) * ffn(p, mstau(2), 0., q));
  electroweak(1, 1) = electroweak(1, 1) +
    sqr(g) / sqr(costhDrbar) * sqr(geL) *
    (sqr(ctau) * ffn(p, mstau(1), mz, q) +
     sqr(stau) * ffn(p, mstau(2), mz, q)) +
    sqr(g) * 0.5 * ffn(p, msnu(3), mw, q);
  electroweak(1, 1) = electroweak(1, 1) +
    sqr(g) * 0.25 * (sqr(ctau) * a0(mstau(1), q) +
                     sqr(stau) * a0(mstau(2), q) +
                     2.0 * a0(msnu(3), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) *
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(yeL) *
    (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * yeL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));


  electroweak(2, 2) = electroweak(2, 2) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geR) * a0(mz, q) +
    sqr(e) *
    (sqr(stau) * ffn(p, mstau(1), 0., q) +
     sqr(ctau) * ffn(p, mstau(2), 0., q));
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(g) / sqr(costhDrbar) * sqr(geR) *
    (sqr(stau) * ffn(p, mstau(1), mz, q) +
     sqr(ctau) * ffn(p, mstau(2), mz, q));
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(gp) * 0.25 * sqr(yeR) *
    (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q));
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * yeR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  electroweak(1, 2) = electroweak(1, 2) +
    sqr(gp) * stau * ctau * 0.25 * yeL * yeR *
    (a0(mstau(1), q) - a0(mstau(2), q)) +
    sqr(e) * stau * ctau *
    (ffn(p, mstau(1), 0., q) - ffn(p, mstau(2), 0., q)) -
    sqr(g) / sqr(costhDrbar) * geL * geR * stau * ctau *
    (ffn(p, mstau(1), mz, q) - ffn(p, mstau(2), mz, q));

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), mtau, q);
    chargino(1, 1) = chargino(1, 1) + fChNuStauLL(i) * one;
    chargino(1, 2) = chargino(1, 2) + fChNuStauLR(i) * one;
    chargino(2, 2) = chargino(2, 2) + fChNuStauRR(i) * one;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p, mneut(i), mtau, q);
    neutralino(1, 1) = neutralino(1, 1) + fChi0TauStauLL(i) * one;
    neutralino(2, 2) = neutralino(2, 2) + fChi0TauStauRR(i) * one;
    neutralino(1, 2) = neutralino(1, 2) + fChi0TauStauLR(i) * one;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (stop + sbottom + higgs + chargino + neutralino + electroweak);

  piSq(2, 1) = piSq(1, 2);

  mass = mass - piSq;
}

void MssmSoftsusy::addSdownCorrection(DoubleMatrix & mass, int family) {

// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) {
    if (PRINTOUT > 2) cout << " sdown tachyonic ";
    flagTachyon(true);
    return;
  }

  // one-loop correction matrix
  DoubleMatrix piSq(2, 2); // Self-energy matrix

  // brevity
  double    mz      = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    e       = g * sinthDrbar;
  double    beta    = atan(displayTanb());

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  DoubleVector msd(2);
  msd(1)  = forLoops.md(1, family);
  msd(2)  = forLoops.md(2, family);
  DoubleVector msup(2);
  msup(1) = forLoops.mu(1, family);
  msup(2) = forLoops.mu(2, family);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double    ctau    = cos(thetatau);
  double    stau    = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)          = forLoops.msnu(1);
  msnu(2)          = forLoops.msnu(2);
  msnu(3)          = forLoops.msnu(3);

  double    mg      = forLoops.mGluino;
  double q = displayMu(), g3sq = sqr(displayGaugeCoupling(3)),
    sinb = sin(beta), cosb = cos(beta);

  double    p1 = msd(1), p2 = msd(2);
  p1 = p2 = sqrt(msd(1) * msd(2)); // DEBUG

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleMatrix lsSbotLSbotLR(4, 2), lsSbotLSbot12(4, 2);
  DoubleMatrix lsSbotRSbotLR(4, 2), lsSbotRSbot12(4, 2);
  // Order (s1 s2 G A, L R)
  lsSbotLSbotLR(1, 1) = g * mz * gdL * cosb / costhDrbar;
  lsSbotLSbotLR(2, 1) = -g * mz * gdL * sinb / costhDrbar;

  lsSbotRSbotLR(1, 1) = lsSbotLSbotLR(1, 2);
  lsSbotRSbotLR(1, 2) = g * mz * gdR * cosb / costhDrbar;
  lsSbotRSbotLR(2, 1) = lsSbotLSbotLR(2, 2);
  lsSbotRSbotLR(2, 2) = -g * mz * gdR * sinb / costhDrbar;
  lsSbotRSbotLR(3, 1) = -lsSbotLSbotLR(3, 2);
  lsSbotRSbotLR(4, 1) = -lsSbotLSbotLR(4, 2);

  // Mix sbots up
  int i;
  DoubleVector temp(2), temp2(2);
  DoubleMatrix lHSbotLSbot12(lsSbotLSbotLR), lHSbotRSbot12(lsSbotRSbotLR);
  // Mix CP-even Higgses up
  for (i=1; i<=2; i++) { // i is the L/R label
    temp(1) = lsSbotLSbot12(1, i);
    temp(2) = lsSbotLSbot12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotLSbot12(1, i) = temp2(1);
    lHSbotLSbot12(2, i) = temp2(2);
    temp(1) = lsSbotRSbot12(1, i);
    temp(2) = lsSbotRSbot12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotRSbot12(1, i) = temp2(1);
    lHSbotRSbot12(2, i) = temp2(2);
  }

  // Charged Higgs Feynman rules
  DoubleMatrix lChHsbotLstopLR(2, 2); // (G+ H+, L R) basis
  lChHsbotLstopLR(1, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHsbotLstopLR(2, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;
  DoubleMatrix lChHsbotLstop12(2, 2);
  temp(1) = lChHsbotLstopLR(1, 1);
  temp(2) = lChHsbotLstopLR(1, 2);
  temp2 = temp;
  lChHsbotLstop12(1, 1) = temp2(1);
  lChHsbotLstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotLstopLR(2, 1);
  temp(2) = lChHsbotLstopLR(2, 2);
  temp2 = temp;
  lChHsbotLstop12(2, 1) = temp2(1);
  lChHsbotLstop12(2, 2) = temp2(2);
  // (none for sdownL since they are all Yukawa suppressed)

  // Neutralino Feynman rules
  DoubleVector aPsi0BSbotr(4), bPsi0BSbotr(4), aPsi0BSbotl(4),
    bPsi0BSbotl(4);
  aPsi0BSbotr(1) = gp * ydR / root2;
  bPsi0BSbotl(1) = gp * ydL / root2;
  bPsi0BSbotl(2) = -0.5 * g * root2;

  ComplexVector aChi0BSbotl(4), bChi0BSbotl(4), aChi0BSbotr(4),
    bChi0BSbotr(4);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0BSbotl = n.complexConjugate() * aPsi0BSbotl;
  bChi0BSbotl = n * bPsi0BSbotl;
  aChi0BSbotr = n.complexConjugate() * aPsi0BSbotr;
  bChi0BSbotr = n * bPsi0BSbotr;

  DoubleVector gChi0BotSbotLL(4), fChi0BotSbotLL(4);
  DoubleVector gChi0BotSbotLR(4), fChi0BotSbotLR(4);
  DoubleVector gChi0BotSbotRR(4), fChi0BotSbotRR(4);
  for (i=1; i<=4; i++) {
    fChi0BotSbotLL(i) = (aChi0BSbotl(i) * aChi0BSbotl(i).conj() +
      bChi0BSbotl(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLL(i) = (bChi0BSbotl(i).conj() * aChi0BSbotl(i) +
      bChi0BSbotl(i) * aChi0BSbotl(i).conj()).real();
    fChi0BotSbotRR(i) = (aChi0BSbotr(i) * aChi0BSbotr(i).conj() +
      bChi0BSbotr(i) * bChi0BSbotr(i).conj()).real();
    gChi0BotSbotRR(i) = (bChi0BSbotr(i).conj() * aChi0BSbotr(i) +
      bChi0BSbotr(i) * aChi0BSbotr(i).conj()).real();
    fChi0BotSbotLR(i) = (aChi0BSbotr(i) * aChi0BSbotl(i).conj() +
      bChi0BSbotr(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLR(i) = (bChi0BSbotl(i).conj() * aChi0BSbotr(i) +
      bChi0BSbotr(i) * aChi0BSbotl(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicTSbotl(2), bPsicTSbotr(2), aPsicTSbotl(2),
    aPsicTSbotr(2);

  bPsicTSbotl(1) = g;

  ComplexVector aChicTSbotr(2), aChicTSbotl(2), bChicTSbotl(2),
      bChicTSbotr(2);
  ComplexMatrix aChTSbot(2, 2), bChTSbot(2, 2);

  aChicTSbotl = v.complexConjugate() * aPsicTSbotl;
  bChicTSbotl = u * bPsicTSbotl;
  aChicTSbotr = v.complexConjugate() * aPsicTSbotr;
  bChicTSbotr = u * bPsicTSbotr;

  DoubleVector fChTSbotLL(2), gChTSbotLL(2) ;
  DoubleVector fChTSbotLR(2), gChTSbotLR(2);
  DoubleVector fChTSbotRR(2), gChTSbotRR(2);
  for (i=1; i<=2; i++) {
    fChTSbotLL(i) = (aChicTSbotl(i).conj() * aChicTSbotl(i) +
                      bChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    gChTSbotLL(i) = (bChicTSbotl(i).conj() * aChicTSbotl(i) +
                      aChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    fChTSbotLR(i) = (aChicTSbotl(i).conj() * aChicTSbotr(i) +
                      bChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    gChTSbotLR(i) = (bChicTSbotl(i).conj() * aChicTSbotr(i) +
                      aChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    fChTSbotRR(i) = (aChicTSbotr(i).conj() * aChicTSbotr(i) +
                      bChicTSbotr(i).conj() * bChicTSbotr(i)).real();
    gChTSbotRR(i) = (bChicTSbotr(i).conj() * aChicTSbotr(i) +
                      aChicTSbotr(i).conj() * bChicTSbotr(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2),
    higgs(2, 2), chargino(2, 2), neutralino(2, 2);
  double a0t1 = a0(msd(1), q), a0t2 = a0(msd(2), q);
  double ft1 = ffn(p1, msd(1), 0.0, q), ft2 = ffn(p2, msd(2), 0.0, q);
  double ggt1 = gfn(p1, mg, 0., q), ggt2 = gfn(p2, mg, 0., q);
  strong(1, 1) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt1 + ft1 + a0t1);
  strong(2, 2) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt2 + ft2 + a0t2);

  for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) +
      0.5 * (- sqr(g) * gdL * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) +
      0.5 * (- sqr(g) * gdR * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) + a0p *
          (sqr(g) * (gdL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) = higgs(2, 2) + a0p *
      (sqr(g) * gdR / (2.0 * sqr(costhDrbar)) * cn(i));
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsm(i), msd(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHSbotLSbot12(i, j)) * b0p;
      b0p = b0(p2, higgsm(i), msd(j), q);
      higgs(2, 2) = higgs(2, 2) + sqr(lHSbotRSbot12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsc(i), msup(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lChHsbotLstop12(i, j)) * b0p;
      b0p = b0(p2, higgsc(i), msup(j), q);
    }

  DoubleMatrix electroweak(2, 2);
  // line by line....
  electroweak(1, 1) = electroweak(1, 1) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e / 3.0) *
    ffn(p1, msd(1), 0., q);
  electroweak(1, 1) = electroweak(1, 1) +
    sqr(g) / sqr(costhDrbar) * sqr(gdL) * ffn(p1, msd(1), mz, q) +
    sqr(g) * 0.5 * ffn(p1, msup(1), mw, q);
  electroweak(1, 1) = electroweak(1, 1) +
    sqr(g) * 0.25 * (a0(msd(1), q) + 2.0 * (a0(msup(1), q)));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) *
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(ydL) * a0(msd(1), q);
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * ydL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  electroweak(2, 2) = electroweak(2, 2) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdR) * a0(mz, q) +
    sqr(e / 3.0) *
    (ffn(p2, msd(2), 0., q));
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(g) / sqr(costhDrbar) * sqr(gdR) *
    (ffn(p2, msd(2), mz, q));
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(gp) * 0.25 * sqr(ydR) *
    (a0(msd(2), q));
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * ydR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  for (i=1; i<=2; i++) {
    double one = gfn(p1, mch(i), 0., q);
    chargino(1, 1) = chargino(1, 1) + fChTSbotLL(i) * one;
    one = gfn(p2, mch(i), 0., q);
    chargino(2, 2) = chargino(2, 2) + fChTSbotRR(i) * one;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p1, mneut(i), 0., q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0BotSbotLL(i) * one;
    one = gfn(p2, mneut(i), 0., q);
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0BotSbotRR(i) * one;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (strong + higgs + chargino + neutralino + electroweak);

  mass = mass - piSq;
}

// Checked 16.08.2005
void MssmSoftsusy::addSbotCorrection(double p,
                                     DoubleMatrix & mass, double mt) {

// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) {
    if (PRINTOUT > 2) cout << " sbottom tachyonic ";
    flagTachyon(true);
    return;
  }

  // one-loop correction matrix
  DoubleMatrix piSq(2, 2); // Self-energy matrix

  // brevity
  double    mz      = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    e       = g * sinthDrbar;
  double    beta    = atan(displayTanb());

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double    thetatau= forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double    ctau    = cos(thetatau);
  double    stau    = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)          = forLoops.msnu(1);
  msnu(2)          = forLoops.msnu(2);
  msnu(3)          = forLoops.msnu(3);

  double    mg      = forLoops.mGluino;
  double    smu     = -displaySusyMu();
  double q = displayMu(), g3sq = sqr(displayGaugeCoupling(3)),
    ht = displayYukawaElement(YU, 3, 3),
    hb = displayYukawaElement(YD, 3, 3),
    mb = dataSet.displayMass(mBottom),
    htsq = sqr(displayYukawaElement(YU, 3, 3)),
    sinb = sin(beta), cosb = cos(beta),
    hbsq = sqr(displayYukawaElement(YD, 3, 3));

  //  double    p       = sqrt(msbot(1) * msbot(2));

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleMatrix lsSbotLSbotLR(4, 2), lsSbotLSbot12(4, 2);
  DoubleMatrix lsSbotRSbotLR(4, 2), lsSbotRSbot12(4, 2);
  // Order (s1 s2 G A, L R)
  lsSbotLSbotLR(1, 1) = g * mz * gdL * cosb / costhDrbar +
    root2 * hb * mb;
  lsSbotLSbotLR(1, 2) = displayTrilinear(DA, 3, 3) / root2;
  lsSbotLSbotLR(3, 2) = -1.0 / root2 *
    (smu * sinb * hb + displayTrilinear(DA, 3, 3) * cosb);
  lsSbotLSbotLR(4, 2) = -1.0 / root2 *
    (smu * cosb * hb - displayTrilinear(DA, 3, 3) * sinb);
  lsSbotLSbotLR(2, 1) = -g * mz * gdL * sinb / costhDrbar;
  lsSbotLSbotLR(2, 2) = hb * smu / root2;

  lsSbotRSbotLR(1, 1) = lsSbotLSbotLR(1, 2);
  lsSbotRSbotLR(1, 2) = g * mz * gdR * cosb / costhDrbar +
    root2 * hb * mb;
  lsSbotRSbotLR(2, 1) = lsSbotLSbotLR(2, 2);
  lsSbotRSbotLR(2, 2) = -g * mz * gdR * sinb / costhDrbar;
  lsSbotRSbotLR(3, 1) = -lsSbotLSbotLR(3, 2);
  lsSbotRSbotLR(4, 1) = -lsSbotLSbotLR(4, 2);

  // Mix sbots up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsSbotLSbotLR(i, 1);
    temp(2) = lsSbotLSbotLR(i, 2);
    temp2 = rot2d(thetab) * temp;
    lsSbotLSbot12(i, 1) = temp2(1);
    lsSbotLSbot12(i, 2) = temp2(2);
    temp(1) = lsSbotRSbotLR(i, 1);
    temp(2) = lsSbotRSbotLR(i, 2);
    temp2 = rot2d(thetab) * temp;
    lsSbotRSbot12(i, 1) = temp2(1);
    lsSbotRSbot12(i, 2) = temp2(2);
 }

  DoubleMatrix lHSbotLSbot12(lsSbotLSbot12), lHSbotRSbot12(lsSbotRSbot12);
  // Mix CP-even Higgses up
  for (i=1; i<=2; i++) { // i is the s1/s2 label
    temp(1) = lsSbotLSbot12(1, i);
    temp(2) = lsSbotLSbot12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotLSbot12(1, i) = temp2(1);
    lHSbotLSbot12(2, i) = temp2(2);
    temp(1) = lsSbotRSbot12(1, i);
    temp(2) = lsSbotRSbot12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotRSbot12(1, i) = temp2(1);
    lHSbotRSbot12(2, i) = temp2(2);
  }


  // Charged Higgs Feynman rules
  DoubleMatrix lChHsbotLstopLR(2, 2); // (G+ H+, L R) basis
  lChHsbotLstopLR(1, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2
    - ht * mt * sinb + hb * mb * cosb;
  lChHsbotLstopLR(1, 2) = -ht * smu * cosb - displayTrilinear(UA, 3, 3) * sinb;
  lChHsbotLstopLR(2, 1) = g * displayMwRun() * sin(2.0 * beta) / root2
    - ht * mt * cosb - hb * mb * sinb;
  lChHsbotLstopLR(2, 2) = ht * smu * sinb - displayTrilinear(UA, 3, 3) * cosb;

  DoubleMatrix lChHsbotLstop12(2, 2);
  temp(1) = lChHsbotLstopLR(1, 1);
  temp(2) = lChHsbotLstopLR(1, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotLstop12(1, 1) = temp2(1);
  lChHsbotLstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotLstopLR(2, 1);
  temp(2) = lChHsbotLstopLR(2, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotLstop12(2, 1) = temp2(1);
  lChHsbotLstop12(2, 2) = temp2(2);

  DoubleMatrix lChHsbotRstopLR(2, 2); // (G+ H+, L R) basis
  lChHsbotRstopLR(2, 1) = hb * smu * cosb - displayTrilinear(DA, 3, 3) * sinb;
  lChHsbotRstopLR(2, 2) = -ht * mb * cosb - hb * mt * sinb;
  lChHsbotRstopLR(1, 1) = hb * smu * sinb + displayTrilinear(DA, 3, 3) * cosb;
  DoubleMatrix lChHsbotRstop12(2, 2);
  temp(1) = lChHsbotRstopLR(1, 1);
  temp(2) = lChHsbotRstopLR(1, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotRstop12(1, 1) = temp2(1);
  lChHsbotRstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotRstopLR(2, 1);
  temp(2) = lChHsbotRstopLR(2, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotRstop12(2, 1) = temp2(1);
  lChHsbotRstop12(2, 2) = temp2(2);

  // Neutralino Feynman rules
  DoubleVector aPsi0BSbotr(4), bPsi0BSbotr(4), aPsi0BSbotl(4),
    bPsi0BSbotl(4);
  aPsi0BSbotr(1) = ydR * gp / root2;
  bPsi0BSbotl(1) = gp / (3.0 * root2);
  bPsi0BSbotl(2) = -g / root2;
  aPsi0BSbotl(3) = hb;
  bPsi0BSbotr(3) = hb;

  ComplexVector aChi0BSbotl(4), bChi0BSbotl(4), aChi0BSbotr(4),
    bChi0BSbotr(4);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0BSbotl = n.complexConjugate() * aPsi0BSbotl;
  bChi0BSbotl = n * bPsi0BSbotl;
  aChi0BSbotr = n.complexConjugate() * aPsi0BSbotr;
  bChi0BSbotr = n * bPsi0BSbotr;

  DoubleVector gChi0BotSbotLL(4), fChi0BotSbotLL(4);
  DoubleVector gChi0BotSbotLR(4), fChi0BotSbotLR(4);
  DoubleVector gChi0BotSbotRR(4), fChi0BotSbotRR(4);
  for (i=1; i<=4; i++) {
    fChi0BotSbotLL(i) = (aChi0BSbotl(i) * aChi0BSbotl(i).conj() +
      bChi0BSbotl(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLL(i) = (bChi0BSbotl(i).conj() * aChi0BSbotl(i) +
      bChi0BSbotl(i) * aChi0BSbotl(i).conj()).real();
    fChi0BotSbotRR(i) = (aChi0BSbotr(i) * aChi0BSbotr(i).conj() +
      bChi0BSbotr(i) * bChi0BSbotr(i).conj()).real();
    gChi0BotSbotRR(i) = (bChi0BSbotr(i).conj() * aChi0BSbotr(i) +
      bChi0BSbotr(i) * aChi0BSbotr(i).conj()).real();
    fChi0BotSbotLR(i) = (aChi0BSbotr(i) * aChi0BSbotl(i).conj() +
      bChi0BSbotr(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLR(i) = (bChi0BSbotl(i).conj() * aChi0BSbotr(i) +
      bChi0BSbotr(i) * aChi0BSbotl(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicTSbotl(2), bPsicTSbotr(2), aPsicTSbotl(2),
    aPsicTSbotr(2);

  bPsicTSbotl(1) = g;
  bPsicTSbotr(2) = -hb;
  aPsicTSbotl(2) = -ht;

  ComplexVector aChicTSbotr(2), aChicTSbotl(2), bChicTSbotl(2),
      bChicTSbotr(2);
  ComplexMatrix aChTSbot(2, 2), bChTSbot(2, 2);

  aChicTSbotl = v.complexConjugate() * aPsicTSbotl;
  bChicTSbotl = u * bPsicTSbotl;
  aChicTSbotr = v.complexConjugate() * aPsicTSbotr;
  bChicTSbotr = u * bPsicTSbotr;

  DoubleVector fChTSbotLL(2), gChTSbotLL(2) ;
  DoubleVector fChTSbotLR(2), gChTSbotLR(2);
  DoubleVector fChTSbotRR(2), gChTSbotRR(2);
  for (i=1; i<=2; i++) {
    fChTSbotLL(i) = (aChicTSbotl(i).conj() * aChicTSbotl(i) +
                      bChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    gChTSbotLL(i) = (bChicTSbotl(i).conj() * aChicTSbotl(i) +
                      aChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    fChTSbotLR(i) = (aChicTSbotl(i).conj() * aChicTSbotr(i) +
                      bChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    gChTSbotLR(i) = (bChicTSbotl(i).conj() * aChicTSbotr(i) +
                      aChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    fChTSbotRR(i) = (aChicTSbotr(i).conj() * aChicTSbotr(i) +
                      bChicTSbotr(i).conj() * bChicTSbotr(i)).real();
    gChTSbotRR(i) = (bChicTSbotr(i).conj() * aChicTSbotr(i) +
                      aChicTSbotr(i).conj() * bChicTSbotr(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2),
    higgs(2, 2), chargino(2, 2), neutralino(2, 2);
  double a0t1 = a0(msbot(1), q), a0t2 = a0(msbot(2), q);
  double ft1 = ffn(p, msbot(1), 0.0, q), ft2 = ffn(p, msbot(2), 0.0, q);
  double ggt = gfn(p, mg, mb, q);
  strong(1, 1) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(cb) * (ft1 + a0t1) + sqr(sb) * (ft2 + a0t2));
  strong(2, 2) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(sb) * (ft1 + a0t1) + sqr(cb) * (ft2 + a0t2));
  strong(1, 2) = 4.0 * g3sq / 3.0 *
    (4.0 * mg * mb * b0(p, mg, mb, q) + sb * cb * (ft1 - a0t1 - ft2 + a0t2));

  stop(1, 1) = hbsq * (sqr(sb) * a0t1 + sqr(cb) * a0t2);
  stop(2, 2) = hbsq * (sqr(cb) * a0t1 + sqr(sb) * a0t2);
  stop(1, 2) = hbsq * cb * sb * 3.0 * (a0t1 - a0t2);

  sbottom(1, 1) =
    htsq * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q));
  sbottom(2, 2) =
    hbsq * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q));

  for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) +
      0.5 * (hbsq * dnd(i) - sqr(g) * gdL * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) +
      0.5 * (hbsq * dnd(i) - sqr(g) * gdR * cn(i)
             / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) + a0p *
          (htsq * dnd(i) +
           sqr(g) * (gdL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) = higgs(2, 2) + a0p *
      (hbsq * dnu(i) + sqr(g) * gdR / (2.0 * sqr(costhDrbar)) * cn(i));
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsm(i), msbot(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHSbotLSbot12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) +
        lHSbotLSbot12(i, j) * lHSbotRSbot12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lHSbotRSbot12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), mstop(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lChHsbotLstop12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) +
        lChHsbotLstop12(i, j) * lChHsbotRstop12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lChHsbotRstop12(i, j)) * b0p;
    }

  DoubleMatrix electroweak(2, 2);
  // line by line....
  electroweak(1, 1) = electroweak(1, 1) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e / 3.0) *
    (sqr(cb) * ffn(p, msbot(1), 0., q) +
     sqr(sb) * ffn(p, msbot(2), 0., q));
  electroweak(1, 1) = electroweak(1, 1) +
    sqr(g) / sqr(costhDrbar) * sqr(gdL) *
    (sqr(cb) * ffn(p, msbot(1), mz, q) +
     sqr(sb) * ffn(p, msbot(2), mz, q)) +
    sqr(g) * 0.5 * (sqr(ct) * ffn(p, mstop(1), mw, q) +
                    sqr(st) * ffn(p, mstop(2), mw, q));
  electroweak(1, 1) = electroweak(1, 1) +
    sqr(g) * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q) +
                     2.0 * (sqr(ct) * a0(mstop(1), q) +
                            sqr(st) * a0(mstop(2), q)));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) *
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(ydL) *
    (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * ydL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  electroweak(2, 2) = electroweak(2, 2) +
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdR) * a0(mz, q) +
    sqr(e / 3.0) *
    (sqr(sb) * ffn(p, msbot(1), 0., q) +
     sqr(cb) * ffn(p, msbot(2), 0., q));
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(g) / sqr(costhDrbar) * sqr(gdR) *
    (sqr(sb) * ffn(p, msbot(1), mz, q) +
     sqr(cb) * ffn(p, msbot(2), mz, q));
  electroweak(2, 2) = electroweak(2, 2) +
    sqr(gp) * 0.25 * sqr(ydR) *
    (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q));
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * ydR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(2), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));


  electroweak(1, 2) = electroweak(1, 2) +
    sqr(gp) * sb * cb * 0.25 * ydL * ydR *
    (a0(msbot(1), q) - a0(msbot(2), q)) +
    sqr(e / 3.0) * sb * cb *
    (ffn(p, msbot(1), 0., q) - ffn(p, msbot(2), 0., q)) -
    sqr(g) / sqr(costhDrbar) * gdL * gdR * sb * cb *
    (ffn(p, msbot(1), mz, q) - ffn(p, msbot(2), mz, q));

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), mt, q);
    double two = mch(i) * mt * b0(p, mch(i), mt, q) * 2.0;
    chargino(1, 1) = chargino(1, 1) + fChTSbotLL(i) * one -
      gChTSbotLL(i) * two;
    chargino(1, 2) = chargino(1, 2) + fChTSbotLR(i) * one -
      gChTSbotLR(i) * two;
    chargino(2, 2) = chargino(2, 2) + fChTSbotRR(i) * one -
      gChTSbotRR(i) * two;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p, mneut(i), mb, q);
    double two = 2.0 * mneut(i) * mb * b0(p, mneut(i), mb, q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0BotSbotLL(i) * one - gChi0BotSbotLL(i) * two;
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0BotSbotRR(i) * one - gChi0BotSbotRR(i) * two;
    neutralino(1, 2) = neutralino(1, 2) +
      fChi0BotSbotLR(i) * one - gChi0BotSbotLR(i) * two;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (strong + stop + sbottom + higgs + chargino + neutralino + electroweak);

  piSq(2, 1) = piSq(1, 2);

  mass = mass - piSq;
}

// As in BPMZ appendix, INCLUDING weak boson loops.
void MssmSoftsusy::addSupCorrection(DoubleMatrix & mass, int family) {

// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) {
    if (PRINTOUT > 2) cout << " sup tachyonic ";
    flagTachyon(true);
    return;
  }

  // one-loop correction matrix
  DoubleMatrix piSq(2, 2); // Self-energy matrix

  // brevity
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha   = forLoops.thetaH;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());
  double    e       = g * sinthDrbar; // DRbar value of e

  DoubleVector msbot(2);
  msbot(1) = forLoops.md(1, 3);
  msbot(2) = forLoops.md(2, 3);
  double    thetat  = forLoops.thetat;
  double    thetab  = forLoops.thetab;
  double   thetatau = forLoops.thetatau;
  double    ct      = cos(thetat) ;
  double    st      = sin(thetat) ;
  double    cb      = cos(thetab) ;
  double    sb      = sin(thetab) ;
  double  ctau      = cos(thetatau);
  double  stau      = sin(thetatau);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector msup(2);
  msup(1)          = forLoops.mu(1, family);
  msup(2)          = forLoops.mu(2, family);
  DoubleVector msd(2);
  msd(1)          = forLoops.mu(1, family);
  msd(2)          = forLoops.mu(2, family);

  double    mg      = forLoops.mGluino;
  double q = displayMu(), g3sq = sqr(displayGaugeCoupling(3)),
    sinb = sin(beta), cosb = cos(beta),
    mz = displayMzRun();

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  // top now stands for up
  DoubleMatrix lsStopLStopLR(4, 2), lsStopLStop12(4, 2);
  DoubleMatrix lsStopRStopLR(4, 2), lsStopRStop12(4, 2);
  // Order (s1 s2 G A, L R)
  lsStopLStopLR(1, 1) = g * mz * guL * cosb / costhDrbar;
  lsStopLStopLR(2, 1) = -g * mz * guL * sinb / costhDrbar;

  lsStopRStopLR(1, 1) = lsStopLStopLR(1, 2);
  lsStopRStopLR(1, 2) = g * mz * guR * cosb / costhDrbar;
  lsStopRStopLR(2, 1) = lsStopLStopLR(2, 2);
  lsStopRStopLR(2, 2) = -g * mz * guR * sinb / costhDrbar;
  lsStopRStopLR(3, 1) = -lsStopLStopLR(3, 2);
  lsStopRStopLR(4, 1) = -lsStopLStopLR(4, 2);

  // Mix stops up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsStopLStopLR(i, 1);
    temp(2) = lsStopLStopLR(i, 2);
    temp2 = temp;
    lsStopLStop12(i, 1) = temp2(1);
    lsStopLStop12(i, 2) = temp2(2);
    temp(1) = lsStopRStopLR(i, 1);
    temp(2) = lsStopRStopLR(i, 2);
    temp2 = temp;
    lsStopRStop12(i, 1) = temp2(1);
    lsStopRStop12(i, 2) = temp2(2);
 }

  DoubleMatrix lHStopLStop12(lsStopLStop12), lHStopRStop12(lsStopRStop12);
  // Mix CP-even Higgses up
  for (i=1; i<=2; i++) { // i is the s1/s2 label
    temp(1) = lsStopLStop12(1, i);
    temp(2) = lsStopLStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopLStop12(1, i) = temp2(1);
    lHStopLStop12(2, i) = temp2(2);
    temp(1) = lsStopRStop12(1, i);
    temp(2) = lsStopRStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopRStop12(1, i) = temp2(1);
    lHStopRStop12(2, i) = temp2(2);
  }

  // Charged Higgs Feynman rules
  DoubleMatrix lChHstopLsbotLR(2, 2); // (G+ H+, L R) basis
  lChHstopLsbotLR(2, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;
  lChHstopLsbotLR(1, 1) =-g * displayMwRun() * cos(2.0 * beta) / root2;

  DoubleMatrix lChHstopLsbot12(2, 2);
  temp(1) = lChHstopLsbotLR(1, 1);
  temp(2) = lChHstopLsbotLR(1, 2);
  temp2 = temp;
  lChHstopLsbot12(1, 1) = temp2(1);
  lChHstopLsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopLsbotLR(2, 1);
  temp(2) = lChHstopLsbotLR(2, 2);
  temp2 = temp;
  lChHstopLsbot12(2, 1) = temp2(1);
  lChHstopLsbot12(2, 2) = temp2(2);

  DoubleMatrix lChHstopRsbotLR(2, 2); // (G+ H+, L R) basis
  DoubleMatrix lChHstopRsbot12(2, 2);
  temp(1) = lChHstopRsbotLR(1, 1);
  temp(2) = lChHstopRsbotLR(1, 2);
  temp2 = temp;
  lChHstopRsbot12(1, 1) = temp2(1);
  lChHstopRsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopRsbotLR(2, 1);
  temp(2) = lChHstopRsbotLR(2, 2);
  temp2 = temp;
  lChHstopRsbot12(2, 1) = temp2(1);
  lChHstopRsbot12(2, 2) = temp2(2);

  // Neutralino Feynman rules
  DoubleVector aPsi0TStopr(4), bPsi0TStopr(4), aPsi0TStopl(4),
    bPsi0TStopl(4);
  aPsi0TStopr(1) = - 4.0 * gp / (3.0 * root2);
  bPsi0TStopl(1) = gp / (3.0 * root2);
  bPsi0TStopl(2) = g / root2;

  ComplexVector aChi0TStopl(4), bChi0TStopl(4), aChi0TStopr(4),
    bChi0TStopr(4);
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz);
  DoubleVector mch(forLoops.mchBpmz);

  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  DoubleVector gChi0TopStopLL(4), fChi0TopStopLL(4);
  DoubleVector gChi0TopStopLR(4), fChi0TopStopLR(4);
  DoubleVector gChi0TopStopRR(4), fChi0TopStopRR(4);
  for (i=1; i<=4; i++) {
    fChi0TopStopLL(i) = (aChi0TStopl(i) * aChi0TStopl(i).conj() +
      bChi0TStopl(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLL(i) = (bChi0TStopl(i).conj() * aChi0TStopl(i) +
      bChi0TStopl(i) * aChi0TStopl(i).conj()).real();
    fChi0TopStopRR(i) = (aChi0TStopr(i) * aChi0TStopr(i).conj() +
      bChi0TStopr(i) * bChi0TStopr(i).conj()).real();
    gChi0TopStopRR(i) = (bChi0TStopr(i).conj() * aChi0TStopr(i) +
      bChi0TStopr(i) * aChi0TStopr(i).conj()).real();
    fChi0TopStopLR(i) = (aChi0TStopr(i) * aChi0TStopl(i).conj() +
      bChi0TStopr(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLR(i) = (bChi0TStopl(i).conj() * aChi0TStopr(i) +
      bChi0TStopr(i) * aChi0TStopl(i).conj()).real();
  }

  // Chargino Feynman Rules
  DoubleVector bPsicBStopl(2), bPsicBStopr(2), aPsicBStopl(2),
    aPsicBStopr(2);

  aPsicBStopl(1) = g;

  DoubleVector aPsicCSbotl(2);
  ComplexVector aChicBStopr(2), aChicBStopl(2), bChicBStopl(2),
      bChicBStopr(2);
  ComplexMatrix aChBStop(2, 2), bChBStop(2, 2);

  aChicBStopl = v.complexConjugate() * aPsicBStopl;
  bChicBStopl = u * bPsicBStopl;
  aChicBStopr = v.complexConjugate() * aPsicBStopr;
  bChicBStopr = u * bPsicBStopr;

  DoubleVector fChBStopLL(2), gChBStopLL(2) ;
  DoubleVector fChBStopLR(2), gChBStopLR(2);
  DoubleVector fChBStopRR(2), gChBStopRR(2);
  for (i=1; i<=2; i++) {
    fChBStopLL(i) = (aChicBStopl(i).conj() * aChicBStopl(i) +
                      bChicBStopl(i).conj() * bChicBStopl(i)).real();
    gChBStopLL(i) = (bChicBStopl(i).conj() * aChicBStopl(i) +
                      aChicBStopl(i).conj() * bChicBStopl(i)).real();
    fChBStopLR(i) = (aChicBStopl(i).conj() * aChicBStopr(i) +
                      bChicBStopl(i).conj() * bChicBStopr(i)).real();
    gChBStopLR(i) = (bChicBStopl(i).conj() * aChicBStopr(i) +
                      aChicBStopl(i).conj() * bChicBStopr(i)).real();
    fChBStopRR(i) = (aChicBStopr(i).conj() * aChicBStopr(i) +
                      bChicBStopr(i).conj() * bChicBStopr(i)).real();
    gChBStopRR(i) = (bChicBStopr(i).conj() * aChicBStopr(i) +
                      aChicBStopr(i).conj() * bChicBStopr(i)).real();
  }

  // Corrections themselves start here
  // Full corrections from BPMZ w/ g=g'=e=0
  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2),
    higgs(2, 2), chargino(2, 2), neutralino(2, 2);
  double a0t1 = a0(msup(1), q), a0t2 = a0(msup(2), q);
  double p1 = msup(1), p2 = msup(2);
  p1 = p2 = sqrt(msup(1) * msup(2)); // DEBUG

  double ft1 = ffn(p1, msup(1), 0.0, q), ft2 = ffn(p2, msup(2), 0.0, q);
  double ggt1 = gfn(p1, mg, 0., q), ggt2 = gfn(p2, mg, 0., q);
  strong(1, 1) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt1 + (ft1 + a0t1));
  strong(2, 2) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt2 + (ft2 + a0t2));

  for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) +
      0.5 * (- sqr(g) * guL * 0.5 / sqr(costhDrbar) * cn(i))
      * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) +
      0.5 * (- sqr(g) * guR * 0.5 / sqr(costhDrbar) * cn(i))
      * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) +
      (sqr(g) * (guL * 0.5 / costhDrbar2 - 0.5) * cn(i))* a0p;
    higgs(2, 2) = higgs(2, 2) +
      (sqr(g) * guR * 0.5 / costhDrbar2 * cn(i))* a0p;
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsm(i), msup(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHStopLStop12(i, j)) * b0p;
      b0p = b0(p2, higgsm(i), msup(j), q);
      higgs(2, 2) = higgs(2, 2) + sqr(lHStopRStop12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsc(i), msbot(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lChHstopLsbot12(i, j)) * b0p;
      b0p = b0(p2, higgsc(i), msbot(j), q);
      higgs(2, 2) = higgs(2, 2) + sqr(lChHstopRsbot12(i, j)) * b0p;
    }

  // EW bosons
  higgs(1, 1) = higgs(1, 1) +
    4.0 * sqr(g) / costhDrbar2 * sqr(guL) * a0(mz, q) +
    2.0 * sqr(g) * a0(displayMwRun(), q) + sqr(2.0 / 3.0 * e) *
    ffn(p1, msup(1), 0.0, q)
    + sqr(g * guL / costhDrbar) * ffn(p1, msup(1), mz, q) +
    sqr(g) * 0.5 * ffn(p1, msd(1), displayMwRun(), q) +
    sqr(g) * 0.25 *
    a0(msup(1), q) + 2.0 * a0(msd(1), q) +
    sqr(g) * 0.5 *
    (1.5 * a0(forLoops.mu(1, 1), q) +
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) +
            sqr(st) * a0(forLoops.mu(2, 2), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
            sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) +
            a0(forLoops.msnu(2), q) +
            a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(yuL) * a0(msup(1), q) +
    sqr(gp) * 0.25 * yuL *
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) +
                  a0(forLoops.mu(1, 2), q) +
                  sqr(ct) * a0(forLoops.mu(1, 3), q) +
                  sqr(st) * a0(forLoops.mu(2, 2), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) +
                  a0(forLoops.mu(2, 2), q) +
                  sqr(st) * a0(forLoops.mu(1, 3), q) +
                  sqr(ct) * a0(forLoops.mu(2, 2), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) +
                  a0(forLoops.md(1, 2), q) +
                  sqr(cb) * a0(forLoops.md(1, 3), q) +
                  sqr(sb) * a0(forLoops.md(2, 2), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) +
                  a0(forLoops.md(2, 2), q) +
                  sqr(sb) * a0(forLoops.md(1, 3), q) +
                  sqr(cb) * a0(forLoops.md(2, 2), q)) +
     yeL * (a0(forLoops.me(1, 1), q) +
            a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 2), q)) +
     yeR * (a0(forLoops.me(2, 1), q) +
            a0(forLoops.me(2, 2), q) +
            sqr(stau) * a0(forLoops.me(1, 3), q) +
            sqr(ctau) * a0(forLoops.me(2, 2), q)) +
     ynuL * (a0(forLoops.msnu(1), q) +
             a0(forLoops.msnu(2), q) +
             a0(forLoops.msnu(3), q)));

  higgs(2, 2) = higgs(2, 2) +
    4.0 * sqr(g) / costhDrbar2 * sqr(guR) * a0(mz, q) +
    sqr(2.0 / 3.0 * e) *
    (ffn(p2, msup(2), 0.0, q))
    + sqr(g * guR / costhDrbar) *
    (ffn(p2, msup(2), mz, q)) +
    sqr(gp) * 0.25 * sqr(yuR) *
    (a0(msup(2), q)) +
    sqr(gp) * 0.25 * yuR *
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) +
                  a0(forLoops.mu(1, 2), q) +
                  sqr(ct) * a0(forLoops.mu(1, 3), q) +
                  sqr(st) * a0(forLoops.mu(2, 2), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) +
                  a0(forLoops.mu(2, 2), q) +
                  sqr(st) * a0(forLoops.mu(1, 3), q) +
                  sqr(ct) * a0(forLoops.mu(2, 2), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) +
                  a0(forLoops.md(1, 2), q) +
                  sqr(cb) * a0(forLoops.md(1, 3), q) +
                  sqr(sb) * a0(forLoops.md(2, 2), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) +
                  a0(forLoops.md(2, 2), q) +
                  sqr(sb) * a0(forLoops.md(1, 3), q) +
                  sqr(cb) * a0(forLoops.md(2, 2), q)) +
     yeL * (a0(forLoops.me(1, 1), q) +
            a0(forLoops.me(1, 2), q) +
            sqr(ctau) * a0(forLoops.me(1, 3), q) +
            sqr(stau) * a0(forLoops.me(2, 2), q)) +
     yeR * (a0(forLoops.me(2, 1), q) +
            a0(forLoops.me(2, 2), q) +
            sqr(stau) * a0(forLoops.me(1, 3), q) +
            sqr(ctau) * a0(forLoops.me(2, 2), q)) +
     ynuL * (a0(forLoops.msnu(1), q) +
             a0(forLoops.msnu(2), q) +
             a0(forLoops.msnu(3), q)));

  for (i=1; i<=2; i++) {
    double one = gfn(p1, mch(i), 0., q);
    double two = 0.;
    chargino(1, 1) = chargino(1, 1) + fChBStopLL(i) * one -
      gChBStopLL(i) * two;
    one = gfn(p2, mch(i), 0., q);
    chargino(2, 2) = chargino(2, 2) + fChBStopRR(i) * one -
      gChBStopRR(i) * two;
  }

  for (i=1; i<=4; i++) {
    double one = gfn(p1, mneut(i), 0., q);
    double two = 0.;
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0TopStopLL(i) * one - gChi0TopStopLL(i) * two;
    one = gfn(p2, mneut(i), 0., q);
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0TopStopRR(i) * one - gChi0TopStopRR(i) * two;
  }

  piSq = 1.0 / (16.0 * sqr(PI)) *
    (strong + stop + sbottom + higgs + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);

  mass = mass - piSq;
}

void MssmSoftsusy::doUpSquarks(double mt, double pizztMS, double sinthDRbarMS,
                               int accuracy) {

  // first two families are simpler
  int family; for (family = 1; family <= 2; family++) {

    DoubleMatrix mStopSquared(2, 2);
    treeUpSquark(mStopSquared, mt, pizztMS, sinthDRbarMS, family);

    if (accuracy > 0) addSupCorrection(mStopSquared, family);

    double theta;
    DoubleVector
      physicalStopMassesSquared(mStopSquared.sym2by2(theta));
    if (physicalStopMassesSquared(1) < 0.0 ||
        physicalStopMassesSquared(2) < 0.0) {
      flagTachyon(true);
      if (PRINTOUT > 2) cout << " tachyonic up squarks ";
    }

    DoubleVector physicalStopMasses(physicalStopMassesSquared.apply(ccbSqrt));

    physpars.mu(1, family) = physicalStopMasses(1);
    physpars.mu(2, family) = physicalStopMasses(2);
  }

  family = 3;

  DoubleMatrix mStopSquared(2, 2), mStopSquared2(2, 2);
  treeUpSquark(mStopSquared, mt, pizztMS, sinthDRbarMS, family);
  mStopSquared2 = mStopSquared; // StopSquared2 is now tree-level

  // one loop corrections
  if (accuracy > 0) {
    double pLight = minimum(forLoops.mu(1, 3), forLoops.mu(2, 3));
    double pHeavy = maximum(forLoops.mu(1, 3), forLoops.mu(2, 3));

    addStopCorrection(pLight, mStopSquared, mt);
    addStopCorrection(pHeavy, mStopSquared2, mt);
  }

  double theta;
  DoubleVector physicalStopMassesSquared2(mStopSquared2.sym2by2(theta));
  DoubleVector physicalStopMassesSquared(mStopSquared.sym2by2(theta));
  physpars.thetat = theta; // theta is calculated at p=mstop1

  if (physicalStopMassesSquared(1) < 0.0 || physicalStopMassesSquared(2) < 0.0
      ||  physicalStopMassesSquared2(1) < 0.0 ||
      physicalStopMassesSquared2(2) < 0.0) {
    flagTachyon(true);
    if (PRINTOUT > 2) cout << " tachyonic stops ";
  }

  DoubleVector physicalStopMasses(physicalStopMassesSquared.apply(ccbSqrt));
  DoubleVector
      physicalStopMasses2(physicalStopMassesSquared2.apply(ccbSqrt));

  double lightStopMass = minimum(physicalStopMasses(1), physicalStopMasses(2));
  double heavyStopMass = maximum(physicalStopMasses2(1),
                                 physicalStopMasses2(2));
  // twisted measures the ordering of the sbottom masses. If mstop1 > mstop2,
  // twisted is defined to be true (mstop2 > mstop1 is defined "untwisted").
  bool twisted = false;
  if (physicalStopMassesSquared(1) > physicalStopMassesSquared(2))
    twisted = true;

  if (twisted) {
    physpars.mu(1, family) = heavyStopMass;
    physpars.mu(2, family) = lightStopMass;
  } else {
    physpars.mu(1, family) = lightStopMass;
    physpars.mu(2, family) = heavyStopMass;
  }
}

void MssmSoftsusy::treeDownSquark(DoubleMatrix & mass, double mbrun,
                                  double pizztMS, double sinthDRbarMS,
                                  int family) {
  const double cd = 1.0 / 3.0;
  double mz2 = sqr(displayMzRun()), mb2 = sqr(mbrun);
  double beta = atan(displayTanb()), mu = displaySusyMu(),
    c2b = cos(2 * beta), tanb = displayTanb(), sinth2 = sqr(sinthDRbarMS);

  mass(1, 1) = displaySoftMassSquared(mQl, family, family) +
      mz2 * (-0.5 + cd * sinth2) * c2b;
  mass(2, 2) = displaySoftMassSquared(mDr, family, family) -
      mz2 * cd * sinth2 * c2b;

  if (family != 3) mass(1, 2) = 0.0;
  else {
    mass(1, 2) = mbrun * (displaySoftA(DA, 3, 3) - mu * tanb);
    mass(1, 1) = mass(1, 1) + mb2;
    mass(2, 2) = mass(2, 2) + mb2;
  }

  mass(2, 1) = mass(1, 2);
}


void MssmSoftsusy::doDownSquarks(double mb, double pizztMS, double
                            sinthDRbarMS, int accuracy, double mt) {
  int family; for (family = 1; family <= 2; family++) {

    DoubleMatrix mSbotSquared(2, 2);
    treeDownSquark(mSbotSquared, mb, pizztMS, sinthDRbarMS, family);

    if (accuracy > 0) addSdownCorrection(mSbotSquared, family);

    double theta;
    DoubleVector
      physicalSbotMassesSquared(mSbotSquared.sym2by2(theta));
    if (physicalSbotMassesSquared(1) < 0.0 ||
        physicalSbotMassesSquared(2) < 0.0) {
      flagTachyon(true);
      if (PRINTOUT > 2) cout << " tachyonic down squarks ";
    }

    DoubleVector physicalSbotMasses(physicalSbotMassesSquared.apply(ccbSqrt));

    physpars.md(1, family) = physicalSbotMasses(1);
    physpars.md(2, family) = physicalSbotMasses(2);
  }

  // start 3rd family computation: more complicated!
    family = 3;
    DoubleMatrix mSbotSquared(2, 2), mSbotSquared2(2, 2);
    treeDownSquark(mSbotSquared, mb, pizztMS, sinthDRbarMS, family);
    mSbotSquared2 = mSbotSquared; // tree-level matrices

  if (accuracy > 0) {
    double pLight = minimum(forLoops.md(1, 3), forLoops.md(2, 3));
    double pHeavy = maximum(forLoops.md(1, 3), forLoops.md(2, 3));

    addSbotCorrection(pLight, mSbotSquared, mb);
    addSbotCorrection(pHeavy, mSbotSquared2, mb);
  }

  double theta;
  DoubleVector physicalSbotMassesSquared2(mSbotSquared2.sym2by2(theta));
  // The ordering and the mixing angle is defined by the following matrix,
  // with the LIGHT sbottom as the external momentum
  DoubleVector physicalSbotMassesSquared(mSbotSquared.sym2by2(theta));
  physpars.thetab = theta;

  // Check for tachyonic sbottoms
  if (minimum(physicalSbotMassesSquared(1), physicalSbotMassesSquared(2))
      < 0.0 || minimum(physicalSbotMassesSquared(1),
                       physicalSbotMassesSquared(2)) < 0.0) {
    flagTachyon(true);
    if (PRINTOUT > 2) cout << " tachyonic sbottom ";
  }

  DoubleVector physicalSbotMasses(physicalSbotMassesSquared.apply(ccbSqrt));
  DoubleVector physicalSbotMasses2(physicalSbotMassesSquared2.apply(ccbSqrt));

  // twisted measures the ordering of the sbottom masses. If msbot1 > msbot2,
  // twisted is defined to be true (msbot2 > msbot1 is defined "untwisted").
  bool twisted = false;
  if (physicalSbotMassesSquared(1) > physicalSbotMassesSquared(2))
    twisted = true;

  double lightSbotMass = minimum(physicalSbotMasses(1), physicalSbotMasses(2));
  double heavySbotMass = maximum(physicalSbotMasses2(1),
                                 physicalSbotMasses2(2));

  if (twisted) {
    physpars.md(1, family) = heavySbotMass;
    physpars.md(2, family) = lightSbotMass;
  } else {
    physpars.md(1, family) = lightSbotMass;
    physpars.md(2, family) = heavySbotMass;
  }
}

void MssmSoftsusy::treeChargedSlepton(DoubleMatrix & mass, double mtaurun,
                                      double pizztMS, double sinthDRbarMS,
                                      int family) {
  double mz2 = sqr(displayMzRun()), mtau2 = sqr(mtaurun);
  double beta = atan(displayTanb()), mu = displaySusyMu(),
    c2b = cos(2 * beta), tanb = displayTanb(), sinth2 = sqr(sinthDRbarMS);

  mass(1, 1) = displaySoftMassSquared(mLl, family, family) -
    mz2 * (0.5 - sinth2) * c2b;
  mass(2, 2) = displaySoftMassSquared(mEr, family, family) -
    mz2 * sinth2 * c2b;

  if (family != 3) mass(1, 2) = 0.0;
  else {
    mass(1, 2) = mtaurun * (displaySoftA(EA, 3, 3) - mu * tanb);
    mass(1, 1) = mass(1, 1) + mtau2;
    mass(2, 2) = mass(2, 2) + mtau2;
  }

  mass(2, 1) = mass(1, 2);
}

void MssmSoftsusy::doChargedSleptons(double mtau, double pizztMS, double
                            sinthDRbarMS, int accuracy) {

  int family; for (family = 1; family <= 3; family++) {

    DoubleMatrix mStauSquared(2, 2);
    treeChargedSlepton(mStauSquared, mtau, pizztMS, sinthDRbarMS, family);
    if (family == 3 && accuracy > 0) addStauCorrection(mStauSquared, mtau);
    else if (accuracy > 0) addSlepCorrection(mStauSquared, family);

    double theta;
    DoubleVector
      physicalStauMassesSquared(mStauSquared.sym2by2(theta));
    if (physicalStauMassesSquared(1) < 0.0 ||
        physicalStauMassesSquared(2) < 0.0) {
      flagTachyon(true);
      if (PRINTOUT > 2) cout << " stau tachyonic ";
    }
    DoubleVector physicalStauMasses(physicalStauMassesSquared.apply(ccbSqrt));

    if (family == 3) physpars.thetatau = theta;

    physpars.me(1, family) = physicalStauMasses(1);
    physpars.me(2, family) = physicalStauMasses(2);
  }
}

void MssmSoftsusy::doSnu(double pizztMS, int accuracy) {
  double mSnuSquared;
  int family; for (family = 1; family <= 3; family++) {
    treeSnu(mSnuSquared, pizztMS, family);
    if (family == 3 && accuracy > 0) addSnuTauCorrection(mSnuSquared);
    else if (accuracy > 0) addSnuCorrection(mSnuSquared, family);

    physpars.msnu(family) = ccbSqrt(mSnuSquared);
  }
}

void MssmSoftsusy::treeSnu(double & mSnuSquared, double pizztMS, int family) {
  double mz2 = sqr(displayMzRun());
  double beta = atan(displayTanb());
  double c2b = cos(2.0 * beta);

  mSnuSquared = displaySoftMassSquared(mLl, family, family) + 0.5 * mz2
      * c2b;
}

// Organises calculation of physical quantities such as sparticle masses etc
// Call AT MSusy
void MssmSoftsusy::physical(int accuracy) {
  double sinthDRbarMS, piwwtMS, pizztMS;

  calcDrBarPars();

  if (accuracy == 0) {
    sinthDRbarMS = 0.0;
    piwwtMS = 0.0;
    pizztMS = 0.0;
  }
  else {
  sinthDRbarMS = calcSinthdrbar();
  piwwtMS = sqr(displayMwRun()) - sqr(displayMw());
  pizztMS = sqr(displayMzRun()) - sqr(displayMz());
  }

  // Running masses at MSUSY
  double beta = atan(displayTanb());
  double vev = displayHvev();
  double mt = fabs(displayYukawaElement(YU, 3, 3) * vev * sin(beta)) / root2;
  double mb = fabs(displayYukawaElement(YD, 3, 3) * vev * cos(beta)) / root2;
  double mtau = fabs(displayYukawaElement(YE, 3, 3) * vev * cos(beta)) / root2;

  // Sfermion masses: all three families in each
  doUpSquarks(mt, pizztMS, sinthDRbarMS, accuracy);
  doDownSquarks(mb, pizztMS, sinthDRbarMS, accuracy, mt);
  doChargedSleptons(mtau, pizztMS, sinthDRbarMS, accuracy);
  doSnu(pizztMS, accuracy);

  // Charginos/neutralinos/higgs
  MssmSoftsusy * ppp;
  ppp = this;
  ppp->higgs(accuracy, piwwtMS, pizztMS);
  if (accuracy > 0) ppp->higgs(accuracy, piwwtMS, pizztMS); // iterate
  if (accuracy > 1) ppp->higgs(accuracy, piwwtMS, pizztMS); // iterate
  physpars.mhiggs = ppp->displayPhys().mhiggs;
  gluino(accuracy);
  charginos(accuracy, piwwtMS);
  neutralinos(accuracy, piwwtMS, pizztMS);
}

// For a given trial value of the log of field H2, gives the value of the
// potential at the minimum. The following global variables must be set before
// it is called:
static double unificationScale, minTol;
extern double minimufb3(double lnH2) {

  // Save initial parameters
  double initialMu = tempSoft1->displayMu();
  DoubleVector savePars(tempSoft1->display());

  MssmSoftsusy temp(*tempSoft1);

  double mx = unificationScale;
  double mu = -temp.displaySusyMu();

  int family = 1; // family could be 2 as well if they're very different

  // field VEVS
  double htau = temp.displayYukawaElement(YE, 3, 3);
  double h2 =  - exp(lnH2),
    eR = sqrt(fabs(mu * h2 / htau)),
    Lisq = -4.0 * tempSoft1->displaySoftMassSquared(mLl, family, family) /
    (0.6 * sqr(tempSoft1->displayGaugeCoupling(1)) +
     sqr(tempSoft1->displayGaugeCoupling(2))) + sqr(h2) + sqr(eR);

  // qhat is scale at which one-loop potential corrections are smallest:
  // iterate
  double qhat = getQhat(minTol,eR, h2, Lisq, mx, temp);
  if (qhat < 0.0) {
    ostringstream ii;
    ii << "Error: qhat returned negative in ufb3min" << endl;
    throw ii.str();
  }

  double vufb3 = ufb3fn(mu, htau, h2, family, temp);
  if (PRINTOUT > 1) cout << vufb3 << endl;

  // Restore initial parameters
  tempSoft1->setMu(initialMu);
  tempSoft1->set(savePars);

  return vufb3;
}

// from hep-ph/9507294 -- debugged 19/11/04
double ufb3fn(double mu, double htau, double h2, int family, const MssmSoftsusy
              & temp) {
  double vufb3 = 0.0;
  // potential value for these VEVs
  if (fabs(h2) >
      sqrt(sqr(mu) / (4.0 * sqr(htau)) +
           4.0 * temp.displaySoftMassSquared(mLl, family, family) /
           (0.6 * sqr(temp.displayGaugeCoupling(1)) +
            sqr(temp.displayGaugeCoupling(2)))) - fabs(mu) /
      temp.displayYukawaElement(YE, 3, 3) * 0.5)
    vufb3 =
      sqr(h2) * (temp.displayMh2Squared() +
                 temp.displaySoftMassSquared(mLl, family, family)) +
      fabs(mu * h2) / htau *
      (temp.displaySoftMassSquared(mLl, 3, 3) +
       temp.displaySoftMassSquared(mEr, 3, 3)
       + temp.displaySoftMassSquared(mLl, family, family)) -
      2.0 * sqr(temp.displaySoftMassSquared(mLl, family, family)) /
      (0.6 * sqr(temp.displayGaugeCoupling(1)) +
       sqr(temp.displayGaugeCoupling(2)));
  else
    vufb3 =
      sqr(h2) * temp.displayMh2Squared() +
      fabs(mu * h2) / htau *
      (temp.displaySoftMassSquared(mLl, 3, 3) +
       temp.displaySoftMassSquared(mEr, 3, 3)) +
      1.0 / 8.0 * (0.6 * sqr(temp.displayGaugeCoupling(1)) +
                   sqr(temp.displayGaugeCoupling(2))) *
      sqr(sqr(h2) + fabs(mu * h2) / htau);

  if (PRINTOUT > 1) cout << vufb3 << endl;
  return vufb3;
}

// For ufb3direction, returns scale at which one-loop corrections are smallest
double getQhat(double inminTol,double eR, double h2, double Lisq, double mx,
                MssmSoftsusy & temp) {
  double oldQhat = -1.0e16;
  int maxNum = 40;

  int d; for (d = 1; d <= maxNum; d++)     {
    double qhat =
      maximum(maximum(maximum(temp.displayGaugeCoupling(2) * eR,
                     temp.displayGaugeCoupling(2) * fabs(h2)),
                temp.displayGaugeCoupling(2) * sqrt(fabs(Lisq))),
           temp.displayYukawaElement(YU, 3, 3) * fabs(h2));
    // Run all paramaters to that scale
    if (qhat < mx) temp.runto(qhat);
    else temp.runto(mx);
    if (PRINTOUT > 1) cout << qhat << " ";

    if (fabs((qhat - oldQhat) / qhat) < inminTol) return qhat;
    oldQhat = qhat;
  }
  // Return NOB if no convergence on qhat
  return -6.66e66;
}

// Input mx the scale up to which you search for minima
// Returns minimum value of potential along that direction
// Does ufbs truly properly but takes ages.
double MssmSoftsusy::ufb3sl(double mx) {

  // Save initial parameters
  double initialMu = displayMu();
  DoubleVector parSave(display());

  unificationScale = mx;
  // Running should be faster from high scales on the whole
  runto(mx);

  double ax = log(1.0e3), bx = 0.5*(log(mx) + log(1.0e3)), cx = log(mx);

  double tol = 1.0e-3, lnx; minTol = tol * 10.0;
  // Numerical recipes routine to determine minimum of potential specified in
  // minimufb3
  double Vmin = findMinimum(ax, bx, cx, minimufb3, tol, &lnx);

  // Restore initial parameters
  setMu(initialMu);
  set(parSave);

  return Vmin;
}

// Does SUSY (and other) threshold corrections to alphaS
// Input alphas in MSbar and it returns it in DRbar scheme.
// From hep-ph/9606211
double MssmSoftsusy::qcdSusythresh(double alphasMSbar, double q) const {
  drBarPars tree(displayDrBarPars());
  double mt = tree.mt;
  double deltaAlphas = alphasMSbar / (2.0 * PI) *
    (0.5 - 2.0 / 3.0 * log(mt / q) -
     2.0 * log(fabs(tree.mGluino) / q));

  int i,j;
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++)
      deltaAlphas = deltaAlphas - alphasMSbar / (12.0 * PI) *
        (log(tree.mu(i, j) / q) +
         log(tree.md(i, j) / q));
  return alphasMSbar / (1.0 - deltaAlphas);
}

/* BPMZ version of alpha calculation
// Does SUSY (and other) threshold corrections to alphaEm - returns alpha in
// DRbar scheme at scale Q. From hep-ph/9606211.....seeing if it affects sin^2
// theta_eff calculation
double MssmSoftsusy::qedSusythresh(double alphaEm, double q) const {

  drBarPars tree(displayDrBarPars());

  if (tree.mhiggs(4) < TOLERANCE) return 0.0;

  double alpha = 1.0 / 137.036;


  double deltaAlpha = 0.0682 -  alpha / (2.0 * PI) *
    (-7.0 * log(MW / q)
     + 16.0 / 9.0 * log(tree.mt / q)
     + log(tree.mhiggs(4) / q) / 3.0 + 4.0 / 9.0 *
     (log(tree.mu(1,1) / q) +
      log(tree.mu(2,1) / q) +
      log(tree.mu(1,2) / q) +
      log(tree.mu(2,2) / q) +
      log(tree.mu(1,3) / q) +
      log(tree.mu(2,3) / q)) +
     (log(tree.md(1,1) / q) +
      log(tree.md(2,1) / q) +
      log(tree.md(1,2) / q) +
      log(tree.md(2,2) / q) +
      log(tree.md(1,3) / q) +
      log(tree.md(2,3) / q))
     / 9.0 +
     (log(tree.me(1,1) / q) +
      log(tree.me(2,1) / q) +
      log(tree.me(1,2) / q) +
      log(tree.me(2,2) / q) +
      log(tree.me(1,3) / q) +
      log(tree.me(2,3) / q)) / 3.0
     + (log(fabs(tree.mch(1)) / q)
        + log(fabs(tree.mch(2)) / q))
     * 4.0 / 3.0);

  return alpha / (1.0 - deltaAlpha);
}
*/

// Does SUSY (and other) threshold corrections to alphaEm - returns alpha in
// DRbar scheme at scale Q. From hep-ph/9606211. Input empirical value of
// alpha at MZ external momentum....
double MssmSoftsusy::qedSusythresh(double alphaEm, double q) const {

  drBarPars tree(displayDrBarPars());

  if (tree.mhiggs(4) < TOLERANCE) return 0.0;

  double deltaAlpha = -  alphaEm / (2.0 * PI) *
    // commented out since alpha(MZ) includes it!
    (//-7.0 * log(tree.mw / q)
     -1.0 / 3.0 + 16.0 / 9.0 * log(tree.mt / q)
     + log(tree.mhiggs(4) / q) / 3.0 + 4.0 / 9.0 *
     (log(tree.mu(1,1) / q) +
      log(tree.mu(2,1) / q) +
      log(tree.mu(1,2) / q) +
      log(tree.mu(2,2) / q) +
      log(tree.mu(1,3) / q) +
      log(tree.mu(2,3) / q)) +
     (log(tree.md(1,1) / q) +
      log(tree.md(2,1) / q) +
      log(tree.md(1,2) / q) +
      log(tree.md(2,2) / q) +
      log(tree.md(1,3) / q) +
      log(tree.md(2,3) / q))
     / 9.0 +
     (log(tree.me(1,1) / q) +
      log(tree.me(2,1) / q) +
      log(tree.me(1,2) / q) +
      log(tree.me(2,2) / q) +
      log(tree.me(1,3) / q) +
      log(tree.me(2,3) / q)) / 3.0
     + (log(fabs(tree.mch(1)) / q)
        + log(fabs(tree.mch(2)) / q))
     * 4.0 / 3.0);

  return alphaEm / (1.0 - deltaAlpha);
}


// Prints out what the lsp is
void recogLsp(int temp, int posj) {
  switch(temp) {
  case 0: cout << "neutralino"; break;
  case 1:
    switch(posj)
      {
      case 3: cout << "stop"; break;
      case 2: cout << "scharm"; break;
      case 1: cout << "sup"; break;
      } break;
  case 2:
    switch(posj)
      {
      case 3: cout << "sbottom"; break;
      case 2: cout << "sstange"; break;
      case 1: cout << "sdown"; break;
      } break;
  case 3:
    switch(posj)
      {
      case 3: cout << "stau"; break;
      case 2: cout << "smu"; break;
      case 1: cout << "selectron"; break;
      } break;
  case 4: cout << "chargino"; break;
  case 5: cout << "sneutrino"; break;
  case 6: cout << "gluino"; break;
  default:
    ostringstream ii;
    ii << "Wrong input to lsp printing routine\n";
    throw ii.str(); break;
  }

}

// Returns lsp mass in mass and function return labels which particle is lsp:
// 0 is neutralino posi = #, posj = 0
int MssmSoftsusy::lsp(double & mass, int & posi, int & posj) const {
  int temp = 0, pos1 = 0, pos2 = 0;
  sPhysical s(displayPhys());

  double minmass = fabs(s.mneut(1)); posi = pos1; posj = 0;

  // up squarks 1
  double lightest = s.mu.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) {
    minmass = lightest; posi = pos1; posj = pos2; temp = 1; }

  // down squarks 2
  lightest = s.md.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) {
    minmass = lightest; posi = pos1; posj = pos2; temp = 2; }

  // sleptons 3
  lightest = s.me.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) {
    minmass = lightest; posi = pos1; posj = pos2; temp = 3; }

  // charginos 4
  lightest = fabs(s.mch(1));
  if (lightest < minmass) {
    minmass = lightest; posi = pos1; posj = 0; temp = 4; }

  // sneutrinos 5
  lightest = s.msnu.apply(fabs).min(pos1);
  if (lightest < minmass) {
    minmass = lightest; posi = pos1; posj = 0; temp = 5; }

  // gluino 6
  lightest = s.mGluino;
  if (lightest < minmass) {
    minmass = lightest; posi = 0; posj = 0; temp = 6; }

  mass = minmass;
  return temp;
}

double MssmSoftsusy::maxMass() const {
  int pos1, pos2;
  sPhysical s(displayPhys());

  double maxmass = s.mneut.apply(fabs).min(pos1);

  // up squarks 1
  double heaviest = s.mu.apply(fabs).max(pos1, pos2);
  if (heaviest > maxmass) maxmass = heaviest;

  // down squarks 2
  heaviest = s.md.apply(fabs).max(pos1, pos2);
  if (heaviest > maxmass) maxmass = heaviest;

  // sleptons 3
  heaviest = s.me.apply(fabs).max(pos1, pos2);
  if (heaviest > maxmass) maxmass = heaviest;

  // charginos 4
  heaviest = s.mch.apply(fabs).max();
  if (heaviest > maxmass) maxmass = heaviest;

  // sneutrinos 5
  heaviest = s.msnu.apply(fabs).max();
  if (heaviest > maxmass) maxmass = heaviest;

  // gluino 6
  heaviest = s.mGluino;
  if (heaviest > maxmass) maxmass = heaviest;

  return maxmass;
}

double MssmSoftsusy::calcMs() const {

  drBarPars tree(displayDrBarPars());

  if (QEWSB < MZ)
    return maximum(QEWSB * sqrt(tree.mu(2, 3) * tree.mu(1, 3)), displayMz());
  else return QEWSB;
}

// Provides the first guess at a SUSY object at mt, inputting tanb and oneset
// (should be at MZ) - it's very crude, doesn't take radiative corrections
// into account etc.
MssmSusy guessAtSusyMt(double tanb, const QedQcd & oneset) {
  // This bit gives a guess at a SUSY object
  QedQcd leAtMt(oneset);

  DoubleVector a(3), g(3);
  double sinth2 = 1.0 - sqr(MW / MZ);
  // Gauge couplings at mt
  a = leAtMt.getGaugeMu(oneset.displayPoleMt(), sinth2);

  MssmSusy t;
  t.setTanb(tanb);
  double beta = atan(tanb);

  // Yukawa couplings
  double vev = 246.22;
  double ht = (leAtMt.displayMass(mTop) - 30.0) * root2
    / (vev * sin(beta));
  double hb =  ht * tanb * leAtMt.displayMass(mBottom) /
    leAtMt.displayMass(mTop);
  double htau =  hb * leAtMt.displayMass(mTau) /
    leAtMt.displayMass(mBottom);

  const double pi = 3.141592653589793;
  int i; for(i=1; i<=3; i++) g(i) = sqrt(4.0 * pi * a(i));
  t.setAllGauge(g);

  t.setYukawaElement(YU, 3, 3, ht);
  t.setYukawaElement(YD, 3, 3, hb);
  t.setYukawaElement(YE, 3, 3, htau);
  t.setHvev(vev);

  t.setMu(oneset.displayPoleMt());

  return t;
}


// Returns low energy softsusy object consistent with BC's m0 etc at MGUT.
// oneset should be at MZ and contains the SM data to fit the model to.
// If the running comes into difficulty, eg if a Landau pole is reached, it
// returns a ZERO object: no result is possible!
// Boundary condition is the theoretical condition on parameters at the high
// energy scale mx: the parameters themselves are contained within the vector.
double MssmSoftsusy::lowOrg
(void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &),
 double mxGuess,
 const DoubleVector & pars, int sgnMu, double tanb, const QedQcd &
 oneset, bool gaugeUnification, bool ewsbBCscale) {

  double mx = 0.0;

  try {

    initialise();

    double muFirst = displaySusyMu(); // Remember initial value

    const static MssmSoftsusy empty;
    setSoftsusy(empty); // Always starts from an empty object
    setMw(MW);
    double mz = displayMz();

    // Here all was same
    if (mxGuess > 0.0)
      mx = mxGuess;
    else {
      string ii("Trying to use negative mx in MssmSoftsusy::lowOrg.\n");
      ii = ii + "Now illegal! Use positive mx for first guess of mx.\n";
      throw ii;
    }

    if (oneset.displayMu() != mz) {
      cout << "WARNING: lowOrg in softsusy.cpp called at " << displayMu() <<
        "with oneset at scale\n" <<
        oneset.displayMu() << "\ninstead of " << mz << endl;
    }

    const int maxtries = 50;
    double tol = TOLERANCE;

    MssmSusy t(guessAtSusyMt(tanb, oneset));

    t.setLoops(2); // 2 loops should protect against ht Landau pole
    t.runto(mx);

    setSusy(t);
    setData(oneset);

    // Initial guess: B=0, mu=1st parameter, need better guesses
    boundaryCondition(*this, pars);
    if (sgnMu == 1 || sgnMu == -1) {
      setSusyMu(sgnMu * 1.0);
      setM3Squared(0.);
    }
    else {
      setSusyMu(muFirst);
      setM3Squared(muFirst);
    }

    run(mx, mz);

    if (sgnMu == 1 || sgnMu == -1) rewsbTreeLevel(sgnMu);

    physical(0);

    setThresholds(3); setLoops(2);

    itLowsoft(maxtries, mx, sgnMu, tol, tanb, boundaryCondition, pars,
              gaugeUnification, ewsbBCscale);

    if (displayProblem().nonperturbative
        || displayProblem().higgsUfb || displayProblem().tachyon
        || displayProblem().noRhoConvergence)
      return mx;

    runto(maximum(displayMsusy(), MZ));
    if (ewsbBCscale) boundaryCondition(*this, pars);
    physical(3);
    runto(mz);

    if (PRINTOUT) cout << " end of iteration" << endl;

    // check to see if there are any negative mass squared scalars: if there
    // are, flag scalarMinimum
    int posi, posj, lspCode; double minmass;
    lspCode = lsp(minmass, posi, posj);
    if (lspCode > 0 && lspCode < 6 && lspCode !=4 && minmass < 0.0) {
      flagTachyon(true);
      if (PRINTOUT > 2) recogLsp(posi, posj); cout << " tachyonic ";
    }

  }
  catch(const char *a) {
    cerr << "Problem: " << a << " pars=" << pars << " tanb=" << tanb
         << " oneset=" << oneset <<
      endl;
  }
  catch(const string & a) {
    cerr << "Problem: " << a << " pars=" << pars << " tanb=" << tanb
         << " oneset=" << oneset << endl;
  }

  return mx;
}


// You should evaluate this at a scale MSusy average of stops.
// Depth of electroweak minimum: hep-ph/9507294. Bug-fixed 19/11/04
double MssmSoftsusy::realMinMs() const {
  MssmSusy temp(displaySusy());
  temp.runto(calcMs(), TOLERANCE);

  double beta = atan(temp.displayTanb());
  return - 1.0 / 32.0 *
    (0.6 * sqr(temp.displayGaugeCoupling(1)) +
     sqr(temp.displayGaugeCoupling(2))) *
    sqr(sqr(displayHvev()) * cos(2.0 * beta));
}

// Difference between two SOFTSUSY objects in and out: EWSB terms only
double sumTol(const MssmSoftsusy & in, const MssmSoftsusy & out) {

  drBarPars inforLoops(in.displayDrBarPars()),
    outforLoops(out.displayDrBarPars());

  DoubleVector sT(32);
  int k = 1;
  int i; for (i=1; i<=4; i++) {
    double sTin  = fabs(inforLoops.mhiggs(i));
    double sTout = fabs(outforLoops.mhiggs(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k = k + 1;
  }
  for (i=1; i<=3; i++) {
    double sTin  = fabs(inforLoops.msnu(i));
    double sTout = fabs(outforLoops.msnu(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k = k + 1;
  }
  for (i=1; i<=2; i++) {
    double sTin = fabs(inforLoops.mch(i));
    double sTout = fabs(outforLoops.mch(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k = k + 1;
  }
  for (i=1; i<=4; i++) {
    double sTin = fabs(inforLoops.mneut(i));
    double sTout = fabs(outforLoops.mneut(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k = k + 1;
  }
  double sTin = fabs(inforLoops.mGluino);
  double sTout = fabs(outforLoops.mGluino);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
  k = k + 1;
  int j; for (j=1; j<=3; j++)
    for(i=1; i<=2; i++) {
      sTin = fabs(inforLoops.mu(i, j));
      sTout = fabs(outforLoops.mu(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k = k + 1;
      sTin = fabs(inforLoops.md(i, j));
      sTout = fabs(outforLoops.md(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k = k + 1;
      sTin = fabs(inforLoops.me(i, j));
      sTout = fabs(outforLoops.me(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k = k + 1;
    }
  return sT.max();
}

/*
// Yukawa coupling at any scale - should give it vev at correct scale
double MssmSoftsusy::calcHt(double vev) {

  double mtrun = calcRunningMt();
  double beta = atan(displayTanb());
  double ht = mtrun * root2 / (vev * sin(beta));

  return ht;
}

double MssmSoftsusy::calcHb(double vev) const {
  double mbrun = calcRunningMb();
  double beta = atan(displayTanb());
  double hb = mbrun * root2 / (vev * cos(beta));
  return hb;
}

double MssmSoftsusy::calcHtau(double vev) const {
  double mtaurun = calcRunningMtau();
  double beta = atan(displayTanb());
  double htau = mtaurun * root2 / (vev * cos(beta));
  return htau;
}
*/

// Calculates sin theta at the current scale
double MssmSoftsusy::calcSinthdrbar() const {

  double tanth = displayGaugeCoupling(1) /
    displayGaugeCoupling(2) * sqrt(0.6);
  return sin(atan(tanth));
}

//VEV at current scale, using an input value of Z self-energy
double MssmSoftsusy::getVev(double pizzt) const {

  double vsquared = 4.0 * (sqr(displayMz()) + pizzt) /
    (sqr(displayGaugeCoupling(2)) +
     sqr(displayGaugeCoupling(1)) * 0.6);

  if (vsquared < 200.0 || testNan(vsquared))
    return 246.22;

  return sqrt(vsquared);
}

//VEV at current scale: calculates Z self energy first
double MssmSoftsusy::getVev() const {
  double pizzt = piZZT(displayMz(), displayMu());

  return getVev(pizzt);
}

// It'll set the important SUSY couplings: supposed to be applied at MZ
// You should set up an iteration here since Yuk's depend on top mass which
// depends on Yuk's etc.
void MssmSoftsusy::sparticleThresholdCorrections(double tb) {
  double mz = displayMz();  if (displayMu() != mz) {
    ostringstream ii;
    ii << "Called MssmSoftsusy::sparticleThresholdCorrections "
       << "with scale" << displayMu() << endl;
    throw ii.str();
  }

  setTanb(tb);
  calcDrBarPars(); // for the up-coming loops
  setHvev(getVev());

  double alphaMsbar = dataSet.displayAlpha(ALPHA);
  double alphaDrbar = qedSusythresh(alphaMsbar, displayMu());

  double alphasMZDRbar =
    qcdSusythresh(displayDataSet().displayAlpha(ALPHAS), displayMu());

  // Do gauge couplings
  double outrho = 1.0, outsin = 0.48, tol = TOLERANCE * 1.0e-5;
  int maxTries = 20;
   double pizztMZ = piZZT(displayMz(), displayMu(), true);
   double piwwt0  = piWWT(0., displayMu(), true);
   double piwwtMW = piWWT(displayMw(), displayMu(), true);
   rhohat(outrho, outsin, alphaDrbar, pizztMZ, piwwt0, piwwtMW, tol, maxTries);
   if (problem.noRhoConvergence)
     outsin = sqrt(1.0 - sqr(displayMw() / displayMz()));

  double eDR = sqrt(4.0 * PI * alphaDrbar), costhDR = cos(asin(outsin));

  DoubleVector newGauge(3);
  newGauge(1) = eDR / costhDR * sqrt(5.0 / 3.0);
  newGauge(2) = eDR / outsin;
  newGauge(3) = sqrt(4.0 * PI * alphasMZDRbar);

  double vev = displayHvev();

  setGaugeCoupling(1, newGauge(1));
  setGaugeCoupling(2, newGauge(2));
  setGaugeCoupling(3, newGauge(3));

  setMw(sqrt(0.25 * sqr(newGauge(2)) * sqr(vev) - piWWT(displayMw(), displayMu())));

  if (MIXING < -1 || MIXING > 2) {
    ostringstream ii;
    ii << "In MssmSoftsusy::sparticleThresholdCorrections(double tb) ";
    ii << "\n MIXING=" << MIXING << " is out of range (-1 -> 2)\n";
    throw ii.str();
  }

  double beta = atan(displayTanb());
  DoubleMatrix mUq(3, 3), mDq(3, 3), mLep(3, 3);
  massFermions(displayDataSet(), mDq, mUq, mLep);

  mDq(3, 3) = calcRunningMb();
  mUq(3, 3) = calcRunningMt();
  mLep(3, 3) = calcRunningMtau();

  // 3-family mixed-up Yukawa couplings: From PDG 2000
  static DoubleMatrix ckm(3, 3);
  if (ckm(1, 1) == 0.0) {
    double s12 = 0.2225, s23 = 0.04, s13 = 0.0035;
    double c12 = sqrt(1.0 - sqr(s12)), c23 = sqrt(1.0 - sqr(s23)),
      c13 = sqrt(1.0 - sqr(s13));
    ckm(1, 1) = c12 * c13;
    ckm(1, 2) = s12 * c13;
    ckm(1, 3) = s13;
    ckm(2, 1) = -s12 * c23 - c12 * s23 * s13;
    ckm(2, 2) = c12 * c23 - s12 * s23 * s13;
    ckm(2, 3) = s23 * c13;
    ckm(3, 1) = s12 * s23 - c12 * c23 * s13;
    ckm(3, 2) = -c12 * s23 - s12 * c23 * s13;
    ckm(3, 3) = c23 * c13;
  }
  doQuarkMixing(ckm, mDq, mUq);

  setYukawaMatrix(YU, mUq * (root2 / (vev * sin(beta))));
  setYukawaMatrix(YD, mDq * (root2 / (vev * cos(beta))));
  setYukawaMatrix(YE, mLep * (root2 / (vev * cos(beta))));
}

double MssmSoftsusy::displayMzRun() const {
  return displayHvev() * 0.5 *
    sqrt(sqr(displayGaugeCoupling(2)) + 0.6 * sqr(displayGaugeCoupling(1)));
}

// calculates masses all at tree-level in the DRbar scheme, useful for
// radiative corrections.
void MssmSoftsusy::calcDrBarPars() {
  drBarPars eg;
    // First, must define mstop,sbot,stau and mixing angles in DRbar scheme
  double beta = atan(displayTanb()), mzPole = displayMz();
  double sinthDRbar = calcSinthdrbar();
  double mz = displayMzRun(), mz2 = sqr(mz);
  double pizzt = sqr(mz) - sqr(mzPole);

  sw2 = sqr(sinthDRbar);
  guL = 0.5 - 2.0 * sw2 / 3.0;
  gdL = -0.5 + sw2 / 3.0;
  geL = -0.5 + sw2;
  guR = 2.0 * sw2 / 3.0;
  gdR = -sw2 / 3.0;
  geR = -sw2;

  double vev = displayHvev();

  eg.mt   = displayYukawaElement(YU, 3, 3) * vev * sin(beta) / root2;
  eg.mb   = displayYukawaElement(YD, 3, 3) * vev * cos(beta) / root2;
  eg.mtau = displayYukawaElement(YE, 3, 3) * vev * cos(beta) / root2;

  eg.mGluino = displayGaugino(3);

  if (pizzt == 0.) {
    flagTachyon(true);
    if (PRINTOUT > 1) cout << " MZ tachyonic ";
  }

  DoubleVector mSq(2);
  int family; for(family = 1; family <= 3; family++) {

  DoubleMatrix mSquared(2, 2);
  treeUpSquark(mSquared, eg.mt, pizzt, sinthDRbar, family);
  mSq = mSquared.sym2by2(eg.thetat);
  DoubleVector mstopDRbar(mSq.apply(ccbSqrt));
  if (mSq(1) < 0. || mSq(2) < 0.) {
    flagTachyon(true);
    if (PRINTOUT > 1) cout << " tree sup(" << family << ") tachyon ";
  }

  treeDownSquark(mSquared, eg.mb, pizzt, sinthDRbar, family);
  mSq = mSquared.sym2by2(eg.thetab);
  DoubleVector msbotDRbar(mSq.apply(ccbSqrt));
  if (mSq(1) < 0. || mSq(2) < 0.) {
    flagTachyon(true);
    if (PRINTOUT > 1) cout << " tree sdown(" << family << ") tachyon ";
  }

  treeChargedSlepton(mSquared, eg.mtau, pizzt, sinthDRbar, family);
  mSq = mSquared.sym2by2(eg.thetatau);
  DoubleVector mstauDRbar(mSq.apply(ccbSqrt));
  if (mSq(1) < 0. || mSq(2) < 0.) {
    flagTachyon(true);
    if (PRINTOUT > 1) cout << " tree selectron(" << family << ") tachyon ";
  }

  int i; for (i=1; i<=2; i++) {
    eg.mu(i, family) = mstopDRbar(i);    eg.md(i, family) = msbotDRbar(i);
    eg.me(i, family) = mstauDRbar(i);
    }
  double mSnuSquared;
  treeSnu(mSnuSquared, pizzt, family);
  if (mSnuSquared < 0.) {
    flagTachyon(true);
    if (PRINTOUT > 1) cout << " tree sneutrino(" << family << ") tachyon ";
  }
  eg.msnu(family) = ccbSqrt(mSnuSquared);
  }

  DoubleMatrix mCh(2, 2);
  double mw = displayMwRun();
  double mw2 = sqr(mw);
  mCh(1, 1) = displayGaugino(2);
  mCh(2, 1) = root2 * mw * cos(beta);
  mCh(1, 2) = mCh(2, 1) * displayTanb();
  mCh(2, 2) = displaySusyMu();
  eg.mch = mCh.asy2by2(eg.thetaL, eg.thetaR);

  DoubleMatrix mNeut(4, 4);
  mNeut(1, 1) = displayGaugino(1);
  mNeut(2, 2) = displayGaugino(2);
  mNeut(1, 3) = - mz * cos(beta) * sinthDRbar;
  mNeut(1, 4) = - mNeut(1, 3) * displayTanb();
  mNeut(2, 3) = mw * cos(beta);
  mNeut(2, 4) = - mNeut(2, 3) * displayTanb();
  mNeut(3, 4) = - displaySusyMu();
  mNeut.symmetrise();
  if (mNeut.diagonaliseSym(eg.mixNeut, eg.mneut) > TOLERANCE *
      1.0e-3) {
    ostringstream ii;
    ii << "accuracy bad in neutralino diagonalisation"<< flush;
    throw ii.str();
    }

  double mAsq = displayM3Squared() / (sin(beta) * cos(beta));

  if (mAsq < 0.0) {
    // If it's for the EWSB BC at MZ, we simply use the one-loop corrected
    // mass for mA instead. You could indeed make this option permanent (ie
    // even for mAsq > 0)
    if (displayMu() == MZ) {
      double mApole = physpars.mhiggs(2); // physical value
      double piaa = piAA(mApole, displayMu());
      double t1Ov1 = doCalcTadpole1oneLoop(eg.mt, sinthDRbar),
        t2Ov2 = doCalcTadpole2oneLoop(eg.mt, sinthDRbar);
      double poleMasq =
        (displayMh2Squared() - t2Ov2 -
         displayMh1Squared() + t1Ov1) /
        cos(2.0 * beta) - mz2 - piaa +
        sqr(sin(beta)) * t1Ov1 + sqr(cos(beta)) * t2Ov2;

      mAsq = poleMasq;

      if (mAsq < 0.) { flagTachyon(true); mAsq = fabs(poleMasq); }
    }
    else {
      flagTachyon(true);
      if (PRINTOUT > 1) cout << " mA^2(tree)=" << mAsq << " since m3sq=" <<
                          displayM3Squared() << " @ "<< displayMu() <<
                          " " << endl;
      mAsq = fabs(mAsq);
    }
  }

  DoubleMatrix mH(2, 2);
  mH(1, 1) = mAsq * sqr(sin(beta)) + mz2 * sqr(cos(beta));
  mH(1, 2) = - sin(beta) * cos(beta) * (mAsq + mz2);
  mH(2, 2) = mAsq * sqr(cos(beta)) + mz2 * sqr(sin(beta));
  mH(2, 1) = mH(1 ,2);
  mSq = mH.sym2by2(eg.thetaH);
  if (mSq(1) < 0. || mSq(2) < 0.) {
    flagTachyon(true);
    if (PRINTOUT > 1) cout << " mH/h tachyon ";
  }
  DoubleVector temp(mSq.apply(ccbSqrt));
  if (temp(2) > temp(1)) eg.thetaH = eg.thetaH + PI * 0.5;

  int pos;
  eg.mhiggs(1) = temp.min(pos); eg.mhiggs(2) = temp.max();
  eg.mhiggs(3) = sqrt(mAsq); eg.mhiggs(4) = sqrt(mAsq + mw2);

  eg.mw = mw;
  eg.mz = mz;
  eg.mpzNeutralinos();
  eg.mpzCharginos();

  setDrBarPars(eg);

  return;
}

void MssmSoftsusy::itLowsoft
(int maxTries, double & mx, int sgnMu, double tol, double tanb,
 void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &),
 const DoubleVector & pars, bool gaugeUnification, bool ewsbBCscale) {

  static MssmSoftsusy old;
  static int numTries = 0;
  double mz = displayMz();

  if (numTries != 0 && sqr(displayMu() / mz - 1.0) > TOLERANCE) {
    cout << "WARNING: itLowsoft called at inappropriate";
    cout << " scale:" << displayMu() << endl;
    cout << "whereas it should be " << mz << endl;
  }

  if (numTries - 1 > maxTries) {// Iterating too long: bail out
    flagNoConvergence(true);
    if (PRINTOUT) cout << "itLowsoft reached maxtries\n";
    numTries = 0; //old = MssmSoftsusy();
    return;
  }
  double c = sumTol(*this, old);

  if (numTries !=0 && c < tol) {// Accuracy achieved: bail out
    numTries = 0; // Reset the number of iterations for the next time
    if (PRINTOUT > 0) cout << " sT=" << c << " " << flush;
    if (displayProblem().test() && PRINTOUT > 0)
      cout << " ***problem point***: " << displayProblem() << ".";
    //    old = MssmSoftsusy();
    return;
  }

  if (PRINTOUT > 1) cout << displayProblem();

  // All problems should be reset since only the ones of the final iteration
  // should count (sometimes problems disappear). This can mean that problems
  // only show up as no rho convergence....
  flagAllProblems(false);

  double mtpole, mtrun;

  mtpole = displayDataSet().displayPoleMt();
  // On first iteration, don't bother with finite corrections

  numTries = numTries + 1;

  // Old iteration is 'old': these are the parameters by which convergence is
  // measured.
  old.setDrBarPars(displayDrBarPars());

  try {
    sparticleThresholdCorrections(tanb);

  if (problem.noRhoConvergence) {
    if (PRINTOUT) cout << "No convergence in rhohat\n";
  }

  // precision of running/RGE integration: start off low and increase
  double eps = maximum(exp(double(- numTries) * log(10.0)), tol * 0.01);

  // first stab at MSUSY: root(mstop1(MZ) mstop2(MZ))
  if (numTries == 1) setMsusy(calcMs());

  int err = run(mz, displayMsusy(), eps);

  if (err) {
    // problem with running: bail out
    if (PRINTOUT) cout
      << "itLowsoft gone non-perturbative on way to MX(u)\n";
    if (PRINTOUT > 1) printObj();
    numTries = 0;
    return;
  }

  mtrun = mtpole - 20.0;

  if (numTries != 1) {
    double beta = atan(displayTanb());
    double vev = displayHvev();
    mtrun = displayYukawaElement(YU, 3, 3) * vev * sin(beta) /
      root2;
  }

  double MZp = mz;
  double tb; MZp = sqrt(fabs(predMzsq(tb)));
  // if predicted MZ is ridiculous, bail out to avoid further problems
  if (fabs(MZp) > 1.0e10) {
    if (PRINTOUT > 2) cout << "MZ predicted too high";
    flagTachyon(true); numTries = 0;
  //old = MssmSoftsusy();
  return; }

  if (PRINTOUT > 0) {
    cout << "\nsT=" << c << " mu=" << displaySusyMu() <<  " m3sq=" <<
      displayM3Squared() << " MWp=" << displayMw() << " MZp=" << MZp << flush;
  }

  err = runto(mx, eps);

  // Guard against the top Yukawa fixed point
  if (displayYukawaElement(YU, 3, 3) > 3.0
      || displayYukawaElement(YD, 3, 3) > 3.0
      || displayYukawaElement(YE, 3, 3) > 3.0) {
    setYukawaElement(YU, 3, 3, minimum(3.0, displayYukawaElement(YU, 3, 3)));
    setYukawaElement(YD, 3, 3, minimum(3.0, displayYukawaElement(YD, 3, 3)));
    setYukawaElement(YE, 3, 3, minimum(3.0, displayYukawaElement(YE, 3, 3)));
    flagIrqfp(true); // will do for now, but you need it as an option!
  }

  if (err) {
    // problem with running: bail out
    flagNonperturbative(true);
    if (PRINTOUT)
      cout << "itLowsoft gone non-perturbative approaching mgut\n";
    if (PRINTOUT > 1) printObj();
    numTries = 0; // old = MssmSoftsusy();
    return;
  }

  if (gaugeUnification) {
    sBrevity (dummy);
    MssmSusy a(this -> MssmSusy::beta(dummy));

    // Equal gauge couplings BC -- linear approximation
    mx = mx * exp((displayGaugeCoupling(2) - displayGaugeCoupling(1))
                  / (a.displayGaugeCoupling(1) - a.displayGaugeCoupling(2)));

    // tachyons can cause a problem with mx
    if (mx < 1.0e16) {
      mx = 1.0e16;
      if (PRINTOUT > 2) cout << " mx too low ";
      flagTachyon(true);
     }
    if (mx > 1.0e17) {
      if (PRINTOUT > 2) cout << " mx =" << mx <<" too high ";
      mx = 1.0e17;
      flagTachyon(true);
    }
  }

  boundaryCondition(*this, pars);

  if (!ewsbBCscale) err = runto(displayMsusy(), eps);
  if (err) {
    // problem with running: bail out
    flagNonperturbative(true);
    if (PRINTOUT) cout << "itLowsoft gone non-perturbative on way to MZ\n";
    if (PRINTOUT > 1) printObj();
    //    old = MssmSoftsusy();
    numTries = 0;
    return;
  }

  calcDrBarPars();
  setMsusy(calcMs());
  if (ewsbBCscale) mx = displayMsusy();
  if (PRINTOUT > 0) cout << " mgut=" << mx << flush;

  rewsb(sgnMu, mtrun);

  if (problem.noMuConvergence) {
    if (PRINTOUT)
      cout << "itLowsoft doesn't break EWSB\n";
    if (PRINTOUT > 1) printObj();
    numTries = 0; //old = MssmSoftsusy();
    return;
  }

  err = runto(mz, eps);
  if (err) {
    // problem with running: bail out
    flagNonperturbative(true);
    if (PRINTOUT) cout << "itLowsoft gone non-perturbative on way to MZ\n";
    if (PRINTOUT > 1) printObj();
    //    old = MssmSoftsusy();
    numTries = 0;
    return;
  }

  itLowsoft(maxTries, mx, sgnMu, tol, tanb, boundaryCondition, pars,
            gaugeUnification, ewsbBCscale);
  }
  catch(const char *a) {
    numTries = 0;
    throw a;
  }
  catch(const string &a) {
    numTries = 0;
    throw a;
  }
}

// Transverse part of Z self-energy: has been checked
double MssmSoftsusy::piZZT(double p, double q, bool usePoleMt) const {

  drBarPars tree(displayDrBarPars());

  // fermions: these are valid at MZ
  double    mtop =  tree.mt;
  // We utilise pole mt for these corrections (which the 2-loop Standard Model
  // pieces assume)
  if (usePoleMt) mtop = displayDataSet().displayPoleMt();

  double    mb   =  tree.mb;
  double    mtau =  tree.mtau;
  double    ms   =  displayDataSet().displayMass(mStrange) ;
  double    mc   =  displayDataSet().displayMass(mCharm) ;
  double    mmu  =  displayDataSet().displayMass(mMuon) ;
  double    mE  =  displayDataSet().displayMass(mElectron) ;
  double    mD  =  displayDataSet().displayMass(mDown) ;
  double    mU  =  displayDataSet().displayMass(mUp);
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    g       = displayGaugeCoupling(2);
  double    alpha   = tree.thetaH ;
  double    beta    = atan(displayTanb());
  double    thetat = tree.thetat ;
  double    thetab = tree.thetab;
  double    thetatau= tree.thetatau ;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau) ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    ct2     = sqr(ct) ;
  double    cb2     = sqr(cb) ;
  double    ctau2   = sqr(ctau) ;
  double    st2     = (1.0 - ct2);
  double    sb2     = (1.0 - cb2);
  double    stau2   = (1.0 - ctau2);
  double    mz      = displayMzRun();

  double rhs = 0.0;

  double higgs1 = 0.0, higgs2 = 0.0, higgs3 = 0.0, charginos = 0.0,
    neutralinos = 0.0, squarks = 0.0, thirdFamily = 0.0, sneutrinos = 0.0;

  higgs1 =
    - sqr(sin(alpha - beta)) *
    (b22bar(p, tree.mhiggs(2), tree.mhiggs(3), q) +
     b22bar(p, mz, tree.mhiggs(1), q) -
     sqr(mz) * b0(p, mz, tree.mhiggs(1), q));

  higgs2 =
    - sqr(cos(alpha - beta)) *
    (b22bar(p, mz, tree.mhiggs(3), q) +
     b22bar(p, tree.mhiggs(2), tree.mhiggs(1), q) -
     sqr(mz) * b0(p, mz, tree.mhiggs(3), q));

  higgs3 =
    - 2.0 * sqr(cw2DRbar) * (2 * sqr(p) + sqr(displayMwRun()) - sqr(mz) *
                             sqr(sw2DRbar) / cw2DRbar)
    * b0(p, displayMwRun(), displayMwRun(), q)
    - (8.0 * sqr(cw2DRbar) + sqr(cos(2.0 * thetaWDRbar))) *
    b22bar(p, displayMwRun(), displayMwRun(), q) -
    sqr(cos(2.0 * thetaWDRbar)) * b22bar(p, tree.mhiggs(4),
                                         tree.mhiggs(4), q);

  //static
  DoubleVector vu(2), vd(2), ve(2);
  //static
  DoubleMatrix vt(2, 2), vb(2, 2), vtau(2, 2);

  vu(1) = guL;
  vu(2) = guR;

  vt(1, 1) = guL * ct2 - guR * st2;
  vt(1, 2) = (guL + guR) * ct * st;
  vt(2, 2) = guR * ct2 - guL * st2;
  vt(2, 1) = vt(1, 2);

  vd(1) = gdL;
  vd(2) = gdR;

  vb(1, 1) = gdL * cb2 - gdR * sb2;
  vb(1, 2) = (gdL + gdR) * cb * sb;
  vb(2, 2) = gdR * cb2 - gdL * sb2;
  vb(2, 1) = vb(1, 2);

  ve(1) = geL;
  ve(2) = geR;

  vtau(1, 1) = geL * ctau2 - geR * stau2;
  vtau(1, 2) = (geL + geR) * ctau * stau;
  vtau(2, 2) = geR * ctau2 - geL * stau2;
  vtau(2, 1) = vtau(1, 2);

  int i, j, family;
  for (family = 1; family<=2; family++) {
    // sneutrinos
    sneutrinos = sneutrinos -
      b22bar(p, tree.msnu(family), tree.msnu(family), q);
    for (i=1; i<=2; i++) {
      // up squarks
      squarks = squarks - 12.0 * sqr(vu(i)) *
        b22bar(p, tree.mu(i, family),
               tree.mu(i, family), q);
      // down squarks
      squarks = squarks - 12.0 * sqr(vd(i)) *
        b22bar(p, tree.md(i, family),
               tree.md(i, family), q);
      // selectrons
      squarks = squarks - 4.0 * sqr(ve(i)) *
        b22bar(p, tree.me(i, family),
               tree.me(i, family), q);
    }
  }
  family = 3;

  // THIRD FAMILY
  // sneutrinos
  thirdFamily = thirdFamily -
    b22bar(p, tree.msnu(family),
           tree.msnu(family), q);
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      // up squarks
      thirdFamily = thirdFamily - 12.0 * sqr(vt(i, j)) *
        b22bar(p, tree.mu(i, family),
               tree.mu(j, family), q);

      // down squarks
      thirdFamily = thirdFamily - 12.0 * sqr(vb(i, j)) *
        b22bar(p, tree.md(i, family),
               tree.md(j, family), q);
      // selectrons
       thirdFamily = thirdFamily - 4.0 * sqr(vtau(i, j)) *
        b22bar(p, tree.me(i, family),
               tree.me(j, family), q);
    }

  double quarks = 0.0;

  quarks = quarks + 3.0 * hfn(p, mU, mU, q) * (sqr(guL) + sqr(guR));
  quarks = quarks + 3.0 * hfn(p, mc, mc, q) * (sqr(guL) + sqr(guR));
  quarks = quarks + 3.0 * (hfn(p, mtop, mtop, q) *
                     (sqr(guL) + sqr(guR)) - 4.0 * guL * guR * sqr(mtop) *
                     b0(p, mtop, mtop, q));

  quarks = quarks + 3.0 * hfn(p, mD, mD, q) * (sqr(gdL) + sqr(gdR));
  quarks = quarks + 3.0 * hfn(p, ms, ms, q) * (sqr(gdL) + sqr(gdR));
  quarks = quarks + 3.0 * (hfn(p, mb, mb, q) * (sqr(gdL) + sqr(gdR))
                           - 4.0 * gdL * gdR * sqr(mb) * b0(p, mb, mb, q));

  quarks = quarks + hfn(p, mE, mE, q) * (sqr(geL) + sqr(geR));
  quarks = quarks + hfn(p, mmu, mmu, q) * (sqr(geL) + sqr(geR));
  quarks = quarks + hfn(p, mtau, mtau, q) *
    (sqr(geL) + sqr(geR)) - 4.0 * geL * geR
    * sqr(mtau) * b0(p, mtau, mtau, q);

  quarks = quarks + 3.0 * hfn(p, 0., 0., q) * 0.25;

  // Neutralinos
  //static
  ComplexMatrix aPsi(4, 4), bPsi(4, 4), aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  aPsi(3, 3) = g / (2.0 * cos(thetaWDRbar)); aPsi(4, 4) = 0. - aPsi(3, 3);
  bPsi = 0. - aPsi;

  aChi = n.complexConjugate() * aPsi * n.transpose();
  bChi = n * bPsi * n.hermitianConjugate();

  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      neutralinos = neutralinos + cw2DRbar / (2.0 * sqr(g)) *
        ((sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod())) *
         hfn(p, mneut(i), mneut(j), q)
         + 4.0 * (bChi(i, j).conj() * aChi(i, j)).real() *
         mneut(i) * mneut(j) * b0(p, mneut(i), mneut(j), q));
    }

  // Charginos
  //  static
  ComplexMatrix aPsiCh(2, 2), aCh(2, 2), bCh(2, 2);
  aPsiCh(1, 1) = g * cos(thetaWDRbar);
  aPsiCh(2, 2) = g * cos(2.0 * thetaWDRbar) / (2.0 * cos(thetaWDRbar));

  aCh = v.complexConjugate() * aPsiCh * v.transpose();
  bCh = u * aPsiCh * u.hermitianConjugate();

  for (i=1; i<=2; i++)
    for(j=1; j<=2; j++) {
      charginos = charginos + cw2DRbar / sqr(g) *
        ((sqr(aCh(i, j).mod()) + sqr(bCh(i, j).mod())) *
         hfn(p, mch(i), mch(j), q)
         + 4.0 * (bCh(i, j).conj() * aCh(i, j)).real() * mch(i) *
         mch(j) * b0(p, mch(i), mch(j), q));
    }

  rhs = higgs1 + higgs2 + higgs3 + charginos + neutralinos + squarks +
    thirdFamily + sneutrinos + quarks;

  double pi = rhs * sqr(g) / (cw2DRbar * 16.0 * sqr(PI));


  if (pi + sqr(mz) < 0.0) {
    if (PRINTOUT > 2) cout << " tachyon MZ ";
    return 0.;
  }

  return pi;
}

// W propagator to 1 loop in MSSM
// It's all been checked
double MssmSoftsusy::piWWT(double p, double q, bool usePoleMt) const {
  drBarPars tree(displayDrBarPars());

  double    mtop =  tree.mt;
  // We utilise pole mt for these corrections (which the 2-loop Standard Model
  // pieces assume)
  if (usePoleMt) mtop = displayDataSet().displayPoleMt();
  double    mb   =  tree.mb;
  double    mtau =  tree.mtau;
  double    beta    = atan(displayTanb());
  double    alpha   = tree.thetaH ;
  double    mH = tree.mhiggs(3);
  double    mh0 = tree.mhiggs(1);
  double    mHc = tree.mhiggs(4);
  double    mA = tree.mhiggs(2);
  // fermions: these are valid at MZ
  double    ms   =  displayDataSet().displayMass(mStrange) ;
  double    mc   =  displayDataSet().displayMass(mCharm) ;
  double    mmu  =  displayDataSet().displayMass(mMuon) ;
  double    mE  =  displayDataSet().displayMass(mElectron) ;
  double    mD  =  displayDataSet().displayMass(mDown) ;
  double    mU  =  displayDataSet().displayMass(mUp);
  double    thetat = tree.thetat ;
  double    thetab = tree.thetab;
  double    thetatau= tree.thetatau ;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau) ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    g       = displayGaugeCoupling(2);
  double    mz      = displayMzRun();

  double ans = 0.0;

  double higgs = 0.0;

  higgs = higgs -
    sqr(sin(alpha - beta)) *
    (b22bar(p, mH, mHc, q) + b22bar(p, mh0, displayMwRun(), q)
     - sqr(displayMwRun()) * b0(p, mh0, displayMwRun(), q));

  higgs = higgs -
    sqr(cos(alpha - beta)) *
    (b22bar(p, mh0, mHc, q) + b22bar(p, mH, displayMwRun(), q)
     - sqr(displayMwRun()) * b0(p, mH, displayMwRun(), q));

  higgs = higgs -
    b22bar(p, mA, mHc, q) -
    (1.0 + 8.0 * cw2DRbar) * b22bar(p, mz, displayMwRun(), q)
    - sw2DRbar * (8.0 * b22bar(p, displayMwRun(), 0.0, q) + 4.0 * sqr(p) *
                  b0(p, displayMwRun(), 0.0, q))
    - ((4.0 * sqr(p) + sqr(mz) + sqr(displayMwRun())) * cw2DRbar - sqr(mz)  *
       sqr(sw2DRbar)) * b0(p, mz, displayMwRun(), q);

  ans = ans + higgs;

  double fermions =
    1.5 * (hfn(p, mU, mD, q) + hfn(p, mc, ms, q) + hfn(p, mtop, mb, q)) + 0.5
    * (hfn(p, 0.0, mE, q) + hfn(p, 0.0, mmu, q) + hfn(p, 0.0, mtau, q));
  ans = ans + fermions;

  // sfermions
  double sfermions = - 6.0 *
    (b22bar(p, tree.mu(1, 1), tree.md(1, 1), q) +
     b22bar(p, tree.mu(1, 2), tree.md(1, 2), q)) -
    2.0 * (b22bar(p, tree.msnu(1), tree.me(1, 1), q) +
           b22bar(p, tree.msnu(2), tree.me(1, 2), q));
  ans = ans + sfermions;

  // stop/sbottom
  DoubleMatrix w(2, 2);
  double stopBot = 0.0;
  w(1, 1) = ct * cb; w(1, 2) = ct * sb; w(2, 1) = st * cb; w(2, 2) = st * sb;
  int i, j;
  for(i=1; i<=2; i++)
    for(j=1; j<=2; j++)
      stopBot = stopBot - 6.0 * sqr(w(i, j)) *
        b22bar(p, tree.mu(i, 3), tree.md(j, 3), q);
  ans = ans + stopBot;

  // LH slepton
  double slepton = - 2.0 *
    (sqr(ctau) *
     b22bar(p, tree.msnu(3), tree.me(1, 3), q)
     + sqr(stau) *
     b22bar(p, tree.msnu(3), tree.me(2, 3), q));
  ans = ans + slepton;

  ComplexMatrix aPsi0PsicW(4, 2), bPsi0PsicW(4, 2), aChi0ChicW(4, 2),
    bChi0ChicW(4, 2);
  DoubleMatrix fW(4, 2), gW(4, 2);

  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;
  bPsi0PsicW(3, 2) = -g / root2;

  ComplexMatrix aPsi(4, 4), bPsi(4, 4), aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  // These ought to be in physpars
  aChi0ChicW = n.complexConjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.hermitianConjugate();

  double gauginos = 0.0;

  for(i=1;i<=4;i++)
    for(j=1;j<=2;j++) {
      fW(i, j) = sqr(aChi0ChicW(i, j).mod()) + sqr(bChi0ChicW(i, j).mod());
      gW(i, j) = 2.0 * (bChi0ChicW(i, j).conj() * aChi0ChicW(i, j)).real();
      gauginos = gauginos +
        (fW(i, j) * hfn(p, mneut(i), mch(j), q)
         + 2.0 * gW(i, j) * mneut(i) * mch(j) * b0(p, mneut(i), mch(j), q))
        / sqr(g);
    }

  ans = ans + gauginos;

  double pi = ans * sqr(g) / (16.0 * sqr(PI));

  if (pi + sqr(displayMw()) < 0.0) {
    return 0.0; // tachyonic
  }

  return pi;
}

double MssmSoftsusy::pis1s1(double p, double q) const {
  drBarPars tree(displayDrBarPars());

  double    beta    = atan(displayTanb());
  double    mtau    = tree.mtau;
  double    mb      = tree.mb;
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    msbot1  = tree.md(1, 3);
  double    msbot2  = tree.md(2, 3);
  double    mstau1  = tree.me(1, 3);
  double    mstau2  = tree.me(2, 3);
  double    mstop1  = tree.mu(1, 3);
  double    mstop2  = tree.mu(2, 3);
  double    smu     = -displaySusyMu(); // minus sign taken into acct here!
  double    st      = sin(thetat) ;
  double    mz      = displayMzRun();
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau);
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    g       = displayGaugeCoupling(2);
  double    mHc     = tree.mhiggs(4);
  double    mA      = tree.mhiggs(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    alpha   = tree.thetaH;
  double    calpha2 = sqr(cos(alpha)), salpha2 = sqr(sin(alpha)),
    s2alpha = sin(2.0 * alpha), c2alpha = cos(2.0 * alpha);
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta),
    sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);
  double ht = displayYukawaElement(YU, 3, 3),
    hb = displayYukawaElement(YD, 3, 3),
    htau = displayYukawaElement(YE, 3, 3);

  // fermions: 3rd family only for now
  double fermions = 3.0 * sqr(displayYukawaElement(YD, 3, 3)) *
    ((sqr(p) - 4.0 * sqr(mb)) *
     b0(p, mb, mb, q) - 2.0 * a0(mb, q));

  fermions = fermions + sqr(displayYukawaElement(YE, 3, 3)) *
    ((sqr(p) - 4.0 * sqr(mtau)) * b0(p, mtau, mtau, q) -
     2.0 * a0(mtau, q));

  double sbots = 3.0 * sqr(displayYukawaElement(YD, 3, 3)) *
    (a0(msbot1, q) + a0(msbot2, q));

  double staus =  sqr(displayYukawaElement(YE, 3, 3)) *
    (a0(mstau1, q) + a0(mstau2, q));

  double stops = 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (guL * (sqr(ct) * a0(mstop1, q) + sqr(st) * a0(mstop2, q)) +
     guR * (sqr(st) * a0(mstop1, q) + sqr(ct) * a0(mstop2, q)));

  sbots = sbots + 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (gdL * (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q)) +
     gdR * (sqr(sb) * a0(msbot1, q) + sqr(cb) * a0(msbot2, q)));

  staus = staus + sqr(g) / (2.0 * cw2DRbar) *
    (geL * (sqr(ctau) * a0(mstau1, q) + sqr(stau) * a0(mstau2, q)) +
     geR * (sqr(stau) * a0(mstau1, q) + sqr(ctau) * a0(mstau2, q)));

  double sups = 0.0, sdowns = 0.0, sneutrinos = 0.0, sleps = 0.0;
  int fam; for(fam=1; fam<=2; fam++) {
    sups = sups + 3.0 * sqr(g) / (2.0 * cw2DRbar) *
      (guL * a0(tree.mu(1, fam), q) +
       guR * a0(tree.mu(2, fam), q));
    sdowns = sdowns + 3.0 * sqr(g) / (2.0 * cw2DRbar) *
      (gdL * a0(tree.md(1, fam), q) +
       gdR * a0(tree.md(2, fam), q));
    sleps = sleps + sqr(g) / (2.0 * cw2DRbar) *
      (geL * a0(tree.me(1, fam), q) +
       geR * a0(tree.me(2, fam), q));
  }

  // stop couplings to s1 Higgs state
  DoubleMatrix ls1tt(2, 2);
  ls1tt(1, 1) = g * mz * guL * cosb / costhDrbar;
  ls1tt(1, 2) = ht * smu / root2; ls1tt(2, 1) = ls1tt(1, 2);
  ls1tt(2, 2) = g * mz * guR * cosb / costhDrbar;

  // sbottom couplings to s1 Higgs state
  DoubleMatrix ls1bb(2, 2);
  ls1bb(1, 1) = g * mz * gdL * cosb / costhDrbar + root2 * hb * mb;
  ls1bb(1, 2) = displayTrilinear(DA, 3, 3) / root2;
  ls1bb(2, 1) = ls1bb(1, 2);
  ls1bb(2, 2) = g * mz * gdR * cosb / costhDrbar + root2 * hb * mb;

  // stau couplings to s1 Higgs state
  DoubleMatrix ls1tautau(2, 2);
  ls1tautau(1, 1) = g * mz * geL * cosb / costhDrbar + root2 * htau * mtau;
  ls1tautau(1, 2) = displayTrilinear(EA, 3, 3) / root2;
  ls1tautau(2, 1) = ls1tautau(1, 2);
  ls1tautau(2, 2) = g * mz * geR * cosb / costhDrbar + root2 * htau * mtau;


  // Mix 3rd family up
  ls1tt = rot2d(thetat) * ls1tt * rot2d(-thetat);
  ls1bb = rot2d(thetab) * ls1bb * rot2d(-thetab);
  ls1tautau = rot2d(thetatau) * ls1tautau * rot2d(-thetatau);

  int i, j; for (i=1; i<=2; i++) {
    for (j=1; j<=2; j++) {
      stops = stops + 3.0 * sqr(ls1tt(i, j)) *
        b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      sbots = sbots + 3.0 * sqr(ls1bb(i, j)) *
        b0(p, tree.md(i, 3), tree.md(j, 3), q);
      staus = staus +  sqr(ls1tautau(i, j)) *
        b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }}

  // selectron couplings to s1 Higgs state: neglect Yukawas + mixing
  double ls1eeLL = g * mz * geL * cosb / costhDrbar;
  double ls1eeRR = g * mz * geR * cosb / costhDrbar;
  double ls1ddLL = g * mz * gdL * cosb / costhDrbar;
  double ls1ddRR = g * mz * gdR * cosb / costhDrbar;
  double ls1uuLL = g * mz * guL * cosb / costhDrbar;
  double ls1uuRR = g * mz * guR * cosb / costhDrbar;

  int k;
  for (k=1; k<=2; k++) {
    sups = sups +
      3.0 * sqr(ls1uuLL) *
      b0(p, tree.mu(1, k), tree.mu(1, k), q) +
      + 3.0 * sqr(ls1uuRR) *
      b0(p, tree.mu(2, k), tree.mu(2, k), q);
    sdowns = sdowns +
      3.0 * sqr(ls1ddLL) *
      b0(p, tree.md(1, k), tree.md(1, k), q) +
      + 3.0 * sqr(ls1ddRR) *
      b0(p, tree.md(2, k), tree.md(2, k), q);
    sleps = sleps +
      sqr(ls1eeLL) *
      b0(p, tree.me(1, k), tree.me(1, k), q) +
      + sqr(ls1eeRR) *
      b0(p, tree.me(2, k), tree.me(2, k), q);
    }

  sneutrinos = sneutrinos +
    sqr(g) / (2.0 * cw2DRbar) * gsnuL *
    (a0(tree.msnu(1), q) + a0(tree.msnu(2), q) +
     a0(tree.msnu(3), q)) +
    sqr(g * mz / sqrt(cw2DRbar) * gsnuL * cosb) *
    (b0(p, tree.msnu(1), tree.msnu(1), q) +
     b0(p, tree.msnu(2), tree.msnu(2), q) +
     b0(p, tree.msnu(3), tree.msnu(3), q));

  double higgs = sqr(g) * 0.25 *
    (sinb2 * (2.0 * ffn(p, mHc, displayMwRun(), q) +
              ffn(p, mA, mz, q) / cw2DRbar) +
     cosb2 * (2.0 * ffn(p, displayMwRun(), displayMwRun(), q)  + ffn(p, mz, mz, q) / cw2DRbar)) +
    1.75 * sqr(g) * cosb2 *
    (2.0 * sqr(displayMwRun()) * b0(p, displayMwRun(), displayMwRun(), q) +
     sqr(mz) * b0(p, mz, mz, q) / cw2DRbar) +
    sqr(g) * (2.0 * a0(displayMwRun(), q) + a0(mz, q) / cw2DRbar);

  // Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  // are symmetric (ie hHs1 = Hhs1)
  DoubleMatrix hhs1(4, 4);
  hhs1(1, 1) = cosb * (3.0 * calpha2 - salpha2) - sinb * s2alpha;
  hhs1(2, 2) = cosb * (3.0 * salpha2 - calpha2) + sinb * s2alpha;
  hhs1(1, 2) = -2.0 * cosb * s2alpha - sinb * c2alpha;
  hhs1(2, 1) = hhs1(1, 2);
  hhs1(3, 3) = cos2b * cosb;
  hhs1(4, 4) = -cos2b * cosb;
  hhs1(3, 4) = -sin2b * cosb; hhs1(4, 3) = hhs1(3, 4);
  hhs1 = hhs1 * (g * mz / (2.0 * costhDrbar));

  // Quadrilinear Higgs couplings
  DoubleVector hhs1s1(4);
  hhs1s1(1) = 3.0 * calpha2 - salpha2;
  hhs1s1(2) = 3.0 * salpha2 - calpha2;
  hhs1s1(3) = cos2b; hhs1s1(4) = -cos2b;
  hhs1s1 = hhs1s1 * (sqr(g) * 0.25 / (sqr(costhDrbar)));

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  for (i=1; i<=4; i++) {
    for (j=1; j<=4; j++) {
      higgs = higgs + 0.5 * sqr(hhs1(i, j)) * b0(p, higgsm(i), higgsm(j), q);
      //      cout << "higgs(" << i << "," << j << ")=" << higgs;
    }
    higgs = higgs + 0.5 * hhs1s1(i) * a0(higgsm(i), q);
    //      cout << "higgs(" << i << "," << j << ")=" << higgs;
  }

  //  cout << hhs1 << hhs1s1 << higgsm;

  // Basis (G+ H+, G- H-)
  DoubleMatrix hphps1(2, 2);
  hphps1(1, 1) = cos2b * cosb;
  hphps1(2, 2) = -cos2b * cosb + 2.0 * cw2DRbar * cosb;
  hphps1(1, 2) = -sin2b * cosb + cw2DRbar * sinb;
  hphps1(2, 1) = hphps1(1, 2);
  hphps1 = hphps1 * (g * mz * 0.5 / costhDrbar);

  // (G+ H+)
  DoubleVector hphps1s1(2);
  hphps1s1(1) = cw2DRbar + sw2DRbar * cos2b;
  hphps1s1(2) = cw2DRbar - sw2DRbar * cos2b;
  hphps1s1 = hphps1s1 * (sqr(g) * 0.25 / cw2DRbar);

  for (i=1; i<=2; i++) {
    for (j=1; j<=2; j++)
      higgs = higgs + sqr(hphps1(i, j)) * b0(p, higgsc(i), higgsc(j), q);
    higgs = higgs + hphps1s1(i) * a0(higgsc(i), q);
  }

  // Neutralino contribution
  double neutralinos = 0.0;

  DoubleMatrix aPsi(4, 4);
  ComplexMatrix aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  aPsi(1, 3) = -gp * 0.5;
  aPsi(2, 3) = g * 0.5;
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();

  DoubleMatrix fChiChis1s1(4, 4), gChiChis1s1(4, 4);
  for(i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      fChiChis1s1(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChis1s1(i, j) = (bChi(i, j).conj() * aChi(i, j) +
        aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 *
        (fChiChis1s1(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
         gChiChis1s1(i, j) * mneut(i) * mneut(j) *
         b0(p, mneut(i), mneut(j), q));
    }

  // Chargino contribution
  double chargino = 0.0;
  DoubleMatrix aPsic(2, 2);
  aPsic(1, 2) = g / root2;
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fChiChis1s1(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
      gChiChis1s1(i, j) = (bChic(i, j).conj() * aChic(i, j) +
        aChic(i, j).conj() * bChic(i, j)).real();
      chargino = chargino +
        (fChiChis1s1(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
         gChiChis1s1(i, j) * mch(i) * mch(j) *
         b0(p, mch(i), mch(j), q));
    }

  return
    (sups  + sdowns  + sleps  + stops  + sbots  + staus  + sneutrinos +
     fermions + higgs + neutralinos + chargino) / (16.0 * sqr(PI));
}

double MssmSoftsusy::pis1s2(double p, double q) const {
  drBarPars tree(displayDrBarPars());

  double    beta    = atan(displayTanb());
  double    mt   =  displayYukawaElement(YU, 3, 3) * displayHvev()
    * sin(beta) / root2;
  double    mb   =  displayYukawaElement(YD, 3, 3) * displayHvev()
    * cos(beta) / root2;
  double    mtau =  displayYukawaElement(YE, 3, 3) * displayHvev()
    * cos(beta) / root2;
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    smu     = -displaySusyMu(); // minus sign taken into acct here!
  double    g       = displayGaugeCoupling(2);
  double    mHc     = tree.mhiggs(4);
  double    mA      = tree.mhiggs(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    alpha   = tree.thetaH;
  double    calpha2 = sqr(cos(alpha)), salpha2 = sqr(sin(alpha)),
    s2alpha = sin(2.0 * alpha), c2alpha = cos(2.0 * alpha);
  double cosb = cos(beta), cos2b = cos(2.0 * beta),
    sinb = sin(beta), sin2b = sin(2.0 * beta);
  double ht = displayYukawaElement(YU, 3, 3), hb = displayYukawaElement(YD, 3, 3), htau = displayYukawaElement(YE, 3, 3), mz = displayMzRun();



  // stop couplings to s1 Higgs state
  DoubleMatrix ls1tt(2, 2);
  ls1tt(1, 1) = g * mz * guL * cosb / costhDrbar;
  ls1tt(1, 2) = ht * smu / root2; ls1tt(2, 1) = ls1tt(1, 2);
  ls1tt(2, 2) = g * mz * guR * cosb / costhDrbar;

  // sbottom couplings to s1 Higgs state
  DoubleMatrix ls1bb(2, 2);
  ls1bb(1, 1) = g * mz * gdL * cosb / costhDrbar + root2 * hb * mb;
  ls1bb(1, 2) = displayTrilinear(DA, 3, 3) / root2;
  ls1bb(2, 1) = ls1bb(1, 2);
  ls1bb(2, 2) = g * mz * gdR * cosb / costhDrbar + root2 * hb * mb;

  // stau couplings to s1 Higgs state
  DoubleMatrix ls1tautau(2, 2);
  ls1tautau(1, 1) = g * mz * geL * cosb / costhDrbar + root2 * htau * mtau;
  ls1tautau(1, 2) = displayTrilinear(EA, 3, 3) / root2;
  ls1tautau(2, 1) = ls1tautau(1, 2);
  ls1tautau(2, 2) = g * mz * geR * cosb / costhDrbar + root2 * htau * mtau;

  // Mix 3rd family up
  ls1tt = rot2d(thetat) * ls1tt * rot2d(-thetat);
  ls1bb = rot2d(thetab) * ls1bb * rot2d(-thetab);
  ls1tautau = rot2d(thetatau) * ls1tautau * rot2d(-thetatau);

  // stop couplings to s2 Higgs state
  DoubleMatrix ls2tt(2, 2);
  ls2tt(1, 1) = - g * mz * guL * sinb / costhDrbar + root2 * ht * mt;
  ls2tt(1, 2) = displayTrilinear(UA, 3, 3) / root2;
  ls2tt(2, 1) = ls2tt(1, 2);
  ls2tt(2, 2) = - g * mz * guR * sinb / costhDrbar + root2 * ht * mt;

  // sbottom couplings to s2 Higgs state
  DoubleMatrix ls2bb(2, 2);
  ls2bb(1, 1) = -g * mz * gdL * sinb / costhDrbar;
  ls2bb(1, 2) = displayYukawaElement(YD, 3, 3) / root2 * smu;
  ls2bb(2, 1) = ls2bb(1, 2);
  ls2bb(2, 2) = - g * mz * gdR * sinb / costhDrbar;

  // stau couplings to s2 Higgs state
  DoubleMatrix ls2tautau(2, 2);
  ls2tautau(1, 1) = -g * mz * geL * sinb / costhDrbar;
  ls2tautau(1, 2) = displayYukawaElement(YE, 3, 3) / root2 * smu;
  ls2tautau(2, 1) = ls2tautau(1, 2);
  ls2tautau(2, 2) = -g * mz * geR * sinb / costhDrbar;

  // Mix 3rd family up
  ls2tt = rot2d(thetat) * ls2tt * rot2d(-thetat);
  ls2bb = rot2d(thetab) * ls2bb * rot2d(-thetab);
  ls2tautau = rot2d(thetatau) * ls2tautau * rot2d(-thetatau);

  double sfermions = 0.0;
  int i, j; for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      sfermions = sfermions + 3.0 * ls1tt(i, j) * ls2tt(i, j) *
        b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      sfermions = sfermions + 3.0 * ls1bb(i, j) * ls2bb(i, j) *
        b0(p, tree.md(i, 3), tree.md(j, 3), q);
      sfermions = sfermions + ls1tautau(i, j) * ls2tautau(i, j) *
        b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }

  // sneutrinos
  double l1 = g * mz / costhDrbar * cosb * gsnuL,
    l2 = -g * mz / costhDrbar * gsnuL * sinb;
  sfermions = sfermions +
    l1 * l2 * (b0(p, tree.msnu(1), tree.msnu(1), q) +
               b0(p, tree.msnu(2), tree.msnu(2), q) +
               b0(p, tree.msnu(3), tree.msnu(3), q));

  // selectron couplings to s1 Higgs state: neglect Yukawas + mixing
  double ls1uuLL = g * mz * guL * cosb / costhDrbar;
  double ls1uuRR = g * mz * guR * cosb / costhDrbar;
  double ls1eeLL = g * mz * geL * cosb / costhDrbar;
  double ls1eeRR = g * mz * geR * cosb / costhDrbar;
  double ls1ddLL = g * mz * gdL * cosb / costhDrbar;
  double ls1ddRR = g * mz * gdR * cosb / costhDrbar;
  // couplings to s2 Higgs state: neglect Yukawas + mixing
  double ls2uuLL = -g * mz * guL * sinb / costhDrbar;
  double ls2uuRR = -g * mz * guR * sinb / costhDrbar;
  double ls2eeLL = -g * mz * geL * sinb / costhDrbar;
  double ls2eeRR = -g * mz * geR * sinb / costhDrbar;
  double ls2ddLL = -g * mz * gdL * sinb / costhDrbar;
  double ls2ddRR = -g * mz * gdR * sinb / costhDrbar;

  int k;
  for (k=1; k<=2; k++) {
    sfermions = sfermions + 3.0 * ls1uuLL * ls2uuLL *
      b0(p, tree.mu(1, k), tree.mu(1, k), q) +
      + 3.0 * ls1uuRR * ls2uuRR *
      b0(p, tree.mu(2, k), tree.mu(2, k), q);
    sfermions = sfermions + 3.0 * ls1ddLL * ls2ddLL *
      b0(p, tree.md(1, k), tree.md(1, k), q) +
      + 3.0 * ls1ddRR * ls2ddRR *
      b0(p, tree.md(2, k), tree.md(2, k), q);
    sfermions = sfermions + ls1eeLL * ls2eeLL *
      b0(p, tree.me(1, k), tree.me(1, k), q) +
      + ls1eeRR * ls2eeRR *
      b0(p, tree.me(2, k), tree.me(2, k), q);
    }

  double higgs = sqr(g) * 0.25 * sinb * cosb *
    (2.0 * ffn(p, displayMwRun(), displayMwRun(), q) -2.0 * ffn(p, mHc, displayMwRun(), q) +
     (ffn(p, mz, mz, q) - ffn(p, mA, mz, q)) / cw2DRbar +
     7.0 * (2.0 * sqr(displayMwRun()) * b0(p, displayMwRun(), displayMwRun(), q) +
            sqr(mz) * b0(p, mz, mz, q) / cw2DRbar));

  // Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  // are symmetric (ie hHs1 = Hhs1)
  DoubleMatrix hhs1(4, 4);
  hhs1(1, 1) = cosb * (3.0 * calpha2 - salpha2) - sinb * s2alpha;
  hhs1(2, 2) = cosb * (3.0 * salpha2 - calpha2) + sinb * s2alpha;
  hhs1(1, 2) = -2.0 * cosb * s2alpha - sinb * c2alpha;
  hhs1(2, 1) = hhs1(1, 2);
  hhs1(3, 3) = cos2b * cosb;
  hhs1(4, 4) = -cos2b * cosb;
  hhs1(3, 4) = -sin2b * cosb; hhs1(4, 3) = hhs1(3, 4);
  hhs1 = hhs1 * (g * mz / (2.0 * costhDrbar));
  // Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  // are symmetric (ie hHs2 = Hhs2)
  DoubleMatrix hhs2(4, 4);
  hhs2(1, 1) = sinb * (3.0 * salpha2 - calpha2) - cosb * s2alpha;
  hhs2(2, 2) = sinb * (3.0 * calpha2 - salpha2) + cosb * s2alpha;
  hhs2(1, 2) = 2.0 * sinb * s2alpha - cosb * c2alpha;
  hhs2(2, 1) = hhs2(1, 2);
  hhs2(3, 3) = -cos2b * sinb;
  hhs2(4, 4) = cos2b * sinb;
  hhs2(3, 4) = sin2b * sinb; hhs2(4, 3) = hhs2(3, 4);
  hhs2 = hhs2 * (g * mz / (2.0 * costhDrbar));

  // Quadrilinear Higgs couplings
  DoubleVector hhs1s2(4);
  hhs1s2(1) = -s2alpha;
  hhs1s2(2) = s2alpha;
  hhs1s2 = hhs1s2 * (sqr(g) * 0.25 / (sqr(costhDrbar)));

  //  cout << "alpha=" << alpha << " g=" << g << " " << gp << " cw" << costhDrbar << " " << displayTanb();

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  assignHiggs(higgsm, higgsc);

  double trilinear = 0.;
  double quartic = 0., quarticNeut = 0.;

  for (i=1; i<=4; i++) {
    for (j=1; j<=4; j++) {
      higgs = higgs + 0.5 * hhs1(i, j) * hhs2(i, j) *
        b0(p, higgsm(i), higgsm(j), q);
        trilinear = trilinear + 0.5 * hhs1(i, j) * hhs2(i, j) *
          b0(p, higgsm(i), higgsm(j), q);
    }
    higgs = higgs + 0.5 * hhs1s2(i) * a0(higgsm(i), q);
    quarticNeut = quarticNeut + 0.5 * hhs1s2(i) * a0(higgsm(i), q);
  }

  DoubleMatrix hphps2(2, 2);
  hphps2(1, 1) = -cos2b * sinb;
  hphps2(2, 2) = cos2b * sinb + 2.0 * cw2DRbar * sinb;
  hphps2(1, 2) = sin2b * sinb - cw2DRbar * cosb;
  hphps2(2, 1) = hphps2(1, 2);
  hphps2 = hphps2 * (g * mz * 0.5 / costhDrbar);
  DoubleMatrix hphps1(2, 2);
  hphps1(1, 1) = cos2b * cosb;
  hphps1(2, 2) = -cos2b * cosb + 2.0 * cw2DRbar * cosb;
  hphps1(1, 2) = -sin2b * cosb + cw2DRbar * sinb;
  hphps1(2, 1) = hphps1(1, 2);
  hphps1 = hphps1 * (g * mz * 0.5 / costhDrbar);

  DoubleVector hphps1s2(2);
  hphps1s2(1) = - cw2DRbar * sin2b;
  hphps1s2(2) = cw2DRbar * sin2b;
  hphps1s2 = hphps1s2 * (sqr(g) * 0.25 / cw2DRbar);

  for (i=1; i<=2; i++) {
    for (j=1; j<=2; j++) {
      higgs = higgs + hphps1(i, j) * hphps2(i, j)
        * b0(p, higgsc(i), higgsc(j), q);
        trilinear = trilinear +
      hphps1(i, j) * hphps2(i, j)
          * b0(p, higgsc(i), higgsc(j), q);
    }
    higgs = higgs + hphps1s2(i) * a0(higgsc(i), q);
    quartic = quartic + hphps1s2(i) * a0(higgsc(i), q);
  }

  // Neutralino contribution
  double neutralinos = 0.0;

  DoubleMatrix aPsi1(4, 4);
  ComplexMatrix aChi1(4, 4), bChi1(4, 4);
  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  aPsi1(1, 3) = -gp * 0.5;
  aPsi1(2, 3) = g * 0.5;
  aPsi1.symmetrise();
  aChi1 = n.complexConjugate() * aPsi1 * n.hermitianConjugate();
  bChi1 = n * aPsi1 * n.transpose();
  DoubleMatrix aPsi2(4, 4);
  ComplexMatrix aChi2(4, 4), bChi2(4, 4);
  aPsi2(1, 4) = gp * 0.5;
  aPsi2(2, 4) = -g * 0.5;
  aPsi2.symmetrise();
  aChi2 = n.complexConjugate() * aPsi2 * n.hermitianConjugate();
  bChi2 = n * aPsi2 * n.transpose();

  DoubleMatrix fChiChis1s2(4, 4), gChiChis1s2(4, 4);
  for(i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      fChiChis1s2(i, j) = (aChi1(i, j).conj() * aChi2(i, j) +
        bChi1(i, j).conj() * bChi2(i, j)).real();
      gChiChis1s2(i, j) = (bChi1(i, j).conj() * aChi2(i, j) +
        aChi1(i, j).conj() * bChi2(i, j)).real();
      neutralinos = neutralinos + 0.5 *
        (fChiChis1s2(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
         gChiChis1s2(i, j) * mneut(i) * mneut(j) *
         b0(p, mneut(i), mneut(j), q));
    }

  // Chargino contribution
  double chargino = 0.0;
  DoubleMatrix aPsic1(2, 2), aPsic2(2, 2);
  aPsic1(1, 2) = g / root2;
  ComplexMatrix aChic1(2, 2), bChic1(2, 2);
  ComplexMatrix aChic2(2, 2), bChic2(2, 2);
  aChic1 = v.complexConjugate() * aPsic1 * u.hermitianConjugate();
  bChic1 = u * aPsic1.transpose() * v.transpose();
  aPsic2(2, 1) = g / root2;
  aChic2 = v.complexConjugate() * aPsic2 * u.hermitianConjugate();
  bChic2 = u * aPsic2.transpose() * v.transpose();

  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fChiChis1s2(i, j) = (aChic1(i, j).conj() * aChic2(i, j) +
        bChic1(i, j).conj() * bChic2(i, j)).real();
      gChiChis1s2(i, j) = (bChic1(i, j).conj() * aChic2(i ,j) +
        aChic1(i, j).conj() * bChic2(i, j)).real();
      chargino = chargino +
        (fChiChis1s2(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
         gChiChis1s2(i, j) * mch(i) * mch(j) *
         b0(p, mch(i), mch(j), q));
    }

  return (sfermions + higgs + neutralinos + chargino)
    / (16.0 * sqr(PI));
}

// checked 28.10.02
double MssmSoftsusy::pis2s2(double p, double q) const {
  drBarPars tree(displayDrBarPars());

  double    beta    = atan(displayTanb());
  double    mt   =  displayYukawaElement(YU, 3, 3) * displayHvev()
    * sin(beta) / root2;
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    msbot1  = tree.md(1, 3);
  double    msbot2  = tree.md(2, 3);
  double    mstau1  = tree.me(1, 3);
  double    mstau2  = tree.me(2, 3);
  double    mstop1  = tree.mu(1, 3);
  double    mstop2  = tree.mu(2, 3);
  double    smu     = -displaySusyMu(); // minus sign taken into acct here!
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau);
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    g       = displayGaugeCoupling(2);
  double    mHc     = tree.mhiggs(4);
  double    mA      = tree.mhiggs(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    alpha   = tree.thetaH;
  double    calpha2 = sqr(cos(alpha)), salpha2 = sqr(sin(alpha)),
    s2alpha = sin(2.0 * alpha), c2alpha = cos(2.0 * alpha);
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta),
    sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);
  double ht = displayYukawaElement(YU, 3, 3), mz = displayMzRun();

  // fermions: 3rd family only for now, others make negligible difference
  double fermions = 3.0 * sqr(ht) *
    ((sqr(p) - 4.0 * sqr(mt)) *
     b0(p, mt, mt, q) - 2.0 * a0(mt, q));

  // stop contribution
  double sfermions = 3.0 * sqr(ht) * (a0(mstop1, q) + a0(mstop2, q));
  sfermions = sfermions - 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (guL * (sqr(ct) * a0(mstop1, q) + sqr(st) * a0(mstop2, q)) +
     guR * (sqr(st) * a0(mstop1, q) + sqr(ct) * a0(mstop2, q)));

  // sbottom contribution
  sfermions = sfermions - 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (gdL * (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q)) +
     gdR * (sqr(sb) * a0(msbot1, q) + sqr(cb) * a0(msbot2, q)));

  //stau
  sfermions = sfermions - sqr(g) / (2.0 * cw2DRbar) *
    (geL * (sqr(ctau) * a0(mstau1, q) + sqr(stau) * a0(mstau2, q)) +
     geR * (sqr(stau) * a0(mstau1, q) + sqr(ctau) * a0(mstau2, q)));

  // first two families of sparticles
  int fam; for(fam=1; fam<=2; fam++)
    sfermions = sfermions - 3.0 * sqr(g) / (2.0 * cw2DRbar) *
      (guL * a0(tree.mu(1, fam), q) +
       guR * a0(tree.mu(2, fam), q) +
       gdL * a0(tree.md(1, fam), q) +
       gdR * a0(tree.md(2, fam), q)) - sqr(g) / (2.0 * cw2DRbar) *
      (geL * a0(tree.me(1, fam), q) +
       geR * a0(tree.me(2, fam), q)) -
      sqr(g)  / (2.0 * cw2DRbar) * gsnuL *
      (a0(tree.msnu(fam), q));
  sfermions = sfermions -
     sqr(g)  / (2.0 * cw2DRbar) * gsnuL *
      (a0(tree.msnu(3), q));

  // stop couplings to s2 Higgs state
  DoubleMatrix ls2tt(2, 2);
  ls2tt(1, 1) = - g * mz * guL * sinb / costhDrbar + root2 * ht * mt;
  ls2tt(1, 2) = displayTrilinear(UA, 3, 3) / root2;
  ls2tt(2, 1) = ls2tt(1, 2);
  ls2tt(2, 2) = - g * mz * guR * sinb / costhDrbar + root2 * ht * mt;

  // sbottom couplings to s2 Higgs state
  DoubleMatrix ls2bb(2, 2);
  ls2bb(1, 1) = -g * mz * gdL * sinb / costhDrbar;
  ls2bb(1, 2) = displayYukawaElement(YD, 3, 3) / root2 * smu;
  ls2bb(2, 1) = ls2bb(1, 2);
  ls2bb(2, 2) = - g * mz * gdR * sinb / costhDrbar;

  // stau couplings to s2 Higgs state
  DoubleMatrix ls2tautau(2, 2);
  ls2tautau(1, 1) = -g * mz * geL * sinb / costhDrbar;
  ls2tautau(1, 2) = displayYukawaElement(YE, 3, 3) / root2 * smu;
  ls2tautau(2, 1) = ls2tautau(1, 2);
  ls2tautau(2, 2) = -g * mz * geR * sinb / costhDrbar;

  // Mix 3rd family up
  ls2tt = rot2d(thetat) * ls2tt * rot2d(-thetat);
  ls2bb = rot2d(thetab) * ls2bb * rot2d(-thetab);
  ls2tautau = rot2d(thetatau) * ls2tautau * rot2d(-thetatau);

  int i, j; for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      // stop
      sfermions = sfermions + 3.0 * sqr(ls2tt(i, j)) *
        b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      // sbottom
      sfermions = sfermions + 3.0 * sqr(ls2bb(i, j)) *
        b0(p, tree.md(i, 3), tree.md(j, 3), q);
      // stay
      sfermions = sfermions + sqr(ls2tautau(i, j)) *
        b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }

  // couplings to s2 Higgs state: neglect Yukawas + mixing
  double ls2uuLL = -g * mz * guL * sinb / costhDrbar;
  double ls2uuRR = -g * mz * guR * sinb / costhDrbar;
  double ls2eeLL = -g * mz * geL * sinb / costhDrbar;
  double ls2eeRR = -g * mz * geR * sinb / costhDrbar;
  double ls2ddLL = -g * mz * gdL * sinb / costhDrbar;
  double ls2ddRR = -g * mz * gdR * sinb / costhDrbar;

  int k;
  for (k=1; k<=2; k++) {
    sfermions = sfermions + 3.0 * sqr(ls2uuLL) *
      b0(p, tree.mu(1, k), tree.mu(1, k), q) +
      + 3.0 * sqr(ls2uuRR) *
      b0(p, tree.mu(2, k), tree.mu(2, k), q);
    sfermions = sfermions + 3.0 * sqr(ls2ddLL) *
      b0(p, tree.md(1, k), tree.md(1, k), q) +
      + 3.0 * sqr(ls2ddRR) *
      b0(p, tree.md(2, k), tree.md(2, k), q);
    sfermions = sfermions + sqr(ls2eeLL) *
      b0(p, tree.me(1, k), tree.me(1, k), q) +
      + sqr(ls2eeRR) *
      b0(p, tree.me(2, k), tree.me(2, k), q);
    }

  sfermions = sfermions +
    sqr(g * mz * gsnuL * sinb) / cw2DRbar *
    (b0(p, tree.msnu(1), tree.msnu(1), q) +
     b0(p, tree.msnu(2), tree.msnu(2), q) +
     b0(p, tree.msnu(3), tree.msnu(3), q));

  // gauge boson/Higgs contributions
  double higgs = sqr(g) * 0.25 *
    (cosb2 * (2.0 * ffn(p, mHc, displayMwRun(), q) +
              ffn(p, mA, mz, q) / cw2DRbar) +
     sinb2 * (2.0 * ffn(p, displayMwRun(), displayMwRun(), q) +
              ffn(p, mz, mz, q) / cw2DRbar)) +
    1.75 * sqr(g) * sinb2 *
    (2.0 * sqr(displayMwRun()) * b0(p, displayMwRun(), displayMwRun(), q) +
     sqr(mz) * b0(p, mz, mz, q) / cw2DRbar) +
    sqr(g) * (2.0 * a0(displayMwRun(), q) + a0(mz, q) / cw2DRbar);

  double quartic = 0., trilinear = 0.;

  // Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  // are symmetric (ie hHs2 = Hhs2)
  DoubleMatrix hhs2(4, 4);
  hhs2(1, 1) = sinb * (3.0 * salpha2 - calpha2) - cosb * s2alpha;
  hhs2(2, 2) = sinb * (3.0 * calpha2 - salpha2) + cosb * s2alpha;
  hhs2(1, 2) = 2.0 * sinb * s2alpha - cosb * c2alpha;
  hhs2(2, 1) = hhs2(1, 2);
  hhs2(3, 3) = -cos2b * sinb;
  hhs2(4, 4) = cos2b * sinb;
  hhs2(3, 4) = sin2b * sinb; hhs2(4, 3) = hhs2(3, 4);
  hhs2 = hhs2 * (g * mz / (2.0 * costhDrbar));

  // Quadrilinear Higgs couplings
  DoubleVector hhs2s2(4);
  hhs2s2(1) = 3.0 * salpha2 - calpha2;
  hhs2s2(2) = 3.0 * calpha2 - salpha2;
  hhs2s2(3) = -cos2b; hhs2s2(4) = cos2b;
  hhs2s2 = hhs2s2 * (sqr(g) * 0.25 / (sqr(costhDrbar)));

  // define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  assignHiggs(higgsm, higgsc);

  for (i=1; i<=4; i++) {
    for (j=1; j<=4; j++) {
      higgs = higgs + 0.5 * sqr(hhs2(i, j)) * b0(p, higgsm(i), higgsm(j), q);
      trilinear = trilinear +
        0.5 * sqr(hhs2(i, j)) * b0(p, higgsm(i), higgsm(j), q);
    }
    higgs = higgs + 0.5 * hhs2s2(i) * a0(higgsm(i), q);
    quartic = quartic + 0.5 * hhs2s2(i) * a0(higgsm(i), q);
  }

  DoubleMatrix hphps2(2, 2);
  hphps2(1, 1) = -cos2b * sinb;
  hphps2(2, 2) = cos2b * sinb + 2.0 * cw2DRbar * sinb;
  hphps2(1, 2) = sin2b * sinb - cw2DRbar * cosb;
  hphps2(2, 1) = hphps2(1, 2);
  hphps2 = hphps2 * (g * mz * 0.5 / costhDrbar);

  DoubleVector hphps2s2(2);
  hphps2s2(1) = cw2DRbar - sw2DRbar * cos2b;
  hphps2s2(2) = cw2DRbar + sw2DRbar * cos2b;
  hphps2s2 = hphps2s2 * (sqr(g) * 0.25 / cw2DRbar);

  for (i=1; i<=2; i++) {
    for (j=1; j<=2; j++) {
      higgs = higgs + sqr(hphps2(i, j)) * b0(p, higgsc(i), higgsc(j), q);
      trilinear = trilinear + sqr(hphps2(i, j)) *
        b0(p, higgsc(i), higgsc(j), q);
    }
    higgs = higgs + hphps2s2(i) * a0(higgsc(i), q);
    quartic = quartic + hphps2s2(i) * a0(higgsc(i), q);
  }

  // Neutralino contribution
  double neutralinos = 0.0;

  DoubleMatrix aPsi(4, 4);
  ComplexMatrix aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  aPsi(1, 4) = gp * 0.5;
  aPsi(2, 4) = -g * 0.5;
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();

  DoubleMatrix fChiChis2s2(4, 4), gChiChis2s2(4, 4);
  for(i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      fChiChis2s2(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChis2s2(i, j) = (bChi(i, j).conj() * aChi(i, j) +
        aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 *
        (fChiChis2s2(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
         gChiChis2s2(i, j) * mneut(i) * mneut(j) *
         b0(p, mneut(i), mneut(j), q));
    }

  // Chargino contribution
  double chargino = 0.0;
  DoubleMatrix aPsic(2, 2);
  aPsic(2, 1) = g / root2;
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fChiChis2s2(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
      gChiChis2s2(i, j) = (bChic(i, j).conj() * aChic(i, j) +
        aChic(i, j).conj() * bChic(i, j)).real();
      chargino = chargino +
        (fChiChis2s2(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
         gChiChis2s2(i, j) * mch(i) * mch(j) *
         b0(p, mch(i), mch(j), q));
    }

  return (fermions + sfermions + higgs + neutralinos + chargino)
    / (16.0 * sqr(PI));
}

// New routine
double MssmSoftsusy::piHpHm(double p, double q) const {

  drBarPars tree(displayDrBarPars());

  double    ht      = displayYukawaElement(YU, 3, 3);
  double    hb      = displayYukawaElement(YD, 3, 3);
  double    htau    = displayYukawaElement(YE, 3, 3);
  double    beta    = atan(displayTanb());
  double    mt     =  tree.mt;
  double    mb    =  tree.mb;
  double    mtau =    tree.mtau;
  double    alpha = tree.thetaH;
  double    mz = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cwDRbar    = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cwDRbar);
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    stau    = sin(thetatau);
  double    ctau    = cos(thetatau);
  double    g       = displayGaugeCoupling(2);
  double    mh0     = tree.mhiggs(1);
  double    mHc     = tree.mhiggs(4);
  double    mH     = tree.mhiggs(3);
  double    mA      = tree.mhiggs(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    mwRun   = displayMwRun();
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta),
    sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);
  double   mu = - displaySusyMu();

  double fermions = 3.0 * (sqr(ht) * cosb2 + sqr(hb) * sinb2) *
    gfn(p, mt, mb, q) + sinb2 * sqr(htau) * gfn(p, 0.0, mtau, q) -
    6.0 * hb * ht * mt * mb * sin2b * b0(p, mt, mb, q);

  // first two generations: forget lighter Yukawas
  DoubleMatrix lHpud(2, 2), lHpen(2, 2);
  lHpud(1, 1) = g * mwRun * sin2b / root2;

  double sfermions = 3.0 * sqr(lHpud(1, 1)) *
    (b0(p, tree.mu(1, 1), tree.md(1, 1), q) +
     b0(p, tree.mu(1, 2), tree.md(1, 2), 1)) + sqr(lHpud(1, 1)) *
    (b0(p, tree.msnu(1), tree.me(1, 1), q) +
     b0(p, tree.msnu(2), tree.me(1, 2), q));

  // 3rd family
  lHpud(1, 1) = g * mwRun * sin2b / root2
    - ht * mt * cosb - hb * mb * sinb;
  lHpud(1, 2) = hb * mu * cosb - displayTrilinear(DA, 3, 3) * sinb;
  lHpud(2, 1) = ht * mu * sinb - displayTrilinear(UA, 3, 3) * cosb;
  lHpud(2, 2) = -ht * mb * cosb - hb * mt * sinb;
  lHpud = rot2d(tree.thetat) * lHpud * rot2d(-tree.thetab);

  lHpen(1, 1) = g * mwRun * sin2b / root2
     - htau * mtau * sinb;
  lHpen(1, 2) = htau * mu * cosb - displayTrilinear(EA, 3, 3) * sinb;
  lHpen(2, 1) = 0.;
  lHpen(2, 2) = 0.;
  lHpen = lHpen * rot2d(-tree.thetatau);

  double thirdFamSfermions = 0.;
  int i, j; for(i=1; i<=2; i++)
    for (j=1; j<=2; j++)
      thirdFamSfermions = thirdFamSfermions +
        3.0 * sqr(lHpud(i, j)) * b0(p, tree.mu(i, 3), tree.md(j, 3), q);

  for (j=1; j<=2; j++)
    thirdFamSfermions = thirdFamSfermions +
      sqr(lHpen(1, j)) * b0(p, tree.msnu(3), tree.me(j, 3), q);

  // ups
  sfermions = sfermions +
    3.0 * sqr(g) * cos2b * (-0.5 * guL / cw2DRbar + 0.5) *
    (a0(tree.mu(1, 1), q) + a0(tree.mu(1, 2), q)) +
    3.0 * sqr(g) * cos2b * (-0.5 * guR / cw2DRbar) *
    (a0(tree.mu(2, 1), q) + a0(tree.mu(2, 2), q));
  // downs
  sfermions = sfermions +
    3.0 * sqr(g) * cos2b * (-0.5 * gdL / cw2DRbar - 0.5) *
    (a0(tree.md(1, 1), q) + a0(tree.md(1, 2), q)) +
    3.0 * sqr(g) * cos2b * (-0.5 * gdR / cw2DRbar) *
    (a0(tree.md(2, 1), q) + a0(tree.md(2, 2), q));
  // sneutrinos/selectrons
  sfermions = sfermions +
    sqr(g) * cos2b * (-0.5 * gsnuL / cw2DRbar + 0.5) *
    (a0(tree.msnu(1), q) + a0(tree.msnu(2), q));
  // downs
  sfermions = sfermions +
    sqr(g) * cos2b * (-0.5 * geL / cw2DRbar - 0.5) *
    (a0(tree.me(1, 1), q) + a0(tree.me(1, 2), q)) +
    sqr(g) * cos2b * (-0.5 * geR / cw2DRbar) *
    (a0(tree.me(2, 1), q) + a0(tree.me(2, 2), q));
  // 3rd fam ups
  thirdFamSfermions = thirdFamSfermions +
    3.0 * (sqr(hb) * sinb2 - sqr(g) * cos2b * 0.5 * guL / cw2DRbar +
           0.5 * cos2b * sqr(g)) *
    (sqr(ct) * a0(tree.mu(1, 3), q) + sqr(st) * a0(tree.mu(2, 3), q)) +
    3.0 * (sqr(ht) * cosb2 - sqr(g) * cos2b * 0.5 * guR / cw2DRbar) *
    (sqr(st) * a0(tree.mu(1, 3), q) + sqr(ct) * a0(tree.mu(2, 3), q));
  // 3rd fam downs
  thirdFamSfermions = thirdFamSfermions +
    3.0 * (sqr(ht) * cosb2 - sqr(g) * cos2b * 0.5 * gdL / cw2DRbar -
           0.5 * cos2b * sqr(g)) *
    (sqr(cb) * a0(tree.md(1, 3), q) + sqr(sb) * a0(tree.md(2, 3), q)) +
    3.0 * (sqr(hb) * sinb2 - sqr(g) * cos2b * 0.5 * gdR / cw2DRbar) *
    (sqr(sb) * a0(tree.md(1, 3), q) + sqr(cb) * a0(tree.md(2, 3), q));
  // 3rd fam snus
  thirdFamSfermions = thirdFamSfermions +
    (sqr(htau) * sinb2 - sqr(g) * cos2b * 0.5 * gsnuL / cw2DRbar +
     0.5 * cos2b * sqr(g)) * a0(tree.msnu(3), q);
  // 3rd fam es
  thirdFamSfermions = thirdFamSfermions +
    (-sqr(g) * cos2b * 0.5 * geL / cw2DRbar - 0.5 * cos2b * sqr(g)) *
    (sqr(ctau) * a0(tree.me(1, 3), q) + sqr(stau) * a0(tree.me(2, 3), q)) +
    (sqr(htau) * sinb2 - sqr(g) * cos2b * 0.5 * geR / cw2DRbar) *
    (sqr(stau) * a0(tree.me(1, 3), q) + sqr(ctau) * a0(tree.me(2, 3), q));

  double weak = sqr(g) * 0.25 *
    (sqr(sin(alpha - beta)) * ffn(p, mH, mwRun, q) +
     sqr(cos(alpha - beta)) * ffn(p, mh0, mwRun, q) +
     ffn(p, mA, mwRun, q) + sqr(cos(2.0 * thetaWDRbar)) / cw2DRbar *
     ffn(p, mHc, mz, q)) +
    sqr(g) * sw2 * ffn(p, mHc, 0., q) + 2.0 * sqr(g) * a0(mwRun, q) +
    sqr(g) * sqr(cos(2.0 * thetaWDRbar)) / cw2DRbar * a0(mz, q) +
    sqr(g) * sqr(mwRun) * 0.25 * b0(p, mwRun, mA, q);

  double higgs = 0.;
  DoubleMatrix lHpH0Hm(2, 2);
  // BPMZ not so clear: assuming (D.67) that s_i G+ H- coupling is same as
  // s_i G- H+ (by C conservation)
  lHpH0Hm(1, 1) = -sin2b * cosb + cw2DRbar * sinb;
  lHpH0Hm(2, 1) =  sin2b * sinb - cw2DRbar * cosb;
  lHpH0Hm(1, 2) = -cos2b * cosb + 2.0 * cw2DRbar * cosb;
  lHpH0Hm(2, 2) =  cos2b * sinb + 2.0 * cw2DRbar * sinb;
  // in (H h) basis
  lHpH0Hm = rot2d(alpha) * lHpH0Hm * (g * mz * 0.5 / cwDRbar);
  DoubleVector spech(2), specHm(2);
  spech(1) = mH; spech(2) = mh0;
  specHm(1) = mwRun; specHm(2) = mHc;
  for(i=1; i<=2; i++)
    for(j=1; j<=2; j++)
      higgs = higgs + sqr(lHpH0Hm(i, j)) * b0(p, spech(i), specHm(j), q);

  double l1 = (2.0 * sqr(sin2b) - 1.0) * sqr(g) * 0.25 / cw2DRbar;
  double l2 = 2.0 * sqr(cos2b) * sqr(g) * 0.25 / cw2DRbar;
  higgs = higgs + l1 * a0(specHm(1), q) + l2 * a0(specHm(2), q);
  DoubleMatrix lHpHmH0H0(2, 2);
  lHpHmH0H0(1, 1) = cw2DRbar - sw2DRbar * cos2b;
  lHpHmH0H0(1, 2) = cw2DRbar * sin2b;
  lHpHmH0H0(2, 1) = lHpHmH0H0(1, 2);
  lHpHmH0H0(2, 2) = cw2DRbar + sw2DRbar * cos2b;
  lHpHmH0H0 = rot2d(alpha) * lHpHmH0H0 * rot2d(-alpha) *
    (0.25 * sqr(g) / cw2DRbar);
  double lHpHmG0G0 = (cw2DRbar * (1.0 + sqr(sin2b)) - sw2DRbar * sqr(cos2b)) *
    0.25 * sqr(g) / cw2DRbar;
  double lHpHmA0A0 = sqr(cos2b) * 0.25 * sqr(g) / cw2DRbar;
  higgs = higgs + 0.5 * (lHpHmH0H0(1, 1) * a0(spech(1), q) +
    lHpHmH0H0(2, 2) * a0(spech(2), q) + lHpHmG0G0 * a0(mz, q) +
    lHpHmA0A0 * a0(mA, q));

  double gauginos = 0.;
  ComplexMatrix apph1(4, 2), apph2(4, 2), bpph1(4, 2), bpph2(4, 2);
  apph1(1, 2) = gp / root2;
  bpph2(1, 2) = gp / root2;
  apph1(2, 2) = g / root2;
  bpph2(2, 2) = g / root2;
  apph1(3, 1) = -g;
  bpph2(4, 1) = g;

  // Get to physical gaugino estates, and G+H+
  apph1 = -sinb * tree.nBpmz.complexConjugate() *
    apph1 * tree.uBpmz.transpose();
  bpph2 = cosb * tree.nBpmz * bpph2 * tree.vBpmz.transpose();

  for (i=1; i<=4; i++)
    for(j=1; j<=2; j++) {
      double f = sqr(apph1(i, j).mod()) + sqr(bpph2(i, j).mod());
      double ga = 2.0 * (bpph2(i, j).conj() * apph1(i, j)).real();
      gauginos = gauginos + f * gfn(p, tree.mchBpmz(j), tree.mnBpmz(i), q) -
        2.0 * ga * tree.mchBpmz(j) * tree.mnBpmz(i) *
        b0(p, tree.mchBpmz(j), tree.mnBpmz(i), q);
    }

  const static double loopFac = 1.0 / (16 * sqr(PI));

  double pihh = fermions + sfermions + thirdFamSfermions + weak + higgs +
    gauginos;

  return loopFac * pihh;
}


double MssmSoftsusy::piAA(double p, double q) const {// checked 30.07.03
  drBarPars tree(displayDrBarPars());

 if (tree.mu(1, 3) == 0.0 || tree.mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate piAA without having first calculated"
         << " the DRbar masses.\n";
   return 0.0;
  }

  double    beta    = atan(displayTanb());
  double    mt     =  tree.mt;
  double    mb    =  tree.mb;
  double    mtau =    tree.mtau;
  double    mmu =  displayDataSet().displayMass(mMuon);
  double    mE =  displayDataSet().displayMass(mElectron);
  double    mz = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    stau    = sin(thetatau);
  double    ctau    = cos(thetatau);
  double    g       = displayGaugeCoupling(2);
  double    mh0     = tree.mhiggs(1);
  double    mHc     = tree.mhiggs(4);
  double    mH     = tree.mhiggs(3);
  double    mA      = tree.mhiggs(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta),
    sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);

  // Up fermions/sfermions
  double fermions = 0., upsfermions = 0.0;
  fermions = fermions + cosb2 *
    (3.0 * sqr(displayYukawaElement(YU, 3, 3)) *
     (sqr(p) * b0(p, mt, mt, q) - 2.0 * a0(mt, q))); // ignore 1st 2 families

  // LH gens 1-2
  upsfermions = upsfermions +
    3.0 * (-sqr(g) / (cw2DRbar) * cos2b * 0.5 * guL) *
    (a0(tree.mu(1, 1), q) + a0(tree.mu(1, 2), q));

  // RH gen 1-2
  upsfermions = upsfermions +
    3.0 * (-sqr(g) / (cw2DRbar) * cos2b * 0.5 * guR) *
    (a0(tree.mu(2, 1), q) + a0(tree.mu(2, 2), q));

  upsfermions = upsfermions +
    3.0 * (sqr(displayYukawaElement(YU, 1, 1) *
               (-displaySusyMu() * sinb -
                displaySoftA(UA, 1, 1) * cosb))) * 0.5 *
    (b0(p, tree.mu(1, 1), tree.mu(2, 1), q) +
     b0(p, tree.mu(2, 1), tree.mu(1, 1), q)) +
    3.0 * (sqr(displayYukawaElement(YU, 2, 2) *
               (-displaySusyMu() * sinb -
                displaySoftA(UA, 2, 2) * cosb))) * 0.5 *
    (b0(p, tree.mu(1, 2), tree.mu(2, 2), q) +
     b0(p, tree.mu(2, 2), tree.mu(1, 2), q));

  double stops = 0.0;
  stops = stops +
    3.0 * (sqr(displayYukawaElement(YU, 3, 3)) * cosb2 -
           sqr(g) / (cw2DRbar) * cos2b * 0.5 * guL) *
    (sqr(ct) * a0(tree.mu(1, 3), q) +
     sqr(st) * a0(tree.mu(2, 3), q));

  stops = stops +
    3.0 * (sqr(displayYukawaElement(YU, 3, 3)) * cosb2 -
           sqr(g) / (cw2DRbar) * cos2b * 0.5 * guR) *
    (sqr(st) * a0(tree.mu(1, 3), q) +
     sqr(ct) * a0(tree.mu(2, 3), q));

  double auu12 =
    (displayYukawaElement(YU, 3, 3) *
     displaySusyMu() * sinb + displayTrilinear(UA, 3, 3) * cosb) / root2;

  int i, j;
  stops = stops +
    3.0 * sqr(auu12) *
    (b0(p, tree.mu(1, 3), tree.mu(2, 3), q) +
     b0(p, tree.mu(2, 3), tree.mu(1, 3), q));

  // Other Higgs'
  double higgs = 0.;
  higgs = higgs +
    sqr(g) * 0.25 *
    (2.0 * ffn(p, tree.mhiggs(4), displayMwRun(), q) +
     sqr(sin(tree.thetaH - beta)) / cw2DRbar *
     ffn(p, tree.mhiggs(3), mz, q) +
     sqr(cos(tree.thetaH - beta)) / cw2DRbar *
     ffn(p, tree.mhiggs(1), mz, q));

  // trilinear Higgs coupling Feynman rules
  DoubleVector lAas(2), lAah(2);
  lAas(1) = -cos2b * cosb; lAas(2) = cos2b * sinb;
  lAas = lAas * (g * mz / (2.0 * cos(thetaWDRbar)));
  lAah = rot2d(tree.thetaH) * lAas;

  // quartic Higgs coupling Ferynman rules: check
  DoubleMatrix lAa12(2, 2), lAahh(2, 2);
  lAa12(1, 1) = -cos(2.0 * beta);
  lAa12(2, 2) = cos(2.0 * beta);
  lAahh = rot2d(tree.thetaH) * lAa12 *
    rot2d(-tree.thetaH);
  lAahh = lAahh * sqr(g) / (4.0 * cw2DRbar);
  double lAaaa = sqr(g) / (4.0 * cw2DRbar) * 3.0 * sqr(cos2b),
    lAagg = (3.0 * sqr(sin2b) - 1.0) * sqr(g) / (4.0 * cw2DRbar);

  higgs = higgs +
    0.5 * (sqr(lAah(2)) * b0(p, mA, mh0, q) +
           sqr(lAah(1)) * b0(p, mA, mH, q) +
           lAahh(1, 1) * a0(mH, q) + lAahh(2, 2) * a0(mh0, q) +
           lAaaa * a0(mA, q) + lAagg * a0(mz, q)) +
    sqr(g) * sqr(displayMwRun()) * 0.5 * b0(p, displayMwRun(), mHc, q) +
    sqr(g) / (4.0 * cw2DRbar)
    * ((cw2DRbar * (1.0 + sqr(sin2b)) - sw2DRbar * sqr(cos2b)) *
       a0(displayMwRun(), q) + sqr(cos2b) * a0(mHc, q));

  higgs = higgs + sqr(g) *
    (2.0 * a0(displayMwRun(), q) + a0(mz, q) / cw2DRbar);

  DoubleMatrix aPsiP1(4, 4), aPsiP2(4, 4), bPsiP1(4, 4),
    bPsiP2(4, 4);
  aPsiP1(1, 3) = - gp * 0.5;   aPsiP1(3, 1) = - gp * 0.5;
  aPsiP1(2, 3) = g * 0.5;      aPsiP1(3, 2) = g * 0.5;
  aPsiP2(1, 4) = - gp * 0.5;   aPsiP2(4, 1) = - gp * 0.5;
  aPsiP2(2, 4) = g * 0.5;      aPsiP2(4, 2) = g * 0.5;

  bPsiP1 = -1.0 * aPsiP1; bPsiP2 = -1.0 * aPsiP2;

  ComplexMatrix aPsi(4, 4), bPsi(4, 4), aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  // Rotate to A neutralino basis
  ComplexMatrix aChiChiA(4, 4), bChiChiA(4, 4);
  aChiChiA = -sinb * n.complexConjugate() * aPsiP1 * n.hermitianConjugate() +
    cosb * n.complexConjugate() * aPsiP2 * n.hermitianConjugate();
  bChiChiA = -sinb * n * bPsiP1 * n.transpose() + cosb * n * bPsiP2 *
    n.transpose();

  double neutralinos = 0.0;
  DoubleMatrix fChiChiA(4, 4), gChiChiA(4, 4);
  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      fChiChiA(i, j) = sqr(aChiChiA(i, j).mod()) + sqr(bChiChiA(i, j).mod());
      gChiChiA(i, j) = 2.0 * (bChiChiA(i, j).conj() * aChiChiA(i, j)).real();
      neutralinos = neutralinos +
        0.5 * fChiChiA(i, j) *
        gfn(p, mneut(i), mneut(j), q) -
        gChiChiA(i, j) * mneut(i) * mneut(j) *
        b0(p, mneut(i), mneut(j), q);
    }

  ComplexMatrix aChPsiP1(2, 2), aChPsiP2(2, 2), bChPsiP1(2, 2),
    bChPsiP2(2, 2), aChPsiA(2, 2), bChPsiA(2, 2);
  aChPsiP1(1, 2) = g / root2;   aChPsiP2(2, 1) = - g / root2;
  bChPsiP1 = - 1.0 * aChPsiP1.transpose();
  bChPsiP2 = - 1.0 * aChPsiP2.transpose();

  // To A-mass basis
  aChPsiA = -sinb * aChPsiP1 + cosb * aChPsiP2;
  bChPsiA = -sinb * bChPsiP1 + cosb * bChPsiP2;

  ComplexMatrix aChChiA(2, 2), bChChiA(2, 2);
  DoubleMatrix  fChChiA(2, 2), gChChiA(2, 2);
  aChChiA = v.complexConjugate() * aChPsiA * u.hermitianConjugate();
  bChChiA = u * bChPsiA * v.transpose();

  double charginos = 0.0;
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fChChiA(i, j) = sqr(aChChiA(i, j).mod()) + sqr(bChChiA(i, j).mod());
      gChChiA(i, j) = 2.0 * (bChChiA(i, j).conj() * aChChiA(i, j)).real();

      charginos = charginos +
        fChChiA(i, j) * gfn(p, mch(i), mch(j), q) -
        2.0 * gChChiA(i, j) * mch(i) * mch(j) * b0(p, mch(i), mch(j), q);
    }
  //  cout << "f=" << fChChiA << " g=" << gChChiA;

  fermions = fermions + sinb2 *
    (sqr(displayYukawaElement(YE, 3, 3)) *
     (sqr(p) * b0(p, mtau, mtau, q) - 2.0 * a0(mtau, q)) +
     sqr(displayYukawaElement(YE, 2, 2)) *
     (sqr(p) * b0(p, mmu, mmu, q) - 2.0 * a0(mmu, q)) +
     sqr(displayYukawaElement(YE, 1, 1)) *
     (sqr(p) * b0(p, mE, mE, q) - 2.0 * a0(mE, q)));

  DoubleMatrix aee(2, 2);
  aee(1, 2) = - displayYukawaElement(YE, 3, 3) / root2 *
    (-displaySusyMu() * sinb - displaySoftA(EA, 3, 3) * cosb);
  aee(2, 1) = - aee(1, 2);
  aee = rot2d(thetatau) * aee * rot2d(-thetatau);

  // Down fermions
  fermions = fermions + sinb2 *
    (3.0 * sqr(displayYukawaElement(YD, 3, 3)) *
     (sqr(p) * b0(p, mb, mb, q) - 2.0 * a0(mb, q)));

  DoubleMatrix add(2, 2);
  add(1, 2) = - displayYukawaElement(YD, 3, 3) / root2 *
    (-displaySusyMu() * sinb - displaySoftA(DA, 3, 3) * cosb);
  add(2, 1) = - add(1, 2);
  add = rot2d(thetab) * add * rot2d(-thetab);

  // sneutrinos
  double sneutrinos = -
    sqr(g) / cw2DRbar * 0.5 * gsnuL * cos2b *
    (a0(tree.msnu(1), q) + a0(tree.msnu(2), q) +
     a0(tree.msnu(3), q));

  // LH gens 1-2
  double sleptons = 0.;
  sleptons = -(sqr(g) / (cw2DRbar) * cos2b * 0.5 * geL) *
    (a0(tree.me(1, 1), q) + a0(tree.me(1, 2), q));

  // RH gens 1-2
  sleptons = sleptons -
    (sqr(g) / (cw2DRbar) * cos2b * 0.5 * geR) *
    (a0(tree.me(2, 1), q) + a0(tree.me(2, 2), q));

  double staus =
    (sqr(displayYukawaElement(YE, 3, 3)) * sinb2 -
           sqr(g) / (cw2DRbar) * cos2b * 0.5 * geL) *
    (sqr(ctau) * a0(tree.me(1, 3), q) +
     sqr(stau) * a0(tree.me(2, 3), q));

  staus = staus +
    (sqr(displayYukawaElement(YE, 3, 3)) * sinb2 -
           sqr(g) / (cw2DRbar) * cos2b * 0.5 * geR) *
    (sqr(stau) * a0(tree.me(1, 3), q) +
     sqr(ctau) * a0(tree.me(2, 3), q));

  sleptons = sleptons +
    (sqr(displayYukawaElement(YE, 1, 1) *
               (-displaySusyMu() * cosb -
                displaySoftA(EA, 1, 1) * sinb))) * 0.5 *
    (b0(p, tree.me(1, 1), tree.me(2, 1), q) +
     b0(p, tree.me(2, 1), tree.me(1, 1), q)) +
    (sqr(displayYukawaElement(YE, 2, 2) *
               (-displaySusyMu() * cosb -
                displaySoftA(EA, 2, 2) * sinb))) * 0.5 *
    (b0(p, tree.me(1, 2), tree.me(2, 2), q) +
     b0(p, tree.me(2, 2), tree.me(1, 2), q));


  double aee12 = (displayYukawaElement(YE, 3, 3) *
    displaySusyMu() * cosb + displayTrilinear(EA, 3, 3) * sinb) / root2;

  staus = staus +
    sqr(aee12) *
    (b0(p, tree.me(1, 3), tree.me(2, 3), q) +
     b0(p, tree.me(2, 3), tree.me(1, 3), q));

  // Sign of certain contributions bug-fixed 30-07-03
  double dsfermions = 0.0;
  dsfermions = dsfermions -
    3.0 * (sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdL) *
    (a0(tree.md(1, 1), q) + a0(tree.md(1, 2), q));

  double sbots =
    3.0 * (sqr(displayYukawaElement(YD, 3, 3)) * sinb2 -
           sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdL) *
    (sqr(cb) * a0(tree.md(1, 3), q) +
     sqr(sb) * a0(tree.md(2, 3), q));

  dsfermions = dsfermions -
    3.0 * (sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdR) *
    (a0(tree.md(2, 1), q) + a0(tree.md(2, 2), q));

  sbots = sbots +
    3.0 * (sqr(displayYukawaElement(YD, 3, 3)) * sinb2 -
           sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdR) *
    (sqr(sb) * a0(tree.md(1, 3), q) +
     sqr(cb) * a0(tree.md(2, 3), q));

  dsfermions = dsfermions +
    3.0 * (sqr(displayYukawaElement(YD, 1, 1) *
               (-displaySusyMu() * cosb -
                displaySoftA(DA, 1, 1) * sinb))) * 0.5 *
    (b0(p, tree.md(1, 1), tree.md(2, 1), q) +
     b0(p, tree.md(2, 1), tree.md(1, 1), q)) +
    3.0 * (sqr(displayYukawaElement(YD, 2, 2) *
               (-displaySusyMu() * cosb -
                displaySoftA(DA, 2, 2) * sinb))) * 0.5 *
    (b0(p, tree.md(1, 2), tree.md(2, 2), q) +
     b0(p, tree.md(2, 2), tree.md(1, 2), q));

  double add12 = (displayYukawaElement(YD, 3, 3) *
    displaySusyMu() * cosb + displayTrilinear(DA, 3, 3) * sinb) / root2;

  sbots = sbots +
        3.0 * sqr(add12) *
        (b0(p, tree.md(1, 3), tree.md(2, 3), q) +
         b0(p, tree.md(2, 3), tree.md(1, 3), q));

  double sfermions = sleptons + staus + upsfermions + dsfermions +
    sneutrinos + stops + sbots;

  return (fermions + sfermions + higgs + neutralinos + charginos)
    / (16.0 * sqr(PI));
}


double MssmSoftsusy::piZGT(double p, double q) const {
  drBarPars tree(displayDrBarPars());
  double alphaMsbar = dataSet.displayAlpha(ALPHA);
  double alphaDrbar = qedSusythresh(alphaMsbar, displayMu());

  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    mHc = tree.mhiggs(4);
  double    mtop = displayDataSet().displayPoleMt();
  double    mb   =  tree.mb;
  double    mtau =  tree.mtau;
  double    ms   =  displayDataSet().displayMass(mStrange) ;
  double    mc   =  displayDataSet().displayMass(mCharm) ;
  double    mmu  =  displayDataSet().displayMass(mMuon) ;
  double    mE  =  displayDataSet().displayMass(mElectron) ;
  double    mD  =  displayDataSet().displayMass(mDown) ;
  double    mU  =  displayDataSet().displayMass(mUp);
  double    thetat = tree.thetat ;
  double    thetab = tree.thetab;
  double    thetatau= tree.thetatau ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    ct2     = sqr(ct) ;
  double    cb2     = sqr(cb) ;
  double    ctau2   = sqr(ctau) ;
  double    st2     = (1.0 - ct2);
  double    sb2     = (1.0 - cb2);
  double    stau2   = (1.0 - ctau2);
  double    msbot1  = tree.md(1, 3);
  double    msbot2  = tree.md(2, 3);
  double    mstau1  = tree.me(1, 3);
  double    mstau2  = tree.me(2, 3);
  double    mstop1  = tree.mu(1, 3);
  double    mstop2  = tree.mu(2, 3);
  double    g       = displayGaugeCoupling(2);

  double ans = 0.0;

  ans = ans +
    (12.0 * sw2DRbar - 10.0) * b22bar(p, displayMwRun(), displayMwRun(), q) -
    2.0 * (sqr(displayMwRun()) +
           2.0 * cw2DRbar * sqr(p)) * b0(p, displayMwRun(), displayMwRun(), q);

  ans = ans + 2.0 * (guL - guR) * (4.0 * b22bar(p, mU, mU, q) +
                                   sqr(p) * b0(p, mU, mU, q));
  ans = ans + 2.0 * (guL - guR) * (4.0 * b22bar(p, mc, mc, q) +
                                   sqr(p) * b0(p, mc, mc, q));
  ans = ans + 2.0 * (guL - guR) * (4.0 * b22bar(p, mtop, mtop, q) +
                                   sqr(p) * b0(p, mtop, mtop, q));
  ans = ans - (gdL - gdR) * (4.0 * b22bar(p, mD, mD, q) +
                             sqr(p) * b0(p, mD, mD, q));
  ans = ans - (gdL - gdR) * (4.0 * b22bar(p, ms, ms, q) +
                             sqr(p) * b0(p, ms, ms, q));
  ans = ans - (gdL - gdR) * (4.0 * b22bar(p, mb, mb, q) +
                             sqr(p) * b0(p, mb, mb, q));
  ans = ans - (geL - geR) * (4.0 * b22bar(p, mE, mE, q) +
                             sqr(p) * b0(p, mE, mE, q));
  ans = ans - (geL - geR) * (4.0 * b22bar(p, mmu, mmu, q) +
                             sqr(p) * b0(p, mmu, mmu, q));
  ans = ans - (geL - geR) * (4.0 * b22bar(p, mtau, mtau, q) +
                             sqr(p) * b0(p, mtau, mtau, q));

  ans = ans - 2.0 * cos(2.0 * thetaWDRbar) * b22bar(p, mHc, mHc, q);

  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  int i; for (i=1; i<=2; i++) {
    ans = ans + 0.5 * (sqr(v(i, 1).mod()) + sqr(u(i, 1).mod())
                       + 2.0 * cos(2.0 * thetaWDRbar))
      * (4.0 * b22bar(p, mch(i), mch(i), q) +
         sqr(p) * b0(p, mch(i), mch(i), q));

    ans = ans -
      8.0 * (guL * b22bar(p, tree.mu(1, i), tree.mu(1, i), q) -
             guR * b22bar(p, tree.mu(2, i), tree.mu(2, i), q)) +
      4.0 * (gdL * b22bar(p, tree.md(1, i), tree.md(1, i), q) -
             gdR * b22bar(p, tree.md(2, i), tree.md(2, i), q)) +
      4.0 * (geL * b22bar(p, tree.me(1, i), tree.me(1, i), q) -
             geR * b22bar(p, tree.me(2, i), tree.me(2, i), q));
  }

  ans = ans -  8.0 *
    ((guL * ct2 - guR * st2) * b22bar(p, mstop1, mstop1, q) +
     (guL * st2 - guR * ct2) * b22bar(p, mstop2, mstop2, q));
  ans = ans +  4.0 *
    ((gdL * cb2 - gdR * sb2) * b22bar(p, msbot1, msbot1, q) +
     (gdL * sb2 - gdR * cb2) * b22bar(p, msbot2, msbot2, q));
  ans = ans +  4.0 *
    ((geL * ctau2 - geR * stau2) * b22bar(p, mstau1, mstau1, q) +
     (geL * stau2 - geR * ctau2) * b22bar(p, mstau2, mstau2, q));

  double e = sqrt(alphaDrbar * 4.0 * PI);

  return ans * e * g / (16.0 * sqr(PI) * cos(thetaWDRbar));
}

double rho2(double r) {// checked
  if (r <= 1.9)
    return 19.0 - 16.5 * r + 43.0 * sqr(r) / 12.0 + 7.0 / 120.0 * sqr(r) * r -
      PI * sqrt(r) * (4.0 - 1.5 * r + 3.0 / 32.0 * sqr(r) + sqr(r) * r /
                      256.0) - sqr(PI) * (2.0 - 2.0 * r + 0.5 * sqr(r)) -
      log(r) * (3.0 * r - 0.5 * sqr(r));
  else {
    double rm1 = 1.0 / r, rm2 = sqr(rm1), rm3 = rm2 * rm1, rm4 = rm3 * rm1,
      rm5 = rm4 * rm1;
    return sqr(log(r)) * (1.5 - 9.0 * rm1 - 15.0 * rm2 - 48.0 * rm3 - 168.0
                          * rm4 - 612.0 * rm5) -
      log(r) * (13.5 + 4.0 * rm1 - 125.0 / 4.0 * rm2 - 558.0 / 5.0 * rm3 -
                8307.0 / 20.0 * rm4 - 109321.0 / 70.0 * rm5)
      + sqr(PI) * (1.0 - 4.0 * rm1 - 5.0 * rm2 - 16.0 * rm3 -
                   56.0 * rm4 - 204.0 * rm5)
      + 49.0 / 4.0 + 2.0 / 3.0 * rm1 + 1613.0 / 48.0 * rm2 + 87.57 * rm3 +
      341959.0 / 1200.0 * rm4 + 9737663.0 / 9800.0 * rm5;
  }
}

double fEff(double x) {
  double arg = 1.0 / (1.0 + x);
  return 2.0 / x + 3.5 - (3.0 + 2.0 / x) * log(x) +
    sqr(1.0 + 1.0 / x) *
    (2.0 * dilog(arg) - sqr(PI) / 3.0 + sqr(log(1.0 + x)));
}

double gEff(double x) {
  double y = sqrt(x / (4.0 - x));
  return (1.0 / x + 0.5) * (atan(y) / y - 1.0) + 9.0 / 8.0 + 0.5 / x -
    (1.0 + 0.5 / x) * 4.0 / x * sqr(atan(y));
}


double MssmSoftsusy::sinSqThetaEff() {
  if (displayMu() != MZ) {
    throw("Should call MssmSoftsusy::sinSqThetaEff() at MZ only\n");
  }
  double kl = 0.;

  calcDrBarPars();

  double alphaMsbar = dataSet.displayAlpha(ALPHA);
  double alphaDrbar = qedSusythresh(alphaMsbar, displayMu());
  double sinthDrbar = calcSinthdrbar();
  double costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double costhW     = MW / MZ;

  double vLmzSq = 0.5 * fEff(1.0 / sqr(costhW)) +
    4.0 * sqr(costhDrbar) * gEff(1.0 / sqr(costhW)) -
    (1.0 - 6.0 * sqr(sinthDrbar) + 8.0 * sqr(sinthDrbar) * sqr(sinthDrbar)) /
    (4.0 * sqr(costhDrbar)) * fEff(1.0);

  kl = 1.0 + costhDrbar / sinthDrbar *
    (piZGT(MZ, displayMu()) - piZGT(0., displayMu())) / sqr(MZ) +
     alphaDrbar / PI * sqr(costhDrbar) / sqr(sinthDrbar) * 2.0 * log(costhW)
     - alphaDrbar / (4.0 * PI * sqr(sinthDrbar)) * vLmzSq;

  return sqr(sinthDrbar) * kl;
}

// outrho, outsin represent the DRbar values
double MssmSoftsusy::deltaVb(double outrho, double outsin,
                            double alphaDRbar, double pizztMZ) const {
  drBarPars tree(displayDrBarPars());

  double g       = displayGaugeCoupling(2);
  double gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double costh   = (displayMw() / displayMz());
  double cw2     = sqr(costh) ;
  double sw2     = (1.0 - cw2);
  double outcos  = sqrt(1.0 - sqr(outsin));
  double q       = displayMu();

  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz);
  DoubleVector mch(tree.mchBpmz);

  double deltaVbSm = outrho * alphaDRbar / (4.0 * PI * sqr(outsin)) *
    (6.0 + log(cw2) / sw2 *
     (3.5 - 2.5 * sw2 - sqr(outsin) * (5.0 - 1.5 * cw2 / sqr(outcos))));

  DoubleVector bPsi0NuNul(4), bPsicNuSell(2);
  DoubleVector bPsi0ESell(4), aPsicESnul(2);
  ComplexVector bChi0NuNul(4), bChicNuSell(2);
  ComplexVector bChi0ESell(4), aChicESnul(2);

  bPsicNuSell(1) = g;
  bPsi0NuNul(2) = root2 * g * 0.5;
  bPsi0NuNul(1) = - gp / root2;
  aPsicESnul(1) = g;
  bPsi0ESell(1) = -gp / root2;
  bPsi0ESell(2) = -g * root2 * 0.5;

  bChicNuSell = u * bPsicNuSell;
  bChi0ESell =  n * bPsi0ESell;
  bChi0NuNul = n * bPsi0NuNul;

  aChicESnul = v.complexConjugate() * aPsicESnul;

  double deltaZnue = 0.0, deltaZe = 0.0;
  double mselL = tree.me(1, 1),
    msnue = tree.msnu(1);
  int i; for(i=1; i<=4; i++) {
   if (i < 3) {
      deltaZnue = deltaZnue -
        sqr(bChicNuSell(i).mod()) * b1(0.0, mch(i), mselL, q);
      deltaZe = deltaZe -
        sqr(aChicESnul(i).mod()) * b1(0.0, mch(i), msnue, q);
    }
    deltaZnue = deltaZnue -
      sqr(bChi0NuNul(i).mod()) * b1(0.0, mneut(i), msnue, q);
    deltaZe = deltaZe -
      sqr(bChi0ESell(i).mod()) * b1(0.0, mneut(i), mselL, q);
  }

  DoubleVector bPsicNuSmul(2);
  DoubleVector bPsi0MuSmul(4), aPsicMuSnul(2);
  ComplexVector bChicNuSmul(2);
  ComplexVector bChi0MuSmul(4), aChicMuSnul(2);

  double hmu = displayYukawaElement(YE, 2, 2);
  bPsicNuSmul(1) = g;
  bPsicNuSmul(2) = -hmu;
  aPsicMuSnul(1) = g;
  aPsicMuSnul(2) = -hmu;
  bPsi0MuSmul(1) = -gp / root2;
  bPsi0MuSmul(2) = -g * root2 * 0.5;

  bChicNuSmul = u * bPsicNuSmul;
  bChi0MuSmul =  n * bPsi0MuSmul;
  bChi0NuNul = n * bPsi0NuNul;
  aChicMuSnul = v.complexConjugate() * aPsicMuSnul;

  double deltaZnumu = 0.0, deltaZmu = 0.0;
  double msnumu = tree.msnu(2),
    msmuL = tree.me(1, 2);
  for(i=1; i<=4; i++) {
    if (i < 3) {
      deltaZnumu = deltaZnumu -
        sqr(bChicNuSmul(i).mod()) * b1(0.0, mch(i), msmuL, q);
      deltaZmu = deltaZmu -
        sqr(aChicMuSnul(i).mod()) * b1(0.0, mch(i), msnumu, q);
    }
    deltaZnumu = deltaZnumu -
      sqr(bChi0NuNul(i).mod()) * b1(0.0, mneut(i), msnumu, q);
    deltaZmu = deltaZmu -
      sqr(bChi0MuSmul(i).mod()) * b1(0.0, mneut(i), msmuL, q);
  }

  DoubleMatrix aPsi0PsicW(4, 2), bPsi0PsicW(4, 2), fW(4, 2), gW(4, 2);
  ComplexMatrix aChi0ChicW(4, 2), bChi0ChicW(4, 2);

  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;
  bPsi0PsicW(3, 2) = -g / root2;

  // These ought to be in physpars
  aChi0ChicW = n.complexConjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.hermitianConjugate();

  Complex deltaVE = 0.0;
  int j; for(i=1; i<=2; i++)
    for(j=1; j<=4; j++) {
      deltaVE = deltaVE + bChicNuSell(i) * bChi0ESell(j).conj() *
        (- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
         c0(mselL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
         bChi0ChicW(j, i) *
         (b0(0.0, mch(i), mneut(j), q) + sqr(mselL) *
          c0(mselL, mch(i), mneut(j)) - 0.5));
      deltaVE = deltaVE - aChicESnul(i) * bChi0NuNul(j) *
        (- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
         c0(msnue, mch(i), mneut(j)) + 1.0 / (root2 * g) *
         aChi0ChicW(j, i) *
         (b0(0.0, mch(i), mneut(j), q) + sqr(msnue) *
          c0(msnue, mch(i), mneut(j)) - 0.5));
      if (i == 1)
        deltaVE = deltaVE +
          0.5 * bChi0ESell(j).conj() * bChi0NuNul(j) *
          (b0(0.0, mselL, msnue, q) + sqr(mneut(j)) * c0(mneut(j), mselL,
                                                         msnue) + 0.5);
    }

  Complex deltaVMu = 0.0;
  for(i=1; i<=2; i++)
    for(j=1; j<=4; j++) {
      deltaVMu = deltaVMu + bChicNuSmul(i) * bChi0MuSmul(j).conj() *
        (- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
         c0(msmuL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
         bChi0ChicW(j, i) *
         (b0(0.0, mch(i), mneut(j), q) + sqr(msmuL) *
          c0(msmuL, mch(i), mneut(j)) - 0.5));
      deltaVMu = deltaVMu - aChicMuSnul(i) * bChi0NuNul(j) *
        (- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
         c0(msnumu, mch(i), mneut(j)) + 1.0 / (root2 * g) *
         aChi0ChicW(j, i) *
         (b0(0.0, mch(i), mneut(j), q) + sqr(msnumu) *
          c0(msnumu, mch(i), mneut(j)) - 0.5));
      if (i == 1)
        deltaVMu = deltaVMu +
          0.5 * bChi0MuSmul(j).conj() * bChi0NuNul(j) *
          (b0(0.0, msmuL, msnumu, q) + sqr(mneut(j)) * c0(mneut(j), msmuL,
                                                          msnumu) + 0.5);
    }

  Complex a1(0.0, 0.0);
  for(i=1; i<=2; i++)
    for(j=1; j<=4; j++) {
      a1 = a1 + 0.5 * aChicMuSnul(i) * bChicNuSell(i).conj() *
        bChi0NuNul(j) * bChi0ESell(j) * mch(i) * mneut(j) *
        d0(mselL, msnumu, mch(i), mneut(j));
      a1 = a1 + 0.5 * aChicESnul(i).conj() * bChicNuSmul(i) *
        bChi0NuNul(j).conj() * bChi0MuSmul(j).conj() * mch(i) * mneut(j) *
        d0(msmuL, msnue, mch(i), mneut(j));
      a1 = a1 + bChicNuSmul(i) * bChicNuSell(i).conj() *
        bChi0MuSmul(j).conj() * bChi0ESell(j) *
        d27(msmuL, mselL, mch(i), mneut(j));
      a1 = a1 + aChicMuSnul(i).conj() * aChicESnul(i) *
        bChi0NuNul(j) * bChi0NuNul(j).conj() *
        d27(msnumu, msnue, mch(i), mneut(j));
    }

  double deltaVbSusy =
    (-sqr(outsin) * sqr(outcos) / (2.0 * PI * alphaDRbar) * sqr(displayMz())
     * a1.real() + deltaVE.real() + deltaVMu.real() +
     0.5 * (deltaZe + deltaZnue + deltaZmu + deltaZnumu) ) /
    (sqr(PI) * 16.0);

  double deltaVb = deltaVbSusy + deltaVbSm;

  return deltaVb;
}

double MssmSoftsusy::dRho(double outrho, double outsin, double alphaDRbar,
                             double pizztMZ, double piwwtMW) {
  double mz = displayMz();
  drBarPars tree(displayDrBarPars());

  // 2 loop SM contribution
  double mt   = dataSet.displayPoleMt();
  double sinb = sin(atan(displayTanb()));

  double xt = 3.0 * GMU * sqr(mt) / (8.0 * sqr(PI) * root2);

  double deltaRho2LoopSm = alphaDRbar * sqr(displayGaugeCoupling(3)) /
    (16.0 * PI * sqr(PI) * sqr(outsin)) * // bug-fixed 24.08.2002
    (-2.145 * sqr(mt) / sqr(displayMw()) + 1.262 * log(mt / mz) - 2.24
     - 0.85 * sqr(mz)
     / sqr(mt)) + sqr(xt) * sqr(cos(tree.thetaH)) / sqr(sinb) *
    rho2(tree.mhiggs(1) / mt) / 3.0;

  double deltaRhoOneLoop = pizztMZ / (outrho * sqr(mz))
    - piwwtMW / sqr(displayMw());

  double deltaRho = deltaRhoOneLoop + deltaRho2LoopSm;

  return deltaRho;
}

double MssmSoftsusy::dR(double outrho, double outsin, double alphaDRbar,
                           double pizztMZ, double piwwt0) {
  drBarPars tree(displayDrBarPars());

  double outcos = cos(asin(outsin));
  // 2 loop SM contribution
  double mt   = dataSet.displayPoleMt();
  double sinb = sin(atan(displayTanb()));

  double xt = 3.0 * GMU * sqr(mt) / (8.0 * sqr(PI) * root2);

  double dvb = deltaVb(outrho, outsin, alphaDRbar, pizztMZ);

  double mz = displayMz();

  double deltaR =  outrho * piwwt0 / sqr(displayMw()) -
    pizztMZ / sqr(mz) + dvb;

  // Dominant two-loop SM term
  double deltaR2LoopSm = alphaDRbar * sqr(displayGaugeCoupling(3)) /
    (16.0 * sqr(PI) * PI * sqr(outsin) * sqr(outcos)) *
    (2.145 * sqr(mt) / sqr(mz) + 0.575 * log(mt / mz) - 0.224
     - 0.144 * sqr(mz) / sqr(mt)) -
    sqr(xt) * sqr(cos(tree.thetaH)) / sqr(sinb) *
    rho2(tree.mhiggs(1) / mt) * (1.0 - deltaR) * outrho / 3.0;

  deltaR = deltaR + deltaR2LoopSm;

  return deltaR;
}

// Checked 20.11.00
// Flags noconvergence if there's trouble...then don't believe outrho and
// outsin produced - they are fudged!
void MssmSoftsusy::rhohat(double & outrho, double & outsin, double alphaDRbar,
                          double pizztMZ, double piwwt0, double piwwtMW,
                          double tol, int maxTries) {

  static double oldrho = 0.23, oldsin = 0.8;

  double mz = displayMz();
  if (displayMu() != mz) {
    ostringstream ii;
    ii << "Called MssmSoftsusy::rhohat "
       << "with scale" << displayMu() << endl;
    throw ii.str();
  }

  static int numTries = 0;

  if ((outrho < TOLERANCE || outsin < TOLERANCE) ||
      (numTries - 1 > maxTries)) {
    oldrho = 0.23, oldsin = 0.8;
    numTries = 0;
    flagNoRhoConvergence(true);
    if (PRINTOUT) cout << flush << "rhohat reached maxtries\n";
    return;
  }

  // Difference to last iteration
  double sumTol;
  sumTol = fabs(oldrho / outrho - 1.0) + fabs(oldsin / outsin - 1.0);

  if (numTries != 0 && sumTol < tol) {
    numTries = 0;
    oldrho = 0.23, oldsin = 0.8;
    if (PRINTOUT > 2) cout << "sin rho converged\n";
    return;
  }

  numTries = numTries + 1;

  oldrho = outrho; oldsin = outsin;

  if (PRINTOUT > 2) cout << " st2=" << sumTol;

  double deltaR = dR(outrho, outsin, alphaDRbar, pizztMZ, piwwt0);

  double sin2thetasqO4 = PI * alphaDRbar /
    (root2 * sqr(mz) * GMU * (1.0 - deltaR));

  if (sin2thetasqO4 >= 0.25) sin2thetasqO4 = 0.25;
  if (sin2thetasqO4 < 0.0) sin2thetasqO4 = 0.0;

  double sin2theta = sqrt(4.0 * sin2thetasqO4);

  double theta = 0.5 * asin(sin2theta);

  outsin = sin(theta);

  double deltaRho = dRho(outrho, outsin, alphaDRbar, pizztMZ, piwwtMW);

  if (fabs(deltaRho) < 1.0) outrho = 1.0 / (1.0 - deltaRho);
  else outrho = 1.0;

  if (PRINTOUT > 2) cout << " drho=" << deltaRho << " sw=" << outsin << endl;

  rhohat(outrho, outsin, alphaDRbar, pizztMZ, piwwt0, piwwtMW, tol, maxTries);
}

void MssmSoftsusy::methodBoundaryCondition(const DoubleVector & pars) {
  ostringstream ii;
  ii << "Should only use MssmSoftsusy::methodBoundaryCondition in derived"
     << " objects.\n";
  throw ii.str();
}

// isajet routine for getting alpha_s - only used here for getting numbers
// from isawig interface
double sualfs(double qsq, double alam4, double tmass) {
  double anf = 6.0, bmass = 5.0, alam = 0.1, alam5 = 0.1, b0 = 0.;
  double alamsq = 0., b1 = 0., b2 = 0., x = 0., tt = 0., t = 0., sualfs = 0.;
  if (qsq < 4.0 * sqr(bmass)) {
    anf   = 4.0;
    alam  = alam4;
  }
  else if (qsq < 4.0 * sqr(tmass)) {
    anf   = 5.0;
    alam  = alam4 * pow(alam4 / (2.0*bmass), 2.0 / 23.0)
      * pow(log(4.0 * sqr(bmass) / sqr(alam4)), -963.0 / 13225.0);
  }
  else {
    anf   = 6.0;
        alam5 = alam4 * pow(alam4/(2.0*bmass), 2.0/23.0)
          * pow(log(4.0*sqr(bmass)/sqr(alam4)), -963.0/13225.0);
        alam  = alam5 * pow(alam5/(2.0*tmass), 2.0/21.0)
          * pow(log(4.0*sqr(tmass)/sqr(alam5)), -107.0/1127.0);
  }
  b0       = 11.0-2.0/3.0*anf;
  alamsq   = sqr(alam);
  t        = log(qsq/alamsq);
  if (t <= 1.0) t = log(4.0/alamsq);
  double alphas   = 4*PI/b0/t;
  b1 = 102.0-38.0/3.0*anf;
  b2 = 0.5*(2857.0-5033.0/9.0*anf+325.0/27.0*sqr(anf));
  x  = b1/(sqr(b0)*t);
  tt = log(t);
  sualfs = alphas*(1.0-x*tt+sqr(x)*(sqr(tt-0.5)
                                  +b2*b0/sqr(b1)-1.25));
  return sualfs;
}

// another function which mimics Isajet's handling of quark masses
double ssmqcd(double dm, double dq, double mtopPole) {
  double dlam4=0.177;
  double dqbt=10.;
  double dqtp=2*mtopPole;
  double ssmqcd=0., renorm = 0.;
  int dneff=4;
  double po=12.0/(33.0-2.*dneff);
  if(dq < dqbt) {
    renorm=pow(log(2*dm/dlam4)/log(dq/dlam4), po);
    ssmqcd=renorm*dm;
    return ssmqcd;
  }
  else
    renorm=pow(log(2*dm/dlam4)/log(dqbt/dlam4), po);

  dneff=5;
  po=12.0/(33.0-2.*dneff);
  double dlam5=exp((25.0*log(dlam4)-log(sqr(dqbt)))/23.0);
  if(dq >= dqbt && dq < dqtp) {
    renorm=renorm *pow(log(dqbt/dlam5)/log(dq/dlam5), po);
    ssmqcd=renorm*dm;
    return ssmqcd;
  }
  else {
    renorm=renorm
      *pow(log(dqbt/dlam5)/log(dqtp/dlam5), po);
  }

  dneff=6;
  po=12.0/(33.0-2.*dneff);
  double dlam6=exp((25.0*log(dlam4)-log(sqr(dqbt))
              -log(4*sqr(mtopPole)))/21.0);
  renorm=renorm
    *pow(log(dqtp/dlam6)/log(dq/dlam6), po);
  ssmqcd=renorm*dm;
  return ssmqcd;
}

// Works out how best to fit the isajet numbers to the spectrum.
// There are problems with the Higgs and sbottoms because ISAJET assumes
// certain tree-level relations between masses that are broken by SOFTSUSY's
// higher accuracy. The differences get large for high tan beta around 50, at
// around 10 they're typically only a percent.
void MssmSoftsusy::isajetNumbers764
(double & mtopPole, double & mGPole, double & smu, double & mA, double & tanb,
 double & mq1l, double & mdr, double & mur, double & meL, double & meR,
 double & mql3, double & mdr3, double & mur3, double &  mtauL,
 double & mtauR, double & at, double & ab, double & atau, double & mq2l,
 double & msr, double & mcr, double & mmuL, double & mmuR, double & m1,
 double & m2)
  const {

  // Store a copy of the current object
  MssmSoftsusy store(*this);

  // Run to MSUSY
  store.runto(displayMsusy(), TOLERANCE * 0.01);

  DoubleMatrix mNeut(4, 4);

  // Store current object in vector
  DoubleVector storeObject(store.display());
  double muInitial = store.displayMu();

  m1 = store.displayGaugino(1);
  m2 = store.displayGaugino(2);
  smu = store.displaySusyMu();

  // restore object at MSUSY
  store.setMu(muInitial);
  store.set(storeObject);

  // tree level
  mNeut(1, 1) = m1;
  mNeut(2, 2) = m2;
  mNeut(3, 4) = - smu;
  store.addNeutralinoLoop(mNeut);

  m1 = fabs(mNeut(1, 1));
  m2 = fabs(mNeut(2, 2));
  smu = -mNeut(3, 4);

  double beta = atan(store.displayTanb());

  // define ISAJET constants
  double amdn=0.0099;
  double  amup=0.0056;
  double  amst=0.199;
  double  amch=1.35;
  double  ambt=5.0;
  double  ame=0.511e-3;
  double  ammu=0.105;
  double  amz=91.17;
  double  amw = 80.0;
  double  sn2thw=0.232;

  mq1l = sqrt(sqr(displayPhys().mu(1, 1)) - sqr(amup) -
              (0.5 - 2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));
  mq2l = sqrt(sqr(displayPhys().mu(1, 2)) - sqr(amch) -
              (0.5 - 2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));

  mur = sqrt(sqr(displayPhys().mu(2, 1)) - sqr(amup) -
              (2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));
  mcr = sqrt(sqr(displayPhys().mu(2, 2)) - sqr(amch) -
              (2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));

  mdr = sqrt(sqr(displayPhys().md(2, 1)) - sqr(amdn) -
              (- 1./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));
  msr = sqrt(sqr(displayPhys().md(2, 2)) - sqr(amst) -
              (- 1./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));

  meL = sqrt(sqr(displayPhys().me(1, 1)) - sqr(ame) +
              (0.5 - sn2thw) * sqr(amz) * cos(2.0 * beta));
  mmuL = sqrt(sqr(displayPhys().me(1, 2)) - sqr(ammu) +
              (0.5 - sn2thw) * sqr(amz) * cos(2.0 * beta));

  meR = sqrt(sqr(displayPhys().me(2, 1)) - sqr(ame) -
              (- sn2thw) * sqr(amz) * cos(2.0 * beta));
  mmuR = sqrt(sqr(displayPhys().me(2, 2)) - sqr(ammu) -
              (- sn2thw) * sqr(amz) * cos(2.0 * beta));

  mtopPole = store.displayDataSet().displayPoleMt();

  double asmt = sualfs(sqr(mtopPole), 0.36, mtopPole);
  double mtmt = mtopPole/(1.+4*asmt/3./PI+(16.11-1.04*(5.-6.63/mtopPole))*
                          sqr(asmt/PI));

  double qsusy = displayMsusy();
  DoubleMatrix mStopSquared(2, 2);
  DoubleMatrix m(2, 2);
  m(1, 1) = sqr(displayPhys().mu(1, 3));
  m(2, 2) = sqr(displayPhys().mu(2, 3));
  mStopSquared = rot2d(-displayPhys().thetat) * m *
    rot2d(displayPhys().thetat);

  double mtrun;

  // Loop is for iterative solution
  int i; for (i=1; i<=3; i++) {
    mtrun = ssmqcd(mtmt, qsusy, mtopPole);

    mql3 = sqrt(mStopSquared(1, 1) - sqr(mtrun) -
                (0.5 - 2.0 / 3.0 * sn2thw) * sqr(amz) * cos(2.0 * beta));
    mur3 = sqrt(mStopSquared(2, 2) - sqr(mtrun) -
                2.0 / 3.0 * sn2thw * sqr(amz) * cos(2.0 * beta));
    qsusy = sqrt(mql3 * mur3);
  }

  DoubleMatrix mSbotSquared(2, 2);
  m(1, 1) = sqr(displayPhys().md(1, 3));
  m(2, 2) = sqr(displayPhys().md(2, 3));
  mSbotSquared = rot2d(-displayPhys().thetab) * m *
    rot2d(displayPhys().thetab);

  double asmb = sualfs(sqr(ambt), 0.36, mtopPole);
  double mbmb = ambt * (1.0 - 4.0 * asmb / (3.0 * PI));
  double mbq = ssmqcd(mbmb, qsusy, mtopPole);

  // We're fiddling these variables to reproduce inside isajet the squark mass
  // matrices that we've found in softsusy
  mdr3 = sqrt(mSbotSquared(2, 2) - sqr(mbq)
              + 1.0 / 3.0 * sn2thw * sqr(amz) * cos(2.0 * beta));

  DoubleMatrix mStauSquared(2, 2);
  m(1, 1) = sqr(displayPhys().me(1, 3));
  m(2, 2) = sqr(displayPhys().me(2, 3));
  mStauSquared = rot2d(-displayPhys().thetatau) * m *
    rot2d(displayPhys().thetatau);

  double mlq = 1.7463;

  mtauL = sqrt(mStauSquared(1, 1) - sqr(mlq) + (sqr(amw) -  sqr(amz) * 0.5) *
               cos(2.0 * beta));
  mtauR = sqrt(mStauSquared(2, 2) - sqr(mlq) - (sqr(amw) - sqr(amz)) *
               cos(2.0 * beta));

  at  = mStopSquared(1, 2) / mtrun + smu / tan(beta);
  ab  = mSbotSquared(1, 2) / mbq + smu * tan(beta);
  atau = mStauSquared(1, 2) / mlq + smu * tan(beta);

  mGPole   = store.displayPhys().mGluino;
  mA       = store.displayPhys().mhiggs(2);

  tanb = store.displayTanb();
}

// First name input is the name of an OUTPUT file from ssrun, the second name
// is the name of the interface file for INPUT to ssrun
void MssmSoftsusy::ssrunInterface764(char fname [80],
                                           char softfname [80])
  const {
  fstream softOutput(softfname, ios::out);

  ssrunInterface764Inside(fname, softOutput);

  softOutput.close();
}

void MssmSoftsusy::ssrunInterface764Inside(char fname [80],
                                           fstream & softOutput)
  const {

  double mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur,
    meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau,
    mq2l, msr, mcr, mmuL, mmuR, m1, m2;

  isajetNumbers764(mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur,
                   meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau,
                   mq2l, msr, mcr, mmuL, mmuR, m1, m2);

  softOutput << "'" << fname << "'" << endl;

  softOutput << mtopPole << endl;

  softOutput << mGPole  << "," << smu << "," << mA << "," << tanb
       << endl;

  softOutput << mq1l    << "," << mdr << "," << mur << ","
       << meL     << "," << meR << endl;

  softOutput << mql3    << "," << mdr3  << "," << mur3 << ","
       << mtauL   << "," << mtauR << ","
       << at      << "," << ab    << "," << atau
       << endl;

  softOutput << mq2l << "," << msr << "," << mcr << ","
       << mmuL << "," << mmuR
       << endl;

  softOutput << m1 << "," << m2 << endl << "/" << endl;

}

void MssmSoftsusy::isawigInterface764(char herwigInputFile [80],
                                      char isajetOutputFile [80],
                                      char softOutputFile [80])
  const {

  fstream softOutput(softOutputFile, ios::out);

  ssrunInterface764Inside(isajetOutputFile, softOutput);

  softOutput << "n" << endl << "n" << endl;

  softOutput << "'" << herwigInputFile << "'" << endl;

  softOutput.close();
}

// For input into isajet parameter file called fname
void MssmSoftsusy::isajetInterface764(char fname[80]) const {

  fstream softOutput(fname, ios::out);

  double mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur,
    meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau,
    mq2l, msr, mcr, mmuL, mmuR, m1, m2;

  isajetNumbers764(mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur,
                   meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau,
                   mq2l, msr, mcr, mmuL, mmuR, m1, m2);

  // PRINT out the fiddle in ISAJET PARAMETER file format
  softOutput << "TMASS"  << endl;
  softOutput << mtopPole << "/" << endl;

  softOutput << "MSSMA" << endl;
  softOutput << mGPole  << "," << smu << "," << mA << "," << tanb
       << "/"     << endl;

  softOutput << "MSSMB" << endl;
  softOutput << mq1l    << "," << mdr << "," << mur << ","
       << meL     << "," << meR << "/" << endl;

  softOutput << "MSSMC" << endl;
  softOutput << mql3    << "," << mdr3  << "," << mur3 << ","
       << mtauL   << "," << mtauR << ","
       << at      << "," << ab    << "," << atau
       << "/" << endl;

  softOutput << "MSSMD" << endl;
  softOutput << mq2l << "," << msr << "," << mcr << ","
       << mmuL << "," << mmuR
       << "/" << endl;

  softOutput << "MSSME" << endl;
  softOutput << m1 << "," << m2 << "/" << endl;

  softOutput.close();
}

void generalBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  MssmSusy s; SoftParsMssm r;
  double m3sq = m.displayM3Squared();
  s = m.displaySusy();
  r.set(inputParameters);
  r.setM3Squared(m3sq);
  m.setSoftPars(r);
  m.setSusy(s);

  return;
}

// This one doesn't overwrite mh1sq or mh2sq at the high scale
void generalBcs2(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  MssmSusy s; SoftParsMssm r;
  double mh1sq = m.displayMh1Squared();
  double mh2sq = m.displayMh2Squared();
  double m3sq = m.displayM3Squared();
  s = m.displaySusy();
  r.set(inputParameters);
  r.setMh1Squared(mh1sq);
  r.setMh2Squared(mh2sq);
  r.setM3Squared(m3sq);
  m.setSoftPars(r);
  m.setSusy(s);

  return;
}

void MssmSoftsusy::microMegasIn(char fname[80]) {
  fstream softInput(fname, ios::in);
  char c[80];

  double m1, m2, m3, ssmu;
  softInput >> c >> ssmu;
  softInput >> c >> m1;
  softInput >> c >> m2;
  softInput >> c >> m3;
  setSusyMu(ssmu);
  setGauginoMass(1, m1);  setGauginoMass(2, m2);  setGauginoMass(3, m3);

  double ml1, ml2, ml3, me1, me2, me3, mq1, mq2, mq3, mu1, mu2, mu3, md1,
    md2, md3, ht, hb, htau;
  softInput >> c >> ml1;
  softInput >> c >> ml2;
  softInput >> c >> ml3;
  softInput >> c >> me1;
  softInput >> c >> me2;
  softInput >> c >> me3;
  softInput >> c >> mq1;
  softInput >> c >> mq2;
  softInput >> c >> mq3;
  softInput >> c >> mu1;
  softInput >> c >> mu2;
  softInput >> c >> mu3;
  softInput >> c >> md1;
  softInput >> c >> md2;
  softInput >> c >> md3;
  softInput >> c >> ht;
  softInput >> c >> hb;
  softInput >> c >> htau;
  setSoftMassElement(mLl, 1, 1, ml1);
  setSoftMassElement(mLl, 2, 2, ml2);
  setSoftMassElement(mLl, 3, 3, ml3);
  setSoftMassElement(mEr, 1, 1, me1);
  setSoftMassElement(mEr, 2, 2, me2);
  setSoftMassElement(mEr, 3, 3, me3);
  setSoftMassElement(mQl, 1, 1, mq1);
  setSoftMassElement(mQl, 2, 2, mq2);
  setSoftMassElement(mQl, 3, 3, mq3);
  setSoftMassElement(mUr, 1, 1, mu1);
  setSoftMassElement(mUr, 2, 2, mu2);
  setSoftMassElement(mUr, 3, 3, mu3);
  setSoftMassElement(mDr, 1, 1, md1);
  setSoftMassElement(mDr, 2, 2, md2);
  setSoftMassElement(mDr, 3, 3, md3);

  double mh3, tb;
  softInput >> c >> mh3;
  softInput >> c >> tb; setTanb(tb);
  double beta = atan(displayTanb());
  setM3Squared(mh3 * sin(beta) * cos(beta));

  double hmu;
  softInput >> c >> hmu; // needs converting

  double mtq, mbq, mtauq;
  softInput >> c >> mtq;
  softInput >> c >> mbq;
  softInput >> c >> mtauq;

  softInput >> c >> physpars.mhiggs(1) >> c >> physpars.mhiggs(2) >> c
            >> physpars.mhiggs(3) >> c >> physpars.mhiggs(4) >> c
            >> physpars.thetaH;
  softInput >> c >> physpars.msnu(1) >> c >> physpars.msnu(2)
            >> c >> physpars.msnu(3);
  softInput >> c >> physpars.mu >> c >> physpars.md >> c >> physpars.me;
  softInput >> c >> physpars.mGluino;
  softInput >> c >> physpars.mch;
  softInput >> c >> physpars.mneut;
  softInput >> c >> c >> c;
  double q;
  softInput >> c >> q; setMu(q);
  // set Yukawas
  setYukawaElement(YU, 3, 3, mtq * root2 / (displayHvev() * sin(beta)));
  setYukawaElement(YD, 3, 3, mbq * root2 / (displayHvev() * cos(beta)));
  setYukawaElement(YE, 3, 3, mtauq * root2 / (displayHvev() * cos(beta)));
  DoubleVector g1(3);
  softInput >> c >> c >> g1; setAllGauge(g1);
  double mmu;
  softInput >> c >> mmu;
  setYukawaElement(YE, 2, 2, mmu * root2 / (displayHvev() * cos(beta)));

  setTrilinearElement(UA, 3, 3, ht * displayYukawaElement(YU, 3, 3));
  setTrilinearElement(DA, 3, 3, hb * displayYukawaElement(YD, 3, 3));
  setTrilinearElement(EA, 3, 3, htau * displayYukawaElement(YE, 3, 3));
  setTrilinearElement(EA, 2, 2, hmu * displayYukawaElement(YE, 2, 2));

  softInput.close();
}

// outputs to softoutput a micromegas file (the name of which is input)
void MssmSoftsusy::microMegasInterface(char fname[80]) const {
  fstream softOutput(fname, ios::out);

  MssmSoftsusy r(*this);
  double beta = atan(r.displayTanb());
  softOutput << "mu " << r.displaySusyMu() << endl;
  softOutput << "MG1 " << r.displayGaugino(1) << endl;
  softOutput << "MG2 " << r.displayGaugino(2) << endl;
  softOutput << "MG3 " << r.displayGaugino(3) << endl;
  softOutput << "Ml1 " << sqrt(r.displaySoftMassSquared(mLl, 1, 1)) << endl;
  softOutput << "Ml2 " << sqrt(r.displaySoftMassSquared(mLl, 2, 2)) << endl;
  softOutput << "Ml3 " << sqrt(r.displaySoftMassSquared(mLl, 3, 3)) << endl;
  softOutput << "Mr1 " << sqrt(r.displaySoftMassSquared(mEr, 1, 1)) << endl;
  softOutput << "Mr2 " << sqrt(r.displaySoftMassSquared(mEr, 2, 2)) << endl;
  softOutput << "Mr3 " << sqrt(r.displaySoftMassSquared(mEr, 3, 3)) << endl;
  softOutput << "Mq1 " << sqrt(r.displaySoftMassSquared(mQl, 1, 1)) << endl;
  softOutput << "Mq2 " << sqrt(r.displaySoftMassSquared(mQl, 2, 2)) << endl;
  softOutput << "Mq3 " << sqrt(r.displaySoftMassSquared(mQl, 3, 3)) << endl;
  softOutput << "Mu1 " << sqrt(r.displaySoftMassSquared(mUr, 1, 1)) << endl;
  softOutput << "Mu2 " << sqrt(r.displaySoftMassSquared(mUr, 2, 2)) << endl;
  softOutput << "Mu3 " << sqrt(r.displaySoftMassSquared(mUr, 3, 3)) << endl;
  softOutput << "Md1 " << sqrt(r.displaySoftMassSquared(mDr, 1, 1)) << endl;
  softOutput << "Md2 " << sqrt(r.displaySoftMassSquared(mDr, 2, 2)) << endl;
  softOutput << "Md3 " << sqrt(r.displaySoftMassSquared(mDr, 3, 3)) << endl;
  softOutput << "Atop " << r.displaySoftA(UA, 3, 3) << endl;
  softOutput << "Ab " << r.displaySoftA(DA, 3, 3) << endl;
  softOutput << "Atau " << r.displaySoftA(EA, 3, 3) << endl;
  softOutput << "MH3 " << sqrt(r.displayM3Squared() /
             (sin(beta) * cos(beta))) << endl;
  softOutput << "tb " << r.displayTanb() << endl;
  softOutput << "Am " << r.displaySoftA(EA, 2, 2) << endl;

  softOutput << "mt(Q) " << r.displayYukawaElement(YU, 3, 3) *
    displayHvev() * sin(beta) / root2 << endl;
  softOutput << "mb(Q) " << r.displayYukawaElement(YD, 3, 3) *
    displayHvev() * cos(beta) / root2 << endl;
  softOutput << "mtau(Q) " << r.displayYukawaElement(YE, 3, 3) *
    displayHvev() * cos(beta) / root2 << endl;

  softOutput << "mh^0 " << r.displayPhys().mhiggs.display(1) << endl
             << "mA^0 " << r.displayPhys().mhiggs.display(2) << endl
             << "mH^0 " << r.displayPhys().mhiggs.display(3) << endl
             << "mH^+- " << r.displayPhys().mhiggs.display(4) << endl
             << "alpha " << r.displayPhys().thetaH << "\n";

  softOutput << "Msnu1 " << r.displayPhys().msnu(1) << "\nMsnu2 "
             << r.displayPhys().msnu(2)
             << "\nMsnu3 " << r.displayPhys().msnu(3);
  softOutput << "\nmU~" << r.displayPhys().mu << "mD~" << r.displayPhys().md << "mE~" << r.displayPhys().me;
  softOutput << "mGluino  " << r.displayPhys().mGluino << "\n";
  softOutput << "charginos" << r.displayPhys().mch;

  softOutput << "neutralinos" << r.displayPhys().mneut;
  softOutput << "higgs VEV(Q) " << r.displayHvev() << endl;
  softOutput << "Q " << r.displayMu() << endl;
  softOutput << "gauge couplings" << r.displayGauge();
  softOutput << "mmu(Q) " << r.displayYukawaElement(YE, 2, 2) *
    displayHvev() * cos(beta) / root2 << endl;

  softOutput.close();
}


// SUSY Les Houches accord for interfacing to Monte-Carlos, decay programs etc.
void MssmSoftsusy::lesHouchesAccordOutput(char model[],
                                          const DoubleVector & pars,
                                          int sgnMu, double tanb,
                                          double qMax,
                                          int numPoints, double mb,
                                          double mtau, double mgut,
                                          bool altEwsb) {
  outputCharacteristics(8);
  cout << "# SOFTSUSY" << VERSION << endl;
  cout << "# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331,";
  cout << " hep-ph/0104145\n";

  cout << "Block SPINFO         # Program information\n"
       << "     1   SOFTSUSY    # spectrum calculator\n";
  cout << "     2   " << VERSION << "         # version number\n";
  if (displayProblem().noConvergence)
    cout << "     3   Possible problem: Not achieved desired accuracy to factor of 10\n";
  if (displayProblem().testSeriousProblem())
    cout << "     4   Point invalid: " << displayProblem() << endl;

  cout << "Block MODSEL  # Select model\n";
  int modsel = 0;
  if (!strcmp(model, "sugra")) modsel = 1;
  if (!strcmp(model, "gmsb")) modsel = 2;
  if (!strcmp(model, "amsb")) modsel = 3;
  if (!strcmp(model, "mirage")) modsel = 4;

  cout << "     1    " << modsel << "   # " << model << "\n"; // Les Houches
                                                            // accord codes

  QedQcd d(displayDataSet());
  cout << "Block SMINPUTS   # Standard Model inputs\n";
  cout << "     1   "; printRow(1.0 / d.displayAlpha(ALPHA));
  cout << "   # alpha_em^(-1)(MZ) SM MSbar\n";
  cout << "     2   "; printRow(GMU); cout << "   # G_Fermi\n";
  cout << "     3   "; printRow(d.displayAlpha(ALPHAS));
  cout << "   # alpha_s(MZ)MSbar\n";
  cout << "     4   "; printRow(displayMz()); cout << "   # MZ(pole)\n";
  cout << "     5   "; printRow(mb); cout << "   # mb(mb)\n";
  cout << "     6   "; printRow(d.displayPoleMt()); cout << "   # Mtop(pole)\n";
  cout << "     7   "; printRow(mtau); cout << "   # Mtau(pole)\n";

  cout << "Block MINPAR  # SUSY breaking input parameters\n";
  cout << "     3   "; printRow(tanb)            ; cout << "   # tanb" << endl;
  if (!altEwsb) {
    cout << "     4   ";
    printRow(double(sgnMu));
    cout << "   # sign(mu)"<< endl;
  }
  if (!strcmp(model, "sugra")) {
    cout << "     1   "; printRow(pars.display(1)); cout << "   # m0" << endl;
    cout << "     2   "; printRow(pars.display(2)) ; cout << "   # m12" << endl;
    cout << "     5   "; printRow(pars.display(3)) ; cout << "   # A0" << endl;
  }
  else if (!strcmp(model, "gmsb")) {
    cout << "     1   "; printRow(pars.display(3)); cout << "   # lambda" << endl;
    cout << "     2   "; printRow(pars.display(2)) ; cout << "   # M_mess"
         << endl;
    cout << "     5   "; printRow(pars.display(1)) ; cout << "   # N5" << endl;
    cout << "     6   "; printRow(pars.display(4)) ; cout << "   # cgrav" << endl;
  }
  else if (!strcmp(model, "amsb")) {
    cout << "     1   "; printRow(pars.display(1)); cout << "   # m3/2" << endl;
    cout << "     2   "; printRow(pars.display(2)) ; cout << "   # m0"
         << endl;
  }
  else if (!strcmp(model, "mirage")) {
    cout << "     1   "; printRow(pars.display(1)); cout << "   # alpha" << endl;
    cout << "     2   "; printRow(pars.display(2)); cout << "   # M0" <<endl;
    cout << "    51   "; printRow(pars.display(3)); cout << "   # aq" <<endl;
    cout << "    52   "; printRow(pars.display(4)); cout << "   # al" <<endl;
    cout << "    61   "; printRow(pars.display(5)); cout << "   # ahd" <<endl;
    cout << "    62   "; printRow(pars.display(6)); cout << "   # ahu" <<endl;
    cout << "    71   "; printRow(pars.display(7)); cout << "   # cq" <<endl;
    cout << "    72   "; printRow(pars.display(8)); cout << "   # cl" <<endl;
    cout << "    81   "; printRow(pars.display(9)); cout << "   # chd" <<endl;
    cout << "    82   "; printRow(pars.display(10)); cout << "   # chu" <<endl;
  }
  else if (!strcmp(model, "nonUniversal")) {
    cout << "Block EXTPAR          # non-universal SUSY breaking parameters\n";
    int i;
    for (i=1; i<=3; i++) {
      cout << "     " << i << "   ";
      printRow(pars.display(i));
      cout << "         # M_" << i << "(MX)" << endl;
    }
    cout << "     11   "; printRow(pars.display(11)) ;
    cout << "         # At(MX)" << endl;
    cout << "     12   "; printRow(pars.display(12)) ;
    cout << "         # Ab(MX)" << endl;
    cout << "     13   "; printRow(pars.display(13)) ;
    cout << "         # Atau(MX)" << endl;
    if (!altEwsb) {
      cout << "     21   "; printRow(pars.display(21)) ;
      cout << "         # mHd^2(MX)" << endl;
      cout << "     22   "; printRow(pars.display(22)) ;
      cout << "         # mHu^2(MX)" << endl;
    } else {
      cout << "     23   "; printRow(pars.display(23)) ;
      cout << "         # mu(MX)" << endl;
      cout << "     26   "; printRow(pars.display(26)) ;
      cout << "         # mA(pole)" << endl;
    }
    cout << "     31   "; printRow(pars.display(31)) ;
    cout << "         # meL(MX)" << endl;
    cout << "     32   "; printRow(pars.display(32)) ;
    cout << "         # mmuL(MX)" << endl;
    cout << "     33   "; printRow(pars.display(33)) ;
    cout << "         # mtauL(MX)" << endl;
    cout << "     34   "; printRow(pars.display(34)) ;
    cout << "         # meR(MX)" << endl;
    cout << "     35   "; printRow(pars.display(35)) ;
    cout << "         # mmuR(MX)" << endl;
    cout << "     36   "; printRow(pars.display(36)) ;
    cout << "         # mtauR(MX)" << endl;
    cout << "     41   "; printRow(pars.display(41)) ;
    cout << "         # mqL1(MX)" << endl;
    cout << "     42   "; printRow(pars.display(42)) ;
    cout << "         # mqL2(MX)" << endl;
    cout << "     43   "; printRow(pars.display(43)) ;
    cout << "         # mqL3(MX)" << endl;
    cout << "     44   "; printRow(pars.display(44)) ;
    cout << "         # muR(MX)" << endl;
    cout << "     45   "; printRow(pars.display(45)) ;
    cout << "         # mcR(MX)" << endl;
    cout << "     46   "; printRow(pars.display(46)) ;
    cout << "         # mtR(MX)" << endl;
    cout << "     47   "; printRow(pars.display(47)) ;
    cout << "         # mdR(MX)" << endl;
    cout << "     48   "; printRow(pars.display(48)) ;
    cout << "         # msR(MX)" << endl;
    cout << "     49   "; printRow(pars.display(49)) ;
    cout << "         # mbR(MX)" << endl;
  }
  else {
    ostringstream ii;
    ii << "Attempting to use SUSY Les Houches Accord for model "
       << model << " - cannot do at present\n";
    throw ii.str();
  }

  cout << "# Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE="
       << TOLERANCE << endl;
  cout << "# mgut=" << mgut << " GeV\n";

  sPhysical s(displayPhys());

  cout << "Block MASS   # Mass spectrum\n";
  cout << "#PDG code      mass              particle\n";
  // cout << "          6   "; printRow(displayDataSet().displayPoleMt());
  // cout << "   # top\n";
  cout << "        24    "; printRow(displayMw()); cout << "   # MW\n";
  cout << "        25    "; printRow(s.mhiggs(1)); cout << "   # h0\n";
  // cout << "         32   "; printRow(displayMz()); cout << "   # MZ\n";
  cout << "        35    "; printRow(s.mhiggs(3)); cout << "   # H0\n";
  cout << "        36    "; printRow(s.mhiggs(2)); cout << "   # A0\n";
  cout << "        37    "; printRow(s.mhiggs(4)); cout << "   # H+\n";
  cout << "   1000001    "; printRow(s.md(1, 1)); cout << "   # ~d_L\n";
  cout << "   1000002    "; printRow(s.mu(1, 1)); cout << "   # ~u_L\n";
  cout << "   1000003    "; printRow(s.md(1, 2)); cout << "   # ~s_L\n";
  cout << "   1000004    "; printRow(s.mu(1, 2)); cout << "   # ~c_L\n";
  cout << "   1000005    "; printRow(minimum(s.md(1, 3), s.md(2, 3)));
  cout << "   # ~b_1\n";   cout << "   1000006    ";
  printRow(minimum(s.mu(1, 3), s.mu(2, 3)));
  cout << "   # ~t_1\n";
  cout << "   1000011    "; printRow(s.me(1, 1)); cout << "   # ~e_L\n";
  cout << "   1000012    "; printRow(s.msnu(1)); cout << "   # ~nue_L\n";
  cout << "   1000013    "; printRow(s.me(1, 2)); cout << "   # ~mu_L\n";
  cout << "   1000014    "; printRow(s.msnu(2)); cout << "   # ~numu_L\n";
  cout << "   1000015    "; printRow(minimum(s.me(1, 3), s.me(2, 3)));
  cout << "   # ~stau_1\n";
  cout << "   1000016    "; printRow(s.msnu(3)); cout << "   # ~nu_tau_L\n";
  cout << "   1000021    "; printRow(s.mGluino); cout << "   # ~g\n";
  cout << "   1000022    "; printRow(s.mneut(1));
  cout << "   # ~neutralino(1)\n";
  cout << "   1000023    "; printRow(s.mneut(2));
  cout << "   # ~neutralino(2)\n";
  cout << "   1000024    "; printRow(s.mch(1)); cout << "   # ~chargino(1)\n";
  cout << "   1000025    "; printRow(s.mneut(3));
  cout << "   # ~neutralino(3)\n";
  cout << "   1000035    "; printRow(s.mneut(4));
  cout << "   # ~neutralino(4)\n";
  cout << "   1000037    "; printRow(s.mch(2)); cout << "   # ~chargino(2)\n";
  cout << "   2000001    "; printRow(s.md(2, 1)); cout << "   # ~d_R\n";
  cout << "   2000002    "; printRow(s.mu(2, 1)); cout << "   # ~u_R\n";
  cout << "   2000003    "; printRow(s.md(2, 2)); cout << "   # ~s_R\n";
  cout << "   2000004    "; printRow(s.mu(2, 2)); cout << "   # ~c_R\n";
  cout << "   2000005    "; printRow(maximum(s.md(1, 3), s.md(2, 3)));
  cout << "   # ~b_2\n";
  cout << "   2000006    "; printRow(maximum(s.mu(1, 3), s.mu(2, 3)));
  cout << "   # ~t_2\n";
  cout << "   2000011    "; printRow(s.me(2, 1)); cout << "   # ~e_R\n";
  cout << "   2000013    "; printRow(s.me(2, 2)); cout << "   # ~mu_R\n";
  cout << "   2000015    "; printRow(maximum(s.me(1, 3), s.me(2, 3)));
  cout << "   # ~stau_2\n";

  const double underflow = 1.0e-120;
  if (fabs(displayGravitino()) > underflow)
    cout << "   1000039     " << displayGravitino() << "   # ~gravitino\n";
  cout << "# Higgs mixing\n";
  double ms = displayMsusy();
  cout << "Block alpha " <<
    "  # Effective Higgs mixing parameter\n";
  cout << "          "; printRow(s.thetaH);        cout << "   # alpha\n";

  DoubleMatrix m(2, 2);

  cout << "Block stopmix  # stop mixing matrix\n";
  if (s.mu(1, 3) < s.mu(2, 3)) m = rot2d(s.thetat);
  else m = rot2dTwist(s.thetat);
  int i, j;
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      cout << "  " << i << "  " << j << "    ";  printRow(m(i, j));
      cout << "   # O_{" << i << j << "}" << endl;
    }

  cout << "Block sbotmix  # sbottom mixing matrix\n";
  if (s.md(1, 3) < s.md(2, 3)) m = rot2d(s.thetab);
  else m = rot2dTwist(s.thetab);
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      cout << "  " << i << "  " << j << "    "; printRow(m(i, j));
      cout << "   # O_{" << i << j << "}" << endl;
    }

  cout << "Block staumix  # stau mixing matrix\n";
  if (s.me(1, 3) < s.me(2, 3)) m = rot2d(s.thetatau);
  else m = rot2dTwist(s.thetatau);
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      cout << "  " << i << "  " << j << "    "; printRow(m(i, j));
      cout << "   # O_{" << i << j << "}" << endl;
   }

  cout << "Block nmix  # neutralino mixing matrix\n";
  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      cout << "  " << i << "  " << j << "    ";
      printRow(s.mixNeut(j, i));
      cout << "   # N_{" << i << "," << j << "}\n";
    }

  DoubleMatrix u(rot2d(s.thetaL)), v(rot2d(s.thetaR));

  i = 1;
  cout << "Block Umix  # chargino U mixing matrix \n";
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      cout << "  " << i << "  " << j << "    "; printRow(u(i, j));
      cout << "   # U_{" << i << "," << j << "}\n";
    }

  cout << "Block Vmix  # chargino V mixing matrix \n";
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      cout << "  " << i << "  " << j << "    "; printRow(v(i, j));
      cout << "   # V_{" << i << "," << j << "}\n";
    }

  int n = 0;

  while (n < numPoints) {
    n++;
    // Starting non-essential information. The following decides what scale to
    // output the running parameters at. It depends upon what qMax is and how
    // many points the user has requested.
    // For qMax = 0 and 1 point (defaults), Q=MSUSY is printed out.
    // For qMax = 0 and n points, points are spaced logarithmically between MZ
    // and MSUSY.
    // For qMax != 0 and 1 point, Q=qMax is printed.
    // For qMax != 0 and n points, points are log spaced between MZ and qMax.
    double q = ms;
    if (numPoints == 1 && qMax < EPSTOL) q = ms;
    else if (numPoints == 1 && qMax > EPSTOL) q = qMax;
    else if (numPoints > 1 && qMax < EPSTOL) qMax = ms;

    if (numPoints > 1) {
      if (n > 1) {
        double logq = (log(qMax) - log(MZ)) * double(n-1) /
          double(numPoints-1) + log(MZ);
        q = exp(logq);
      }
      else q = MZ;
    }

    runto(q);
    double gp = displayGaugeCoupling(1) * sqrt(0.6);
    cout << "Block gauge Q= " << q << "  \n";
    cout << "     1     " << gp << "   # g'(Q)MSSM DRbar"
         << endl;
    cout << "     2     " << displayGaugeCoupling(2)
         << "   # g(Q)MSSM DRbar" << endl;
    cout << "     3     " << displayGaugeCoupling(3) << "   # g3(Q)MSSM DRbar"
         << endl;

    cout << "Block yu Q= " << q << "  \n"
         << "  3  3     " << displayYukawaElement(YU, 3, 3)
         << "   # Yt(Q)MSSM DRbar" << endl;
    cout << "Block yd Q= " << q << "  \n"
         << "  3  3     " << displayYukawaElement(YD, 3, 3)
         << "   # Yb(Q)MSSM DRbar" << endl;
    cout << "Block ye Q= " << q << "  \n"
         << "  3  3     " << displayYukawaElement(YE, 3, 3)
         << "   # Ytau(Q)MSSM DRbar" << endl;

    cout << "Block hmix Q= " << q <<
      "  # Higgs mixing parameters\n";
    cout << "     1    "; printRow(displaySusyMu());
    cout << "   # mu(Q)MSSM DRbar\n";
    cout << "     2    "; printRow(displayTanb());
    cout << "   # tan beta(Q)MSSM DRbar\n";
    cout << "     3    "; printRow(displayHvev());
    cout << "   # higgs vev(Q)MSSM DRbar\n";
    cout << "     4    ";
    printRow(displayM3Squared() /
             (sin(atan(displayTanb())) * cos(atan(displayTanb()))));
    cout << "   # mA^2(Q)MSSM DRbar\n";

    cout << "Block msoft Q=" << q
         << " # MSSM DRbar SUSY breaking parameters\n";
    int i;
    for (i=1; i<=3; i++) {
      cout << "     " << i << "    ";
      printRow(displayGaugino(i));
      cout << "   # M_" << i << "(Q)" << endl;
    }

    cout << "    21    "; printRow(displayMh1Squared());
    cout << "   # mH1^2(Q)" << endl;
    cout << "    22    "; printRow(displayMh2Squared());
    cout << "   # mH2^2(Q)" << endl;

    cout << "    31    "; printRow(ccbSqrt(displaySoftMassSquared(mLl, 1, 1)));
    cout << "   # meL(Q)" << endl;
    cout << "    32    "; printRow(ccbSqrt(displaySoftMassSquared(mLl, 2, 2)));
    cout << "   # mmuL(Q)" << endl;
    cout << "    33    "; printRow(ccbSqrt(displaySoftMassSquared(mLl, 3, 3)));
    cout << "   # mtauL(Q)" << endl;
    cout << "    34    "; printRow(ccbSqrt(displaySoftMassSquared(mEr, 1, 1)));
    cout << "   # meR(Q)" << endl;
    cout << "    35    "; printRow(ccbSqrt(displaySoftMassSquared(mEr, 2, 2)));
    cout << "   # mmuR(Q)" << endl;
    cout << "    36    "; printRow(ccbSqrt(displaySoftMassSquared(mEr, 3, 3)));
    cout << "   # mtauR(Q)" << endl;
    cout << "    41    "; printRow(ccbSqrt(displaySoftMassSquared(mQl, 1, 1)));
    cout << "   # mqL1(Q)" << endl;
    cout << "    42    "; printRow(ccbSqrt(displaySoftMassSquared(mQl, 2, 2)));
    cout << "   # mqL2(Q)" << endl;
    cout << "    43    "; printRow(ccbSqrt(displaySoftMassSquared(mQl, 3, 3)));
    cout << "   # mqL3(Q)" << endl;
    cout << "    44    "; printRow(ccbSqrt(displaySoftMassSquared(mUr, 1, 1)));
    cout << "   # muR(Q)" << endl;
    cout << "    45    "; printRow(ccbSqrt(displaySoftMassSquared(mUr, 2, 2)));
    cout << "   # mcR(Q)" << endl;
    cout << "    46    "; printRow(ccbSqrt(displaySoftMassSquared(mUr, 3, 3)));
    cout << "   # mtR(Q)" << endl;
    cout << "    47    "; printRow(ccbSqrt(displaySoftMassSquared(mDr, 1, 1)));
    cout << "   # mdR(Q)" << endl;
    cout << "    48    "; printRow(ccbSqrt(displaySoftMassSquared(mDr, 2, 2)));
    cout << "   # msR(Q)" << endl;
    cout << "    49    "; printRow(ccbSqrt(displaySoftMassSquared(mDr, 3, 3)));
    cout << "   # mbR(Q)" << endl;

    cout << "Block au Q= " << q << "  \n"
         << "  1  1    "; printRow(displaySoftA(UA, 1, 1));
    cout  << "   # Au(Q)MSSM DRbar" << endl
         << "  2  2    "; printRow(displaySoftA(UA, 2, 2));
    cout  << "   # Ac(Q)MSSM DRbar" << endl
         << "  3  3    "; printRow(displaySoftA(UA, 3, 3));
    cout  << "   # At(Q)MSSM DRbar" << endl;
    cout << "Block ad Q= " << q << "  \n"
         << "  1  1    "; printRow(displaySoftA(DA, 1, 1));
    cout  << "   # Ad(Q)MSSM DRbar" << endl
         << "  2  2    "; printRow(displaySoftA(DA, 2, 2));
    cout  << "   # As(Q)MSSM DRbar" << endl
         << "  3  3    "; printRow(displaySoftA(DA, 3, 3));
    cout << "   # Ab(Q)MSSM DRbar" << endl;
    cout << "Block ae Q= " << q << "  \n"
         << "  1  1    "; printRow(displaySoftA(EA, 1, 1));
    cout  << "   # Ae(Q)MSSM DRbar" << endl
         << "  2  2    "; printRow(displaySoftA(EA, 2, 2));
    cout  << "   # Amu(Q)MSSM DRbar" << endl
         << "  3  3    ";   printRow(displaySoftA(EA, 3, 3));
    cout << "   # Atau(Q)MSSM DRbar" << endl;
  }
}


void extendedSugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
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
  m.setMh1Squared(inputParameters.display(21));
  m.setMh2Squared(inputParameters.display(22));
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

void MssmSoftsusy::lesHouchesAccordOutput(char model[],
                                          const DoubleVector & pars,
                                          int sgnMu, double tanb,
                                          double qMax, char fileName[],
                                          int numPoints, double mb,
                                          double mtau, double mgut,
                                          bool altEwsb) {
  fstream fout(fileName, ios::out);
  fout.setf(ios::scientific, ios::floatfield);
  fout.precision(9);

  fout << "# SOFTSUSY" << VERSION << endl;
  fout << "# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331,";
  fout << " hep-ph/0104145\n";

  fout << "Block SPINFO         # Program information\n"
       << "     1   SOFTSUSY    # spectrum calculator\n";
  fout << "     2   " << VERSION << "         # version number\n";
  if (displayProblem().noConvergence)
    fout << "     3   Possible problem: Not achieved desired accuracy to factor of 10\n";
  if (displayProblem().testSeriousProblem())
    fout << "     4   Point invalid: " << displayProblem() << endl;

  fout << "Block MODSEL  # Select model\n";
  int modsel = 0;
  if (!strcmp(model, "sugra")) modsel = 1;
  if (!strcmp(model, "gmsb")) modsel = 2;
  if (!strcmp(model, "amsb")) modsel = 3;
  if (!strcmp(model, "mirage")) modsel = 4;

  fout << "     1    " << modsel << "   # " << model << "\n"; // Les Houches
                                                            // accord codes

  QedQcd d(displayDataSet());
  fout << "Block SMINPUTS   # Standard Model inputs\n";
  fout << "     1   "; printRow(fout, 1.0 / d.displayAlpha(ALPHA));
  fout << "   # alpha_em^(-1)(MZ) SM MSbar\n";
  fout << "     2   "; printRow(fout, GMU); fout << "   # G_Fermi\n";
  fout << "     3   "; printRow(fout, d.displayAlpha(ALPHAS));
  fout << "   # alpha_s(MZ)MSbar\n";
  fout << "     4   "; printRow(fout, displayMz()); fout << "   # MZ(pole)\n";
  fout << "     5   "; printRow(fout, mb); fout << "   # mb(mb)\n";
  fout << "     6   "; printRow(fout, d.displayPoleMt()); fout << "   # Mtop(pole)\n";
  fout << "     7   "; printRow(fout, mtau); fout << "   # Mtau(pole)\n";

  fout << "Block MINPAR  # SUSY breaking input parameters\n";
  fout << "     3   "; printRow(fout, tanb)            ; fout << "   # tanb" << endl;
  if (!altEwsb) {
    fout << "     4   ";
    printRow(fout, double(sgnMu));
    fout << "   # sign(mu)"<< endl;
  }
  if (!strcmp(model, "sugra")) {
    fout << "     1   "; printRow(fout, pars.display(1)); fout << "   # m0" << endl;
    fout << "     2   "; printRow(fout, pars.display(2)) ; fout << "   # m12" << endl;
    fout << "     5   "; printRow(fout, pars.display(3)) ; fout << "   # A0" << endl;
  }
  else if (!strcmp(model, "gmsb")) {
    fout << "     1   "; printRow(fout, pars.display(3)); fout << "   # lambda" << endl;
    fout << "     2   "; printRow(fout, pars.display(2)) ; fout << "   # M_mess"
         << endl;
    fout << "     5   "; printRow(fout, pars.display(1)) ; fout << "   # N5" << endl;
    fout << "     6   "; printRow(fout, pars.display(4)) ; fout << "   # cgrav" << endl;
  }
  else if (!strcmp(model, "amsb")) {
    fout << "     1   "; printRow(fout, pars.display(1)); fout << "   # m3/2" << endl;
    fout << "     2   "; printRow(fout, pars.display(2)) ; fout << "   # m0"
         << endl;
  }
  else if (!strcmp(model, "mirage")) {
    fout << "     1   "; printRow(fout, pars.display(1)); fout << "   # alphac" << endl;
    fout << "     2   "; printRow(fout, pars.display(2)); fout << "   # M0 " << endl;
    fout << "    51   "; printRow(fout, pars.display(3)); fout << "   # aq" << endl;
    fout << "    52   "; printRow(fout, pars.display(4)); fout << "   # al" << endl;
    fout << "    61   "; printRow(fout, pars.display(5)); fout << "   # ahd" << endl;
    fout << "    62   "; printRow(fout, pars.display(6)); fout << "   # ahu" << endl;
    fout << "    71   "; printRow(fout, pars.display(7)); fout << "   # cq" << endl;
    fout << "    72   "; printRow(fout, pars.display(8)); fout << "   # cl" << endl;
    fout << "    81   "; printRow(fout, pars.display(9)); fout << "   # chd" << endl;
    fout << "    82   "; printRow(fout, pars.display(10)); fout << "   # chu" << endl;
  }
  else if (!strcmp(model, "nonUniversal")) {
    fout << "Block EXTPAR         # non-universal SUSY breaking parameters\n";
    int i;
    for (i=1; i<=3; i++) {
      fout << "     " << i << "   ";
      printRow(fout, pars.display(i));
      fout << "         # M_" << i << "(MX)" << endl;
    }
    fout << "     11   "; printRow(fout, pars.display(11)) ;
    fout << "         # At(MX)" << endl;
    fout << "     12   "; printRow(fout, pars.display(12)) ;
    fout << "         # Ab(MX)" << endl;
    fout << "     13   "; printRow(fout, pars.display(13)) ;
    fout << "         # Atau(MX)" << endl;
    if (!altEwsb) {
      fout << "     21   "; printRow(fout, pars.display(21)) ;
      fout << "         # mHd^2(MX)" << endl;
      fout << "     22   "; printRow(fout, pars.display(22)) ;
      fout << "         # mHu^2(MX)" << endl;
    } else {
      fout << "     23   "; printRow(fout, pars.display(23)) ;
      fout << "         # mu(MSUSY)" << endl;
      fout << "     26   "; printRow(fout, pars.display(26)) ;
      fout << "         # mA(pole)" << endl;
    }
    fout << "     31   "; printRow(fout, pars.display(31)) ;
    fout << "         # meL(MX)" << endl;
    fout << "     32   "; printRow(fout, pars.display(32)) ;
    fout << "         # mmuL(MX)" << endl;
    fout << "     33   "; printRow(fout, pars.display(33)) ;
    fout << "         # mtauL(MX)" << endl;
    fout << "     34   "; printRow(fout, pars.display(34)) ;
    fout << "         # meR(MX)" << endl;
    fout << "     35   "; printRow(fout, pars.display(35)) ;
    fout << "         # mmuR(MX)" << endl;
    fout << "     36   "; printRow(fout, pars.display(36)) ;
    fout << "         # mtauR(MX)" << endl;
    fout << "     41   "; printRow(fout, pars.display(41)) ;
    fout << "         # mqL1(MX)" << endl;
    fout << "     42   "; printRow(fout, pars.display(42)) ;
    fout << "         # mqL2(MX)" << endl;
    fout << "     43   "; printRow(fout, pars.display(43)) ;
    fout << "         # mqL3(MX)" << endl;
    fout << "     44   "; printRow(fout, pars.display(44)) ;
    fout << "         # muR(MX)" << endl;
    fout << "     45   "; printRow(fout, pars.display(45)) ;
    fout << "         # mcR(MX)" << endl;
    fout << "     46   "; printRow(fout, pars.display(46)) ;
    fout << "         # mtR(MX)" << endl;
    fout << "     47   "; printRow(fout, pars.display(47)) ;
    fout << "         # mdR(MX)" << endl;
    fout << "     48   "; printRow(fout, pars.display(48)) ;
    fout << "         # msR(MX)" << endl;
    fout << "     49   "; printRow(fout, pars.display(49)) ;
    fout << "         # mbR(MX)" << endl;
  }
  else {
    ostringstream ii;
    ii << "Attempting to use SUSY Les Houches Accord for model "
       << model << " - cannot do at present\n";
    throw ii.str();
  }

  fout << "# Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE="
       << TOLERANCE << endl;
  fout << "# mgut=" << mgut << " GeV\n";

  sPhysical s(displayPhys());

  fout << "Block MASS   # Mass spectrum\n";
  fout << "#PDG code      mass              particle\n";
  // fout << "          6   "; printRow(fout, displayDataSet().displayPoleMt());
  // fout << "   # top\n";
  fout << "        24    "; printRow(fout, displayMw()); fout << "   # MW\n";
  fout << "        25    "; printRow(fout, s.mhiggs(1)); fout << "   # h0\n";
  // fout << "         32   "; printRow(fout, displayMz()); fout << "   # MZ\n";
  fout << "        35    "; printRow(fout, s.mhiggs(3)); fout << "   # H0\n";
  fout << "        36    "; printRow(fout, s.mhiggs(2)); fout << "   # A0\n";
  fout << "        37    "; printRow(fout, s.mhiggs(4)); fout << "   # H+\n";
  fout << "   1000001    "; printRow(fout, s.md(1, 1)); fout << "   # ~d_L\n";
  fout << "   1000002    "; printRow(fout, s.mu(1, 1)); fout << "   # ~u_L\n";
  fout << "   1000003    "; printRow(fout, s.md(1, 2)); fout << "   # ~s_L\n";
  fout << "   1000004    "; printRow(fout, s.mu(1, 2)); fout << "   # ~c_L\n";
  fout << "   1000005    "; printRow(fout, minimum(s.md(1, 3), s.md(2, 3)));
  fout << "   # ~b_1\n";   fout << "   1000006    ";
  printRow(fout, minimum(s.mu(1, 3), s.mu(2, 3)));
  fout << "   # ~t_1\n";
  fout << "   1000011    "; printRow(fout, s.me(1, 1)); fout << "   # ~e_L\n";
  fout << "   1000012    "; printRow(fout, s.msnu(1)); fout << "   # ~nue_L\n";
  fout << "   1000013    "; printRow(fout, s.me(1, 2)); fout << "   # ~mu_L\n";
  fout << "   1000014    "; printRow(fout, s.msnu(2)); fout << "   # ~numu_L\n";
  fout << "   1000015    "; printRow(fout, minimum(s.me(1, 3), s.me(2, 3)));
  fout << "   # ~stau_1\n";
  fout << "   1000016    "; printRow(fout, s.msnu(3)); fout << "   # ~nu_tau_L\n";
  fout << "   1000021    "; printRow(fout, s.mGluino); fout << "   # ~g\n";
  fout << "   1000022    "; printRow(fout, s.mneut(1));
  fout << "   # ~neutralino(1)\n";
  fout << "   1000023    "; printRow(fout, s.mneut(2));
  fout << "   # ~neutralino(2)\n";
  fout << "   1000024    "; printRow(fout, s.mch(1)); fout << "   # ~chargino(1)\n";
  fout << "   1000025    "; printRow(fout, s.mneut(3));
  fout << "   # ~neutralino(3)\n";
  fout << "   1000035    "; printRow(fout, s.mneut(4));
  fout << "   # ~neutralino(4)\n";
  fout << "   1000037    "; printRow(fout, s.mch(2)); fout << "   # ~chargino(2)\n";
  fout << "   2000001    "; printRow(fout, s.md(2, 1)); fout << "   # ~d_R\n";
  fout << "   2000002    "; printRow(fout, s.mu(2, 1)); fout << "   # ~u_R\n";
  fout << "   2000003    "; printRow(fout, s.md(2, 2)); fout << "   # ~s_R\n";
  fout << "   2000004    "; printRow(fout, s.mu(2, 2)); fout << "   # ~c_R\n";
  fout << "   2000005    "; printRow(fout, maximum(s.md(1, 3), s.md(2, 3)));
  fout << "   # ~b_2\n";
  fout << "   2000006    "; printRow(fout, maximum(s.mu(1, 3), s.mu(2, 3)));
  fout << "   # ~t_2\n";
  fout << "   2000011    "; printRow(fout, s.me(2, 1)); fout << "   # ~e_R\n";
  fout << "   2000013    "; printRow(fout, s.me(2, 2)); fout << "   # ~mu_R\n";
  fout << "   2000015    "; printRow(fout, maximum(s.me(1, 3), s.me(2, 3)));
  fout << "   # ~stau_2\n";

  const double underflow = 1.0e-120;
  if (fabs(displayGravitino()) > underflow)
    fout << "   1000039     " << displayGravitino() << "   # ~gravitino\n";
  fout << "# Higgs mixing\n";
  double ms = displayMsusy();
  fout << "Block alpha " <<
    "  # Effective Higgs mixing parameter\n";
  fout << "          "; printRow(fout, s.thetaH);        fout << "   # alpha\n";

  DoubleMatrix m(2, 2);

  fout << "Block stopmix  # stop mixing matrix\n";
  if (s.mu(1, 3) < s.mu(2, 3)) m = rot2d(s.thetat);
  else m = rot2dTwist(s.thetat);
  int i, j;
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fout << "  " << i << "  " << j << "    ";  printRow(fout, m(i, j));
      fout << "   # O_{" << i << j << "}" << endl;
    }

  fout << "Block sbotmix  # sbottom mixing matrix\n";
  if (s.md(1, 3) < s.md(2, 3)) m = rot2d(s.thetab);
  else m = rot2dTwist(s.thetab);
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fout << "  " << i << "  " << j << "    "; printRow(fout, m(i, j));
      fout << "   # O_{" << i << j << "}" << endl;
    }

  fout << "Block staumix  # stau mixing matrix\n";
  if (s.me(1, 3) < s.me(2, 3)) m = rot2d(s.thetatau);
  else m = rot2dTwist(s.thetatau);
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fout << "  " << i << "  " << j << "    "; printRow(fout, m(i, j));
      fout << "   # O_{" << i << j << "}" << endl;
   }

  fout << "Block nmix  # neutralino mixing matrix\n";
  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      fout << "  " << i << "  " << j << "    ";
      printRow(fout, s.mixNeut(j, i));
      fout << "   # N_{" << i << "," << j << "}\n";
    }

  DoubleMatrix u(rot2d(s.thetaL)), v(rot2d(s.thetaR));

  i = 1;
  fout << "Block Umix  # chargino U mixing matrix \n";
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fout << "  " << i << "  " << j << "    "; printRow(fout, u(i, j));
      fout << "   # U_{" << i << "," << j << "}\n";
    }

  fout << "Block Vmix  # chargino V mixing matrix \n";
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fout << "  " << i << "  " << j << "    "; printRow(fout, v(i, j));
      fout << "   # V_{" << i << "," << j << "}\n";
    }

  int n = 0;

  while (n < numPoints) {
    n++;
    // Starting non-essential information. The following decides what scale to
    // output the running parameters at. It depends upon what qMax is and how
    // many points the user has requested.
    // For qMax = 0 and 1 point (defaults), Q=MSUSY is printed out.
    // For qMax = 0 and n points, points are spaced logarithmically between MZ
    // and MSUSY.
    // For qMax != 0 and 1 point, Q=qMax is printed.
    // For qMax != 0 and n points, points are log spaced between MZ and qMax.
    double q = ms;
    if (numPoints == 1 && qMax < EPSTOL) q = ms;
    else if (numPoints == 1 && qMax > EPSTOL) q = qMax;
    else if (numPoints > 1 && qMax < EPSTOL) qMax = ms;

    if (numPoints > 1) {
      if (n > 1) {
        double logq = (log(qMax) - log(MZ)) * double(n-1) /
          double(numPoints-1) + log(MZ);
        q = exp(logq);
      }
      else q = MZ;
    }

    runto(q);
    double gp = displayGaugeCoupling(1) * sqrt(0.6);
    fout << "Block gauge Q= " << q << "  \n";
    fout << "     1     " << gp << "   # g'(Q)MSSM DRbar"
         << endl;
    fout << "     2     " << displayGaugeCoupling(2)
         << "   # g(Q)MSSM DRbar" << endl;
    fout << "     3     " << displayGaugeCoupling(3) << "   # g3(Q)MSSM DRbar"
         << endl;

    fout << "Block yu Q= " << q << "  \n"
         << "  3  3     " << displayYukawaElement(YU, 3, 3)
         << "   # Yt(Q)MSSM DRbar" << endl;
    fout << "Block yd Q= " << q << "  \n"
         << "  3  3     " << displayYukawaElement(YD, 3, 3)
         << "   # Yb(Q)MSSM DRbar" << endl;
    fout << "Block ye Q= " << q << "  \n"
         << "  3  3     " << displayYukawaElement(YE, 3, 3)
         << "   # Ytau(Q)MSSM DRbar" << endl;

    fout << "Block hmix Q= " << q <<
      "  # Higgs mixing parameters\n";
    fout << "     1    "; printRow(fout, displaySusyMu());
    fout << "   # mu(Q)MSSM DRbar\n";
    fout << "     2    "; printRow(fout, displayTanb());
    fout << "   # tan beta(Q)MSSM DRbar\n";
    fout << "     3    "; printRow(fout, displayHvev());
    fout << "   # higgs vev(Q)MSSM DRbar\n";
    fout << "     4    "; printRow(fout, displayM3Squared() /
             (sin(atan(displayTanb())) * cos(atan(displayTanb()))));
    fout << "   # mA^2(Q)MSSM DRbar\n";

    fout << "Block msoft Q=" << q
         << " # MSSM DRbar SUSY breaking parameters\n";
    int i;
    for (i=1; i<=3; i++) {
      fout << "     " << i << "    ";
      printRow(fout, displayGaugino(i));
      fout << "   # M_" << i << "(Q)" << endl;
    }

    fout << "    21    "; printRow(fout, displayMh1Squared());
    fout << "   # mH1^2(Q)" << endl;
    fout << "    22    "; printRow(fout, displayMh2Squared());
    fout << "   # mH2^2(Q)" << endl;

    fout << "    31    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mLl, 1, 1)));
    fout << "   # meL(Q)" << endl;
    fout << "    32    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mLl, 2, 2)));
    fout << "   # mmuL(Q)" << endl;
    fout << "    33    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mLl, 3, 3)));
    fout << "   # mtauL(Q)" << endl;
    fout << "    34    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mEr, 1, 1)));
    fout << "   # meR(Q)" << endl;
    fout << "    35    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mEr, 2, 2)));
    fout << "   # mmuR(Q)" << endl;
    fout << "    36    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mEr, 3, 3)));
    fout << "   # mtauR(Q)" << endl;
    fout << "    41    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mQl, 1, 1)));
    fout << "   # mqL1(Q)" << endl;
    fout << "    42    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mQl, 2, 2)));
    fout << "   # mqL2(Q)" << endl;
    fout << "    43    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mQl, 3, 3)));
    fout << "   # mqL3(Q)" << endl;
    fout << "    44    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mUr, 1, 1)));
    fout << "   # muR(Q)" << endl;
    fout << "    45    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mUr, 2, 2)));
    fout << "   # mcR(Q)" << endl;
    fout << "    46    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mUr, 3, 3)));
    fout << "   # mtR(Q)" << endl;
    fout << "    47    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mDr, 1, 1)));
    fout << "   # mdR(Q)" << endl;
    fout << "    48    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mDr, 2, 2)));
    fout << "   # msR(Q)" << endl;
    fout << "    49    "; printRow(fout, ccbSqrt(displaySoftMassSquared(mDr, 3, 3)));
    fout << "   # mbR(Q)" << endl;

    fout << "Block au Q= " << q << "  \n"
         << "  1  1    "; printRow(fout, displaySoftA(UA, 1, 1));
    fout  << "   # Au(Q)MSSM DRbar" << endl
         << "  2  2    "; printRow(fout, displaySoftA(UA, 2, 2));
    fout  << "   # Ac(Q)MSSM DRbar" << endl
         << "  3  3    "; printRow(fout, displaySoftA(UA, 3, 3));
    fout  << "   # At(Q)MSSM DRbar" << endl;
    fout << "Block ad Q= " << q << "  \n"
         << "  1  1    "; printRow(fout, displaySoftA(DA, 1, 1));
    fout  << "   # Ad(Q)MSSM DRbar" << endl
         << "  2  2    "; printRow(fout, displaySoftA(DA, 2, 2));
    fout  << "   # As(Q)MSSM DRbar" << endl
         << "  3  3    "; printRow(fout, displaySoftA(DA, 3, 3));
    fout << "   # Ab(Q)MSSM DRbar" << endl;
    fout << "Block ae Q= " << q << "  \n"
         << "  1  1    "; printRow(fout, displaySoftA(EA, 1, 1));
    fout  << "   # Ae(Q)MSSM DRbar" << endl
         << "  2  2    "; printRow(fout, displaySoftA(EA, 2, 2));
    fout  << "   # Amu(Q)MSSM DRbar" << endl
         << "  3  3    ";   printRow(fout, displaySoftA(EA, 3, 3));
    fout << "   # Atau(Q)MSSM DRbar" << endl;
  }
}


// User supplied routine. Inputs m at the unification scale, and uses
// inputParameters vector to output m with high energy soft boundary
// conditions.
void sugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  // Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  // trilinear scalar coupling to be a0
  m.standardSugra(m0, m12, a0);

  return;
}

// Other types of boundary condition
void amsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters)
{
  double m32 = inputParameters.display(1);
  double m0 = inputParameters.display(2);

  m.standardSugra(m0, 0., 0.);
  m.addAmsb(m32);

  return;
}

void gmsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters)
{
  int n5     = int(inputParameters.display(1));
  double mMess  = inputParameters.display(2);
  double lambda = inputParameters.display(3);
  double cgrav = inputParameters.display(4);

  m.minimalGmsb(n5, lambda, mMess, cgrav);

  return;
}

void userDefinedBcs(MssmSoftsusy & m, const DoubleVector & inputParameters) {
  double alphac = inputParameters.display(1);
  double M0 = inputParameters.display(2);
  double aq = inputParameters.display(3);
  double al = inputParameters.display(4);
  double ahd = inputParameters.display(5);
  double ahu = inputParameters.display(6);
  double cq = inputParameters.display(7);
  double cl = inputParameters.display(8);
  double chd = inputParameters.display(9);
  double chu = inputParameters.display(10);

// Sets the Mirage Mediation model parameters which will be related with the
// soft scalar masses ,gaugino masses, trilinear scalar coupings.
// Here the method-mirage is defined in SoftParMssm/softpars.h & softpars.cpp
  m.mirage(alphac,M0,aq,al,ahu,ahd,cq,cl,chu,chd);

  return;
}

// all matrices in inputs should be 3x3: they're output as the quark
// FCNC matrices.
void MssmSoftsusy::sCkmRotation
(DoubleMatrix & deltaULL, DoubleMatrix & deltaURR, DoubleMatrix & deltaULR,
 DoubleMatrix & deltaDLL, DoubleMatrix & deltaDRR, DoubleMatrix & deltaDLR)
  const {

  double v = displayHvev();
  double beta = atan(displayTanb());
  double v1 = v / root2 * cos(beta);
  double v2 = v / root2 * sin(beta);

  DoubleMatrix vdl(3, 3), vdr(3, 3), vul(3, 3), vur(3, 3);
  diagQuarkBasis(vdl, vdr, vul, vur);

  DoubleMatrix   muLL(vul * displaySoftMassSquared(mQl) * vul.transpose());
  DoubleMatrix   muRR(vur * displaySoftMassSquared(mUr) * vur.transpose());
  DoubleMatrix   muLR(v2 * vul * displayTrilinear(UA) * vur.transpose());
  DoubleMatrix   mdLL(vdl * displaySoftMassSquared(mQl) * vdl.transpose());
  DoubleMatrix   mdRR(vdr * displaySoftMassSquared(mDr) * vdr.transpose());
  DoubleMatrix   mdLR(v1 * vdl * displayTrilinear(DA) * vdr.transpose());

  int i, j;
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      deltaULL(i, j) = muLL(i, j) / (sqrt(fabs(muLL(i, i) * muLL(j, j))));
      deltaURR(i, j) = muRR(i, j) / (sqrt(fabs(muRR(i, i) * muRR(j, j))));
      deltaULR(i, j) = muLR(i, j) / (sqrt(fabs(muLL(i, i) * muRR(j, j))));
      deltaDLL(i, j) = mdLL(i, j) / (sqrt(fabs(mdLL(i, i) * mdLL(j, j))));
      deltaDRR(i, j) = mdRR(i, j) / (sqrt(fabs(mdRR(i, i) * mdRR(j, j))));
      deltaDLR(i, j) = mdLR(i, j) / (sqrt(fabs(mdLL(i, i) * mdRR(j, j))));
    }
}

void MssmSoftsusy::outputFcncs() const {
  DoubleMatrix deltaULL(3, 3), deltaURR(3, 3), deltaULR(3, 3);
  DoubleMatrix deltaDLL(3, 3), deltaDRR(3, 3), deltaDLR(3, 3);
  sCkmRotation(deltaULL, deltaURR, deltaULR, deltaDLL, deltaDRR, deltaDLR);

  cout << "FCNC delta parameters\n";
  cout << "deltaULL" << deltaULL << "deltaURR" << deltaURR << "deltaULR"
       << deltaULR << "deltaDLL" << deltaDLL << "deltaDRR" << deltaDRR
       << "deltaDLR" << deltaDLR;
}

// Returns nlsp mass in mass and function return labels which particle is nlsp:
// 0 is neutralino posi = #, posj = 0
int MssmSoftsusy::nlsp(double & mass, int & posi, int & posj) const {
  int temp = 0, ntemp, lsppos1, lsppos2, nlsppos1, nlsppos2, pos1, pos2;
  sPhysical s(displayPhys());

  double minmass = s.mneut.apply(fabs).min(pos1); lsppos1 = pos1; lsppos2 = 0;
  double nminmass;

  // up squarks 1
  double lightest = s.mu.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) {
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=0;
    minmass = lightest; lsppos1 = pos1; lsppos2 = pos2; temp = 1; }
  else { nminmass = lightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 1;}

  // down squarks 2
  lightest = s.md.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) {
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = pos2; temp = 2;}
  else if (lightest < nminmass){
     nminmass = lightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 2;}

  // sleptons 3
  lightest = s.me.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) {
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = pos2; temp = 3;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 3;}

  // charginos 4
  lightest = s.mch.apply(fabs).min(pos1);
  if (lightest < minmass) {
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = 0; temp = 4;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 4;}

  // sneutrinos 5
  lightest = s.msnu.apply(fabs).min(pos1);
  if (lightest < minmass) {
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = 0; temp = 5;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 5;}

  // gluino 6
  lightest = s.mGluino;
  if (lightest < minmass) {
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = 0; temp = 6;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 6;}


  // check that nlsp is not the next-to-lightest in the group of the lsp
  double nlightest;

  switch (temp){
    // neutralinos
  case 0: {
    nlightest = s.mneut.apply(fabs).nmin(pos1);
    if (nlightest < nminmass) {
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 0;}
  }
    // up squarks 1
  case 1: {
    nlightest = s.mu.apply(fabs).nmin(pos1, pos2);
    if (nlightest < nminmass) {
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 1;}
  }
    // down squarks 2
  case 2: {
    nlightest = s.md.apply(fabs).nmin(pos1, pos2);
    if (nlightest < nminmass) {
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 2;}
  }
    // sleptons 3
  case 3: {
    nlightest = s.me.apply(fabs).nmin(pos1, pos2);
    if (nlightest < nminmass){
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 3;}
  }
    // charginos 4
  case 4: {
    nlightest = s.mch.apply(fabs).nmin(pos1);
    if (nlightest < nminmass){
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 4;}
  }
    // sneutrinos 5
  case 5: {
    nlightest = s.msnu.apply(fabs).nmin(pos1);
    if (nlightest < nminmass){
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 5;}
  }
    // gluino 6
    // there is only one gluino mass
  }

  mass = nminmass;
  posi = nlsppos1; posj = nlsppos2;
  return ntemp;
}
