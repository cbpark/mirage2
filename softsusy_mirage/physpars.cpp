
/** \file physpars.cpp
   - Project:     SOFTSUSY 
   - File:        physpars.cpp
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: physpars.cpp,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.11  2004/03/19 21:03:04  allanach
   Added 1-loop tadpole part

   Revision 1.10  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.9  2003/07/22 09:16:04  allanach
   Added running MW, MZ to definition of drbar parameters

   Revision 1.8  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.7  2003/02/24 14:26:39  allanach
   Implementing DRbar parameters in loop corrections. Half way though:
   drbarpars now changes to MPZ notation, but need to get rid of pole masses
   in loop corrections (and re-calculate tadpoles).

   Revision 1.6  2003/02/21 17:59:36  allanach
   Added drbar parameter class and calculation, starting to move to DRbar
   parameters in the 1-loop corrections

   Revision 1.5  2002/10/22 13:12:10  allanach
   Introduced new problem flag for infra-red quasi fixed points

   Revision 1.4  2002/08/30 12:31:55  allanach
   Test RPV version that works!

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#include "physpars.h"

const drBarPars & drBarPars::operator=(const drBarPars &s) {
  if (this == &s) return *this;
  mz = s.mz; mw = s.mw;
  mt = s.mt; mb = s.mb; mtau = s.mtau;
  nBpmz = s.nBpmz; uBpmz = s.uBpmz; vBpmz = s.vBpmz; 
  mnBpmz = s.mnBpmz; mchBpmz = s.mchBpmz;
  setsPhysical(s.displaysPhysical());
  
  return *this;
}

const sPhysical & sPhysical::operator=(const sPhysical &s) {
  if (this == &s) return *this;
  mhiggs = s.mhiggs; msnu = s.msnu; 
  mch = s.mch; mneut = s.mneut; mixNeut = s.mixNeut;
  thetaL = s.thetaL; thetaR = s.thetaR; mGluino = s.mGluino;
  thetat = s.thetat; thetab = s.thetab; thetatau = s.thetatau;
  mu = s.mu; md = s.md; me = s.me; thetaH = s.thetaH;
  t1OV1Ms = s.t1OV1Ms; t2OV2Ms = s.t2OV2Ms;
  t1OV1Ms1loop = s.t1OV1Ms1loop; t2OV2Ms1loop = s.t2OV2Ms1loop;
  aChi0ChicW = s.aChi0ChicW; bChi0ChicW = s.bChi0ChicW;
  return *this;
}

// a should be in C convention ie start from index zero
void sPhysical::display(double *a) const {
  a[0] = mhiggs.display(1); a[1] = mhiggs.display(2); 
  a[2] = mhiggs.display(3); a[3] = mhiggs.display(4);

  a[4] = msnu.display(1); a[5] = msnu.display(2); a[6] = msnu.display(3);

  a[7] = mch.display(1); a[8] = mch.display(2);

  a[9] = mneut.display(1); a[10] = mneut.display(2); 
  a[11] = mneut.display(3); a[12] = mneut.display(4);

  a[13] = mGluino;

  int i, j, k = 13; 
  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      k++;
      a[k] = mixNeut.display(i, j);  
    }

  a[30] = thetaL; a[31] = thetaR; 
  a[32] = thetat; a[33] = thetab; a[34] = thetatau;

  k = 34;
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++) {
      k++;
      a[k] = mu.display(i, j);
      a[k+6] = md.display(i, j);
      a[k+12] = me.display(i, j);
    }  

  a[53] = thetaH;
}

#define HR "---------------------------------------------------------------\n"

ostream & operator <<(ostream & left, const drBarPars &s) {
  left << s.displaysPhysical();
  left << "BPMZ conventions, N" << s.nBpmz << "U" << s.uBpmz << "V" 
       << s.vBpmz;
  left << "mt: "  << s.mt << " mb: " << s.mb << " mtau: " << s.mtau << endl;
  left << "mz: "  << s.mz << " mw: " << s.mw << endl;

  return left;
}

ostream & operator <<(ostream & left, const sPhysical &s) {
  left << "mh^0: " << s.mhiggs.display(1) << " mA^0: " << s.mhiggs.display(2)
       << " mH^0: " << 
    s.mhiggs.display(3) << " mH^+-: " << s.mhiggs.display(4) << "\n";
  left << "alpha: " << s.thetaH << "\n";
  left << "sneutrinos" << s.msnu; 
  left << "mU~" << s.mu << "mD~" << s.md << "mE~" << s.me;
  left << "thetat: " << s.thetat << " thetab: " << s.thetab << 
    " thetatau: " << s.thetatau << "\n";
  left << "mGluino:  " << s.mGluino << "\n";
  left << "charginos" << s.mch;
  left << "thetaL: " << s.thetaL << " thetaR: " << s.thetaR << "\n";
  left << "neutralinos" << s.mneut;
  left << "neutralino mixing matrix " << s.mixNeut;
  return left;
}

istream & operator >>(istream & left, sPhysical &s) {
  char c[70];
  left >> c >> c >> c >> c;
  DoubleVector mh(4);
  left >> c >> s.mhiggs(1) >> c >> s.mhiggs(2)
       >> c >> s.mhiggs(3) >> c >> s.mhiggs(4);
  left >> c >> s.thetaH;
  left >> s.msnu; 
  left >> c >> s.mu >> c >> s.md >> c >> s.me;
  left >> c >> s.thetat >> c >> s.thetab >> 
    c >> s.thetatau;
  left >> c >> s.mGluino;
  left >> s.mch;
  left >> c >> s.thetaL >> c >> s.thetaR;
  left >> s.mneut;
  left >> c >> c >> c >> c >> s.mixNeut;
  return left;
}

#undef HR

ostream & operator <<(ostream &st, const sProblem & p) {
  if (!p.test()) return st;
  st << "[ ";
  if (p.badConvergence) st << "No acceptable solution found";
  if (p.irqfp) st << "Quasi-fixed point breached ";
  if (p.noMuConvergence) st << "No mu convergence ";
  if (p.noRhoConvergence) st << "No rho convergence ";
  if (p.nonperturbative) st << "Non-perturbative ";
  if (p.noConvergence) st << "No convergence ";
  if (p.tachyon) st << "Tachyon ";
  if (p.muSqWrongSign) st << "MuSqWrongsign ";
  if (p.b) st << "B-problem ";
  if (p.higgsUfb) st << "Higgs potential ufb ";
  st << "]";
  return st;
}

const sProblem & sProblem::operator=(const sProblem &s) {
  if (this == &s) return *this;
  irqfp = s.irqfp;
  badConvergence = s.badConvergence;
  noMuConvergence = s.noMuConvergence;
  noRhoConvergence = s.noRhoConvergence;
  nonperturbative = s.nonperturbative;
  noConvergence = s.noConvergence;
  tachyon = s.tachyon;
  muSqWrongSign = s.muSqWrongSign;
  higgsUfb = s.higgsUfb;
  b = s.b;
  return *this;
}

// Returns mixing matrix o and neutralino masses mn in the MPZ convention
// (hep-ph/9606211), n is 4 by 4 and mneut is 1->4.
void drBarPars::mpzNeutralinos() { 
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang

  DoubleVector temp(mneut);
  
  ComplexMatrix K(4, 4);
  int i; for (i=1; i<=4; i++) 
    if (mneut.display(i) < 0.0) K(i, i) = Complex(0.0, 1.0);
    else
      K(i, i) = Complex(1.0, 0.0);
  
  mnBpmz = temp.apply(fabs);
  nBpmz = K.hermitianConjugate() * mixNeut.transpose();
}

// Returns mixing matrices u,v and neutralino masses mneut in the MPZ
// convention (hep-ph/9606211),  u+v are (2,2) and mch is 1->2.
void drBarPars::mpzCharginos() {
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  ComplexMatrix u(2, 2), v(2, 2);
  positivise(thetaL, thetaR, mch, u, v);
  uBpmz = u; vBpmz = v;
  mchBpmz = mch.apply(fabs); 
}
