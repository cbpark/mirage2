
/** \file rgeroutines.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: odds and sods for my own machievellian needs

*/

#include <rgeroutines.h>

const static double factor = 1.0;
const static int highLuminosity = 1;

double spc(const MssmSoftsusy & r) {
  double chi = sqr(r.displayPhys().mneut(1)), 
    xi = sqr(r.displayPhys().mneut(2)), 
    l = sqr(r.displayPhys().me(2, 1)),
    q  = sqr((r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1)) * 0.5);
  return 6.66e66;// DEBUG
  return sqrt(fabs((q - xi) * (l - chi) / (2.0 * l - chi)));
}

double fn1(double a, double b) {
  if (a+b > 500.) return -6.66;
  else return a;
} 

double fn2(double a) { 
  if (a > 250.) return -6.66;
  else return a;
}

void spsPoints(const QedQcd & oneset) {
  double m0, m12, tanb, lambda, a0, mMess, mgut = 1.9e16;
  int sgnMu = 1;

  MssmSoftsusy r[10];
  DoubleVector pars(3);

  bool gaugeUnification = true;

  // SPS 1a
  m0 = 100.; m12 = 250.; a0 = -100.; tanb = 10.; 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  r[0].lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, gaugeUnification);
  // SPS 1b
  m0 = 200.; m12 = 400.; a0 = 0.; tanb = 30.; 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  r[1].lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
  // SPS 2
  m0 = 1450.; m12 = 300.; a0 = 0.; tanb = 10.; 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  r[2].lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
  // SPS 3
  m0 = 90.; m12 = 400.; a0 = 0.; tanb = 10.; 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  r[3].lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
  // SPS 4
  m0 = 400.; m12 = 300.; a0 = 0.; tanb = 50.; 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  r[4].lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
  // SPS 5
  m0 = 150.; m12 = 300.; a0 = -1000.; tanb = 5.; 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  r[5].lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
  // SPS 6
  m0 = 150.; m12 = 300.; a0 = 0.; tanb = 10.; 
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
  r[6].lowOrg(nonUniGauginos, mgut, pars, sgnMu, tanb, oneset, true);
  // SPS 7
  lambda = 4.0e4; mMess = 8.0e4; double n5 = 3; tanb = 15.; mgut = mMess;
  pars(1) = n5; pars(2) = mMess; pars(3) = lambda;
  r[7].lowOrg(gmsbBcs, mgut, pars, sgnMu, tanb, oneset, false);
  // SPS 8
  lambda = 1.0e5; mMess = 2.0e5; n5 = 1; tanb = 15.; mgut = mMess;
  pars(1) = n5; pars(2) = mMess; pars(3) = lambda;
  r[8].lowOrg(gmsbBcs, mgut, pars, sgnMu, tanb, oneset, false);
  // SPS 9 
  m0 = 450.; m12 = 6.0e4; tanb = 10.; 
  pars(1) = m12; pars(2) = m0;
  r[9].lowOrg(amsbBcs, mgut, pars, sgnMu, tanb, oneset, true);

  fstream fout("sps_masses.dat", ios::out);  
  fout.setf(ios::fixed, ios::floatfield);
  fout.precision(1);
  fout.width(7);

  fout << "SOFTSUSY\n\n";
  fout << "{";
  int i = 0;
  for (i=0; i<=9; i++) 
    fout << "{" << r[i].displayPhys().mGluino 
	 << ", " << minimum(r[i].displayPhys().mu(1, 1), 
			    r[i].displayPhys().mu(2, 1)) 
	 << ", " << maximum(r[i].displayPhys().mu(1, 1),
			    r[i].displayPhys().mu(2, 1)) 
	 << ", " << minimum(r[i].displayPhys().mu(1, 3), 
			    r[i].displayPhys().mu(2, 3)) 
	 << ", " << maximum(r[i].displayPhys().mu(1, 3),
			    r[i].displayPhys().mu(2, 3)) 
	 << ", " << minimum(r[i].displayPhys().md(1, 3), 
			    r[i].displayPhys().md(2, 3)) 
	 << ", " << maximum(r[i].displayPhys().md(1, 3),
			    r[i].displayPhys().md(2, 3)) << "}" << endl; 
  fout << "}\n";
  fout.close();

  fstream fout2("lhc_masses.dat", ios::out);  
  fout2.setf(ios::fixed, ios::floatfield);
  fout2.precision(1);
  fout2.width(7);
    
    fout2 << "{\n ";
    for (i=0; i<=9; i++) {
      double mn1 = r[i].displayPhys().mneut(1), 
	me1 = r[i].displayPhys().me(2, 1),
	me3 = r[i].displayPhys().me(2, 3);
      
      fout2 << "{" << fn1(r[i].displayPhys().mneut(1), mn1) << ", " 
	   << fn1(r[i].displayPhys().mneut(2), mn1) << ", "      
	   << fn2(r[i].displayPhys().mch(1)) << ", "
	   << fn2(r[i].displayPhys().msnu(1)) << ", "
	   << fn2(r[i].displayPhys().msnu(3)) << ", "
	   << fn2(r[i].displayPhys().me(2, 1)) << ", "
	   << fn1(r[i].displayPhys().me(1, 2), me1) << ", "
	   << fn1(r[i].displayPhys().me(2, 3), me3) << ", "
	   << fn1(r[i].displayPhys().me(1, 3), me3) 
	   << "} ";
      if (i != 9) fout2 << ",\n";
    }
    fout2 << "};"; fout2.close();
    
}

// Calculates Gaussian Errors from input spectra
void gaussErrors() {
  const int numProgs = 4;
  double output[numProgs+1][10][10] = 
    { 
      { // ISASUSY
	{96.6, 176.8, 176.4, 185.8, 184.9, 142.9, 201.9, 133.3, 206.0},
	{161.2, 299.2, -6., -6., -6., -6., -6., 195.2, -6.},
	{120.2, 221.7, 221.5, -6., -6., -6., -6., -6., -6.},
	{160.5, 297.0, -6., -6., -6., 178.3, -6., 170.6, -6.},
	{119.3, 219.6, 219.6, -6., -6., -6., -6., -6., -6.},
	{119.8, 225.9, 225.9, 244.3, 242.3, 191.5, 256.3, 180.9, 257.7},
	{190.0, 218.1, 215.7, -6., -6., 236.8, -6.6, 227.8, -6.},
	{162.4, 268.0, -6., 248.7, 248.3, 127.3, 261.1, 119.9, 263.3},
	{137.4, 254.6, -6., -6., -6., 175.7, -6., 168.9, -6.},
	{174.8, -6., 175.0, -6., -6., -6., -6., -6., -6.}},
      // softsusy 1.7.1
      {{96.4, 178.2, 177.6, 188.7, 187.9, 144.9, 204.3, 136.3, 207.8} ,
       {161.3, 303.5, -6.7, -6.7, -6.7, -6.7, -6.7, 204.6, -6.7} ,
       {118.6, 232.0, 231.8, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7} ,
       {160.8, 300.4, -6.7, -6.7, -6.7, 182.0, 290.4, 175.3, 292.3} ,
       {118.6, 223.6, 223.4, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7} ,
       {118.7, 229.4, 229.2, 248.2, 246.1, 193.5, 259.4, 183.7, 260.7} ,
       {190.2, 220.6, 217.9, -6.7, -6.7, 240.9, -6.7, 232.8, -6.7} ,
       {163.6, 263.5, -6.7, 247.2, 246.9, 126.4, 259.3, 120.5, 261.2} ,
       {138.4, 261.2, -6.7, -6.7, -6.7, 175.4, -6.7, 169.9, -6.7} ,
       {196.7, -6.7, 196.7, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7}},
      { // spheno1.0
	{97.63, 182.87, 181.28, 190.46, 189.56, 143.87, 206.63, 134.58, 210.33 },
       { 163.68, 310.17,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66, 195.06,  -6.66 },
       { 123.84, 233.62, 232.06,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66 },
       { 162.67, 307.24,  -6.66,  -6.66,  -6.66, 180.17, 290.32, 172.44, 292.40 },
       { 121.42, 228.86, 227.38,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66 },
       { 121.11, 230.03, 229.83, 247.97, 245.88, 191.97, 259.44, 181.66, 260.88 },
       { 191.19, 225.33, 222.30,  -6.66,  -6.66, 237.58,  -6.66, 229.05,  -6.66 },
       { 163.43, 271.07,  -6.66,  251.6,  251.3, 131.03, 265.26, 123.78, 267.34 },
       { 139.19, 266.09,  -6.66,  -6.66,  -6.66, 180.31,  -6.66, 173.47,  -6.66 },
       { 167.98,  -6.66, 168.42,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66 }},
      {// SUSPECT
	{97.4235712,  179.766409,  179.113673,  188.464702,  187.501867,  144.936659, 
	 204.413939,  135.740046,  208.05465},
	{162.565654,  305.495483,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66,  199.664616, 
	 -6.66},
	{123.104075,  232.972834,  232.610638,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66, 
	 -6.66},
	{161.738427,  302.379932,  -6.66,  -6.66,  -6.66,  182.084159,  290.714381, 
	 174.800729,  292.607019},
	{121.078702,  225.898117,  225.677587,  -6.66,  -6.66,  -6.66,  -6.66,  -6.66, 
	 -6.66},
	{120.689454,  230.554202,  230.380472,  247.832625,  245.495339,  193.505402, 
	 259.46528,  182.760299,  260.600115},
	{190.233548,  222.183697,  219.613826,  -6.66,  -6.66,  240.920867,  -6.66, 
	 232.198286,  -6.66},
	{163.593105,  262.153395,  -6.66,  246.918512,  246.558705,  127.82699, 
	 259.431956,  121.584763,  261.41081},
	{139.982733,  263.537761,  -6.66,  -6.66,  -6.66,  177.632754,  -6.66, 
	 171.83456,  -6.66},
	{167.3, -6.6,167.3, -6.6,-6.6, -6.6,-6.6, -6.6,-6.6}}
};

  int i, j, k; 
    for (j=0; j<10; j++)
      for (k=0; k<10; k++) {
	double sumSquares = 0., mean = 0.;
	double effNum = numProgs;
	for (i=0; i<numProgs; i++) {
	  if (output[i][j][k] < 0) effNum = effNum - 1; 
	  else mean = mean + output[i][j][k]; 
	}
	mean = mean / maximum(effNum, 1.0);

	for (i=0; i<numProgs; i++) {
	  if (output[i][j][k] > 0) 
	    sumSquares = sumSquares + sqr(output[i][j][k] - mean); 
	}
	if (effNum == 0) output[4][j][k] = -6.6;
	else output[4][j][k] = sqrt(sumSquares) / 
	       maximum(double(effNum - 1), 1.0);
      }

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(1);
    cout.width(7);

    char * progName[5] =
      {"ISASUSY","SOFTSUSY","SPHENO","SUSPECT","{\\bf error}"};

    for (j=0; j<10; j++)
      for (i=0; i<numProgs+1; i++) {
	if (i == 0) 
	  if (j == 0) cout << "1a";
	  else if (j == 1) cout << "1b";
	  else cout << j;
	
	cout << " & " << progName[i]; 
	for (k=0; k<9; k++) { 
	  cout << " & ";
	  if (i == 4) cout << "{\\bf "; 
	  if (output[i][j][k] < 0.0) cout << "--";
	  else cout << output[i][j][k];
	  if (i == 4) cout << "}"; 
	}
	cout << "\\\\ ";
	if (i == 4) cout << " \\hline";
	cout << endl;
      }
}

void scanPars(double startM0, double endM0, double startM12, double endM12, 
	      double tanb, int numPoints, int sgnMu, int numRoutine) {
  double a0 = 0.;
  int i, j; for (i=0; i<=numPoints; i++) 
    for (j=0; j<=numPoints; j++) {
      double m0 = (endM0 - startM0) * double(i) / numPoints + startM0;
      double m12 = (endM12 - startM12) * double(j) / numPoints + startM12;
      cout << m0 << " " << m12 << " " 
	   << findChiSqSugra(m0, m12, a0, tanb, sgnMu, numRoutine);
    }
  cout << endl << endl; 
}

void scanPars2(double starttb, double endtb, double startM0, double endM0, 
	      double m12, int numPoints, int sgnMu, int numRoutine) {
  double a0 = 0.;
  int i, j; for (i=0; i<=numPoints; i++) 
    for (j=0; j<=numPoints; j++) {
      double tanb = (endtb - starttb) * double(i) / numPoints + starttb;
      double m0 = (endM0 - startM0) * double(j) / numPoints + startM0;
      cout << tanb << " " << m0 << " " 
	   << findChiSqSugra(m0, m12, a0, tanb, sgnMu, numRoutine);
    }
  cout << endl << endl; 
}

double findChiSqSugra(double m0, double m12, double a0, double tanb, 
		      int sgnMu, int numRoutine) {
  // set parameters in vector form
  DoubleVector pars(3);
  pars(1) = m0; pars(2) = m12; pars(3) = a0;

  QedQcd oneset;
  readIn(oneset, "massIn");
  oneset.toMz();
  
  MssmSoftsusy r; 

  double mgut = 1.9e16;
  r.lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
  
  double chiSqPoint;

  switch(numRoutine) {
  case 1: chiSqPoint = chiSqP1(r); break;
  case 2: chiSqPoint = chiSqP2(r); break;
  case 3: chiSqPoint = chiSqP3(r); break;
  case 4: chiSqPoint = chiSqP4(r); break;
  case 5: chiSqPoint = chiSqP5(r); break;
  case 6: chiSqPoint = chiSqP6(r); break;
  case 7: chiSqPoint = chiSqSPS1a(r); break;
  default: 
    throw "\nIllegal number of routine called from findChiSqSugra\n";
    break;
  }
  
  // Trap nan chi squared values
  if (chiSqPoint * 0.0 != 0.0) {
    cout << "Found Nan in chi squared. Bailing out \n"; exit(-1);
  }
    
  return chiSqPoint;
}

// decides which things to switch on
static int LHCexp = 1, LHCth = 1, LCexp = 0, LCth = 0;
double chiSqSPS1a(const MssmSoftsusy & r) {
  // LC errors etc -- checked
  double mmuRExp = 142.9, mmuRErr = 0.2, mmuRTh = 0.5;
  double mmuLExp = 201.9, mmuLErr = 0.5, mmuLTh = 3.2;
  double meLExp = 201.9, meLErr = 0.2, meLTh = 3.2;
  double meRExp = 142.9, meRErr = 0.05, meRTh = 0.5;
  double mstau1Exp = 133.3, mstau1Err = 0.3, mstau1Th = 0.6;
  double mstau2Exp = 206.0, mstau2Err = 1.1, mstau2Th = 3.0;
  double msnueExp = 185.8, msnueErr = 0.7, msnueTh = 3.1;
  double mstop1Exp = 379.4, mstop1Err = 2.0, mstop1Th = 11.9;
  double mchi1pmExp = 176.4, mchi1pmErr = 0.55, mchi1pmTh = 2.6;
  double mchi2pmExp = 379.8, mchi2pmErr = 3.0, mchi2pmTh = 1.2;
  double mchi10Exp = 96.6, mchi10Err = 0.05, mchi10Th = 0.6;
  double mchi20Exp = 176.8, mchi20Err = 1.2, mchi20Th = 3.2;
  double mA0Exp = 394.9, mA0Err = 1.5, mA0Th = 2.8;
  double mHpExp = 403.0, mHpErr = 1.5, mHpTh = 3.1;
  double mH0Exp = 395.4, mH0Err = 1.5, mH0Th = 6.2;

  double chimmuR = 0., chimmuL = 0., chimeL = 0., 
    chimeR = 0., chimstau1 = 0., chimstau2 = 0., 
    chimchi1pm = 0., chimchi10 = 0., chimchi20 = 0., 
    chimA0 = 0., chimH0 = 0., chimHp = 0., chimstop1 = 0.,
    chimsnue = 0., chimchi2pm = 0.;
  
  double mmuRPred = r.displayPhys().me(2, 2),
    mmuLPred = r.displayPhys().me(1, 2),
    msnuePred = r.displayPhys().msnu(1),
    meLPred = r.displayPhys().me(1, 1),
    meRPred = r.displayPhys().me(2, 1),
    mstau1Pred = minimum(r.displayPhys().me(1, 3), r.displayPhys().me(2, 3)),
    mstop1Pred = minimum(r.displayPhys().mu(1, 3), r.displayPhys().mu(2, 3)),
    mstau2Pred = maximum(r.displayPhys().me(1, 3), r.displayPhys().me(2, 3)),
    mH0Pred = r.displayPhys().mhiggs(3),
    mA0Pred = r.displayPhys().mhiggs(2),
    mHpPred = r.displayPhys().mhiggs(4),
    mchi1pmPred = fabs(r.displayPhys().mch(1)),    
    mchi2pmPred = fabs(r.displayPhys().mch(2)),    
    mchi10Pred = fabs(r.displayPhys().mneut(1)),    
    mchi20Pred = fabs(r.displayPhys().mneut(2));
  
  // Include LC theory errors
  if (LCth == 1) {
    mmuRErr = sqrt(sqr(mmuRErr) + sqr(mmuRTh));
    mmuLErr = sqrt(sqr(mmuLErr) + sqr(mmuLTh));
    msnueErr = sqrt(sqr(msnueErr) + sqr(msnueTh));
    meRErr = sqrt(sqr(meRErr) + sqr(meRTh));
    meLErr = sqrt(sqr(meLErr) + sqr(meLTh));
    mstau1Err = sqrt(sqr(mstau1Err) + sqr(mstau1Th));
    mstop1Err = sqrt(sqr(mstop1Err) + sqr(mstop1Th));
    mstau2Err = sqrt(sqr(mstau2Err) + sqr(mstau2Th));
    mchi1pmErr = sqrt(sqr(mchi1pmErr) + sqr(mchi1pmTh));  
    mchi2pmErr = sqrt(sqr(mchi2pmErr) + sqr(mchi2pmTh));  
    mchi10Err = sqrt(sqr(mchi10Err) + sqr(mchi10Th));  
    mchi20Err = sqrt(sqr(mchi20Err) + sqr(mchi20Th));  
    mH0Err = sqrt(sqr(mH0Err) + sqr(mH0Th));  
    mA0Err = sqrt(sqr(mA0Err) + sqr(mA0Th));  
    mHpErr = sqrt(sqr(mHpErr) + sqr(mHpTh));  
  }

  chimmuR = sqr((mmuRExp - mmuRPred) / mmuRErr);
  chimmuL = sqr((mmuLExp - mmuLPred) / mmuLErr);
  chimeR = sqr((meRExp - meRPred) / meRErr);
  chimeL = sqr((meLExp - meLPred) / meLErr);
  chimsnue = sqr((msnueExp - msnuePred) / msnueErr);
  chimstau1 = sqr((mstau1Exp - mstau1Pred) / mstau1Err);
  chimstop1 = sqr((mstop1Exp - mstop1Pred) / mstop1Err);
  chimstau2 = sqr((mstau2Exp - mstau2Pred) / mstau2Err);
  chimchi1pm = sqr((mchi1pmExp - mchi1pmPred) / mchi1pmErr);
  chimchi2pm = sqr((mchi2pmExp - mchi2pmPred) / mchi2pmErr);
  chimchi10 = sqr((mchi10Exp - mchi10Pred) / mchi10Err);
  chimchi20 = sqr((mchi20Exp - mchi20Pred) / mchi20Err);
  chimH0 = sqr((mH0Exp - mH0Pred) / mH0Err);
  chimA0 = sqr((mA0Exp - mA0Pred) / mA0Err);
  chimHp = sqr((mHpExp - mHpPred) / mHpErr);

  /*  cout << mmuRPred << " " << mmuLPred << " " << msnumuPred << " " 
       << meLPred << " " << meRPred << " " << mstau1Pred << " " 
       << mstau2Pred << " " << msnutauPred << " " << mchi1pmPred 
       << " " << mchi10Pred << " " << mchi20Pred << endl;

  cout << chimmuR << " " << chimmuL << " " << chimeR << " " << chimeL 
       << " " << chimsnumu << " " << chimsnutau << " " << chimstau1 
       << " " << chimstau2 << " " << chimchi1pm << " " << chimchi10 
       << " " << chimchi20 << endl; */

  /*  double chisqLC = chimmuR + chimmuL + chimeR + chimeL + chimsnue + 
     chimstau1 + chimstop1 + chimstau2 + chimchi1pm + chimchi2pm + chimchi10 + 
     chimchi20 + chimH0 + chimA0 + chimHp; */
  double chisqLC = chimmuR + chimmuL + chimeR + chimeL + chimsnue + 
     chimstau1 + chimstau2 + chimchi1pm + chimchi10 + 
    chimchi20;

  // LHC errors from Giacomo's talk: ATLAS only, and central values from
  // ISAJET7.64 -- 300 fb^-1
  double mllmaxExp = 76.698, mllmaxErr = 0.03, mllmaxTh = 3.2;
  double mttmaxExp = 80.04, mttmaxErr = 5.06, mttmaxTh = 2.7;
  double mll4maxExp = 282.07, mll4maxErr = 2.3, mll4maxTh = 1.0;
  double mllqmaxExp = 426.12, mllqmaxErr = 1.4, mllqmaxTh = 15.4;
  double mllqminExp = 200.87, mllqminErr = 1.6, mllqminTh = 9.87;
  double mlqmaxLowExp = 299.4, mlqmaxLowErr = 0.87, mlqmaxLowTh = 10.6;
  double mlqmaxHighExp = 375.0, mlqmaxHighErr = 1.0, mlqmaxHighTh = 11.32;
  double mllbminExp = 182.9, mllbminErr = 3.6, mllbminTh = 8.43;
  double mSbot1gExp = 114.33, mSbot1gErr = 1.5, mSbot1gTh = 18.2;
  double mSbot2gExp = 85.73, mSbot2gErr = 2.5, mSbot2gTh = 19.7;
  double mgExp = 511.02, mgErr = 2.3, mgTh = 8.98;
  double mqrExp = 423.2, mqrErr = 10.0, mqrTh = 14.7;
  double mhExp = 113.6, mhErr = 1.0, mhTh = 1.59;
  double mlLchiExp = 105.34, mlLchiErr = 1.4, mlLchiTh = 2.57;
  
  if (LCexp == 1) mhErr = 0.05; // for LC fits

  // including LHC systematic errors
  mllmaxErr = 0.08;
  mll4maxErr = 2.32;
  mllqmaxErr = 4.49;
  mllqminErr = 2.57;
  mlqmaxLowErr = 3.12;
  mlqmaxHighErr = 3.88;
  mllbminErr = 4.04;
  mSbot1gErr = 1.89;
  mSbot2gErr = 2.64;
  mgErr = 5.6;
  mlLchiErr = 1.4;

  // include LHC theory errors
  if (LHCth == 1) {
    mllmaxErr = sqrt(sqr(mllmaxErr) + sqr(mllmaxTh));
    mttmaxErr = sqrt(sqr(mttmaxErr) + sqr(mttmaxTh));
    mll4maxErr = sqrt(sqr(mll4maxErr) + sqr(mll4maxTh));
    mllqmaxErr = sqrt(sqr(mllqmaxErr) + sqr(mllqmaxTh));
    mllqminErr = sqrt(sqr(mllqminErr) + sqr(mllqminTh));
    mlqmaxLowErr = sqrt(sqr(mlqmaxLowErr) + sqr(mlqmaxLowTh));
    mlqmaxHighErr = sqrt(sqr(mlqmaxHighErr) + sqr(mlqmaxHighTh));
    mllbminErr = sqrt(sqr(mllbminErr) + sqr(mllbminTh));
    mSbot1gErr = sqrt(sqr(mSbot1gErr) + sqr(mSbot1gTh));
    mSbot2gErr = sqrt(sqr(mSbot2gErr) + sqr(mSbot2gTh));
    mgErr = sqrt(sqr(mgErr) + sqr(mgTh));
    mqrErr = sqrt(sqr(mqrErr) + sqr(mqrTh));
    mhErr = sqrt(sqr(mhErr) + sqr(mhTh)); 
    mlLchiErr = sqrt(sqr(mlLchiErr) + sqr(mlLchiTh)); 
  }

  double chimllmax = 0.0, chimll4max = 0., chimllqmax = 0.0, chimllqmin = 0.0, 
    chimlqmaxLow = 0.0, chimlqmaxHigh = 0.0, chimllbmin = 0.0, 
    chimsbot1g = 0.0, chimsbot2g = 0., chimg = 0.0, chimqr = 0.0, chimh = 0.0,
    chimlLchi = 0., chimttmax;

  double mllMaxPred = mllMax(r), mll4MaxPred = mll4max(r), 
    mttMaxPred = mttMax(r), 
    mllqMaxPred = mllqMax(r), 
    mllqMinPred = mllqMin(r), mlqMaxPred = mlqMax(r),
    mlqMaxFarPred = mlqMaxFar(r), 
    mlqMaxHighPred = maximum(mlqMaxPred, mlqMaxFarPred), 
    mlqMaxLowPred = minimum(mlqMaxPred, spc(r)), 
    mllbMinPred = mllbMin(r),
    msbot1gPred = r.displayPhys().mGluino - 
    minimum(r.displayPhys().md(1, 3), r.displayPhys().md(2, 3)),
    msbot2gPred = r.displayPhys().mGluino - 
    maximum(r.displayPhys().md(1, 3), r.displayPhys().md(2, 3)),
    mgPred = r.displayPhys().mGluino - 0.99 * fabs(r.displayPhys().mneut(1)), 
    mqrPred = r.displayPhys().mu(2, 1) - fabs(r.displayPhys().mneut(1)),
    mhPred = r.displayPhys().mhiggs(1),
    mlLchiPred = r.displayPhys().me(1, 1) - fabs(r.displayPhys().mneut(1));

  chimllmax = sqr((mllmaxExp - mllMaxPred) / mllmaxErr);
  chimttmax = sqr((mttmaxExp - mttMaxPred) / mttmaxErr);
  chimll4max = sqr((mll4maxExp - mll4MaxPred) / mll4maxErr);
  chimllqmax = sqr((mllqmaxExp - mllqMaxPred) / mllqmaxErr);
  chimllqmin = sqr((mllqminExp - mllqMinPred) / mllqminErr);
  chimlqmaxLow = sqr((mlqmaxLowExp - mlqMaxLowPred) / mlqmaxLowErr);
  chimlqmaxHigh = sqr((mlqmaxHighExp - mlqMaxHighPred) / mlqmaxHighErr);
  chimllbmin = sqr((mllbminExp - mllbMinPred) / mllbminErr);
  chimsbot1g = sqr((mSbot1gExp - msbot1gPred) / mSbot1gErr);
  chimsbot2g = sqr((mSbot2gExp - msbot2gPred) / mSbot2gErr);
  chimg = sqr((mgExp - mgPred) / mgErr);
  chimqr = sqr((mqrExp - mqrPred) / mqrErr);
  chimh = sqr((mhExp - mhPred) / mhErr);
  chimlLchi = sqr((mlLchiExp - mlLchiPred) / mlLchiErr);
  
  double chisqLHC = chimllmax + chimll4max + chimllqmax + 
    chimllqmin + chimlqmaxLow + 
    chimlqmaxHigh + chimllbmin + chimsbot1g + chimsbot2g + chimg + 
    chimqr + chimh + chimlLchi + chimttmax; 
    	       
  /*  cout << mllMaxPred << " " << mttMaxPred << " " << mll4MaxPred << " " << 
    mllqMaxPred << " " << mllqMinPred << " " 
       << mlqMaxLowPred << " " << mlqMaxHighPred << " " << mllbMinPred << " " 
       <<    msbot1gPred << " "<< msbot2gPred << " " <<  mgPred << " " << mqrPred << " " << mhPred << " " << mlLchiPred<<  endl; // DEBUG

    cout << chimllmax << " " << chimttmax << " " << chimll4max << " " 
	 << chimllqmax << " " << chimllqmin << " " 
       << chimlqmaxLow << " " << chimlqmaxHigh << " " << chimllbmin << " " 
	 << chimsbot1g << " " << chimsbot2g << " " << chimg << " " << chimqr << " " << chimh << " " << chimlLchi; // DEBUG
  
         exit(1);
  */
  double chisq = 0.0;
  if (LHCexp == 1) chisq = chisq + chisqLHC;
  if (LCexp == 1)  chisq = chisq + chisqLC;  

  if (r.displayProblem().test()) { 
    //    chisq = chisq * 2.0;
  }

  return chisq;
}


double mll4max(const MssmSoftsusy & r) {

  double mchi10 = fabs(r.displayPhys().mneut(1)), 
    mchi40 = fabs(r.displayPhys().mneut(4)), 
    mer = r.displayPhys().me(1, 1);

  double sqrtArg = (sqr(mchi40) - sqr(mer)) *
    (sqr(mer) - sqr(mchi10));

  double answer = fabs(sqrtArg) /  sqr(mer);

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double mlqMaxFar(const MssmSoftsusy & r) {

  double mchi20 = r.displayPhys().mneut(2), mchi10 = r.displayPhys().mneut(1),
    msqL = 0.5 * (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1)), 
    mer = r.displayPhys().me(2, 1);

  double sqrtArg = (sqr(msqL) - sqr(mchi20)) * (sqr(mer) - sqr(mchi10)) 
    / (sqr(mer));

  double answer = sqrtArg;

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mllqMax(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut(1), 
    mchi20 = r.displayPhys().mneut(2), 
    msqL = 0.5 * (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1)), 
    mer = r.displayPhys().me(2, 1);

  double q = sqr(msqL), xi = sqr(mchi20), chi = sqr(mchi10), sl = sqr(mer);
  double m1 = (q - xi) * (xi - chi) / xi;
  double m2 = (q - sl) * (sl - chi) / sl;
  double m3 = (q * sl - xi * chi) * (xi - sl) / (xi * sl);

  if (sqr(sl) < q * chi && q * chi < sqr(xi) && sqr(xi) * chi < q * sqr(sl)) 
    return sqr(msqL - mchi10);

  double answer = maximum(m3, maximum(m1, m2));

  if (answer > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mllbMin(const MssmSoftsusy & r) {

  double mchi10 = fabs(r.displayPhys().mneut(1)), 
    mchi20 = fabs(r.displayPhys().mneut(2)), 
    msqL = minimum(r.displayPhys().md(1, 3), r.displayPhys().md(2, 3)), 
    mer = r.displayPhys().me(2, 1); 
 
  double sqrtArg = (sqr(mchi10) * sqr(mchi10) + sqr(mer) * sqr(mer)) * 
	     sqr(sqr(mchi20) + sqr(mer)) + 2.0 * sqr(mchi10) * sqr(mer) *
	     (sqr(mchi20) * sqr(mchi20) - 6.0 * sqr(mchi20) * sqr(mer) + 
	      sqr(mer) * sqr(mer));

  double answer = 
    ( - sqr(mchi10) * sqr(mchi20) * sqr(mchi20) + 
      3.0 * sqr(mchi10) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mer) * sqr(mer) - 
      sqr(mchi10) * sqr(mchi20) * sqr(msqL) - 
      sqr(mchi10) * sqr(mer) * sqr(msqL) + 
      3.0 * sqr(mchi20) * sqr(mer) * sqr(msqL) - 
      sqr(mer) * sqr(mer) * sqr(msqL) + (sqr(mchi20) - sqr(msqL)) 
      * sqrt(fabs(sqrtArg)))
    / (4.0 * sqr(mchi20) * sqr(mer));

  if (sqrtArg > 0.0 && answer > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double chiSqP1(const MssmSoftsusy & r) {

  // From ATLAS TDR: low lumi
   double mhExp = 95.4, mhErr = 1.0;
   double mhqMaxExp = 758.3, mhqMaxErr = 25.0;
   double msqrExp = 959., msqrErr = 40.;
   double mgExp = 1004., mgErr = 25.;
   double mstExp = 647., mstErr = 100.0;

  if (highLuminosity == 1) {
    mhErr = 1.0; mhqMaxErr = 13.0; msqrErr = 15.0; mgErr = 12.0; 
  }

  double chiMh, chiMhqMax, chiMsqr, chiMg, chiMst;
  chiMh = sqr((mhExp - r.displayPhys().mhiggs(1)) / mhErr);
  chiMhqMax = sqr((mhqMaxExp - mhqMax(r)) / mhqMaxErr);
  chiMsqr = sqr((msqrExp - 0.5 * (r.displayPhys().mu(2, 1) +
				  r.displayPhys().md(2, 1))) / msqrErr);
  chiMg = sqr((mgExp - r.displayPhys().mGluino) / mgErr);
  if (highLuminosity == 1) 
    chiMst = sqr((mstExp -
		  minimum(r.displayPhys().mu(1, 3), r.displayPhys().mu(2, 3)))
		 / mstErr);
  else chiMst = 0.0;

  double chisq;
  chisq = chiMh + chiMhqMax + chiMsqr + chiMg + chiMst;

  return chisq;
}

double chiSqP2(const MssmSoftsusy & r) {

  // From ATLAS TDR: low lumi
   double mhExp = 115.3, mhErr = 1.0;
   double mhqMaxExp = 758.3, mhqMaxErr = 25.0;
   double msqrExp = 959., msqrErr = 40.;
   double mgExp = 1004., mgErr = 25.;
   double mstExp = 713., mstErr = 100.0;

  if (highLuminosity == 1) {
    mhErr = 1.0; mhqMaxErr = 13.0; msqrErr = 15.0; mgErr = 12.0; 
  }

  double chiMh, chiMhqMax, chiMsqr, chiMg, chiMst;
  chiMh = sqr((mhExp - r.displayPhys().mhiggs(1)) / mhErr);
  chiMhqMax = sqr((mhqMaxExp - mhqMax(r)) / mhqMaxErr);
  chiMsqr = sqr((msqrExp - 0.5 * (r.displayPhys().mu(2, 1) +
				  r.displayPhys().md(2, 1))) / msqrErr);
  chiMg = sqr((mgExp - r.displayPhys().mGluino) / mgErr);

  if (highLuminosity == 1) 
    chiMst = sqr((mstExp -
		  minimum(r.displayPhys().mu(1, 3), r.displayPhys().mu(2, 3)))
		 / mstErr);
  else chiMst = 0.0;

  double chisq;
  chisq = chiMh + chiMhqMax + chiMsqr + chiMg + chiMst;

  return chisq;
}

double chiSqP3(const MssmSoftsusy & r) {

  // From ATLAS TDR: low lumi
  double mhExp = 68.5, mhErr = 3.0;
  double mneutExp = 52.42, mneutErr = 0.05;
  double mb1Exp = 278.1, mb1Err = 3.0;
  double msqExp = 320.5, msqErr = 20.0;
  double mgExp = 20.3, mgErr = 2.0;

  if (highLuminosity == 1) msqErr = 10.0;

  double sqAv;
  sqAv = 0.0;
  int j; 
  for (j=1; j<=2; j++) {
    sqAv = sqAv + r.displayPhys().mu(1, j) + r.displayPhys().md(1, j);
  }
  sqAv = sqAv / 4.0;

  double chiMh, chiMsq, chiMg, chiMneut, chib1;
  chiMh = sqr((mhExp - r.displayPhys().mhiggs(1)) / mhErr);
  chiMneut =  sqr((mneutExp - r.displayPhys().mneut(2) + 
		   r.displayPhys().mneut(1)) / mneutErr);
  chib1 = sqr((mb1Exp - minimum(r.displayPhys().md(1,3),
			    r.displayPhys().md(2,3))) / mb1Err);
  chiMsq = sqr((msqExp - sqAv) / msqErr);
  chiMg = sqr((mgExp - (r.displayPhys().mGluino - 
			minimum(r.displayPhys().md(1,3),
			    r.displayPhys().md(2,3)))) / mgErr);

  double chisq;
  chisq = chiMh + chiMsq + chiMg + chiMneut + chib1;

  return chisq;
}

double chiSqP4(const MssmSoftsusy & r) {

  // From ATLAS TDR: low lumi
  double mhExp = 118.1, mhErr = 1.0;
  double mneutExp = 68.7, mneutErr = 0.8;
  double mchExp = 315.0, mchErr = 20.0;
  double msqExp = 915.0, msqErr = 25.;
  double mgExp = 434.0, mgErr = 12.;
  
  if (highLuminosity == 1) {
    mneutErr = 0.25; mgErr = 6.0; mchErr = 7.0; 
  }

  double sqAv;
  sqAv = 0.0;
  int i, j; for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      sqAv = sqAv + r.displayPhys().mu(i, j) + r.displayPhys().md(i, j);
    }
  sqAv = sqAv / 8.0;

  double chiMh, chiMsq, chiMg, chiMneut, chiCh;
  chiMh = sqr((mhExp - r.displayPhys().mhiggs(1)) / mhErr);
  chiMneut =  sqr((mneutExp - r.displayPhys().mneut(2) + 
		   r.displayPhys().mneut(1)) / mneutErr);
  chiCh = sqr((mchExp - r.displayPhys().mch(2)) / mchErr);
  chiMsq = sqr((msqExp - sqAv) / msqErr);
  chiMg = sqr((mgExp - (r.displayPhys().mGluino - r.displayPhys().mneut(2))) 
	      / mgErr);

  double chisq;
  chisq = chiMh + chiMsq + chiMg + chiMneut + chiCh;

  return chisq;
}


double chiSqP5(const MssmSoftsusy & r) {
  double chisq, chiMllMax, chiR, chiMlqMax, chiMllqMin, chiMhqMin, chiMh, 
    chiMhqMax;

  // exptl data taken from ATLAS TDR: low lumi
  double mllMaxExp = 108.92, mllMaxErr = 0.5;
  double rExp = 0.865, rErr = 0.06;
  double mlqMaxExp = 478.1, mlqMaxErr = 11.5;
  double mllqMinExp = 271.8, mllqMinErr = 14.0;
  double mhqMinExp = 346.5, mhqMinErr = 17.0;
  double mhExp = 92.9, mhErr = 1.0;
  double mhqMaxExp = 552.5, mhqMaxErr = 10.0;

  if (highLuminosity == 1) {
    mllMaxErr = 0.1; mlqMaxErr = 5.0; rErr = 0.02; mllqMinErr = 5.4; 
  }

  // Warning - in general, you should check to see if they are available!
  chiMhqMax = sqr((mhqMaxExp - mhqMax(r)) / mhqMaxErr);
  chiMllMax = sqr((mllMaxExp - mllMax(r)) / mllMaxErr);
  chiR = sqr((rExp - edgeRatio(r)) / rErr);
  chiMlqMax = sqr((mlqMaxExp - mlqMax(r)) / mlqMaxErr);
  chiMllqMin = sqr((mllqMinExp - mllqMin(r)) / mllqMinErr);
  chiMhqMin = sqr((mhqMinExp - mhqMin(r)) / mhqMinErr);
  chiMh = sqr((r.displayPhys().mhiggs(1) - mhExp) / mhErr);

  chisq = chiMllMax + chiR + chiMlqMax + chiMllqMin + chiMhqMin + chiMh 
    + chiMhqMax;

  // Try to dissuade any minimisation from finding problem points.
  if (r.displayProblem().test()) chisq = chisq * factor;
  
  return chisq;
}

double chiSqP6(const MssmSoftsusy & r) {

  // From ATLAS TDR: low lumi
  double mhExp = 111.9, mhErr = 1.0;
  double mttExp = 59.6, mttErr = 3.0;
  double mb1Exp = 150.0, mb1Err = 20.0;
  double msqExp = 498.0, msqErr = 50.0;
  double mgExp = 540.0, mgErr = 60.0;
  
  if (highLuminosity == 1) {
    mttErr = 1.2; mgErr = 30.; mb1Err = 10.; msqErr = 25.;
  }

  double sqAv;
  sqAv = 0.0;
  int j; 
  for (j=1; j<=2; j++) {
    sqAv = sqAv + r.displayPhys().mu(2, j) + r.displayPhys().md(2, j);
  }
  sqAv = sqAv / 4.0;

  double chiMh, chiMsq, chiMg, chiMtt, chib1;
  chiMh = sqr((mhExp - r.displayPhys().mhiggs(1)) / mhErr);
  double mttt = mttMax(r);
  chiMtt =  sqr((mttExp - mttt) / mttErr);
  chib1 = sqr((mb1Exp - (r.displayPhys().mGluino - 
			 minimum(r.displayPhys().md(1,3),
				 r.displayPhys().md(2,3)))) / mb1Err);
  chiMsq = sqr((msqExp - sqAv) / msqErr);
  chiMg = sqr((mgExp - r.displayPhys().mGluino) / mgErr);

  double chisq;
  chisq = chiMh + chiMsq + chiMg + chiMtt + chib1;

  return chisq;
}

double mttMax(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut(1), 
    mchi20 = r.displayPhys().mneut(2), 
    mstau1 = minimum(r.displayPhys().me(1, 3),  r.displayPhys().me(2, 3));

  double sqArg1 = sqr(mchi20) - sqr(mstau1);
  double sqArg2 = 1.0 - sqr(mchi10) / sqr(mstau1);

  if (sqArg1 > 0.0 && sqArg2 > 0.0) return sqrt(sqArg1 * sqArg2);
  else 
    return -sqrt(fabs(sqArg1 * sqArg2));
}

// The following edge parameters are defined in Paige, Bachacou, Hinchcliffe
// PRD 62 015009 (2000). They return minus values if the chain doesn't exist
double mhqMax(const MssmSoftsusy & r) {

  double mh0 = r.displayPhys().mhiggs(1), mchi10 = r.displayPhys().mneut(1), 
    mchi20 = r.displayPhys().mneut(2), 
    msqL = 0.5 * (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1));

  double sqrtArg = sqr(sqr(mchi20) - sqr(mh0) - sqr(mchi10))
    - 4.0 * sqr(mh0) * sqr(mchi10);

  double answer = sqr(mh0) + (sqr(msqL) - sqr(mchi20)) * 
      (sqr(mchi20) + sqr(mh0) - sqr(mchi10) + 
       sqrt(fabs(sqrtArg))) / (2.0 * sqr(mchi20));

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double mllMax(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut(1), mchi20 = r.displayPhys().mneut(2), 
    mer = r.displayPhys().me(2, 1);

  double sqrtArg = (sqr(mchi20) - sqr(mer)) *
    (sqr(mer) - sqr(mchi10));

  double answer = fabs(sqrtArg) /  sqr(mer);

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double edgeRatio(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut(1), mchi20 = r.displayPhys().mneut(2),
    mer = r.displayPhys().me(2, 1);

  double sqrtArg = (sqr(mchi20) - sqr(mer)) / (sqr(mchi20) - sqr(mchi10));

  double answer = fabs(sqrtArg);

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(answer);
}

double mlqMax(const MssmSoftsusy & r) {

  double mchi20 = r.displayPhys().mneut(2), 
    msqL = 0.5 * (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1)), 
    mer = r.displayPhys().me(2, 1);

  double sqrtArg = (sqr(msqL) - sqr(mchi20)) * (sqr(mchi20) - sqr(mer)) 
    / (sqr(mchi20));

  double answer = sqrtArg;

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mllqMin(const MssmSoftsusy & r) {

  double mchi10 = r.displayPhys().mneut(1), 
    mchi20 = r.displayPhys().mneut(2), 
    msqL = 0.5 * (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1)), 
    mer = r.displayPhys().me(2, 1);

  double sqrtArg = (sqr(mchi10) * sqr(mchi10) + sqr(mer) * sqr(mer)) * 
	     sqr(sqr(mchi20) + sqr(mer)) + 2.0 * sqr(mchi10) * sqr(mer) *
	     (sqr(mchi20) * sqr(mchi20) - 6.0 * sqr(mchi20) * sqr(mer) + 
	      sqr(mer) * sqr(mer));

  double answer = 
    ( - sqr(mchi10) * sqr(mchi20) * sqr(mchi20) + 
      3.0 * sqr(mchi10) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mchi20) * sqr(mer) - 
      sqr(mchi20) * sqr(mer) * sqr(mer) - 
      sqr(mchi10) * sqr(mchi20) * sqr(msqL) - 
      sqr(mchi10) * sqr(mer) * sqr(msqL) + 
      3.0 * sqr(mchi20) * sqr(mer) * sqr(msqL) - 
      sqr(mer) * sqr(mer) * sqr(msqL) + (sqr(mchi20) - sqr(msqL)) 
      * sqrt(fabs(sqrtArg)))
    / (4.0 * sqr(mchi20) * sqr(mer));

  if (sqrtArg > 0.0 && answer > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

double mhqMin(const MssmSoftsusy & r) {

  double mh0 = r.displayPhys().mhiggs(1), mchi10 = r.displayPhys().mneut(1), 
    mchi20 = r.displayPhys().mneut(2), 
    msqL = 0.5 * (r.displayPhys().mu(1, 1) + r.displayPhys().md(1, 1));

  double sqrtArg = sqr(sqr(mchi20) - sqr(mh0) - sqr(mchi10)) - 
    4.0 * sqr(mchi10) * sqr(mh0);

  double answer = 0.5 / sqr(mchi20) * 
    (sqr(msqL) - sqr(mchi20)) * (sqr(mchi20) + sqr(mh0) - sqr(mchi10) -
     sqrt(fabs(sqrtArg)));

  if (sqrtArg > 0.0) return sqrt(answer);
  else return -factor * sqrt(fabs(answer));
}

//Produces a plot of parameters verses renormalisation scale in mSUGRA
void scaledParsSugra() {
  QedQcd oneset;
  readIn(oneset, "massIn");
  oneset.toMz();

  double mgut = 1.0e16;
  int sgnMu = -1;
  double m0 = 200.0, m12 = 100.0, a0 = 0.0, tanb = 2.0;
  cout << "#m0=" << m0 << " a0=" << a0 << " m12=" << m12 << " tanb=" << tanb
       << " signmu=" << sgnMu << endl;
  DoubleVector pars(3);
  translateSugra(pars, m0, m12, a0);

  MssmSoftsusy r;
  r.lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);  
  MssmSoftsusy s(r);

  double startlnMu = log(MZ), endlnMu = log(1.0e3);
  int numPoints = 50;
  int c; for (c=0; c<=numPoints; c++) {
    double lnmu = startlnMu + double(c) * (endlnMu - startlnMu) / numPoints;
    s.runto(exp(lnmu));
    cout << exp(lnmu) << " " 
	 << (s.displayMh2Squared() + sqr(s.displaySusyMu()))
      * (s.displayMh1Squared() + sqr(s.displaySusyMu())) -
      sqr(s.displayM3Squared())
	 << endl;
  }
  
}

// does a no-scale RPC scan
void contUniversalM0m12(int sgnMu, int accuracy) {

  double a0 = 0., m0 = 0.;
  double startm12 = 100.0, endm12 = 1000.0;
  double startTanb = 3.0, endTanb = 60.0;
  int numPoints = 10; 
      
  int count, count2;

  cout << "%     m12     " << "   tanb    ";
  printShortInitialise(); cout << endl;

  for (count = 0; count <= numPoints; count++) {
    for (count2 = 0; count2 <= numPoints; count2++) { 

      double mgut = 1.9e16;
      double m12 = startm12 + count * (endm12 - startm12) / double(numPoints);
      double tanb = startTanb + count2 * (endTanb - startTanb) / 
	double(numPoints);

      QedQcd oneset;
      readIn(oneset, "massIn");

      oneset.toMz();
      
      DoubleVector pars(3); 
      
      MssmSoftsusy r;   
        
      translateSugra(pars, m0, m12, a0);

      r.lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
      cout << m12 << " " << tanb << " ";
      r.printShort(); 
      double mass; int posi, posj, id;
      cout << " % " << r.displayProblem();
      id = r.lsp(mass, posi, posj);
      recogLsp(id, posj); 
      cout << endl; 
    }
  }
}

// for a given string scale and sign of mu parameter, produces a 50x50 grid of
// soft mass in m32 sintheta vs tan beta plane
// accuracy: 0 no SUSY thresholds except for MSusy
// 1: SUSY particle thresholds
// 2: threshold corrections to SUSY parameters
// susyBreaking: 
// 1: dilaton domination
// 2: moduli domination
void contDomination(double mgut, int sgnMu, int accuracy) {
  // Grid of (numPoints+1)x(numPoints+1)
  int numPoints = 15;

  double startM32 = 50.0, endM32 = 350.0;
  double startTanb = 2.0, endTanb = 47.0;
      
  int count, count2;
  cout << "%    m32     " << "     tanb   ";  
  printShortInitialise();
  cout << "    cMu0   " << "     cB0  " << "     cM32  " 
       << "    cHT    " << endl;

  DoubleVector ft(4);
  MssmSoftsusy r;

  for (count2 = 0; count2 <= numPoints; count2++)
    for (count = 0; count <= numPoints; count++) {
	double m32 = startM32 + count * (endM32 - startM32) / numPoints;
	double tanb = startTanb + count2 * (endTanb - startTanb) /
	  numPoints;
	double m0 = m32, m12 = sqrt(3.0) * m32, a0 = -m12;

	r = doPointUniversalM0m12(m0, m12, tanb, a0, mgut, sgnMu,
				      accuracy, ft);

	cout << m32 << " " << tanb << " ";
	r.printShort(); cout << ft;
	}

  cout.flush(); return;
}

inline void translateSugra(DoubleVector & pars, double m0, double m12,
			    double a0) {
  pars(1) = m0; pars(2) = m12; pars(3) = a0;
}

MssmSoftsusy pointDomination(double tanb, double
			       mgut, int sgnMu, int accuracy, DoubleVector &
			       ft, double m32,
			       double epsilon) {
  QedQcd oneset;
  readIn(oneset, "massIn");
  oneset.toMz();
  
  double m0, m12, a0;

  DoubleVector pars(3); 
  
  m0 = m32; m12 = sqrt(3.0) * m32; a0 = -m12;
  translateSugra(pars, m0, m12, a0);
  ft.setEnd(4);
  MssmSoftsusy r;
  r.lowOrg(sugraBcs, mgut, pars, sgnMu,  tanb, oneset, true); 
  
  return r;
}

void nonUniGauginos(MssmSoftsusy & m, const DoubleVector & inputParameters)
{
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  // Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  // trilinear scalar coupling to be a0
  m.standardSugra(m0, m12, a0);
  m.setGauginoMass(1, 1.6 * m12);
    
  return;
}

// Scans through m12 given other SUGRA parameters, used for making 1d plots.
void m12Scan() {
  
  QedQcd oneset;
  readIn(oneset, "massIn");
  oneset.toMz();

  // Default parameters
  int sgnMu = 1;
  double m0 = 400., M12 = 400., a0 = 0.0, tanb = 10., mgut = 1.9e16;

  double endM12 = 700., startM12 = 100.;
  //  double endLambda = 150.0e3, startLambda = 30.0e3;
  int i, numPts = 50; 
  
  DoubleVector pars(3);

  void (*bCType)(MssmSoftsusy & m, const DoubleVector & inputParameters) =
    sugraBcs;
  
  int option;
  /*
  for (option = 1; option <= 8; option++) {

    cout << "Option: " << option << endl;
    
    for (i=0; i<=numPts; i++)
      {
	
	MssmSoftsusy r;
	
	switch(option) {
	case 1:
	  // Model line A
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = 0.4 * M12;
	  a0 = -0.4 * M12;
	  tanb = 10.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &sugraBcs;
	  break;
	case 2:
	  // Model line B
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = 0.5 * M12;
	  a0 = 0.0;
	  tanb = 10.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &nonUniGauginos;
	  break;
	case 3:
	  // Model line C
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = M12;
	  a0 = 0.;
	  tanb = 35.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &sugraBcs;
	  break;
	case 4:
	  // Model line D
	  { 
	    double lambda = (endLambda - startLambda) / double(numPts) 
	      * double(i) + startLambda; 
	    double mMess = 2.0 * lambda;
	    tanb = 15.; n5 = 3;
	    pars.setEnd(2);
	    pars(1) = mMess; pars(2) = lambda;
	    mgut = mMess;
	    cout << lambda;
	  }
	  bCType = &gmsbBcs;
	  break;
	case 5:
	  // Model line E
	  { 
	    endLambda = 250.0e3;
	    double lambda = (endLambda - startLambda) / double(numPts) 
	      * double(i) + startLambda; 
	    double mMess = 2.0 * lambda;
	    tanb = 15.; n5 = 1;
	    pars.setEnd(2);
	    pars(1) = mMess; pars(2) = lambda;
	    mgut = mMess;
	    cout << lambda;
	  }
	  bCType = &gmsbBcs;
	  break;
	case 6:
	  // Model line F
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = 2.0 * M12 + 800.;
	  a0 = 0.;
	  pars.setEnd(3);
	  tanb = 10.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &sugraBcs;
	  break;
	case 7:
	  // Model line G
	  {
	    double endM32 = 150., startM32 = 400.0;
	    double M32 = 
	      (endM32 - startM32) / double(numPts) * double(i) + startM32;

	    m0 = M32;
	    M32 = 50.0e3;
	    pars.setEnd(2);
	    pars(1) = M32; pars(2) = m0;
	    cout << m0;
	  }
	  tanb = 10.;
	  bCType = &amsbBcs;
	  break;
	case 8:
	  // Model line F
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = 0.25 * M12 - 9.0;
	  a0 = 0.;
	  pars.setEnd(3);
	  tanb = 10.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &sugraBcs;
	  break;
      default:
	cerr << "m12scan called with illegal " << option << "\n";
	exit(1); 
	break;
	}
      
      r.lowOrg(bCType, mgut, pars, sgnMu, tanb, oneset);

      if (option == 1 || option == 2 || option == 3 || option == 6) 
	cout << M12; 

      r.printShort();   

      //      r.runto(r.calcMs());
      //cout << r; exit(1); //DEBUG

      if (r.displayProblem().test()) cout << "%" << r.displayProblem(); 
      cout << endl;      
      }

  cout << endl << endl;
*/

  for (option = 1; option <= 8; option++) {
    for (i=0; i<=numPts; i++) {
	
	MssmSoftsusy r;
	
	switch(option) {
	case 1:
	  // Model line A
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = 0.4 * M12;
	  a0 = -0.4 * M12;
	  tanb = 10.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &sugraBcs;
	  break;
	case 2:
	  // Model line B
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = 0.5 * M12;
	  a0 = 0.0;
	  tanb = 10.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &nonUniGauginos;
	  break;
	case 6:
	  // Model line F
	  M12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
	  m0 = 2.0 * M12 + 800.;
	  a0 = 0.;
	  pars.setEnd(3);
	  tanb = 10.;
	  translateSugra(pars, m0, M12, a0);
	  bCType = &sugraBcs;
	  break;
      }
      
  	if (option==1 || option ==2 || option ==6){
	  r.lowOrg(bCType, mgut, pars, sgnMu, tanb, oneset, true);
	  
	  if (option == 1 || option == 2 || option == 3 || option == 6) 
	    cout << M12; 
	  
	  r.printShort(); 

	  if (r.displayProblem().test()) cout << "%" << r.displayProblem(); 
	  cout << endl;      
	} // option if
    }  // numpts loop
  

      cout << endl << endl;
  }  // option loop
}

// Runs one point in SUGRA space, returning Fine Tuning vector in ft
MssmSoftsusy doPointUniversalM0m12(double m0, double m12, double tanb,
					double a0, double mgut, int sgnMu, 
					int accuracy, DoubleVector & ft) {
  QedQcd oneset;
  readIn(oneset, "massIn");

  oneset.toMz();

  DoubleVector pars(3); translateSugra(pars, m0, m12, a0);

  MssmSoftsusy r;
  r.lowOrg(sugraBcs, mgut, pars, sgnMu, tanb, oneset, true);
  ft = r.fineTune(sugraBcs, pars, mgut);  

  return r;
}

// writes one universal sugra point to standard output
void pointUniversalM0m12(double m0, double m12, double tanb, double a0,
			    double mgut, int sgnMu, int accuracy) {
  DoubleVector ft(6);

  // none=sugra/dil for FT calculation
  MssmSoftsusy r(doPointUniversalM0m12(m0,  m12,  tanb, a0,  mgut,
					    sgnMu, accuracy, ft));

  cout << m0 << " " << m12; cout.flush();  
  r.printShort(); cout << ft;
}

void iterateBound(int maxTries, double tol, double & oldBound, 
		  double m0, 
		  double m12, double a0, double tanb, double &mgut, int sgnMu,
		  int & err, RpvCouplings couplingType, int i, int j, int k,
		  const QedQcd & oneset, double neutrinoBound, 
		  RpvSoftsusy & kw) {
  static int numTries = 0;
  static double newBound;
  const double startBound = 1.0e-4;

  if (numTries == 0) newBound = startBound;

  if (numTries - 1 > maxTries) {     
    numTries = 0; newBound = startBound;
    err = 1; return;
  }
  
  DoubleVector pars(52); 
  RpvSusyPars bc;
  
  // How close to convergence are we?
  double c = 1.0 - 
    minimum(fabs(oldBound), fabs(newBound)) / 
    maximum(fabs(oldBound), fabs(newBound));
    
  if (c < tol) { 
    oldBound = newBound;
    mgut = fabs(mgut);
    numTries = 0; newBound = startBound;

    return; 
  }
	    
  numTries = numTries + 1;
  
  oldBound = newBound;

  bc.setLambda(couplingType, i, j, k, oldBound);

  pars = boundaryCondition(m0, m12, a0, bc);
  const RpvSoftsusy empty;
  kw = empty;
    
  mgut = kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset, true);  

  newBound = minimum(sqrt(neutrinoBound / neutrinoSum(kw)) * oldBound, 1.0);

  /*  cout << "mnu=" << neutrinoSum(kw)
       << " nb=" << newBound << " " << " ob=" << oldBound << " " 
       << kw.displayLambda(couplingType)(i, j, k) << endl; */

  mgut = fabs(mgut);

  iterateBound(maxTries, tol, oldBound, m0, m12, a0, tanb, mgut, 
	       sgnMu, err, couplingType, i, j, k, oneset, neutrinoBound, kw);  
}

// Tachyon bound - adjusts oldStart and oldEnd each time to sandwich the bound
// Initially, startis Ok and end is Ok should be FALSE
void tachyonBound(int maxTries, double tol, double & oldStart, double & oldEnd,
		  bool & startIsOk, bool & endIsOk, double m0, 
		  double m12, double a0, double tanb, double &mgut, int sgnMu,
		  int & err, RpvCouplings couplingType, int i, int j, int k,
		  const QedQcd & oneset, double neutrinoBound, 
		  RpvSoftsusy & kw) {
  static int numTries = 0;
  const RpvSoftsusy empty;

  if (numTries - 1 > maxTries) {     
    numTries = 0; 
    err = 1; return;
  }

  DoubleVector pars(52); 
  RpvSusyPars bc;
  
  // How close to convergence are we?
  double c = 1.0 - 
    minimum(fabs(oldStart), fabs(oldEnd)) / 
    maximum(fabs(oldStart), fabs(oldEnd));

  if (c < tol) { 
    mgut = fabs(mgut);
    numTries = 0; 
    return; 
  }

  numTries = numTries + 1;
  
  bool problemBeginning = false, problemMiddle = false, problemEnd = false;

  bc.setLambda(couplingType, i, j, k, oldStart);
  pars = boundaryCondition(m0, m12, a0, bc);

  if (!startIsOk) {
        kw = empty;
  mgut =  kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset, true);  
  
  if (kw.displayProblem().testSeriousProblem()) 
    problemBeginning = true; // flags a problem
    if (problemBeginning) {
      oldStart = oldStart * 0.5;
      if (PRINTOUT > 0) cout << " Trying start=" << oldStart << endl;
      mgut = 1.0e16;
      tachyonBound(maxTries, tol, oldStart, oldEnd, startIsOk, endIsOk, m0, 
		   m12, a0, tanb, mgut,  
		   sgnMu, err, couplingType, i, j, k, oneset, neutrinoBound, 
		   kw);
    }
    else {
      startIsOk = true;
      if (PRINTOUT > 0) cout << " Found start=" << oldStart << endl;
    }
  }

  if (!endIsOk) {
    bc.setLambda(couplingType, i, j, k, oldEnd);
    pars = boundaryCondition(m0, m12, a0, bc);
        kw = empty;
    mgut = 1.0e16;
    mgut = kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset, true);  
    if (kw.displayProblem().testSeriousProblem()) 
      problemEnd = true; // flags a problem
    if (!problemEnd) {
      oldEnd   = oldEnd * 2.0;
      if (PRINTOUT > 0) cout << " Trying end=" << oldEnd << endl;
      mgut = 1.0e16;
      tachyonBound(maxTries, tol, oldStart, oldEnd, startIsOk, endIsOk, m0, m12, 
		   a0, tanb, mgut, 
		   sgnMu, err, couplingType, i, j, k, oneset, neutrinoBound, 
		   kw);
    }
    else {
      endIsOk = true;
      if (PRINTOUT > 0) cout << " Found end=" << oldEnd << endl;
    } 
  }

  bc.setLambda(couplingType, i, j, k, 0.5 * (oldStart + oldEnd));
  pars = boundaryCondition(m0, m12, a0, bc);
     kw = empty;
  mgut = 1.0e16;
  mgut = kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset, true);  
  if (kw.displayProblem().testSeriousProblem()) 
    problemMiddle = true; // flags a problem

  if (problemMiddle)  {
    oldEnd   = 0.5 * (oldStart + oldEnd);
    if (PRINTOUT > 0) cout << " Mid prob. end=" << oldEnd << endl;
  }
  else {
    oldStart   = 0.5 * (oldStart + oldEnd);
    if (PRINTOUT > 0) cout << " No mid prob. start=" << oldStart << endl;
  }

  mgut = 1.0e16;
  tachyonBound(maxTries, tol, oldStart, oldEnd, startIsOk, endIsOk, m0, m12, 
	       a0, tanb, mgut, sgnMu,
	       err, couplingType, i, j, k, oneset, neutrinoBound, kw);
} // NB a lot of effort is wasted re-doing the bounds. 


void unifyBounds(RpvCouplings couplingType, int i, int j, int k, double m0, 
		 double m12, double a0, double tanb, double & mgut, 
		 int sgnMu, double & tachBound, double & neutBound, 
		 double & tachBoundMZ, double & neutBoundMZ,
		 const QedQcd & oneset, RpvSoftsusy & kw) { 

  int maxTries = 20, err;
  double tol = TOLERANCE * 10.0, oldBound, neutrinoBound = 0.71e-9;
  double start, end;
  
  // First find out if the point is allowed even with zero coupling!
  RpvSusyPars bc;   DoubleVector pars(52); 
  bc.setLambda(couplingType, i, j, k, 0.0);
  pars = boundaryCondition(m0, m12, a0, bc);
  kw = RpvSoftsusy();

  mgut = kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset, true);  
  if (kw.displayProblem().testSeriousProblem()) 
    { // if that's the case, return zero bounds
      tachBound = 0.; tachBoundMZ = 0.; neutBound = 0.; neutBoundMZ = 0.;
      return;
    } 

  start = 0.02; end = 0.2;

  bool startIsOk = false, endIsOk = false; mgut = 1.0e16;

  tachyonBound(maxTries, tol, start, end, startIsOk, endIsOk, m0, m12, 
	       a0, 
	       tanb, mgut, sgnMu, err, couplingType, i, j, k, 
	       oneset, neutrinoBound, kw);  
  tachBound = (start + end) * 0.5; 

  static DoubleMatrix ckm(3, 3);
  if (ckm(1, 1) == 0.0) {
    ckm(1, 1) =  0.9748;  ckm(1, 2) =  0.2229,  ckm(1, 3) = 0.0036;
    ckm(2, 1) = -0.2229;  ckm(2, 2) =  0.9740,  ckm(2, 3) = 0.0412;
    ckm(3, 1) =  0.0057;  ckm(3, 2) = -0.0410,  ckm(3, 3) = 0.9991;
  }
  
  // note MZ bounds should be rotated into the mass eigenbasis, since that's
  // the basis in which people will use the bounds
  if (couplingType == LD && MIXING > 0) {
    if (MIXING == 2)
      kw.setLambda(couplingType, kw.displayLambda(couplingType).product(ckm));
    if (MIXING == 1)
      kw.setLambda(couplingType, kw.displayLambda(couplingType) *
		   ckm.transpose()); 
  }
  tachBoundMZ = kw.displayLambda(couplingType)(i, j, k);
  

  oldBound = 0.01; mgut = 1.0e16;
  iterateBound(maxTries, tol, oldBound, m0, m12, a0, 
	       tanb, mgut, sgnMu, err, couplingType, i, j, k, 
	       oneset, neutrinoBound, kw);  
  
  neutBound = oldBound;
  if (couplingType == LD && MIXING > 0) {
    if (MIXING == 2)
      kw.setLambda(couplingType, kw.displayLambda(couplingType).product(ckm));
    if (MIXING == 1)
      kw.setLambda(couplingType, kw.displayLambda(couplingType) *
		   ckm.transpose()); 
  }
  neutBoundMZ = kw.displayLambda(couplingType)(i, j, k);
}

/*
// Makes a file plottable in the RPVcoupling-M12 plane (no-scale SUGRA)
void rpvScan(RpvCouplings couplingType, int i, int j, int k, double m12Start, 
		 double m12End, double lStart, double lEnd, double tanb, 
		 int sgnMu, const QedQcd & oneset) { 

  double m12, lambda;
  const int numPoints = 10; 

  int ii, jj; for (ii=0; ii<=numPoints; ii++)
    for (jj=0; jj<=numPoints; jj++) {
      double mgut = -1.0e16, m0 = 0.0, a0 = 0.0;
      m12 = (m12End - m12Start) / double(numPoints) * double(ii) + m12Start;
      lambda = (lEnd - lStart) / double(numPoints) * double(jj) + lStart;
      
      RpvSusyPars bc;   DoubleVector pars(52); 
      bc.setLambda(couplingType, i, j, k, lambda);
      pars = boundaryCondition(m0, m12, a0, bc);
      RpvSoftsusy kw;
      kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset);  
      cout << lambda << " " << m12 << " ";
      double mass; int posi, posj, id;
      id = kw.lsp(mass, posi, posj);
      recogLsp(id, posj); cout << " " << mass << " " << kw.displayProblem() <<
      endl;
    }
}
r*/

RpvSoftsusy rpvPoint(RpvCouplings couplingType, int i, int j, int k, 
		     double tanb, 
		     double lambda, double m12, int sgnMu, 
		     const QedQcd & oneset, double & mgut) { 
  double m0 = m12, a0 = 0.0; 
  
  RpvSusyPars bc;   DoubleVector pars(52); 
  bc.setLambda(couplingType, i, j, k, lambda);
  pars = boundaryCondition(m0, m12, a0, bc);
  RpvSoftsusy kw;
  mgut = kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset, true);  
  //  kw.runto(mgut);
  //kw.methodBoundaryCondition(pars);
  //kw.runto(MZ);
  //  kw.runto(mgut); cout << kw;

  return kw;
}


// Makes a file plottable in the RPVcoupling-M12 plane (no-scale SUGRA)
void rpvScan(RpvCouplings couplingType, int i, int j, int k, double tbStart, 
		 double tbEnd, double lStart, double lEnd, double m12, 
		 int sgnMu, const QedQcd & oneset) { 

  const int numPoints = 10; 

  int ii=0, jj=0; 
  for (ii=0; ii<=numPoints; ii++)
    for (jj=0; jj<=numPoints; jj++) { 
      double tanb = (tbEnd - tbStart) / double(numPoints) * double(ii) + 
	tbStart;
      double lambda = (lEnd - lStart) / double(numPoints) * double(jj) + 
	lStart;

      RpvSoftsusy kw; double mgut = 1.0e16;
      kw = rpvPoint(couplingType, i, j, k, tanb, lambda, m12, sgnMu, 
		    oneset, mgut);
      cout << lambda << " " << tanb << " "; kw.printShort();
      cout << neutrinoSum(kw);
      double mass; int posi, posj, id;
      cout << " % " << kw.displayProblem();
      id = kw.lsp(mass, posi, posj);
      recogLsp(id, posj); 
      cout << endl;
    }

}

void organiseBounds() {
  int sgnMu = 1; double mgut = 1.0e16; 
  double m0 = 0., m12 = 500., tanb = 10., a0 = 0.;

  QedQcd oneset;
  readIn(oneset, "massIn");

  oneset.toMz();

  // Set the initialboundary condition for the SUSY parameters:
  RpvCouplings couplingType = LD;
  RpvSoftsusy kw, empty;

  DoubleVector pars(52); RpvSusyPars bc; 
  int i, j=0, k=0;
  double tachBound, neutBound, neutBoundMZ, tachBoundMZ;
    
   cout << "%No-scale scan in m12, tanb, no mixing\n"; MIXING = 0;
   cout << "%Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
        << TOLERANCE << endl << oneset << endl;
   contUniversalM0m12(sgnMu, 3);
   // parameters used for scans
   double tbStart = 3.0, tbEnd = 40., lStart = 0., lEnd = 0.7, m12try = 500.;
   MIXING = 2; // down mixing 
   cout << "%l_231 vs tan beta scan for m12=500 no-scale, down mixing\n";
   cout << "%Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
        << TOLERANCE << endl << oneset << endl;
   rpvScan(LE, 1, 2, 3, tbStart, tbEnd, lStart, lEnd, m12try, sgnMu, oneset);  
   cout << "%l''_323 vs tan beta scan for m12=500 no-scale, down mixing\n";
   cout << "%Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
        << TOLERANCE << endl << oneset << endl;
   rpvScan(LU, 3, 2, 3, tbStart, tbEnd, lStart, lEnd, m12try, sgnMu, oneset);  
   lEnd = 0.2; MIXING = 1; // up mixing
   cout << "%l'_231 vs tan beta scan for m12=500 no-scale, up mixing\n";
   cout << "%Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
        << TOLERANCE << endl << oneset << endl;
	   rpvScan(LD, 1, 2, 3, tbStart, tbEnd, lStart, lEnd, m12try, sgnMu, oneset);


   cout << "%l'333" << endl; MIXING = 2;
   cout << "%Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
        << TOLERANCE << endl << oneset << endl;
   int numPts = 10;
   for (i=0; i<=numPts; i++)
     for (j=0; j<=numPts; j++) { 
       double startTanb = 2.0, endTanb = 50., startM12 = 100., endM12 = 1000.;
       mgut = 1.0e16;
       a0 = 0.; m0 = 0.;
       m12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
       tanb = (endTanb - startTanb) / double(numPts) * double(j) + startTanb;

       unifyBounds(LD, 3, 3, 3, m0, m12, a0, tanb, mgut, sgnMu, 
 		  tachBound, neutBound, tachBoundMZ, neutBoundMZ, oneset, 
 		  kw);
       cout << m12 << " " << tanb << " " << tachBound << " " << neutBound
 	   << " " << tachBoundMZ << " " << neutBoundMZ << " "; 
       kw.printShort(); cout << endl;
     }
   exit(0);
 
   MIXING = 2;
   cout << "%l231: Low energy data in SOFTSUSY: MIXING=" << MIXING 
        << " TOLERANCE=" 
        << TOLERANCE << endl << oneset << endl;

   for (i=0; i<=numPts; i++)
     for (j=0; j<=numPts; j++) { 
       double startTanb = 2.0, endTanb = 50., startM12 = 100., endM12 = 1000.;
       mgut = -1.0e16;
       a0 = 0.; m0 = 0.;
       m12 = (endM12 - startM12) / double(numPts) * double(i) + startM12;
       tanb = (endTanb - startTanb) / double(numPts) * double(j) + startTanb;
       unifyBounds(LE, 1, 2, 3, m0, m12, a0, tanb, mgut, sgnMu, 
 		  tachBound, neutBound, tachBoundMZ, neutBoundMZ, oneset, 
 		  kw);
       cout << m12 << " " << tanb << " " << tachBound << " " << neutBound
        	   << " " << tachBoundMZ << " " << neutBoundMZ << " "; 
       kw.printShort(); cout << endl;
     }
  
  m0 = 100.; m12 = 250.; a0 = -100.; tanb = 10; 
  
  int l;
  for (l=0; l<=2; l++) {
    MIXING = l;
    cout << "%Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
	 << TOLERANCE << endl << oneset << endl;


    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++)
    	for (k=1; k<=3; k++) {

	  couplingType = LD;
	  
	  mgut = 1.0e16;
	  unifyBounds(couplingType, i, j, k, m0, m12, a0, tanb, mgut, sgnMu, 
		      tachBound, neutBound, tachBoundMZ, neutBoundMZ, oneset, 
		      kw);

	  kw.runto(kw.displayMsusy());

         cout << "lambda'_{" << j << "," << k << "," << i << "}(M_GUT) < "
              << tachBound << " " << neutBound << " " << " (MZ) < "
              << tachBoundMZ << " " << neutBoundMZ << endl;
	}
  }
  
  MIXING = 0;
  int m;
  for (i=1; i<=3; i++)
    for (m=1; m<=3; m++) { 
      
      couplingType = LE;
      
      if (m==1) { j = 1; k = 2; }
      if (m==2) { j = 1; k = 3; }
      if (m==3) { j = 2; k = 3; }
      mgut = 1.0e16;	
      unifyBounds(couplingType, i, j, k, m0, m12, a0, tanb, mgut, sgnMu, 
		  tachBound, neutBound, tachBoundMZ, neutBoundMZ, oneset, 
		  kw);
      
      cout << "lambda_{" << j << "," << k << "," << i << "}(M_GUT) < "
	   << tachBound << " " << neutBound << " " << " (MZ) < " 
	   << tachBoundMZ << " " << neutBoundMZ << endl;
    }
  
  return;
}

//  as rpvPoint, but including m0 and a0 as arguments
RpvSoftsusy rpvPointII(RpvCouplings couplingType, int i, int j, int k, 
		       double tanb, 
		       double lambda, double m12, double m0, double a0, int sgnMu, 
		       const QedQcd & oneset, double & mgut) { 
  
  RpvSusyPars bc;   DoubleVector pars(52); 
  bc.setLambda(couplingType, i, j, k, lambda);
  pars = boundaryCondition(m0, m12, a0, bc);
  RpvSoftsusy kw;
  mgut = kw.lowOrg(userDefinedBcs, mgut, pars, sgnMu, tanb, oneset, true);  
  //  kw.runto(mgut);
  //kw.methodBoundaryCondition(pars);
  //kw.runto(MZ);
  //  kw.runto(mgut); cout << kw;

  return kw;
}

