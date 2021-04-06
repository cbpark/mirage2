
/** \File rpvsoft.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: Header file for RP violating MSSM object including all (real)
                soft SUSY breaking parameters and (real) SUSY couplings.

*/

#include "rpvsoft.h"

extern int PRINTOUT, MIXING;
extern double TOLERANCE; 

const RpvSoftsusy & RpvSoftsusy::operator = (const RpvSoftsusy &s) {
  if (this == &s) return *this;
  snuVevs = s.snuVevs;
  setSoftsusy(s.displaySoftsusy());
  setRpvSusyPars(s.displayRpvSusy());
  setRpvSoftPars(s.displayRpvSoft());
  return *this;
}

#define HR "----------------------------------------------------------\n"

ostream & operator <<(ostream &left, const RpvSoftsusy & r) {

  left << HR;
  left << "sneutrino Vevs:" << r.displaySneutrinoVevs();
  left << "RPV soft breaking parameters at Q=" << r.displayMu() << endl;
  left << r.displayRpvSoft();
  left << HR;
  left << "RPV supersymmetric parameters at Q=" << r.displayMu() << endl;
  left << r.displayRpvSusy();
  left << r.displayMssmSoft();
  return left;
}

DoubleVector RpvSoftsusy::display() const {
  DoubleVector parameters(SoftParsMssm::display());
  int k = parameters.displayEnd() + 1;
  parameters.setEnd(numRpvSoftPars);
  RpvSusyPars::display(parameters, k); 
  RpvSoftPars::display(parameters, k); 
  return parameters;
}

void RpvSoftsusy::set(const DoubleVector & v) {
  MssmSoftsusy::set(v);
  int k = numSoftParsMssm + 1;
  RpvSusyPars::set(v, k);
  RpvSoftPars::set(v, k);
}

DoubleVector RpvSoftsusy::beta() const {
  //  dd.setEnd(numRpvSoftPars);
  return (RpvSoftsusy::beta2()).display();
}

// Outputs the anomalous dimensions, of RPV SUSY effects.
// CHECKED: 23/05/02
void RpvSoftsusy::rpvAnomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
				    DoubleMatrix & gQQ, DoubleMatrix & gUU,
				    DoubleMatrix & gDD, 
				    double & gH1H1, double & gH2H2,
					DoubleVector & gH1L) const { 
  // To keep this a const function
  DoubleMatrix d1(displayYukawaMatrix(YD)), e1(displayYukawaMatrix(YE)); 
  // Make life simpler
  Tensor le1(displayLambda(LE)), ld1(displayLambda(LD)),
    lu1(displayLambda(LU));

  // RPV bits of MSSM anomalous dimensions
  if (displayLoops() > 0) {
      const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
      gEE = oneO16Pisq * matrixify(le1.transpose(), le1);
      gLL = oneO16Pisq * (sumProd(le1, le1.transpose()) + 
	 3.0 * sumProd(ld1, ld1.transpose()));
      gQQ = oneO16Pisq * 
	(sumProd(ld1.transpose(), ld1));
      gDD = oneO16Pisq * 2.0 * (sumProd(lu1, lu1.transpose()) 
	     + matrixify(ld1.transpose(), ld1));
      gUU = oneO16Pisq * matrixify(lu1.transpose(), lu1);
      // Convention for this is that it's L on top: 
      gH1L = -oneO16Pisq * (3.0 *(ld1 * d1).trace(2) + 
			    (le1 * e1).trace(2)); // DEBUG +1 orig
    }
}

// RPV bits of anomalous dimensions of superfields: CHECKED 23/05/02
void RpvSoftsusy::rpvAnomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
				    DoubleMatrix & gQQ, DoubleMatrix & gUU,
				    DoubleMatrix & gDD, 
				    double & gH1H1, double & gH2H2,
				    DoubleVector & gH1L) const { 
  Tensor let(displayLambda(LE).transpose()), 
    ldt(displayLambda(LD).transpose()), 
    lut(displayLambda(LU).transpose());
  Tensor her(displayHr(LE)), hdr(displayHr(LD)), hur(displayHr(LU));
  DoubleMatrix he(displayTrilinear(EA)), hd(displayTrilinear(DA));
  
  if (displayLoops() > 0) {
    gH1H1 = 0.0; gH2H2 = 0.0;
    gLL = - sumProd(her, let) - 3.0 * sumProd(hdr, ldt); 
    gEE = - matrixify(let, her);
    gQQ = - sumProd(hdr.transpose(), ldt.transpose());
    gDD = - 2.0 * matrixify(ldt, hdr)
      - 2.0 * sumProd(hur, lut).transpose();
    gUU = - matrixify(lut, hur);
    gH1L = (3.0 * (ldt.transpose() * hd).trace(2) + 
	    (let.transpose() * he).trace(2)) * (1.0); // DEBUG -1

    const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
    gEE  = gEE * oneO16Pisq;    
    gLL  = gLL * oneO16Pisq;
    gQQ  = gQQ * oneO16Pisq;
    gUU  = gUU * oneO16Pisq;
    gDD  = gDD * oneO16Pisq;
    gH1L = gH1L * oneO16Pisq;
  }
}

void RpvSoftsusy::check(const DoubleVector & sneutrinoVevs) const {

  double snuSq, v1, v2;
  usefulVevs(displayHvev(), sneutrinoVevs, snuSq, v1, v2);
  double tb = displayTanb();
  double beta = atan(tb);
  double vSM = 246.22; // DEBUG / root2
  double mz = displayMz();

  // These should be zero if you've truly solved the equations
  double equationOne = displayMh1Squared() * v1 +
    displayMh1lSquared().dot(sneutrinoVevs) +
    sqr(displaySusyMu()) * v1 + displaySusyMu() *
    displayKappa().dot(sneutrinoVevs)
    - displayM3Squared() * v2  + 0.5 * sqr(mz) * cos(2.0 * beta) * v1 +
    v1 * sqr(sin(beta)) * sneutrinoVevs.dot(sneutrinoVevs) * sqr(mz) /
    sqr(vSM);

  double equationTwo = v2 * 
    (displayMh2Squared() + sqr(displaySusyMu()) +
     displayKappa().dot(displayKappa())) - 
    displayM3Squared() * v1 - displayDr().dot(sneutrinoVevs) +
    0.5 * sqr(mz) / sqr(vSM) * (sqr(v2) - sqr(v1) - snuSq) * v2;

  if (PRINTOUT) 
    cout << "    check1=" << equationOne << " " << equationTwo 
	 << endl;
}

// Returns derivatives in form of dr: CHECKED 25/5/02
RpvSoftsusy RpvSoftsusy::beta2() const {

  // Brevity
  const Tensor &le1 = displayLambda(LE), &ld1 = displayLambda(LD),
    &lu1 = displayLambda(LU); 
  const DoubleMatrix &u1 = displayYukawaMatrix(YU), 
    &d1 = displayYukawaMatrix(YD),
    &e1 = displayYukawaMatrix(YE); 
  const DoubleMatrix &hu = displayTrilinear(UA), &hd = displayTrilinear(DA),
    &he = displayTrilinear(EA);
  const Tensor &her = displayHr(LE), &hur = displayHr(LU), hdr = displayHr(LD);

  RpvSoftsusy dR;
  dR.setSoftPars(SoftParsMssm::beta2());

  // Supersymmetric coupling renormalisation first:
  // Extract RPC anomalous dimension part
  DoubleVector dg(3);
  sBrevity a;
  DoubleMatrix gLL(3, 3), gEE(3, 3), gQQ(3, 3), gDD(3, 3),
    gUU(3, 3); 
  double gH1H1, gH2H2;
  gH1H1 = 0.0; gH2H2 = 0.0;
  MssmSusy::anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, a);

  // The RPV bits of anomalous dimensions
  DoubleMatrix gLLrpv(3, 3), gEErpv(3, 3), gQQrpv(3, 3), gDDrpv(3, 3),
    gUUrpv(3, 3);
  DoubleVector gH1lrpv(3);
  double gH1H1rpv, gH2H2rpv;
  gH1H1rpv = 0.0; gH2H2rpv = 0.0;
  rpvAnomalousDimension(gEErpv, gLLrpv, gQQrpv, gUUrpv, gDDrpv, gH1H1rpv,
			gH2H2rpv, gH1lrpv);

  // Add the two together to get TOTAL anomalous dimensions
  DoubleMatrix gLLtot(3, 3), gEEtot(3, 3), gQQtot(3, 3), gDDtot(3, 3),
    gUUtot(3, 3); 
  double gH1H1tot, gH2H2tot;
  gLLtot = gLL + gLLrpv; gEEtot = gEE + gEErpv;
  gQQtot = gQQ + gQQrpv; gUUtot = gUU + gUUrpv;
  gDDtot = gDD + gDDrpv; 
  gH1H1tot = gH1H1 + gH1H1rpv;  gH2H2tot = gH2H2 + gH2H2rpv;  
  //    cout << "gDDtot" << gDDtot << "gQQtot" << gQQtot << "gLLtot" << gLLtot;
      /*  cout << "gUggUtot" << gUUtot << "gEEtot" << gEEtot;
  cout << "gH1H1tot " << gH1H1tot << " gH2H2tot " << gH2H2tot << endl;
  cout << "gH1lrpv " << gH1lrpv;*/
  
   DoubleMatrix g1LLrpv(3, 3), g1EErpv(3, 3), g1QQrpv(3, 3), 
    g1DDrpv(3, 3), g1UUrpv(3, 3);
   DoubleVector g1H1lrpv(3);
   double g1H1H1rpv, g1H2H2rpv;
  g1H1H1rpv = 0.0; g1H2H2rpv = 0.0;
  rpvAnomalousDeriv(g1EErpv, g1LLrpv, g1QQrpv, g1UUrpv, g1DDrpv, g1H1H1rpv,
		    g1H2H2rpv, g1H1lrpv);

  // The RPC bits of derivatives of anomalous dimensions
  DoubleMatrix g1LL(3, 3), g1EE(3, 3), g1QQ(3, 3), 
    g1DD(3, 3), g1UU(3, 3);

  double g1H1H1, g1H2H2;
  g1H1H1 = 0.0; g1H2H2 = 0.0;
  SoftParsMssm::anomalousDeriv(g1EE, g1LL, g1QQ, g1UU, g1DD, g1H1H1,
				   g1H2H2);

  DoubleMatrix g1LLtot(3, 3), g1EEtot(3, 3), g1QQtot(3, 3), 
    g1DDtot(3, 3), g1UUtot(3, 3); 
  double g1H1H1tot, g1H2H2tot;
  g1LLtot = g1LL + g1LLrpv; g1EEtot = g1EE + g1EErpv;
  g1QQtot = g1QQ + g1QQrpv; g1UUtot = g1UU + g1UUrpv;
  g1DDtot = g1DD + g1DDrpv; 
  g1H1H1tot = g1H1H1 + g1H1H1rpv;  g1H2H2tot = g1H2H2 + g1H2H2rpv;  

  // Define the y tilde parameters
  DoubleMatrix yeTildeRpv(3, 3), ydTildeRpv(3, 3), yuTildeRpv(3, 3);
  Tensor leTilde, ldTilde, luTilde;
  rpvyTildes(yeTildeRpv, ydTildeRpv, leTilde, ldTilde, luTilde);
  DoubleMatrix yeTilde(3, 3), ydTilde(3, 3), yuTilde(3, 3);
  DoubleMatrix yeTildeRpc(3, 3), ydTildeRpc(3, 3), yuTildeRpc(3, 3);
  SoftParsMssm::yTildes(yuTildeRpc, ydTildeRpc, yeTildeRpc);
  yeTilde = (yeTildeRpv + yeTildeRpc);
  ydTilde = (ydTildeRpv + ydTildeRpc);

  //cout << "LDTILDE" << ldTilde << "LETILDE" << leTilde << "yeTilde" << yeTilde << "ydTilde" << ydTilde;

  // Differentials of SUSY RPV couplings
  Tensor dLE, dLD, dLU;
  DoubleVector dkappa(3), kappa(displayKappa());
  // Differentials of SUSY RPC couplings
  DoubleMatrix du(3, 3), dd(3, 3), de(3, 3);

  // RGEs of RPV SUSY parameters: checked 29/3/99
  // Later: add onto the already present RPC derivatives
  du = u1 * (gUUrpv + gH2H2rpv) + gQQrpv * u1;
  dd = d1 * (gDDrpv + gH1H1rpv) + gQQrpv * d1 - ld1.dotProd(gH1lrpv, 2);
  de = e1 * (gEErpv + gH1H1rpv) + gLLrpv * e1 - le1.dotProd(gH1lrpv, 2);

  dLU = lu1 * gDDtot + gDDtot.transpose() * lu1 + lu1.product(gUUtot);
  dLD = ld1.product(gDDtot) + ld1 * gQQtot.transpose() + gLLtot * ld1 -
    outerProduct(gH1lrpv, d1, 2);

  dLE = le1.product(gEEtot.transpose()) + le1 * gLLtot.transpose() 
    + outerProduct(gH1lrpv, e1.transpose(), 3)
    + gLLtot * le1 - outerProduct(gH1lrpv, e1, 2);
  double dmu = displaySusyMu() * (gH1H1rpv + gH2H2rpv) + kappa.dot(gH1lrpv);
  dkappa = gH2H2tot * kappa + kappa * gLLtot + displaySusyMu() *
    gH1lrpv;
  // Running tan beta
  double dt = displayTanb() * (gH1H1rpv - gH2H2rpv);

  // SUSY RPC trilinears: RPV contribution
  DoubleMatrix dHerpv(3, 3), dHdrpv(3, 3), dHurpv(3, 3);
  dHerpv = gLLrpv * he + gH1H1rpv * he + (her.dotProd(gH1lrpv, 3)).transpose() 
    + he * gEErpv - 2.0 * g1LLrpv * e1 - 2.0 * g1H1H1rpv * e1 
    - 2.0 * (le1.dotProd(g1H1lrpv, 3)).transpose() 
    - 2.0 * e1 * g1EErpv;

  dHdrpv = gQQrpv * hd + gH1H1rpv * hd - hdr.dotProd(gH1lrpv, 2)
    + hd * gDDrpv - 2.0 * g1QQrpv * d1 - 2.0 * g1H1H1rpv * d1 
    + 2.0 * ld1.dotProd(g1H1lrpv, 2)
    - 2.0 * (d1 * g1DDrpv);

  dHurpv = gQQrpv * hu + gH2H2rpv * hu 
    + hu * gUUrpv - 2.0 * g1QQrpv * u1 - 2.0 * g1H2H2rpv * u1 
    - 2.0 * (u1 * g1UUrpv);

  // derivatives of RPV parts of soft scalar mass squared parameters
  DoubleMatrix dmq(3, 3), dmU(3, 3), dmd(3, 3), dme(3, 3), dml(3, 3);
  double dmH1H1;
  DoubleVector dmH1l(3);
  dme = 2.0 * matrixify(her.transpose(), her) 
    + 2.0 * (e1.transpose() * yeTildeRpv)
    + matrixify(le1.transpose(), leTilde) 
    + 2.0 * (yeTildeRpv.transpose() * e1) 
    + matrixify(leTilde.transpose(), le1);

  dml = 2.0 * sumProd(her, her.transpose()) + 
    6.0 * sumProd(hdr, hdr.transpose()) +
    yeTildeRpv * e1.transpose() + sumProd(leTilde, le1.transpose()) +
    3.0 * sumProd(ldTilde, ld1.transpose()) +
    e1 * yeTildeRpv.transpose() + sumProd(le1, leTilde.transpose()) +
    3.0 * sumProd(ld1, ldTilde.transpose());

  dmq = 2.0 * sumProd(hdr.transpose(), hdr) + ydTildeRpv * d1.transpose() +
    sumProd(ldTilde.transpose(), ld1) + sumProd(ld1.transpose(), ldTilde) +
    d1 * ydTildeRpv.transpose();
    
  dmd =  4.0 * matrixify(hdr.transpose(), hdr) +
    4.0 * sumProd(hur, hur.transpose()) +
    2.0 * matrixify(ld1.transpose(), ldTilde) +
    2.0 * sumProd(lu1, luTilde.transpose()) +
    2.0 * matrixify(ldTilde.transpose(), ld1) +
    2.0 * sumProd(luTilde, lu1.transpose()) +
    2.0 * dt * ydTilde + 
    2.0 * ydTilde.transpose() * d1;

  dmd = 4.0 * matrixify(hdr.transpose(), hdr) + 
    4.0 * sumProd(hur, hur.transpose()) +
    2.0 * d1.transpose() * ydTildeRpv + 
    2.0 * matrixify(ld1.transpose(), ldTilde) +
    2.0 * sumProd(lu1, luTilde.transpose()) +
    2.0 * ydTildeRpv.transpose() * d1 + 
    2.0 * matrixify(ldTilde.transpose(), ld1) +
    2.0 * sumProd(luTilde, lu1.transpose()); 

  dmU = 2.0 * matrixify(hur.transpose(), hur) 
    + matrixify(lu1.transpose(), luTilde) 
    + matrixify(luTilde.transpose(), lu1);

  dmH1H1 = 3.0 * (ydTildeRpv * d1.transpose()).trace() +
    (yeTildeRpv * e1.transpose()).trace() +
    3.0 * (d1 * ydTildeRpv.transpose()).trace() +
    (e1 * yeTildeRpv.transpose()).trace();

  dmH1l = -6.0 * (hdr * hd).trace(2) - 2.0 * (her * he).trace(2)
    -3.0 * (ld1 * ydTilde).trace(2) - (le1 * yeTilde).trace(2) 
    -3.0 * (ldTilde * d1).trace(2) - (leTilde * e1).trace(2);

  //Higgs soft parameters
  double dMh3Squared;
  dMh3Squared = displayM3Squared() * (gH1H1rpv + gH2H2rpv) +
    displayDr().dot(gH1lrpv) - 2.0 * displaySusyMu() * 
    (g1H1H1rpv + g1H2H2rpv) - 2.0 * displayKappa().dot(g1H1lrpv);

  DoubleVector dD(3);
  dD = gLLtot * displayDr() + gH2H2tot * displayDr() + displayM3Squared() *
    gH1lrpv - 2.0 * (g1LLtot * displayKappa() + g1H2H2tot * displayKappa()) -
    2.0 * displaySusyMu() * g1H1lrpv;

  // Finally, we have trilinear RPV soft parameters:
  Tensor dher, dhdr, dhur;
  dher = gLLtot * her - outerProduct(gH1lrpv, he, 2) + her *
    gLLtot.transpose() + outerProduct(gH1lrpv, he, 3).swap(3) +
    her.raise(gEEtot.transpose()) - 2.0 * g1LLtot * le1 + 2.0 * 
    outerProduct(g1H1lrpv, e1, 2) - 2.0 * le1 * g1LLtot.transpose() - 2.0 *
    outerProduct(g1H1lrpv, e1.transpose(), 3) - 2.0 *
    le1.raise(g1EEtot.transpose()); 

  dhdr = gLLtot * hdr - outerProduct(gH1lrpv, hd, 2) + hdr *
    gQQtot.transpose() + 
    hdr.raise(gDDtot.transpose()) - 2.0 * g1LLtot * ld1 + 2.0 * 
    outerProduct(g1H1lrpv, d1, 2) - 2.0 * ld1 * g1QQtot.transpose() - 2.0 *
    ld1.raise(g1DDtot.transpose()); 

  dhur = hur.raise(gUUtot.transpose()) + gDDtot.transpose() * hur +
    hur * gDDtot - 2.0 * lu1.raise(g1UUtot.transpose()) -
    2.0 * g1DDtot.transpose() * lu1 - 2.0 * lu1 * g1DDtot;

  // Normalise by the loop factor
  static const double oneO16pisq = 1.0 / (16.0 * sqr(PI));
  dme = dme * oneO16pisq; dml = dml * oneO16pisq;
  dmq = dmq * oneO16pisq; dmd = dmd * oneO16pisq;
  dmU = dmU * oneO16pisq; dmH1H1 = dmH1H1 * oneO16pisq;
  dmH1l = dmH1l * oneO16pisq;

  // Set the derivatives
  dR.setTanb(dt + dR.displayTanb());
  dR.setSusyMu(dmu + dR.displaySusyMu());
  dR.setYukawaMatrix(YE, de + dR.displayYukawaMatrix(YE));
  dR.setYukawaMatrix(YD, dd + dR.displayYukawaMatrix(YD));
  dR.setYukawaMatrix(YU, du + dR.displayYukawaMatrix(YU));
  dR.setLambda(LU, dLU);
  dR.setLambda(LD, dLD);
  dR.setLambda(LE, dLE);
  dR.setKappa(dkappa);

  dR.setTrilinearMatrix(UA, dR.displayTrilinear(UA) + dHurpv);
  dR.setTrilinearMatrix(DA, dR.displayTrilinear(DA) + dHdrpv);
  dR.setTrilinearMatrix(EA, dR.displayTrilinear(EA) + dHerpv);

  dR.setSoftMassMatrix(mEr, dR.displaySoftMassSquared(mEr) + dme);
  dR.setSoftMassMatrix(mQl, dR.displaySoftMassSquared(mQl) + dmq);
  dR.setSoftMassMatrix(mUr, dR.displaySoftMassSquared(mUr) + dmU);
  dR.setSoftMassMatrix(mDr, dR.displaySoftMassSquared(mDr) + dmd);
  dR.setSoftMassMatrix(mLl, dR.displaySoftMassSquared(mLl) + dml);
  dR.setMh1Squared(dR.displayMh1Squared() + dmH1H1);
  dR.setMh1lSquared(dmH1l);  

  dR.setM3Squared(dR.displayM3Squared() + dMh3Squared);
  dR.setDr(dD);

  dR.setHr(LU, dhur);
  dR.setHr(LD, dhdr);
  dR.setHr(LE, dher);

  return dR;
}

// Returns the RPV parts (only) of the ytildes: CHECKED 23/05/02
void RpvSoftsusy::rpvyTildes(DoubleMatrix & ye, DoubleMatrix & yd, Tensor &
			     letilde, Tensor & ldtilde, Tensor & lutilde)
  const {
  // For brevity
    const Tensor &le = displayLambda(LE), &ld = displayLambda(LD), 
    &lu = displayLambda(LU);
  const DoubleVector &mh1lsq = displayMh1lSquared();
  const DoubleMatrix &mq = displaySoftMassSquared(mQl),
    &mu = displaySoftMassSquared(mUr), &md = displaySoftMassSquared(mDr),
    &ml = displaySoftMassSquared(mLl), &me = displaySoftMassSquared(mEr); 
  const DoubleMatrix &e1 = displayYukawaMatrix(YE), 
  &d1 = displayYukawaMatrix(YD);

  ye = -(le.dotProd(mh1lsq, 3)).transpose(); 
  yd = (ld.dotProd(mh1lsq, 2));
  letilde = ml * le + outerProduct(mh1lsq, e1, 2) + le * ml.transpose() -
    outerProduct(mh1lsq, e1.transpose(), 3) + le.product(me);

  ldtilde = ml * ld + outerProduct(mh1lsq, d1, 2) + ld * mq.transpose() +
    ld.product(md);
  lutilde = lu.product(mu) + md.transpose() * lu + lu * md;
} 

static double mz, sinthDRbarMS;

// mt will only become relevant once 1-loop corrections are incorporated.
// Currently, this is a tree-level algorithm
void RpvSoftsusy::rewsb(int sgnMu, double mt) {
  // Initial values to be iterated
  MssmSoftsusy::rewsbTreeLevel(sgnMu);
  double mu = displaySusyMu(), m3sq = displayM3Squared();

  sinthDRbarMS = calcSinthdrbar();
  calcTadpole1Ms1loop(mt, sinthDRbarMS);  
  calcTadpole2Ms1loop(mt, sinthDRbarMS); 
  double pizztMS = piZZT(displayMz(), displayMu());
  setHvev(getVev(pizztMS)); 
  mz = sqrt(sqr(displayMz()) + piZZT(displayMz(), displayMu()));

  int maxTries = 400; 
  double tol = maximum(TOLERANCE * 1.0e-4, 1.0e-14);

  DoubleVector sneutrinoVevs(3);

  int numTries = 0;
  iterateRewsb(mu, m3sq, sneutrinoVevs, sgnMu, numTries, maxTries, tol, mt);

  if (PRINTOUT > 2) cout << "VEVs" << sneutrinoVevs;
  //  check(sneutrinoVevs);
  // Now do the rotation to get rid of the VEVs
  //  rotateAwayVevs(sneutrinoVevs);
  setSneutrinoVevs(sneutrinoVevs);
}

// Put numTries to zero before using this iterative routine. Does the
// minimisation, returning the VEVs.
void RpvSoftsusy::iterateRewsb(double & mu, double & m3sq, DoubleVector &
			       sneutrinoVevs, int sgnMu, int & numTries, int
			       maxTries, double tol, double mt) {
  static double oldMu, oldM3sq;
  static DoubleVector oldSneutrinoVevs(3);
  oldMu = mu; oldM3sq = m3sq; oldSneutrinoVevs = sneutrinoVevs;

  numTries = numTries + 1;

  // non-convergence
  if (numTries > maxTries) {
    if (PRINTOUT) cout << "REWSB iteration not converged\n";
    flagNoMuConvergence(true);
    return;
  }

  double snuSq, v1, v2;
  if (usefulVevs(displayHvev(), sneutrinoVevs, snuSq, v1, v2)) {
    if (PRINTOUT) cout << "sneutrino VEVs incompatible with MZ, MW\n";
    sneutrinoVevs = DoubleVector(3);
    flagHiggsufb(true);
    return;
  }

  setSusyMu(mu);
  
  // calculate the new one-loop tadpoles with old value of mu
  calcTadpole2Ms1loop(mt, sinthDRbarMS);
  calcTadpole1Ms1loop(mt, sinthDRbarMS);   

  mu = calculateMu(sneutrinoVevs, sgnMu, v1, v2);
  setSusyMu(mu);
  m3sq = calculateM3sq(sneutrinoVevs, snuSq, v1, v2);  
  setM3Squared(m3sq);
  calcDrBarPars();
  sneutrinoVevs = calculateSneutrinoVevs(sneutrinoVevs, tol, snuSq, v1, v2);

  //  check(sneutrinoVevs); DEBUG -- additional check

  DoubleVector sT(5);
  sT(1) = toleranceCheck(mu, oldMu);
  sT(2) = toleranceCheck(m3sq, oldM3sq);
  sT(3) = toleranceCheck(sneutrinoVevs(1), oldSneutrinoVevs(1));
  sT(4) = toleranceCheck(sneutrinoVevs(2), oldSneutrinoVevs(2));
  sT(5) = toleranceCheck(sneutrinoVevs(3), oldSneutrinoVevs(3));
  
  // Fractional difference to last iteration
  double difference = sT.max();

  if (PRINTOUT > 2) 
    cout << numTries << ". D=" << difference << " mu=" << mu 
	 << " m3sq=" << m3sq << "\n    v1=" << sneutrinoVevs(1) << " v2=" 
	 << sneutrinoVevs(2) << " v3=" << sneutrinoVevs(3) << endl; 

  if (difference < tol) return;

  iterateRewsb(mu, m3sq, sneutrinoVevs, sgnMu, numTries, maxTries, tol, mt);
}

// Returns value of mu consistent with sneutrino vevs given
// (which should be a length-3 vector)
double RpvSoftsusy::calculateMu(const DoubleVector & sneutrinoVevs, int sgnMu,
				double v1, double v2) {
  //const double vSM = 246.22 / root2;
  double vSM = displayHvev();
  
  const DoubleVector & kappa = displayKappa();

  // It's a quadratic, so define a,b,c from a x^2 + b x + c = 0:  
  double a = sqr(displayTanb()) - 1.0;

  double b = - displayKappa().dot(sneutrinoVevs) / v1;

  double c = sqr(displayTanb()) * 
    (displayMh2Squared() - displayTadpole2Ms() + 
     kappa.dot(kappa) - displayDr().dot(sneutrinoVevs) / v2 -
     sqr(mz) / sqr(vSM) * sneutrinoVevs.dot(sneutrinoVevs)) - // alter
    (displayMh1Squared() - displayTadpole1Ms() +
     displayMh1lSquared().dot(sneutrinoVevs) / v1) +
    0.5 * sqr(mz) * (sqr(displayTanb()) - 1.0);

  double d = sqr(b) - 4.0 * a * c;

  double mu;
  if (sgnMu > 0)
    mu = (- b + sqrt(fabs(d))) / (2.0 * a);
  else
    mu = (- b - sqrt(fabs(d))) / (2.0 * a);

  if (d < 0) { 
    if (PRINTOUT) cout << " musq<0 "; 
    mu = - mu;
    flagMusqwrongsign(true);
  }

  return mu;
}

// Returns some functions of VEVs, gives 0 if there's a problem with the
// sneutrino VEVs (if they are incompatible with the W and Z masses)
int RpvSoftsusy::usefulVevs(double vSM, 
			    const DoubleVector & sneutrinoVevs, double &
			       snuSq, double & v1, double & v2) const {
  double tb = displayTanb();
  double beta = atan(tb);
  snuSq = sneutrinoVevs.dot(sneutrinoVevs);
  //double vSM = 246.22 / root2;
  double vHiggs = vSM * sqrt(1.0 - snuSq / sqr(vSM));
  v2 = vHiggs * sin(beta);
  v1 = vHiggs * cos(beta);

  if (snuSq < sqr(vSM)) return 0;
  else return 1;
}


// Returns value of m3sq consistent with sneutrino vevs given
// (which should be a length-3 vector)
double RpvSoftsusy::calculateM3sq(const DoubleVector & sneutrinoVevs,
				double snuSq, double v1, double v2) {

  double tb = displayTanb();
  double beta = atan(tb);

  //  if (PRINTOUT > 2) cout << endl << "vd=" << v1 << " vu=" << v2;

  double sqm3 = 0.5 * sin(2.0 * beta) *
    (displayMh1Squared() - displayTadpole1Ms() + 
     displayMh2Squared() - displayTadpole2Ms() + 2.0 * sqr(displaySusyMu()) +
     displayKappa().dot(displayKappa()) - displayDr().dot(sneutrinoVevs) / v2 +
     displayMh1lSquared().dot(sneutrinoVevs) / v1 + displaySusyMu() * 
     displayKappa().dot(sneutrinoVevs) / v1);

  return sqm3;
}

// Input a set of values for sneutrino VEVs and it returns a more accurate set
DoubleVector RpvSoftsusy::calculateSneutrinoVevs
(const DoubleVector & sneutrinoVevs, double tol, double snuSq, double v1,
 double v2) { 
  
  double tb = displayTanb();
  double beta = atan(tb);
  //  const double vSM = 246.22 / root2;
  double vSM = displayHvev();

  DoubleVector n(3), x0(3);
  DoubleMatrix m(3, 3), mInverse(3, 3);

  // Identity matrix
  DoubleMatrix i(3, 3);
  i(1, 1) = i(2, 2) = i(3, 3) = 1.0;

  n = displayDr() * v2 - displaySusyMu() * v1 * displayKappa() 
    - displayMh1lSquared() * v1;
  m = displaySoftMassSquared(mLl).transpose()
    + (0.5 * sqr(mz) * cos(2.0 * beta) 
       + sqr(sin(beta)) * sqr(mz) / sqr(vSM) * snuSq) * i +
    outerProduct(displayKappa(), displayKappa());


  // Must work out how to invert some general matrix m from numerical recipes
  mInverse = m.inverse();

  /*
  cout << "DEBUG:" << mInverse << " " << n << " vd=" << v1 << " vu=" << v2 
       << endl;
  cout << "terms:" << displayDr() * v2 << " " << - displaySusyMu() * v1 * displayKappa()  << " " << 
    - displayMh1lSquared() * v1 << endl << "hvev" << displayHvev() << endl;
  */

  // If this one ever gets flagged, you will have to iterate this bit of the
  // calculation 
  if (fabs(x0.max()) > fabs(n.max()) * tol) // Development flag
    cout << "WARNING IN RPVSOFT: calculation inaccurate\n";

  return mInverse * n;
}

// Transforms all parameters into a new basis where sneutrino vevs are zero:
// needs checking! With our defn of tan beta, it's invariant under this!
// This MUST be checked - it's not working at the moment anyway....
void RpvSoftsusy::rotateAwayVevs(DoubleVector & snVevs) {

  double snuSq, v1, v2;
  usefulVevs(displayHvev(), snVevs, snuSq, v1, v2);

  /*  if (PRINTOUT > 2) cout << endl << " vu=" << v2 << " vd=" << v1 
			 << " v2sq + v1sq + visq=" 
			 << sqrt(sqr(v1) + sqr(v2) + snuSq);*/

  // Form SO(4) rotation matrix
  double c1, s1, c2, s2, c3, s3;
  if (snVevs(2) == 0. && snVevs(3) == 0.) {
    c1 = 1.0; s1 = 0.0;
  } 
  else {
    c1 = snVevs(2) / sqrt(sqr(snVevs(2)) + sqr(snVevs(3)));
    s1 = snVevs(3) / sqrt(sqr(snVevs(2)) + sqr(snVevs(3)));
  }

  double vdv = sqrt(snVevs.dot(snVevs));
  if (vdv == 0.) {
    c2 = 1.; s2 = 0.;
  }
  else {
    c2 = snVevs(1) / vdv; s2 = sqrt(sqr(snVevs(2)) + sqr(snVevs(3))) / vdv;
  }

  double vhl = sqrt(sqr(vdv) + sqr(v1));
  c3 = v1 / vhl; s3 = vdv / vhl;

  static DoubleMatrix O(4, 4);
  O(1, 1) = c3;            O(1, 2) = -s3;

  O(2, 1) = c2 * s3;       O(2, 2) = c2 * c3; O(2, 3) = -s2;

  O(3, 1) = c1 * s2 * s3;  O(3, 2) = c1 * s2 * c3; 
  O(3, 3) = c1 * c2;       O(3, 4) = -s1;

  O(4, 1) = s1 * s2 * s3;  O(4, 2) = s1 * s2 * c3; 
  O(4, 3) = s1 * c2;       O(4, 4) = c1;

  DoubleVector vevs(4);
  vevs(1) = v1;

  int i; for (i=1; i<=3; i++) vevs(i+1) = snVevs(i);
  vevs = O.transpose() * vevs; // rotate to new basis. Note vd is NOT v1!

  // New mu and kappa terms
  static DoubleVector muNew(4), bNew(4);
  muNew(1) = displaySusyMu();
  bNew(1)  = displayM3Squared();
  for (i=1; i<=3; i++) {
    muNew(i+1) = displayKappa()(i);
    bNew(i+1)  = displayDr()(i);
  }
  muNew = O * muNew; // Rotation
  bNew  = O * bNew;

  setSusyMu(muNew(1)); // set new quantities
  for (i=1; i<=3; i++) {
    setKappa(i, muNew(i+1));
    setD(i, bNew(i+1));
  }
  setM3Squared(bNew(1));

  // New down Yukawas and lambda_primes
  DoubleMatrix newYd(3, 3);
  DoubleMatrix newYe(3, 3);
  Tensor newLd, newLe;
  newYd = O(1, 1) * displayYukawaMatrix(YD); // Rotation

  // Do their Rotations:
  int j, k, l, m;
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      for (k=1; k<=3; k++) {
	newYe(i, j) = newYe(i, j) - O(1, i+1) * O(k+1, 1) *
	  displayYukawaMatrix(YE)(k, j) +
	  O(k+1, i+1) * O(1, 1) * displayYukawaMatrix(YE)(k, j);
	newYd(i, j) = newYd(i, j) + O(k+1, 1) * displayLambda(LD)(k, i, j);
	newLd(k, i, j) = O(1, k+1) * displayYukawaMatrix(YD)(i, j);
	for (l=1; l<=3; l++) {
	  newLd(k, i, j) = newLd(k, i, j) + O(l+1, k+1) *
	    displayLambda(LD)(l, i, j);
	  newYe(i, j) = newYe(i, j) + displayLambda(LE)(j, k, l) * O(l+1, 1) *
	    O(k+1, i+1);
	  newLe(k, i, j) = newLe(k, i, j) + displayYukawaMatrix(YE)(l, k) *
	    (O(l+1, i+1) * O(1, j+1) - O(1, i+1) * O(l+1, j+1));
	  for (m=1; m<=3; m++) 
	    newLe(k, i, j) = newLe(k, i, j) + 
	      displayLambda(LE)(k, l, m) * O(l+1, i+1) * O(m+1, j+1);
	}
      }
    }
  setYukawaMatrix(YE, newYe); // set Yukawas again
  setLambda(LE, newLe);
  setYukawaMatrix(YD, newYd);
  setLambda(LD, newLd);

  // New msquared matrix
  DoubleMatrix msq(4, 4);
  msq(1, 1) = displayMh1Squared();
  for (i=1; i<=3; i++) {
    msq(i+1, 1) = displayMh1lSquared()(i);
    msq(1, i+1) = displayMh1lSquared()(i);
    for(j=1; j<=3; j++)
      msq(i+1, j+1) = displaySoftMassSquared(mLl)(i, j);
  }

  msq = O * msq * O.transpose(); // rotations

  setMh1Squared(msq(1, 1)); // set rotated variables
  for (i=1; i<=3; i++) {
    setMh1lSquared(i, msq(1, i+1));
    for(j=1; j<=3; j++) setSoftMassElement(mLl, i, j, msq(i+1, j+1));
  }
  
  // Must also rotate the trilinear terms
}

// You must set SUSY RPV parameters before this
void RpvSoftsusy::standardSugra(double m0,  double m12, double a0) {
  SoftParsMssm::standardSugra(m0, m12, a0);
  setHr(LU, a0 * displayLambda(LU));
  setHr(LD, a0 * displayLambda(LD));
  setHr(LE, a0 * displayLambda(LE));
  DoubleVector n(3); setMh1lSquared(n); // sets them to be zero
  setDr(displayKappa() * displaySusyMu());
}

// Sums up neutrino masses valid for cosmological bound
double neutrinoSum(const RpvSoftsusy & r) {
  /* DoubleMatrix mn(r.neutralinoMassMatrix()), mix(7, 7);
     DoubleVector w(7);
     mn.diagonaliseSym(mix, w);*/

  // in Sakis' notation
  double gp = r.displayGaugeCoupling(1) * sqrt(0.6);
  double g2 = r.displayGaugeCoupling(2);
  double M2 = r.displayGaugino(2);
  double M1 = r.displayGaugino(1);
  double smu = r.displaySusyMu();
  double beta = atan(r.displayTanb()); 
  double vu = r.displayHvev() * sin(beta) / sqrt(2.0);
  double vd = r.displayHvev() * cos(beta) / sqrt(2.0);
  DoubleVector lambda = r.displaySneutrinoVevs() * (1.0 / sqrt(2.0)) - 
    vd / smu * r.displayKappa();
  /*  cout << "lambda=" << lambda << " vevs=" << r.displaySneutrinoVevs() * (1.0/sqrt(2.0)) 
       << " vd=" << vd << " vu=" << vu << " smu=" << smu << " kap=" << r.displayKappa() << " br=" << (M1 * sqr(g2) + M2 * sqr(gp)) << " num=" << smu * (M1 * sqr(g2) + M2 * sqr(gp)) * 
    (sqr(lambda(3)) + sqr(lambda(2)) + sqr(lambda(1))) << 
    " " << (vu * vd * (M1 * sqr(g2) + M2 * sqr(gp)) - smu * M1 * M2) * 2.0;
  double mnuWrong = fabs(smu * (M1 * sqr(g2) + M2 * sqr(gp)) * 
    (sqr(lambda(3)) + sqr(lambda(2)) + sqr(lambda(1))) / 
    (vu * vd * (M1 * sqr(g2) + M2 * sqr(gp)) - 2.0 * smu * M1 * M2));
  */
  double mnuCorrection = fabs(smu * (M1 * sqr(g2) + M2 * sqr(gp)) * 
    (sqr(lambda(3)) + sqr(lambda(2)) + sqr(lambda(1))) / 
    (vu * vd * (M1 * sqr(g2) + M2 * sqr(gp)) - smu * M1 * M2) * 0.5);

  
  return mnuCorrection;
}

// Under construction
void RpvSoftsusy::neutralinos() const {
  DoubleMatrix mn(neutralinoMassMatrix()), mix(7, 7);
  DoubleVector w(7);
  mn.diagonaliseSym(mix, w);
}

// defines 7 by 7 RPV mass matrix
DoubleMatrix RpvSoftsusy::neutralinoMassMatrix() const {
  double v1, v2, snuSq;
  DoubleVector sneutrinoVevs(snuVevs.display());
  if (usefulVevs(displayHvev(), sneutrinoVevs, snuSq, v1, v2)) 
    throw "VEVS in RpvSoftsusy::neutralinoMassMatrix";

  DoubleMatrix mn(7, 7);

  DoubleVector vhiggs1(4), mupar(4);
  vhiggs1(1) = v1; mupar(1) = displaySusyMu();
  int i; for (i=1; i<=3; i++) {
    vhiggs1(i+1) = snuVevs.display(i);
    mupar(i+1)   = displayKappa()(i);
  }

  double gp = sqrt(0.6) * displayGaugeCoupling(1),
    g2 = displayGaugeCoupling(2);

  /* old and possibly wrong 
  for (i=1; i<=3; i++) {
  mn(i, 4) = - gp / root2 * sneutrinoVevs(i);
  mn(i, 5) = g2 / root2 * sneutrinoVevs(i);
  mn(i, 7) = - displayKappa()(i);
  }
  mn(4, 4) = displayGaugino()(1);
  mn(5, 5) = displayGaugino()(2);
  mn(4, 6) = -gp / root2 * v1;
  mn(4, 7) = gp / root2 * v2;
  mn(5, 6) = g2 / root2 * v1;
  mn(5, 7) = -g2 / root2 * v2;
  mn(6, 7) = -displaySusyMu();
  */

  // corrected
  for (i=1; i<=3; i++) {
  mn(i, 4) = - gp * sneutrinoVevs(i) * 0.5;
  mn(i, 5) = g2 * sneutrinoVevs(i) * 0.5;
  mn(i, 7) = - displayKappa()(i);
  }
  mn(4, 4) = displayGaugino()(1);
  mn(5, 5) = displayGaugino()(2);
  mn(4, 6) = -gp * v1 * 0.5;
  mn(4, 7) = gp * v2 * 0.5;
  mn(5, 6) = g2 * v1 * 0.5;
  mn(5, 7) = -g2 * v2 * 0.5;
  mn(6, 7) = -displaySusyMu();

  mn.symmetrise();

  return mn;
}

void RpvSoftsusy::methodBoundaryCondition(const DoubleVector & 
					  inputParameters){
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);
  
  int k; 
  k = 4;
  RpvSusyPars::set(inputParameters, k);

  standardSugra(m0, m12, a0);
}

void RpvSoftsusy::isawigInterface764(char herwigInputFile [80], 
				      char isajetOutputFile [80],
				      char softOutputFile [80])  
  const {

  fstream softOutput(softOutputFile, ios::out);

  ssrunInterface764Inside(isajetOutputFile, softOutput);

  softOutput << "n" << endl << "y" << endl;

  int i=0, j=0, k=0, l=0; bool test = false;
  for (i = 1; i<=3; i++) 
    for (l = 1; l<=3; l++) {
      if (l == 1) { j = 1; k = 2; }
      if (l == 2) { j = 1; k = 3; }
      if (l == 3) { j = 2; k = 3; }
      if (displayLambda(LE)(i, j, k) != 0.0) test = true;
    }

  if (test) softOutput << "y\n" << displayLambda(LE)(1, 1, 2) << " " 
		 << displayLambda(LE)(2, 1, 2) << " " 
		 << displayLambda(LE)(3, 1, 2) << " " 
		 << displayLambda(LE)(1, 1, 3) << " " 
		 << displayLambda(LE)(2, 1, 3) << " " 
		 << displayLambda(LE)(3, 1, 3) << " " 
		 << displayLambda(LE)(1, 2, 3) << " " 
		 << displayLambda(LE)(2, 2, 3) << " " 
		 << displayLambda(LE)(3, 2, 3) << endl; 
  else softOutput << "n\n";
  
  softOutput << "y\n";
  test = false;
  for (i = 1; i<=3; i++) 
    for (j = 1; j<=3; j++) 
      for (k = 1; k<=3; k++) {
	softOutput << displayLambda(LD)(k, i, j) << " ";
    }
  softOutput << endl;

  for (i = 1; i<=3; i++) 
    for (l = 1; l<=3; l++) {
      if (l == 1) { j = 1; k = 2; }
      if (l == 2) { j = 1; k = 3; }
      if (l == 3) { j = 2; k = 3; }
      if (displayLambda(LU)(i, j, k) != 0.0) test = true;
    }

  if (test) softOutput << "y\n" << displayLambda(LU)(1, 1, 2) << " " 
		 << displayLambda(LU)(1, 1, 3) << " " 
		 << displayLambda(LU)(1, 2, 3) << " " 
		 << displayLambda(LU)(2, 1, 2) << " " 
		 << displayLambda(LU)(2, 1, 3) << " " 
		 << displayLambda(LU)(2, 2, 3) << " " 
		 << displayLambda(LU)(3, 1, 2) << " " 
		 << displayLambda(LU)(3, 1, 3) << " " 
		 << displayLambda(LU)(3, 2, 3) << endl; 
  else softOutput << "n\n";

  softOutput << "'" << herwigInputFile << "'" << endl;

  softOutput.close();
}


// It'll set the important SUSY couplings: supposed to be applied at MZ
// You should set up an iteration here since Yuk's depend on top mass which
// depends on Yuk's etc. 
void RpvSoftsusy::sparticleThresholdCorrections(double tb) {
  mz = displayMz();
  if (displayMu() != mz) {
    ostringstream ii;
    ii << "Called MssmSoftsusy::sparticleThresholdCorrections "
	 << "with scale" << displayMu() << endl;
    throw ii.str();
  }
  
  setTanb(tb);
  
  double alphaDrbar = 
    qedSusythresh(displayDataSet().displayAlpha(ALPHA), displayMu());

  double alphasMZDRbar =
    qcdSusythresh(displayDataSet().displayAlpha(ALPHAS), displayMu());
  
  // Do gauge couplings
  double outrho = 1.0, outsin = 0.48, tol = TOLERANCE * 1.0e-3; 
  int maxTries = 20;
  double pizztMZ = piZZT(mz, displayMu());
  double piwwt0 = piWWT(0., displayMu());
  double piwwtMW = piWWT(displayMw(), displayMu());
  rhohat(outrho, outsin, alphaDrbar, pizztMZ, piwwt0, piwwtMW, tol, maxTries);
  if (displayProblem().noRhoConvergence) 
    outsin = sqrt(1.0 - sqr(displayMw() / mz)); 
  
  double eDR = sqrt(4.0 * PI * alphaDrbar), costhDR = cos(asin(outsin));

  DoubleVector newGauge(3);
  newGauge(1) = eDR / costhDR * sqrt(5.0 / 3.0);
  newGauge(2) = eDR / outsin;
  newGauge(3) = sqrt(4.0 * PI * alphasMZDRbar);

  setGaugeCoupling(1, newGauge(1));
  setGaugeCoupling(2, newGauge(2));
  setGaugeCoupling(3, newGauge(3));

  double vev = displayHvev();

  setMw(sqrt(0.25 * sqr(newGauge(2)) * sqr(vev) - piwwtMW));

  if (MIXING < -1 || MIXING > 2) {
    ostringstream ii;
    ii << "In MssmSoftsusy::sparticleThresholdCorrections(double tb) ";
    ii << "\n MIXING=" << MIXING << " is out of range (-1 -> 2)\n";
    throw ii.str();
  }

  DoubleMatrix mUq(3, 3), mDq(3, 3), mLep(3, 3);
  massFermions(displayDataSet(), mDq, mUq, mLep);
  double mtau = calcRunningMtau();
  mDq(3, 3) = calcRunningMb();
  mUq(3, 3) = calcRunningMt();
  mLep(3, 3) = mtau;

  // 3-family mixed-up Yukawa couplings 
  static DoubleMatrix ckm(3, 3);
  if (ckm(1, 1) == 0.0) {
    ckm(1, 1) = 0.9748; ckm(1, 2) = 0.2229, ckm(1, 3) = 0.0036;
    ckm(2, 1) = -0.2229; ckm(2, 2) = 0.9740, ckm(2, 3) = 0.0412;
    ckm(3, 1) = 0.0057; ckm(3, 2) = -0.0410, ckm(3, 3) = 0.9991;
  }
  else doQuarkMixing(ckm, mDq, mUq); 

  double snuSq, v1, v2;
  if (usefulVevs(vev, displaySneutrinoVevs(), snuSq, v1, v2) == 1) 
    flagHiggsufb(true);

  setYukawaMatrix(YU, mUq * root2 / v2);
  setYukawaMatrix(YD, 
		  (mDq - displayLambda(LD).dotProd(displaySneutrinoVevs(), 1)) 
		  * root2 / v1);
  setYukawaMatrix(YE, mLep * root2 / v1);

  int errSignal; DoubleMatrix e(displayYukawaMatrix(YE));
  iterateChargedLeptons(vev, e, tol, maxTries, errSignal, mtau);

  if (errSignal == 1) flagNoRhoConvergence(true);
}

void RpvSoftsusy::iterateChargedLeptons(double vev, DoubleMatrix & yeOld, 
					double tol, int maxTries, int & err, 
					double mtau) { 
  static int numTries = 0;
  static DoubleMatrix yeNew(3, 3);
  const static DoubleMatrix empty(3, 3);

  if (numTries - 1 > maxTries) { 
    if (PRINTOUT > 2) cout << "iterateChargedLeptons reached maxtries\n"; 
    numTries = 0; yeNew = empty;
    err = 1; return;
  }

  double c = yeNew.compare(yeOld) * 1.0e4;

  if (c < tol) { 
    yeOld = yeNew;
    numTries = 0; yeNew = empty;
    if (PRINTOUT > 2) cout << "sTlep=" << c << " leptons converged\n";
    return; 
  }
  
  numTries = numTries + 1;

  yeOld = yeNew;

  // calculate yeNew;
  DoubleMatrix u(5, 5), v(5, 5), approx(5, 5), wMatrix(5, 5);
  DoubleVector w(5);
  approx = chargedLeptons(vev);
  approx.diagonalise(u, v, w);

  if (PRINTOUT > 2) cout << "sTlep=" << c << " me=" << w(1) << " " << w(2) 
			 << " " << w(3) << endl;

  // calculate another matrix with correct lepton masses
  wMatrix(1, 1) = displayDataSet().displayMass(mElectron);
  wMatrix(2, 2) = displayDataSet().displayMass(mMuon);
  wMatrix(3, 3) = mtau;
  wMatrix(4, 4) = w(4);
  wMatrix(5, 5) = w(5);
  approx = u * wMatrix * v.transpose();

  // Now determine new Yukawa matrix corresponding to this matrix
  DoubleMatrix chMass(3, 3);
  int i, j; for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) 
      chMass(i, j) = approx(i, j);
    
  double v1, v2, snuSq;
  usefulVevs(vev, displaySneutrinoVevs(), snuSq, v1, v2);

  yeNew = (chMass - displayLambda(LE).dotProd(displaySneutrinoVevs(), 2)) 
    * root2 / v1;
  setYukawaMatrix(YE, yeNew);
  iterateChargedLeptons(vev, yeOld, tol, maxTries, err, mtau);
}

// Calculates the charged lepton mass matrix once the leptonic Yukawas have
// been set
DoubleMatrix RpvSoftsusy::chargedLeptons(double vev) {
  double beta = atan(displayTanb());
  double v1, v2, snuSq;
  usefulVevs(vev, displaySneutrinoVevs(), snuSq, v1, v2);

  DoubleMatrix temp(3, 3);
  temp = (displayYukawaMatrix(YE) * v1  +
	  displayLambda(LE).dotProd(displaySneutrinoVevs(), 2)) / root2;

  DoubleMatrix m(5, 5);
  m(5, 5) = displayGaugino()(2);
  m(4, 5) = root2 * displayMw() * cos(beta);
  m(5, 4) = root2 * displayMw() * sin(beta);
  m(4, 4) = displaySusyMu();

  int i, j; for (i=1; i<=3; i++) {
    m(i, 4) = displayKappa()(i);
    m(i, 5) = displayGaugeCoupling(2) * displaySneutrinoVevs()(i) / root2;
    for (j=1; j<=3; j++) 
      m(i, j) = temp(i, j);
  }
  return m;
}

void RpvSoftsusy::higgs(int accuracy, double piwwtMS, double pizztMS) {

  //  MssmSoftsusy::higgs(accuracy, piwwtMS, pizztMS); return; // DEBUG

  double mt = displayDataSet().displayPoleMt(), 
    mb = displayDataSet().displayMass(mBottom), 
    alphas = displayDataSet().displayAlpha(ALPHAS); 
  const static double pisq = sqr(PI);
  double tanb = displayTanb();
  double beta = atan(tanb), sinb = sin(beta), cosb = cos(beta);
  double sinb2 = sqr(sinb), cosb2 = sqr(cosb), mz = displayMzRun(),
    mz2 = sqr(mz), mt2 = sqr(mt), mt4 = sqr(mt2);
  double mu = displaySusyMu();
  
  // tree level: there'll be trouble if B has the opp sign to mu
  double mAsq = displayM3Squared() / (sinb * cosb);
  if (mAsq < 0.0) {
    flagTachyon(true);
    if (PRINTOUT) cout << " mA^2(tree)=" << mAsq << " since m3sq=" <<
		    displayM3Squared() << " and mu=" << mu; 
    mAsq = fabs(mAsq);
  }
  
  double costh = sqr(displayMw() / mz), sinth = 1.0 - costh; 
  DoubleMatrix mH(2, 2);
  mH(1, 1) = mAsq * sinb2 + mz2 * cosb2;
  mH(1, 2) = - sinb * cosb * (mAsq + mz2); 
  mH(2, 2) = mAsq * cosb2 + mz2 * sinb2; 
  mH(2, 1) = mH(1 ,2); 

  /*       one loop contributions:       */
  if (accuracy > 0) {
    double msArg = displaySoftMassSquared(mQl, 3, 3) * 
      displaySoftMassSquared(mUr, 3, 3)
      + mt2 * (displaySoftMassSquared(mQl, 3, 3)
	       + displaySoftMassSquared(mUr, 3, 3)) + sqr(mt2);
    
    if (msArg < 0.0) {
      flagTachyon(true);
      if (PRINTOUT > 2) cout << " msArg in mH calc tachyonic";
    }

    // brevity: been checked
    double ms = pow(fabs(msArg), 0.25), 
      ms2 = sqr(ms), ms4 = sqr(ms2), ms6 = ms4 * ms2, ms8 = sqr(ms4);
    double mtOMs2 = sqr(mt / ms), mtOMs4 = sqr(mtOMs2), mtOMs6 =
      mtOMs2 * mtOMs4, mtOMs8 = sqr(mtOMs4);
    double mzOMt2 = mz2 / mt2, mzOMt4 = sqr(mzOMt2);
    double logMtOMssq = log(mtOMs2);
    double mLR = displaySoftA(UA, 3, 3) - mu / tanb, mLR2 = sqr(mLR), mLR4
      = sqr(mLR2), mLR6 = mLR2 * mLR4, mLR8 = sqr(mLR4);
    const static double lambda = 1.0 / 8.0 - sinth / 3.0 + 4.0 / 9.0 *
      sqr(sinth);
    
    // Following are t,stop corrections from hep-ph/9903404: checked, but
    // note that mu and B have opposite signs to my notation from that paper
    DoubleMatrix sigma(2, 2);
    sigma(1, 1) = GMU * root2 * sqr(mz2) * lambda * cosb2 * logMtOMssq 
      / pisq;
    sigma(1, 2) = -GMU * root2 * mz2 / (pisq * tanb) * 
      (-3.0 / 8.0 * mt2 + mz2 * lambda * sinb2) * logMtOMssq;
    sigma(2, 2) = GMU * root2 / (8.0 * pisq * sinb2) * mt4 * 
      (-2.0 * mzOMt2 + 1.1 * mzOMt4 +
       (12.0 - 6.0 * mzOMt2 * sinb2 + 8.0 * mzOMt4 * lambda *
	sqr(sinb2)) * logMtOMssq + mLR2 / ms2 * 
       (-12.0 + 4.0 * mzOMt2 + 6.0 * mtOMs2) + mLR4 / ms4 * 
       (1.0 - 4.0 * mtOMs2 + 3.0 * mtOMs4) + mLR6 / ms6 * 
       (0.6 * mtOMs2 - 2.4 * mtOMs4 + 2.0 * mtOMs6) + mLR8 / ms8 * 
       (3.0 / 7.0 * mtOMs4 - 12.0 / 7.0 * mtOMs6 + 1.5 * mtOMs8)); 
    
    if (accuracy > 1) {
      // Following are non-t/stop corrections from hep-ph/9609331: checked
      double Pb, Pu, Pd, Pl, Pg, Pgp, P2h, P2hp, P1h, Pf;
      double sinth2 = sqr(sinth);
      Pb = 3.0 * (1.0 - 4.0 * sinth / 3.0 + 8.0 / 9.0 * sinth2);
      Pu = 3.0 * (2.0 - 4.0 * sinth  + 64.0 / 9.0 * sinth2);
      Pd = 3.0 * (2.0 - 4.0 * sinth  + 16.0 / 9.0 * sinth2);
      Pl = 3.0 * (2.0 - 4.0 * sinth + 8.0 * sinth2);
      Pf = Pu + Pd + Pl;
      Pg = -44.0 + 106.0 * sinth - 62.0 * sinth2;
      Pgp = 10.0 + 34.0 * sinth - 26.0 * sinth2;
      P2h = -10.0 + 2.0 * sinth - 2.0 * sinth2;
      P2hp = 8.0 - 22.0 * sinth + 10.0 * sinth2;
      double cos2b2 = sqr(cos(2.0 * beta)), cos2b4 = sqr(cos2b2);
      P1h = -9.0 * cos2b4 + (1.0 - 2.0 * sinth + 2.0 * sinth2) * cos2b2;
      
      // non bottom-mixing corrections: checked
      // Following was obtained from FeynHiggsFast.f
      double msusy2 = ms2 - mt2;
      if (msusy2 < 0.0) {
	flagTachyon(true);
	if (PRINTOUT > 2) cout << " msusy2 tachyon in mH calculation ";
      }
      double logMsScale = log(fabs(msusy2) / mz2);
      // Note if B undefined or zero, will have a problem!
      double logA0Scale = log(fabs(mAsq) / mz2);
      if (fabs(mAsq) < EPSTOL) logA0Scale = 0.0;
      double mz4 = sqr(mz2), cosb4 = sqr(cosb2), mb2 = sqr(mb), mb4 =
	sqr(mb2);
      double howie11 = root2 * GMU * cosb2 * mz4 / (24.0 * pisq);
      double howie12 = howie11 * tanb, howie22 = howie12 * tanb;
      
      sigma(1, 1) = sigma(1, 1) - howie11 * 
	(
	 (36.0 * mb4 / (mz4 * cosb4) - 18.0 * mb2 / (mz2 * cosb2) 
	  + Pb + Pf + Pg + P2h) * logMsScale 
	 + (P1h - P2h) * logA0Scale
	 );
      sigma(1, 2) = sigma(1, 2) + howie12 * 
	(
	 (- 9.0 * mb2 / (mz2 * cosb2) + Pb + Pf + Pgp + P2hp) 
	 * logMsScale - (P1h + P2hp) * logA0Scale
	 );
      sigma(2, 2) = sigma(2, 2) - howie22 * 
	(
	 (Pb + Pf + Pg + P2h) * logMsScale + (P1h - P2h) * logA0Scale
	 ); 
      
      // sbottom mixing corrections: checked
      double mblr = displaySoftA(DA, 3, 3) - mu * tanb, 
	mblr2 = sqr(mblr);
      double mu2 = sqr(mu);
      
      sigma(1, 1) = sigma(1, 1) -
	(root2 * GMU / PI * 3.0 * 
	 (-(mb2 / msusy2 * mz2 * (mLR + mu / tanb) * 
	    (mblr + (mLR + mu / tanb) / 3.0)) + 4.0 / cosb2 *
	  mb4 * mblr / msusy2 * (mLR + mu / tanb) * 
	  (1.0 - (mblr / msusy2 * (mLR + mu / tanb)) / 12.0))) / 
	(8.0 * PI);
      sigma(1, 2) = sigma(1, 2)
	+ ((root2 * GMU / PI * 3.0 * 
	    (-4.0 / cosb2 * mb4 * mblr / msusy2
	     * mu * (1.0 - (mblr / msusy2 * (mLR + mu / tanb)) / 6.0) -
	     mb2 / msusy2 * mz2 * (mblr * (mLR + 2.0 * mu / tanb ) + 
				   (mu2 + sqr(mLR + mu / tanb)) / 3.0) *
	     tanb)) / (16.0 * PI));
      sigma(2, 2) = sigma(2, 2) -
	(root2 * GMU / (8.0 * pisq) * 3.0 * 
	 (- ((1.0 / cosb2 * mb4 * mblr2 / sqr(msusy2) * mu2) / 3.0) - mb2
	  / msusy2
	  * mu * mz2 * tanb * (mblr + (mu * tanb) / 3.0))); 
      
      /*         two loop contributions: checked  from hep-ph/9903404   */
      if (accuracy > 2) {
	sigma(2, 2) = sigma(2, 2) - 
	  (root2 * GMU / (PI * pisq) * alphas * mt4 / sinb2 * 
	   ((6.0 * mLR) / ms - (12.0 * mLR4) / (17.0 * ms4) + 
	    (2.0 - (3.0 * mLR2) / ms2) * log(ms2 / mt2) - 
	    3.0 * sqr(log(ms2 / mt2))
	    - 4.0 + 8.0 * mLR2 / ms2));
      }      
    }
    sigma(2, 1) = sigma(1, 2);
    
    mH = mH - sigma;
  }
  
  DoubleVector temp(2);  
  double theta;
  
  temp = mH.sym2by2(theta);

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
  if (temp(2) > temp(1)) { 
    theta = theta + PI * 0.5; 
    if (PRINTOUT > 3) {
      cout << "Changing alpha from " << theta - PI * 0.5;
      cout << " to " << theta << " resulting in:";
      DoubleMatrix O(rot2d(theta));
      cout << O * mH * O.transpose() << flush;
    }
  }
  
  drBarPars forLoops(displayDrBarPars());
  sPhysical physpars(displayPhys());

  double mA = forLoops.mhiggs(2);
  double poleMasq = (displayMh2Squared() - displayMh1Squared() )
    / cos(2.0 * beta) - sqr(mz);
  
  double piaa = 0.0;
  if (accuracy > 0) {
      piaa = piAA(mA, displayMu());
      poleMasq = 
	(displayMh2Squared() - displayTadpole2Ms() - displayMh1Squared() +
	 displayTadpole1Ms()) / cos(2.0 * beta) 
	- sqr(mz) - pizztMS - piaa +
	sqr(sin(beta)) * displayTadpole1Ms() + sqr(cos(beta)) *
	displayTadpole2Ms();
    }

  // FUTURE WORK: could incorporate more accurate calculations on charged
  // Higgs from   hep-ph/9606211: include piH+H-
  double poleMhc = poleMasq + sqr(displayMw()) + piaa + piwwtMS;
  
  physpars.mhiggs(1) = minimum(temp(2), temp(1));
  physpars.mhiggs(2) = ccbSqrt(poleMasq);
  physpars.mhiggs(3) = maximum(temp(1), temp(2)); 
  physpars.mhiggs(4) = ccbSqrt(poleMhc);
  
  // In case there's an mA^2 < 0 problem in an early iteration
  if (displayProblem().tachyon) physpars.mhiggs(2) = physpars.mhiggs(3) = 
				  physpars.mhiggs(4) = ccbSqrt(poleMasq);
  
  // two-loop Yukawa correction to mh^2: checked
  if (accuracy > 2) {
    double mstop1 = forLoops.mu(1, 3);
    double mstop2 = forLoops.mu(2, 3);
    double mstop1sq = sqr(mstop1), mstop2sq = sqr(mstop2);
    double t = 0.5 * log(mstop1sq * mstop2sq / mt4);
    double X = ( sqr( (mstop2sq - mstop1sq) / (4.0 * mt2) * 
		      sqr(sin(2.0 * forLoops.thetat)) ) * 
		 (2.0 - (mstop2sq + mstop1sq) / (mstop2sq - mstop1sq) *
		  log(mstop2sq / mstop1sq)) +
		 (mstop2sq - mstop1sq) / (2.0 * mt2) * 
		 sqr(sin(2.0 * forLoops.thetat)) * log(mstop2sq / mstop1sq)
		 );
    physpars.mhiggs(1) = ccbSqrt(sqr(physpars.mhiggs(1)) + 9.0 / 
			      (16.0 * sqr(pisq)) *
			      sqr(GMU) * mt2 * mt4 * (X * t + sqr(t)));
			      } 
  physpars.thetaH = theta;

  setPhys(physpars);
}
