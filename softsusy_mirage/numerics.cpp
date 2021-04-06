
/** \file numerics.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: numerics.cpp,v $
   Revision 1.4  2006/04/25 11:32:22  allanach
   Fixed B22, B0, B1 - previous expression switched too soon to p=0 
   expression, leading to convergence failures where there shouldn't have been
   any (in particular, at high m0 and m12 for instance).

   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.8  2005/07/15 15:10:47  allanach
   Added analytic dilog routine

   Revision 1.7  2005/06/16 13:57:04  allanach
   Added a Cauchy-distributed random number generator

   Revision 1.6  2005/05/30 14:22:14  allanach
   Fixed stau mixing in ISAWIG interface to 7.64

   Revision 1.5  2005/05/13 16:07:27  allanach
   Edited precision due to double precision

   Revision 1.4  2005/04/13 14:53:54  allanach
   Corrected binning procedure so it does what you expect

   Revision 1.3  2005/04/12 10:44:20  allanach
   Added bin function to calculate the bin number of some data

   Revision 1.2  2005/04/11 14:06:36  allanach
   Added random number routines

   Revision 1.18  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.17  2003/08/19 14:26:22  allanach
   Changing lowOrg to be more sensible about gauge unification. Should now be
   called with POSITIVE mgut and a flag for gauge unification.

   Revision 1.16  2003/07/28 12:11:37  allanach
   More error trapping, and rearranging rpvsoftsusy to use correct Higgs VEV
   (which is sometimes called at MZ)

   Revision 1.15  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.14  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.13  2003/03/25 17:03:05  allanach
   Added extra case for b1(0,m1,0,q)

   Revision 1.12  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.10  2002/12/19 18:26:40  allanach
   Fixed numerical rounding error in d27 function

   Revision 1.7  2002/10/14 17:14:29  allanach
   Bug-fixed c0 function

   Revision 1.6  2002/09/09 10:42:54  allanach
   TOLERANCE replaces EPS as being more user-friendly

   Revision 1.5  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, TOLERANCE out of def.h to make it quicker to
   compile once they are changed.

   Revision 1.4  2002/07/30 12:57:31  allanach
   SOFTSUSY1.5

   Revision 1.3  2002/04/18 14:32:05  allanach
   Changed RGEs and anomalous dimensions to be compatible with new notation;
   started implementation of rewsb in R-parity violation

   Revision 1.6  2001/10/04 19:26:34  allanach
   New version deals with AMSB correctly

   Revision 1.4  2001/09/28 13:50:13  allanach
   More careful treatment of limits in Passarino-Veltman functions b0, b1.

   Revision 1.3  2001/08/08 09:52:33  allanach
   Added dilogarithm function - could be speeded up....

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#include "numerics.h"

// returns >0 if there's a problem:
int integrateOdes(DoubleVector & ystart, double from, double to, double eps,
	      double h1, double hmin, 
	      DoubleVector (*derivs)(double, const DoubleVector &),
	      int (*rkqs)(DoubleVector & y, const DoubleVector & dydx, double
			   *x, double htry, double eps, DoubleVector & yscal,
			   double *hdid, double *hnext, 
			   DoubleVector (*derivs)(double, const DoubleVector
						  &)) ) {  
  int nvar =  ystart.displayEnd();
  int nstp, i;
  double x, hnext, hdid, h;
  DoubleVector yscal(nvar), y(ystart), dydx(nvar);
  
  x = from;
  h = sign(h1, to - from);
  
  const int MAXSTP = 400;
  const double TINY = 1.0e-16;

  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    dydx = (*derivs)(x, y);
    for (i = 1; i <= nvar; i++)
      yscal(i) = fabs(y(i)) + fabs(dydx(i) * h) + TINY;
    if ((x + h - to) * (x + h - from) > 0.0) h = to - x;
    int smallStep = (*rkqs)(y, dydx, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (smallStep) return 1;

    if ((x - to) * (to - from) >= 0.0) {
      for (i = 1; i<= nvar; i++) ystart(i) = y(i);
      return 0;
    }
      
    if (fabs(hnext) <= hmin) {
      nstp = MAXSTP; // bail out
      if (PRINTOUT > 1) {
	cout << "Step size too small in numerics.cpp:integrateOdes\n";
	cout << "**********x = " << x << "*********\n";
	for (i = 1;i<= nvar;i++) 
	  cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	    ") = " << dydx(i) << endl;
	cout.flush();
      }
    }
    
    h = hnext;
  }
  
  if (PRINTOUT > 1) {
    cout << "Bailed out of numerics.cpp:too many steps in integrateOdes\n";
    cout << "**********x = " << x << "*********\n";
    for (i = 1;i<= nvar;i++) 
      cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	") = " << dydx(i) << endl;
    cout.flush();
  }
  
  return 1;
}

int odeStepper(DoubleVector & y, const DoubleVector & dydx, double *x, double
		htry, double eps, DoubleVector & yscal, double *hdid, 
		double *hnext,		
		DoubleVector (*derivs)(double, const DoubleVector &))
{
  const double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;

  int i, n = y.displayEnd();
  double errmax, h, htemp, xnew;
  
  DoubleVector yerr(n), ytemp(n);
  h = htry;
  for (;;) {
    rungeKuttaStep(y, dydx, *x, h, ytemp, yerr, derivs);
    errmax = 0.0;
    for (i = 1; i<= n;i++) errmax = maximum(errmax, fabs(yerr(i) / yscal(i)));
    errmax  /= eps;
    if (errmax <= 1.0) break;
    htemp = SAFETY * h * pow(errmax, PSHRNK);
    h = (h >= 0.0 ? maximum(htemp ,0.1 * h) : minimum(htemp, 0.1 * h));
    xnew = (*x) + h;
    if (xnew == *x) 
      {
	if (PRINTOUT) {
	cout << "At x = " << *x;
	cout << ",stepsize underflow in odeStepper" << flush << endl;
	}
	return 1;
      }
  }
  if (errmax > ERRCON) *hnext = SAFETY * h * pow(errmax,PGROW);
  else *hnext = 5.0 * h;
  *x += (*hdid = h);
  y = ytemp;
  return 0;
}

void rungeKuttaStep(const DoubleVector & y, const DoubleVector & dydx, 
	     double x, double h, DoubleVector & yout, DoubleVector & yerr, 
	     DoubleVector (*derivs)(double, const DoubleVector &)) {
  int i;
  const static double a2 = 0.2,a3 = 0.3,a4 = 0.6,a5 = 1.0,a6 = 0.875,b21 =
    0.2,b31 = 3.0 / 40.0,b32 = 9.0 / 40.0,b41 = 0.3,b42 = -0.9,b43 = 1.2,
    b51 = -11.0 / 54.0, b52 = 2.5,b53 = -70.0 / 27.0,b54 = 35.0 / 27.0,
    b61 = 1631.0 / 55296.0,b62 = 175.0 / 512.0,b63 = 575.0 / 13824.0,
    b64 = 44275.0 / 110592.0,b65 = 253.0 / 4096.0,c1 = 37.0 / 378.0,
    c3 = 250.0 / 621.0,c4 = 125.0 / 594.0,c6 = 512.0 / 1771.0,
    dc5 = -277.00 / 14336.0;
  const double dc1 = c1-2825.0 / 27648.0,dc3 = c3-18575.0 / 48384.0,
    dc4 = c4-13525.0 / 55296.0,dc6 = c6-0.25;
  
  int n = y.displayEnd();
  
  DoubleVector ytemp(y.display() + b21 * h * dydx.display());
  DoubleVector ak2((*derivs)(x + a2 * h, ytemp));

  // Allowing piece-wise calculating of ytemp for speed reasons
  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b31 * dydx.display(i) + b32 * ak2(i));
  DoubleVector ak3((*derivs)(x + a3 * h, ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b41 * dydx.display(i) + b42 * ak2(i) + b43
				   * ak3(i));
  DoubleVector ak4((*derivs)(x+a4*h,ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b51 * dydx.display(i) + b52 * ak2(i) + b53
				   * ak3(i) + b54 * ak4(i));
  DoubleVector ak5((*derivs)(x + a5 * h, ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b61 * dydx.display(i) + b62 * ak2(i) + b63
				   * ak3(i) + b64 * ak4(i) + b65 * ak5(i));
  DoubleVector ak6((*derivs)(x + a6 * h, ytemp));

  for (i = 1; i<= n; i++)
    yout(i) = y.display(i) + h * (c1 * dydx.display(i) + c3 * ak3(i) + c4 *
				  ak4(i) + c6 * ak6(i));
  for (i = 1; i<= n; i++)
    yerr(i) = h * (dc1 * dydx.display(i) + dc3 * ak3(i) + 
		   dc4 * ak4(i) + dc5 * ak5(i) + dc6 * ak6(i));
}

double calcDerivative(double (*func)(double), double x, double h, double
		      *err){
  const double CON = 1.4, CON2 = CON * CON, BIG = 1.0e30, 
    SAFE = 2.0; 
  const int NTAB = 10;
  
  int i, j;
  double errt, fac, hh, ans = 0.0;
  
  if (h == 0.0) throw "h must be nonzero in numerics.cpp:calcDerivative";


  DoubleMatrix a(NTAB, NTAB);
  hh = h;
  a(1, 1) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
  *err = BIG;
  for (i=2; i<=NTAB; i++) {
    hh /= CON;
    a(1, i) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
    fac = CON2;
    for (j=2; j<=i; j++) {
      a(j, i) = (a(j-1, i) * fac - a(j-1, i-1)) / (fac - 1.0);
      fac = CON2 * fac;
      errt = maximum(fabs(a(j, i) - a(j-1, i)), fabs(a(j, i) - a(j-1, i-1)));
      if (errt <= *err) {
	*err = errt;
	ans = a(j, i);
      }
    }
    if (fabs(a(i, i) - a(i-1, i-1)) >= SAFE * (*err)) break;
  }

  return ans;
}

inline void shft2(double & a, double & b, double c) { a = b; b = c; }

inline void shft3(double & a, double & b, double & c, double d) { 
  a = b; b = c; c = d;
}

double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin)
{
  const double R = 0.61803399, C = 1.0 - R;
  double f1, f2, x0, x1, x2, x3;
  
  x0 = ax; 
  x3 = cx; 
  if (fabs(cx - bx) > fabs(bx - ax)) {
    x1 = bx; 
    x2 = bx + C * (cx - bx); 
  } else {
    x2 = bx; 
    x1 = bx - C * (bx - ax); 
  }
  f1 = (*f)(x1); 
  f2 = (*f)(x2); 
  while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2))) {
    if (f2 < f1) {
      shft3(x0, x1, x2, R * x1 + C * x3);
      shft2(f1, f2, (*f)(x2));
    } else {
      shft3(x3, x2, x1, R * x2 + C * x0);
      shft2(f2, f1, (*f)(x1));
	}
  }
  if (f1 < f2) {
    *xmin = x1; 
    return f1; 
  } else {
    *xmin = x2; 
    return f2; 
  }
}



DoubleVector dd(double x, const DoubleVector & y) {
  DoubleVector dydx(1);
  dydx(1) = -integrandThreshbnr(x);
  return dydx;
}

double integrandThreshbnr(double x) {
  return fnfn(x).real();
}

// Integration routine needs these variables
static double mtInt, pInt, m1Int, m2Int;
static int nInt;

Complex fnfn(double x) {
  const static Complex iEpsilon(0.0, TOLERANCE * 1.0e-20);
  
  double xn = 1.0;
  int i; for (i=1; i<=nInt; i++) xn = xn * x;
  return xn * 
    log( ((1 - x) * sqr(m1Int) + x * sqr(m2Int) - x * (1 - x) *
	  sqr(pInt) - iEpsilon)
	 / sqr(mtInt));
}

DoubleVector dilogarg(double t, const DoubleVector & y) {

  const double eps = TOLERANCE * 1.0e-20;

  DoubleVector dydx(1);
  dydx(1) = -log(fabs(1 - t + eps)) / (t + eps);

  return dydx;
}

/*
double dilog(double x) {
  // Set global variables so that integration function can access them
  double from = 0.0, to = x, guess = 0.1, hmin = TOLERANCE * 1.0e-5;

  DoubleVector v(1); 
  double eps = TOLERANCE * 1.0e-5;
  v(1) = 1.0; 

  // Runge-Kutta, f(b) = int^b0 I(x) dx, I is integrand => d f / db = I(b)
  // odeint has a problem at f(0): therefore, define f'(b)=f(b)+1
  integrateOdes(v, from, to, eps, guess, hmin, dilogarg, odeStepper); 
  
  return v(1) - 1.0;
}
*/

// Returns real part of integral
double bIntegral(int n1, double p, double m1, double m2, double mt) {
  // Set global variables so that integration function can access them
  nInt = n1; pInt = p; m1Int = m1; m2Int = m2; mtInt = mt;
  double from = 0.0, to = 1.0, guess = 0.1, hmin = TOLERANCE * 1.0e-5;
  
  DoubleVector v(1); double eps = TOLERANCE * 1.0e-3;
  v(1) = 1.0; 

  // Runge-Kutta, f(b) = int^b0 I(x) dx, I is integrand => d f / db = I(b)
  // odeint has a problem at f(0): therefore, define f'(b)=f(b)+1
  integrateOdes(v, from, to, eps, guess, hmin, dd, odeStepper); 
  
  return v(1) - 1.0;
}

/*
  Analytic expressions follow for above integrals: sometimes useful!
  From hep-ph/9606211
  Note it returns the REAL PART ONLY. 
  */
double b0(double p, double m1, double m2, double q) {

  const double pTolerance = EPSTOL;

  if (sqr(p) > pTolerance * maximum(sqr(m1), sqr(m2))) {
    Complex iEpsilon(0.0, EPSTOL * p);

    double s = sqr(p) - sqr(m2) + sqr(m1); 
    
    Complex 
      xPlus = (s + sqrt(sqr(s) - 4.0 * sqr(p) * (sqr(m1) - iEpsilon))) /
      (2.0 * sqr(p));
    Complex 
      xMinus = (s - sqrt(sqr(s) - 4.0 * sqr(p) * (sqr(m1) - iEpsilon))) /
      (2.0 * sqr(p));
    
    return -log(sqr(p) / sqr(q)) - fB(xPlus) - fB(xMinus);
  }
  else {
    if (close(m1, m2, EPSTOL))
      return - log(sqr(m1 / q));
    else {
      double Mmax2 = maximum(sqr(m1), sqr(m2)), 
	Mmin2 = minimum(sqr(m1), sqr(m2)); 
      if (Mmin2 < sqr(TOLERANCE)) return 1.0 - log(Mmax2 / sqr(q));
      return 
	1.0 - log(Mmax2 / sqr(q)) + Mmin2 * log(Mmax2 / Mmin2) 
	/ (Mmin2 - Mmax2);
    }
  }   
}

double b1(double p, double m1, double m2, double q) {

  const double pTolerance = EPSTOL;

  if (sqr(p) > pTolerance * maximum(sqr(m1), sqr(m2)) ) {
    return (a0(m2, q) - a0(m1, q) + (sqr(p) + sqr(m1) - sqr(m2)) 
	    * b0(p, m1, m2, q)) / (2.0 * sqr(p)); 
  }
  else if (fabs(m1) > EPSTOL && !close(m1, m2, EPSTOL) 
	   && fabs(m2) > EPSTOL) {// checked
    double Mmax2 = maximum(sqr(m1) , sqr(m2)), x = sqr(m2 / m1);
    return 0.5 * (-log(Mmax2 / sqr(q)) + 0.5 + 1.0 / (1.0 - x) + log(x) /
		  sqr(1.0 - x) - theta(1.0 - x) * log(x)); // checked
  }
  
  return bIntegral(1, p, m1, m2, q);
}

double b22(double p,  double m1, double m2, double q) {

  double answer;
  
  const double pTolerance = EPSTOL;

  if (sqr(p) < pTolerance * maximum(sqr(m1), sqr(m2)) ) {
    // m1 == m2 with good accuracy
    if (close(m1, m2, EPSTOL)) 
      answer = -sqr(m1) * log(sqr(m1 / q)) * 0.5 + sqr(m1) * 0.5;
    else
      if (fabs(m1) > EPSTOL && fabs(m2) > EPSTOL)
	answer = 0.375 * (sqr(m1) + sqr(m2)) - 0.25 * 
	  (sqr(sqr(m2)) * log(sqr(m2 / q)) - sqr(sqr(m1)) * 
	   log(sqr(m1 / q))) / (sqr(m2) - sqr(m1)); 
      else
	if (fabs(m1) < EPSTOL)
	  answer = 0.375 * sqr(m2) - 0.25 * sqr(m2) * log(sqr(m2 / q));
	else 
	  answer = 0.375 * sqr(m1) - 0.25 * sqr(m1) * log(sqr(m1 / q));
  }
  else {// checked
    double b0Save = b0(p, m1, m2, q);
    
    answer = 1.0 / 6.0 * 
      (0.5 * (a0(m1, q) + a0(m2, q)) + (sqr(m1) + sqr(m2) - 0.5 * sqr(p))
       * b0Save + (sqr(m2) - sqr(m1)) / (2.0 * sqr(p)) *
       (a0(m2, q) - a0(m1, q) - (sqr(m2) - sqr(m1)) * b0Save) +
       sqr(m1) + sqr(m2) - sqr(p) / 3.0);
  }
  
  return answer;
}

double d0(double m1, double m2, double m3, double m4) {// checked
  if (close(m1, m2, EPSTOL)) {
    return 
      (sqr(m4) / sqr(sqr(m1) - sqr(m4)) * log(sqr(m4) / sqr(m1)) + 
       1.0 / (sqr(m1) * (sqr(m1) - sqr(m4))) -
       sqr(m3) / sqr(sqr(m1) - sqr(m3)) * log(sqr(m3) / sqr(m1)) -
       1.0 / (sqr(m1) * (sqr(m1) - sqr(m3)))) / (sqr(m3) - sqr(m4));
  }
  return (c0(m1, m3, m4) - c0(m2, m3, m4)) / (sqr(m1) - sqr(m2));
}

double d27(double m1, double m2, double m3, double m4) {// checked
  if (close(m1, m2, EPSTOL)) {
    double m1n = m1 + TOLERANCE * 0.01;
    return (sqr(m1n) * c0(m1n, m3, m4) - sqr(m2) * c0(m2, m3, m4)) 
      / (4.0 * (sqr(m1n) - sqr(m2)));
  }
  return (sqr(m1) * c0(m1, m3, m4) - sqr(m2) * c0(m2, m3, m4)) 
    / (4.0 * (sqr(m1) - sqr(m2)));
}

// Bug-fixed 14.10.02 by T. Watari and collaborators - many thanks!
double c0(double m1, double m2, double m3) {

  if (close(m2, m3, EPSTOL)) {
    if (close(m1, m2, EPSTOL)) {
      return ( - 0.5 / sqr(m2) ); // checked 14.10.02
    }
    else {
      return ( sqr(m1) / sqr(sqr(m1)-sqr(m2) ) * log(sqr(m2)/sqr(m1))
               + 1.0 / (sqr(m1) - sqr(m2)) ) ; // checked 14.10.02
    }
  }
  else
    if (close(m1, m2, EPSTOL)) {
      return ( - ( 1.0 + sqr(m3) / (sqr(m2)-sqr(m3)) * log(sqr(m3)/sqr(m2)) )
               / (sqr(m2)-sqr(m3)) ) ; // checked 14.10.02
    }
    else
      if (close(m1, m3, EPSTOL)) {
        return ( - (1.0 + sqr(m2) / (sqr(m3)-sqr(m2)) * log(sqr(m2)/sqr(m3))) 
                 / (sqr(m3)-sqr(m2)) ); // checked 14.10.02
      }
      else return (1.0 / (sqr(m2) - sqr(m3)) * 
		   (sqr(m2) / (sqr(m1) - sqr(m2)) *
		    log(sqr(m2) / sqr(m1)) -
		    sqr(m3) / (sqr(m1) - sqr(m3)) *
		    log(sqr(m3) / sqr(m1))) );
}

double gasdev(long & idum) {
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

double ran1(long & idum) {
  const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32, 
    NDIV = 1+(IM-1)/NTAB;
  const double AM = 1.0 / double(IM), EPS = 1.2e-15, RNMX = 1.0 - EPS;
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (idum <= 0 || !iy) {
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(idum)/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if (idum < 0) idum += IM;
      if (j < NTAB) iv[j] = idum;
    }
    iy=iv[0];
  }
  k=(idum)/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

double cauchyRan(long & idum) {
  double x = ran1(idum) - 0.5;
  double unNormalised = tan(x * PI);
  return unNormalised;
}

int bin(double data, double start, double end, int numBins) {
  double range = end - start;
  double binSize = range / double(numBins);
  return int((data - start) / binSize + 1.);
}

double logOfSum(double a, double b) {
  double max = maximum(a, b);
  double min = minimum(a, b);

  if (max + min < 0. || max == 0.) return asin(1.0);
  double ans = log(max);
  ans = ans + log (1.0 + min / max);

  return ans;
}


double dilog(double x) {
   // The DiLogarithm function
   // Code translated by R.Brun from CERNLIB DILOG function C332

   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};
   
   double T,H,Y,S,A,ALFA,B1,B2,B0;
   
   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

