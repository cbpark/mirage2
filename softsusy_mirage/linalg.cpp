
/** \file linalg.cpp 
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach 
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
    
    $Log: linalg.cpp,v $
    Revision 1.3  2005/11/09 14:12:24  allanach
    Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

    Revision 1.4  2005/07/26 10:51:22  allanach
    new method nmin added to DoubleVector and DoubleMatrix by Bernhardt

    Revision 1.3  2005/05/13 16:06:02  allanach
    Relaxed tolerance a bit in diagonalisation: it was giving me trouble

    Revision 1.1  2004/11/19 16:18:31  allanach
    Initial revision

    Revision 1.25  2004/01/15 13:53:48  allanach
    New time-saving unary operators included.

    Revision 1.24  2003/08/04 16:44:22  allanach
    Added max function for matrix and vector of doubles

    Revision 1.23  2003/07/29 15:26:40  allanach
    Using new and delete rather than old malloc and free C-style routines

    Revision 1.22  2003/07/25 13:39:15  allanach
    Trapped errors properly rather than exiting

    Revision 1.21  2003/06/05 09:17:19  allanach
    Started coding Les Houches Discord
    
    Revision 1.20  2003/05/20 15:18:08  allanach
    Changed to doxygen comments and cout rather than cerr
    
    Revision 1.19  2003/03/09 18:52:19  allanach
    Added check for symmetricity in DoubleMatrix::sym2by2
    
    Revision 1.18  2003/03/04 12:55:02  allanach
    Removed Slavich diagonalisation
    
    Revision 1.17  2003/02/28 18:11:57  allanach
    Added routines to convert to Slavich's stop conventions in order to
    calculate the 2-loop Higgs mass
    
    Revision 1.16  2003/02/21 13:02:06  allanach
    Changed headings to new conventions
    
    Revision 1.15  2002/11/04 15:51:50  allanach
    Added a check in sym2by2 for the case of diagonal matrices
    
    Revision 1.14  2002/10/30 16:39:30  allanach
    Fixed annoying extra spaces bug when outputting vectors
    
    Revision 1.13  2002/10/09 17:18:51  allanach
    Added new softpoint option to specify mu and mA instead of mH1, mH2. Some
    fine-tuning still necessary
    
    Revision 1.12  2002/09/23 18:13:50  allanach
    Eigenvalue order with non symmetric diagonalisation setx
    
    Revision 1.11  2002/09/11 09:52:34  allanach
    TOLERANCE etc now used in preference.
    
    Revision 1.10  2002/08/12 16:29:10  allanach
    Cosmetic change in outputting complex vector
    
    Revision 1.9  2002/07/19 15:54:06  allanach
    SOFTSUSY1.5 version
    
    Revision 1.8  2002/06/14 16:26:30  allanach
    Switches included for 2-loop running of scalar masses, and calulating mt at mt.
    
    Revision 1.7  2002/06/10 16:13:07  allanach
    Speeded things up by adding option for array bounds checking
    
    Revision 1.6  2002/04/30 16:23:59  allanach
    Added matrix inversion routine
    
    Revision 1.5  2002/04/26 15:14:44  allanach
    Deleted all translation routines and defined boundary conditions within
    softsusy.h and softsusy.cpp
    
    Revision 1.2  2002/04/11 12:22:27  allanach
    Traps -0.0 and just prints 0.0
*/

#include "linalg.h"

double * DoubleVector::myDoubleVector(const long nl, const long nh) {
  double *v;
  
  //  v=(double *) malloc((size_t) ((nh - nl + 2) * sizeof(double)));
  v = new double[(nh - nl + 2)];
  if (!v) throw("Allocation failure in vector()"); 
  
  return v - nl + 1;
}

void DoubleVector::setVec(const double *t, long st, long nd) { 
  x = myDoubleVector(st, nd); int i; 
  for(i = st; i <= nd; i++) x[i] = t[i];
}

// a is in old C convention ie index starts at ZERO! (For outputting to
// fortran)
void DoubleVector::display(double * a) const {
  int i; for(i = displayStart(); i <= displayEnd(); i++) a[i-1] = x[i];  
}

DoubleVector::DoubleVector(int e)
  : start(1), end(e) { 
  double * t = myDoubleVector(1, e);
  int i; for(i=1; i<=e; i++) t[i] = double(0);
  setVec(t, 1, e);
  freeMyDoubleVector(t, 1, e);
}

DoubleVector::DoubleVector(int s, int e)
  : start(s), end(e) {
  double * t = myDoubleVector(s, e);
  int i; for(i=s; i<=e; i++) t[i] = double(0);
  setVec(t, s, e);
  freeMyDoubleVector(t, s, e);
}

// Changes the length of a vector - copies as many elements of old one as
// possible, and fills any extra up with zeroes
void DoubleVector::setEnd(int e) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (e < 1) {
  ostringstream ii;
  ii << "DoubleVector.setEnd called with incompatible length " << e;
  ii << *this;
  throw(ii.str());
  }
#endif
  double * t = myDoubleVector(start, e);
  int i; for(i=start; i<=e; i++) {
    if (i <= end) t[i] = x[i]; 
    else
      t[i] = double(0);
  }
  freeMyDoubleVector(x, start, end);
  setVec(t, start, e);
  freeMyDoubleVector(t, start, e);
  
  end = e;
}

const DoubleVector & DoubleVector::operator=(const DoubleVector &v) {
  if (this == &v) return *this;
  // Are the vectors of equivalent size?
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.start || end != v.end) 
    {
      ostringstream ii;      
      ii << "DoubleVector = overload; incompatible lengths\n";
      ii << *this << "=\n" << v;
      throw(ii.str());
    }
#endif
  freeMyDoubleVector(x, start, end);
  setVec(v.x, v.start, v.end);
  
  return *this;
}

DoubleVector DoubleVector::operator+(const DoubleVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.displayStart() || end != v.displayEnd())
    {
      ostringstream ii;
      ii << "DoubleVector + overload; incompatible lengths\n";
      ii << *this << "+\n" << v;
      throw(ii.str());
    }
#endif
  DoubleVector temp(v);
  int i; for(i=start; i<=end; i++) temp(i) = v.display(i) + x[i];
  return temp;
}

DoubleVector DoubleVector::operator-(const DoubleVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.displayStart() || end != v.displayEnd())
    {
      ostringstream ii;
      ii << "DoubleVector - overload; incompatible lengths\n";
      ii << *this << "-\n" << v;
      throw(ii.str());
    }
#endif
  DoubleVector temp(v);
  int i; for(i=start; i<=end; i++) temp(i) = x[i] - v.display(i);
  return temp;
}

// NOT dot product, but product of elements.
DoubleVector DoubleVector::operator*(const DoubleVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.displayStart() || end != v.displayEnd())
    {
      ostringstream ii;
      ii << "DoubleVector * overload; incompatible lengths\n";
      ii << *this << "*\n" << v;
      throw(ii.str());
    }
#endif
  DoubleVector temp(v);
  int i; for(i=start; i<=end; i++) temp(i) = v.display(i) * x[i];
  return temp;
}

DoubleVector DoubleVector::operator*(double f) {
  DoubleVector temp(*this);
  int i; for(i=start; i<=end; i++) temp(i) = temp(i) * f;
  return temp;
}

double DoubleVector::dot(const DoubleVector & v) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (v.displayStart() != start || v.displayEnd() != end)
    {
      ostringstream ii;
      ii << "Scalar product between incompatible vectors " << *this << 
	" and " << v;
      throw ii.str();
    }
#endif
  double f = 0.0;
  int i; for (i=start; i<=end; i++) f += x[i] * v.display(i);
  return f;				   
}

DoubleVector DoubleVector::apply(double (*fn)(double)) {
  DoubleVector temp(displayStart(), displayEnd());
  int i;
  for (i=displayStart(); i<=displayEnd(); i++)
    temp(i) = double(fn(display(i)));
  return temp;
}

double DoubleVector::max() const {
  double m = double(0);
  int i; 
  for (i=displayStart(); i<=displayEnd();  i++)
    if (fabs(display(i)) > m) m = display(i);
  return m;
}

// Most negative / smallest element
double DoubleVector::min(int & p) const {
  double m = display(displayStart());
  p = displayStart();
  int i; 
  for (i=displayStart() + 1; i<=displayEnd();  i++)
    if (display(i) < m) 
      {
	m = display(i);
	p = i;
      }
  return m;
}

void DoubleVector::swap(int i, int j) {
  double m;
  DoubleVector temp(*this);
  m = temp.display(j);
  temp(j) = temp.display(i);
  temp(i) = m;
  *this = temp; 
}

// Next-to-Most negative / smallest element
double DoubleVector::nmin(int & p) const {
  double m1 = display(displayStart()), m2;
  int p1 = displayStart(), p2;
  
  if ( display(displayStart()+1) < m1) {
    m2 = m1;
    p2 = p1;
    
    m1 = display(displayStart()+1);
    p1 = displayStart()+1;
  }
  else {
    m2 = display(displayStart()+1);
    p2 = displayStart()+1;
  }
  
  
  int i; 
  for (i=displayStart() + 2; i<=displayEnd();  i++){
    if (display(i) < m1) {
	m2 = m1;
	p2 = p1;
	
	m1 = display(i);
	p1 = i;
      }
    else if (display(i) < m2)  {
	m2 = display(i);
	p2 = i;
      }
  }
  
  p =  p2;
  return m2;
}


DoubleVector operator*(double f, const DoubleVector & v) {
  DoubleVector temp(v);
  int i;
  for (i=v.displayStart(); i<=v.displayEnd(); i++)
    temp(i) = temp(i) * f;
  return temp;
} 

const int NUMDisplay = 5;
ostream & operator <<(ostream &left, const DoubleVector &v) {
  left << "(" << v.displayStart() << "," << v.displayEnd() << "):\n";
  int i; for(i=v.displayStart(); i<=v.displayEnd(); i++) {
    if (v.display(i) >= 0.0) left << " ";
    left << v.display(i);
    if ((i % NUMDisplay) == 0 || (i == v.displayEnd())) 
      left << '\n';
    else
      left << " ";
  }
  return left;
}

istream & operator>>(istream & left, DoubleVector &v) {
  int i;
  DoubleVector x(v.displayStart(), v.displayEnd());
  char c[70];
  cin >> c; 
  for (i=v.displayStart(); i<=v.displayEnd(); i++) {
    left >> x(i); 
  }
  v = x;
  
  return left;
}


/* -----------------------------------------------------*/
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double ** DoubleMatrix::myDoubleMatrix(long nrl, long nrh, long ncl, long nch)
  const { 
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;
  
  /* allocate pointers to rows */
  //  m=(double **) malloc((size_t)((nrow + 1) * sizeof(double*)));
  m = new double* [(nrow + 1)];
  if (!m) { throw "allocation failure 1 in matrix()";  }
  m += 1;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  //  m[nrl]=(double *) 
  // malloc((size_t)((nrow * ncol + 1) * sizeof(double)));
  m[nrl] = new double [(nrow * ncol + 1)];
  if (!m[nrl]) { throw "allocation failure 2 in matrix()"; }
  m[nrl] += 1;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i] = m[i-1] + ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

double * DoubleMatrix::myDoubleVector(const long nl, const long nh) const {
  double *v;
  
  //  v=(double *) malloc((size_t) ((nh - nl + 1 + 1) * sizeof(double)));
  v = new double [(nh - nl + 1 + 1)];
  if (!v) { throw "allocation failure in vector()"; }
  return v - nl + 1;
}

void DoubleMatrix::setMat(double **t, int r, int c) { 
  x = myDoubleMatrix(1, r, 1, c); 
  int i,j; 
  for(i=1; i<=r; i++)
    for (j=1; j<=c; j++)
      x[i][j] = t[i][j];
}


DoubleMatrix::DoubleMatrix(int r, int c)
  : rows(r), cols(c) {
  double ** t = myDoubleMatrix(1, r, 1, c);
  int i, j; 
  for(i=1; i<=r;i++)
    for(j=1; j<=c; j++)
      t[i][j] = double(0);
  setMat(t, r, c);
  freeMyDoubleMatrix(t, 1, r, 1, c);
}

// Makes diagonal square matrix out of W
DoubleMatrix::DoubleMatrix(const DoubleVector &v)
  : rows(v.displayEnd()), cols(v.displayEnd()) { 
  int e = v.displayEnd();
  double ** t = myDoubleMatrix(1, e, 1, e);
  int i, j;
  for (i=1; i<=e; i++)
    for (j=1; j<=e; j++)
      {	if (i==j) t[i][j] = v.display(i);
      if (i!=j) t[i][j] = double(0);
      }
  setMat(t, e, e); 
  freeMyDoubleMatrix(t, 1, e, 1, e);
}

const DoubleMatrix & DoubleMatrix::operator=(const DoubleMatrix &v) {
  if (this == &v) return *this;
  // Are the vectors of equivalent size?
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != v.rows || cols != v.cols)  {
    ostringstream ii;
    ii << "DoubleMatrix = overload; incompatible sizes\n";
    ii << *this << "=\n" << v;
    throw(ii.str());
  }
#endif
  freeMyDoubleMatrix(x, 1, rows, 1, cols);
  setMat(v.x, v.rows, v.cols);
  return *this;
}

const DoubleMatrix & DoubleMatrix::operator=(const double &v) {
  // Are the vectors of equivalent size?
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "DoubleMatrix = double overload; must be square\n";
    ii << *this << "=\n" << v;
    throw(ii.str());
  }
#endif
  freeMyDoubleMatrix(x, 1, rows, 1, cols);
  DoubleMatrix temp(rows, cols);
  int i; for (i=1; i<=rows; i++) temp(i, i) = v;
  setMat(temp.x, temp.rows, temp.cols);
  return *this;
}

DoubleMatrix DoubleMatrix::operator+(const DoubleMatrix &v) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != v.displayRows() || cols != v.displayCols()) {
    ostringstream ii;
    ii << "DoubleMatrix + overload; incompatible lengths\n";
    ii << *this << "+\n" << v;
    throw(ii.str());
  }
#endif
  DoubleMatrix temp(v);
  int i, j; 
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      temp(i, j) = v.display(i, j) + x[i][j];
  return temp;
}

DoubleMatrix DoubleMatrix::operator-(const DoubleMatrix &v) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != v.displayRows() || cols != v.displayCols()) {
    ostringstream ii;
    ii << "DoubleMatrix - overload; incompatible lengths\n";
    ii << *this << "-\n" << v;
    throw ii.str();
  }
#endif
  DoubleMatrix temp(v);
  int i, j; 
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      temp(i, j) = x[i][j] - v.display(i, j);
  return temp;
}

DoubleMatrix DoubleMatrix::operator+(double f) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols){
    ostringstream ii;
    ii << "DoubleMatrix + double overload; not square\n";
    ii << *this;
    throw ii.str();
  }
#endif
  DoubleMatrix temp(*this);
  int i; 
  for(i=1; i<=rows; i++) 
    temp(i, i) = x[i][i] + f;
  return temp;
}

DoubleMatrix DoubleMatrix::operator-(double f) const {
  DoubleMatrix temp(*this);
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "Subtracting double to non-square matrix\n" << temp;
    throw ii.str();
  }
#endif
  int i; 
  for (i=1; i<=rows; i++)
    temp(i, i) = temp(i, i) - f;
  return temp;  
}

DoubleMatrix DoubleMatrix::operator*(const DoubleMatrix &v) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (cols != v.displayRows()) {
    ostringstream ii;
    ii << "DoubleMatrix * overload; incompatible sizes\n" << *this 
       << "*\n" << v;
    throw ii.str();
  }
#endif
  DoubleMatrix temp(rows, v.displayCols());
  int i, j, k;
  for(i=1; i<=rows; i++) 
    for(j=1; j<=v.displayCols(); j++) {
      for (k=1; k<=cols; k++)
	temp(i, j) = temp(i, j) + x[i][k] * v.display(k, j);
    }
  return temp;
}

DoubleMatrix DoubleMatrix::operator*(const double f) const {
  DoubleMatrix temp(*this);
  int i, j;
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      temp(i, j) *= f;
  return temp;
}

void DoubleMatrix::operator*=(const double f) {

  int i, j;
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      x[i][j] *= f;
  return;
}

void DoubleMatrix::operator+=(DoubleMatrix & f) {

  int i, j;
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      x[i][j] += f.x[i][j];
  return;
}

double DoubleMatrix::trace() const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols)  {
    ostringstream ii;
    ii << "doublerace of non-square matrix\n" << *this;
    throw ii.str();
  }
#endif
  double sum = double(0);
  // this-> after +=
  int i; for(i=1; i<=rows; i++) sum += display(i, i);
  return sum;
}

DoubleMatrix DoubleMatrix::transpose() const {
  DoubleMatrix temp(*this);
  temp.rows = cols;
  temp.cols = rows;
  int i, j; 
  for (i=1; i<=rows; i++)
    for (j=1; j<=cols; j++)
      temp(j, i) = (*this).display(i, j);
  return temp;
}


void DoubleMatrix::swaprows(int i, int j) {
  DoubleMatrix temp(*this);
  for (int k=1; k<= cols; k++)
    {
      temp(i,k) = this->display(j,k);
      temp(j,k) = this->display(i,k);
    };
  *this = temp;
}

void DoubleMatrix::swapcols(int i,int j) {
  DoubleMatrix temp(*this);
  for (int k=1; k<= rows; k++)
    {
      temp(k, i) = this->display(k, j);
      temp(k, j) = this->display(k, i);
    };
  *this = temp;
}

// Perform on a mixing matrix O such that
// diag = O^double M O. It will put the diag evals in abs ascending order amd
// change O accordingly
void DoubleMatrix::associateOrderAbs(DoubleVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if ((displayRows() != v.displayEnd()) ||
      (displayCols() != v.displayEnd()) ||
      (v.displayStart() != 1))
    cout << "Associated ordering incompatibility" << endl;
#endif
  int i,j;
  for (i=v.displayStart(); i<=v.displayEnd(); i++)
    for (j=i+1; j<=v.displayEnd(); j++)
      if (fabs(v.display(i)) > fabs(v.display(j))) {
	v.swap(i,j);
	swapcols(i,j);
      }
}

// Perform such that
// diag = U^T double M V. It will put the diag evals in abs ascending order amd
// change U,V accordingly
void DoubleMatrix::associateOrderAbs(DoubleMatrix & u, DoubleMatrix & v, DoubleVector &w) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if ((u.displayRows() != w.displayEnd()) ||
      (u.displayCols() != w.displayEnd()) ||
      (w.displayStart() != 1))
    cout << "Associated ordering incompatibility" << endl;
#endif
  int i,j;
  for (i=w.displayStart(); i<=w.displayEnd(); i++)
    for (j=i+1; j<=w.displayEnd(); j++)
      if (fabs(w.display(i)) > fabs(w.display(j))) {
	w.swap(i, j);
	v.swapcols(i, j);
	u.swapcols(i, j);
      }
}


// Sum of all element's magnitudes in the matrix
double DoubleMatrix::sumElements() const {
  double temp = 0.0;
  int i, j;
  
  for (i=1; i<=cols; i++)
    for (j=1; j<=rows; j++)
      temp = temp + fabs(display(j, i));
  return temp;
}

void DoubleMatrix::symmetrise() {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "Error: symmetrising rectangular matrix " << *this;
    throw ii.str();
  }
#endif
  int i, j;
  for (i=2; i<=rows; i++)
    for (j=1; j<i; j++)
      x[i][j] = x[j][i];
}

// Gives sum of difference between two matrices
double DoubleMatrix::compare(const DoubleMatrix & a) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != a.displayRows() || displayCols() !=
      a.displayCols()) {
    ostringstream ii;
    ii << "Error: comparing matrices of different sizes" << *this << 
      " and " << a;
    throw ii.str();
  }
#endif
  double c = double(0);
  int i, j; 
  for (i=1; i<=displayRows(); i++)
    for (j=1; j<=displayCols(); j++)
      c = c + double(fabs(a.display(i, j) - display(i, j)));
  return c;
}

// Gives sum of difference between two matrices
double DoubleMatrix::compare(const ComplexMatrix & a) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != a.displayRows() || displayCols() !=
      a.displayCols()) {
    ostringstream ii;
    ii << "Error: comparing matrices of different sizes" << *this << 
      " and " << a;
    throw ii.str();
  }
#endif
  double c = double(0);
  int i, j; 
  for (i=1; i<=displayRows(); i++)
    for (j=1; j<=displayCols(); j++)
      c = c + double(fabs(a.display(i, j).mod() - display(i, j)));
  return c;
}

//Returns diagonal values ONLY of a matrix
DoubleVector DoubleMatrix::diagVals() const {
  DoubleVector temp(1, displayCols()), W(1, displayCols());
  DoubleMatrix U(displayRows(), displayCols()), 
    V(displayRows(), displayCols());
  diagonalise(U, V, temp);
  return temp;
}

// Returns inverse of a matrix
DoubleMatrix DoubleMatrix::inverse() const {
  DoubleMatrix ans(this -> display());
#ifdef ARRAY_BOUNDS_CHECKING
  if (ans.displayCols() != ans.displayRows()) { 
    ostringstream ii;
    ii << "Error in DoubleMatrix::invert, trying to invert non-square";
    ii << " matrix " << *this << endl;
    throw ii.str();
  }
#endif
  
  DoubleMatrix u(ans.displayRows(), ans.displayRows()), 
    v(ans.displayRows(), ans.displayRows());
  DoubleVector w(ans.displayRows());
  ans.diagonalise(u, v, w);
  
  const double underflow = 1.0e-16;
  int i; for (i=1; i<=ans.displayRows(); i++) {
    if (w(i) > underflow) w(i) = 1 / w(i); // Leave zero eigenvalues alone
    else {
      ostringstream ii;
      ii << "Problem in DoubleMatrix::inverse, " << *this 
	 << "non-invertible" << endl;
      throw ii.str();
    }
  }
  DoubleMatrix temp(w);
  
  return (v * temp * u.transpose());
}

// A = U.W.V^T where W is a matrix of the eigenvalues
// therefore W = U^T A V. Uses SVD but is not yet set up for really
// singular matrices
// You should really set up more/faster routines for symmetric case etc
double DoubleMatrix::diagonalise(DoubleMatrix & u, DoubleMatrix & v, 
				 DoubleVector & w) const {
#ifdef ARRAY_BOUNDS_CHECKING
  int c = displayCols();
  if (displayRows() != c || w.displayEnd() != c || v.displayCols() != c
      || v.displayRows() !=c || u.displayCols() !=c || u.displayRows()
      !=c) {
    ostringstream ii;
    ii << "Error: Trying to diagonalise matrix \n" << *this 
       << "with u" << u << "v " << v << "w " << w;
    throw ii.str();
  }
#endif
  
  // Numerical routine replaces argument, so make a copy of elements
  u = *this;
  diagonaliseSvd(u, w, v);
  
  // Phase freedom - tends to give more familiar matrices
  u = -1.0 * u; v = -1.0 * v;
  
  associateOrderAbs(u, v, w);
  // Calculate residue: difference in two matrices, divided by absolute value
  // of largest eigenvalue
  DoubleMatrix temp(w);
  return this->compare(u * temp * v.transpose()) / w.max();
}

// For SYMMETRIC MATRICES ONLY!
// A = V.W.V^T where W is a matrix of the eigenvalue therefore W =
// V^T A V. 
double DoubleMatrix::diagonaliseSym(DoubleMatrix & v, DoubleVector & w) const {
  int c = displayCols();
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != c || w.displayEnd() != c || v.displayCols() != c
      || v.displayRows() !=c || c > 10) {
    ostringstream ii;
    ii << "Error: Trying to diagonalise matrix \n" << *this;
    throw ii.str();
  }
#endif
  // Numerical recipes routine replaces argument, so make a copy of elements
  DoubleMatrix a(*this);
  int nrot;
  diagonaliseJac(a, c, w, v, &nrot);
  
  // So Far, V is such that the eigenvalues are in some random order.
  // We must now re-order the rows of V to get the eigenvalues in the
  // correct order (taken to be increasing fabs(w)).
  v.associateOrderAbs(w);
  
  // Calculate residue: difference in two matrices, divided by absolute
  // value of largest eigenvalue
  DoubleMatrix temp(w); 
  return this->compare(v * temp * v.transpose()) / w.max();
}

// For SYMMETRIC MATRICES ONLY!
// A = V.W.V^double where W is a matrix of the eigenvalue therefore W =
// V^double A V. all diagonal values are made positive!
double DoubleMatrix::diagonaliseSym
(ComplexMatrix & v, DoubleVector & w) const {
  int c = displayCols();
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != c || w.displayEnd() != c || v.displayCols() != c
      || v.displayRows() !=c || c > 10) {
    ostringstream ii;
    ii << "Error: Trying to diagonalise matrix \n" << *this;
    throw ii.str();
  }
#endif
  // Numerical recipes routine replaces argument, so make a copy of elements
  DoubleMatrix a(*this), k(c, c);
  int nrot;
  diagonaliseJac(a, c, w, k, &nrot);
  
  // So Far, V is such that the eigenvalues are in some random order.
  // We must now re-order the rows of V to get the eigenvalues in the
  // correct order (taken to be increasing fabs(w)).
  k.associateOrderAbs(w);
  
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  ComplexMatrix m(c, c);
  
  int i, j; for (i=1; i<=c; i++) 
    if (w(i) < 0.0) {
      for (j=1; j<=c; j++) v(i, j) = Complex(0.0, 1.0) * k(i, j);
      w(i) = - w(i);
    }
    else {
      for (j=1; j<=c; j++) v(i, j) = Complex(k(i, j), 0.0);
    }
  
  // Calculate residue: difference in two matrices, divided by absolute
  // value of largest eigenvalue
  ComplexMatrix temp(c, c); 
  for (i=1; i<=c; i++) temp(i, i) = w(i);
  cout << v * temp * v.transpose();
  return this->compare(v.transpose() * temp * v) / w.max();
}

// Special case that's often used: eigenvalues of 2 by 2 symmetric matrix
// A(2,1) is assumed to be equal to A(1,2). theta is the mixing angle:
// [ cos theta   sin theta   ]    A    [cos theta   -sin theta ]  = [ m1 0 ]
// [ -sin theta  cos theta   ]         [sin theta    cos theta ]    [ 0 m2 ]
// where m1 is NOT NECESSARILY > m2
DoubleVector DoubleMatrix::sym2by2(double & theta) const  {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() !=2 || displayCols() !=2) {
    ostringstream ii;
    ii << "Called sym2by2 with matrix of dimension ("
       << displayRows() << ", " << displayCols() << ")\n";
    throw ii.str();
  }
  if (display(1, 2) != display(2, 1)) {
    ostringstream ii;
    ii << "Called non symmetric sym2by2:" << *this;
    throw ii.str();
  }
#endif
  
  DoubleVector temp(1, 2);
  double sumTol = fabs(fabs(display(1, 1)) - fabs(display(2, 2))) /
    maximum(fabs(display(1, 1)), fabs(display(2, 2)));
  
  double t2t = 2.0 * display(1, 2) / (display(1, 1) - display(2, 2));
  
  if (sumTol > EPSTOL || t2t < EPSTOL) theta = atan(t2t) * 0.5; 
  else {
    if (t2t > 0.0) 
      theta = PI * 0.25;
    else 
      theta = -PI * 0.25;
  }
  
  DoubleMatrix mm(rot2d(theta) * (*this) * rot2d(theta).transpose());
  
  temp(1) = mm(1, 1);
  temp(2) = mm(2, 2);
  
  double maxtol = maximum(fabs(temp(1)), fabs(temp(2))) * TOLERANCE * 1.0e3;
  if (fabs(mm(2, 1)) > maxtol || fabs(mm(1, 2)) > maxtol) {
    ostringstream ii;
    ii << "Inaccurate diagonalisation in LINALG::sym2by2\n";
    ii << "Diagonalising " << * this;
    ii << "Found diagonalised matrix to be " << mm;
    throw ii.str();
  }
  
  return temp; 
}

DoubleMatrix DoubleMatrix::apply(double (*fn)(double)) {
  DoubleMatrix temp(displayRows(), displayCols());
  int i, j;
  for (i=1; i<=displayRows(); i++)
    for (j=1; j<=displayCols(); j++)
      temp(i, j) = fn(display(i, j));
  return temp;
}

// 2 by 2 asymmetric matrices
// [ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ]  = diag
// [ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
DoubleVector DoubleMatrix::asy2by2(double & thetaL, double & thetaR) const {
  DoubleVector temp1(2), temp2(2);
  DoubleMatrix e(*this);
  
  DoubleMatrix h = e.transpose() * e;
  temp1 = (h.sym2by2(thetaR)).apply(sqrt);
  if (temp1(1) > temp1(2))
    thetaR = thetaR + PI * 0.5;
  
  h = e * e.transpose();
  temp2 = (h.sym2by2(thetaL)).apply(sqrt);
  
  // Did the eigenvalues come out in the wrong order? 
  // If so, swap the order of the second case by adding pi/2 to the angle.
  if (temp2(1) > temp2(2)) {
    thetaL = thetaL + PI * 0.5;
  }
  
  DoubleMatrix mm(2, 2);
  mm = rot2d(thetaL) * e * rot2d(-thetaR);
  
  double maxtol = maximum(fabs(temp1(1)), fabs(temp1(2))) * TOLERANCE;
  if (fabs(mm(2, 1)) > maxtol || fabs(mm(1, 2)) > maxtol) {
    outputCharacteristics(16);
    ostringstream ii;
    ii << "Inaccurate diagonalisation in LINALG::asy2by2\n";
    ii << "temp1=" << temp1 << " temp2=" << temp2 << endl;
    ii << "diagonalising " << e << " where thetaR=" << thetaR 
       << " thetaL=" << thetaL << endl;
    ii << "m=" << mm;
    throw ii.str();
  }
  
  temp1(1) = mm(1, 1); temp1(2) = mm(2, 2);
  
  return temp1;
}

double DoubleMatrix::min(int & k, int & l) const {
  double m = display(1, 1);
  int i, j; 
  for (i=1; i<=displayRows();  i++)
    for (j=1; j<=displayCols(); j++)
      if (display(i, j) < m) {
	m = display(i, j);
	k = i; l = j;
      }
  return m;
}

double DoubleMatrix::max(int & k, int & l) const {
  double m = display(1, 1);
  int i, j; 
  for (i=1; i<=displayRows();  i++)
    for (j=1; j<=displayCols(); j++)
      if (display(i, j) > m) {
	m = display(i, j);
	k = i; l = j;
      }
  return m;
}

// You have to free memory for temp before this: dumps matrix in temp
void DoubleMatrix::displayMat(double ** temp) const {
  int i, j; 
  for(i=1; i<=displayRows(); i++)
    for(j=1; j<=displayCols(); j++)
      temp[i][j] = display(i, j);
}

// Assumes the double is multiplied by identity
DoubleMatrix operator*(double f, const DoubleMatrix &m) {
  DoubleMatrix temp(m);
  int i, j; 
  for (i=1; i<=temp.displayRows(); i++)
    for (j=1; j<=temp.displayCols(); j++)
      temp(i, j) = temp(i, j) * double(f);
  return temp;
}

double DoubleMatrix::nmin(int & k, int & l) const {
  double m1 = display(1, 1), m2;
  int k1 = 1, l1 = 1, k2, l2;

  if (display(1, 2) < m1){
    m2 = m1;
    k2 = k1;
    l2 = l1;
    
    m1 = display(1, 2);
    k1 = 1;
    l1 = 2; 
  }
  else{
    m2 = display(1, 2);
    k2 = 1;
    l2 = 2; 
  }
  int i, j;
  for (i=1; i<=displayRows();  i++)
    for (j=1; j<=displayCols(); j++){
      if (display(i, j) < m1) {
	m2 = m1;
	k2 = k1;
	l2 = l1;

	m1 = display(i, j);
	k1 = i; l1 = j;
      }
      else if (display(i, j) < m2){
	m2 = display(i, j);
	k2 = i; l2 = j;
      }
    }
  k = k2; l = l2;  
  return m2;
}

// v_i = M_ij v_j
DoubleVector DoubleMatrix::operator*(const DoubleVector & v) const {
  DoubleMatrix m(*this);
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayCols() != v.displayEnd() || v.displayStart() != 1) {
    ostringstream ii;
    ii << "matrix times vector overload size problem\n" << m << "*" << v;
    throw ii.str();
  }
#endif
  
  DoubleVector temp(m.displayRows());
  int i, j;
  for(j=1; j<=m.displayRows(); j++) {
    temp(j) = double(0);
    for(i=1; i<=m.displayCols(); i++)
      temp(j) += m.display(j, i) * v.display(i);
  }
  return temp;
}

// v_i = v_j M_{ji}
DoubleVector operator*(const DoubleVector & v, const DoubleMatrix & m) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayCols() != v.displayEnd() || v.displayStart() != 1) {
    ostringstream ii;
    ii << "DoubleMatrix times vector overload size problem\n" << m << "*" << v;
    throw ii.str();
  }
#endif
  DoubleVector temp(m.displayRows());
  int i, j;
  for(j=1; j<=m.displayRows(); j++) {
    temp(j) = double(0);
    for(i=1; i<=m.displayCols(); i++)
      temp(j) += v.display(i) * m.display(i, j);
  }
  return temp;
}

DoubleMatrix operator-(double f, const DoubleMatrix &m) {
  DoubleMatrix temp(m);
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayRows() != m.displayCols()) {
    ostringstream ii;
    ii << "Subtracting double from non-square matrix\n" << m;
    throw ii.str();
  }
#endif
  int i, j; 
  for (i=1; i<=temp.displayRows(); i++)
    for (j=1; j<=temp.displayCols(); j++) {
      if (i != j) temp(i, j) = - temp(i, j);
      if (i == j) temp(i, j) = f - temp(i, j);
    }
  return temp;
}

DoubleMatrix operator+(double f, const DoubleMatrix &m) {
  DoubleMatrix temp(m);
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayRows() != m.displayCols()) {
    ostringstream ii;
    ii << "Subtracting double from non-square matrix\n" << m;
    throw ii.str();
  }
#endif
  int i; 
  for (i=1; i<=temp.displayRows(); i++)
    temp(i, i) = temp(i, i) + double(f);
  return temp;
}

// MIj = v1I v2J
DoubleMatrix outerProduct(const DoubleVector & v1, const DoubleVector & v2) {
  DoubleMatrix temp(v1.displayEnd(), v2.displayEnd());
  int i,j;
  for (i=v1.displayStart(); i<=v1.displayEnd(); i++)
    for (j=v2.displayStart(); j<=v2.displayEnd(); j++)
      temp(i, j) = v1.display(i) * v2.display(j);
  return temp;
}


istream & operator>>(istream & left, DoubleMatrix &m) {
  int i, j;
  DoubleMatrix x(m.displayRows(), m.displayCols());
  for (i=1; i<=m.displayRows(); i++)
    for (j=1; j<=m.displayCols(); j++)
      left >> x(i, j);
  m = x;
  
  return left;
}

ostream & operator <<(ostream &left, const DoubleMatrix &v) {
  
  const static double underflow = 1.0e-120;
  double x;
  
  left << "(" << v.displayRows() << "," <<
    v.displayCols() << "):\n";
  int i, j; 
  for(i=1; i<=v.displayRows(); i++) {
    for(j=1; j<=v.displayCols(); j++) {
      x = v.display(i, j);
      if (fabs(x) < underflow) x = 0.0; // Traps -0.0
      if (x >= 0.0) left << " "; 
      left << x << " ";
    }
    left << '\n';
  }
  return left;
}

// 2 dimensional rotation matrix M
// Returns U = [ cos theta    sin theta ]
//             [-sin theta    cos theta ]
DoubleMatrix rot2d(double theta) {
  DoubleMatrix u(2, 2);
  
  u(1, 1) = cos(theta); u(2, 2) = u(1, 1);
  u(1, 2) = sin(theta); u(2, 1) = - u(1, 2);
  return u;
}

// Returns mixing matrix
// [ -sin(theta)  cos(theta) ]
// [  cos(theta)  sin(theta) ] --
// takes into account different conventions for the ordering of a.
DoubleMatrix rot2dTwist(double theta) {
  DoubleMatrix n(2, 2);
  n(1, 1) = -sin(theta);
  n(1, 2) = cos(theta);
  n(2, 1) = cos(theta); 
  n(2, 2) = sin(theta);
  
  return n;
}

// Redefines mixing matrices to be complex such that diagonal values are
// positive for a 2 by 2: if  
// [ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ]  = diag
// [ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
// as given by asy2by2!
// then u^* A v^+ = mdiagpositive
// IT'S BEEN CHECKED!
void positivise(double thetaL, double thetaR, const DoubleVector & diag,
		ComplexMatrix & u, ComplexMatrix & v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (diag.displayStart()!=1 || diag.displayEnd()!=2 ||
      u.displayRows()!=2 || u.displayCols()!=2 ||
      v.displayRows()!=2 || v.displayCols()!=2) {
    ostringstream ii;
    ii << "DoubleMatrix::positivise currentl only available for 2 by 2"
       << "matrices, not " << diag << u << v;
    throw ii.str();
  }
#endif
  
  ComplexMatrix a(2, 2);
  int i; for (i=1; i<=2; i++) if (diag.display(i) < 0.0) 
    a(i, i) = Complex(0., 1.); 
  else a(i, i) = Complex(1.0, 0.0);
  
  u = a * rot2d(thetaL);
  v = a * rot2d(thetaR);
  return;
}

// v_i = M_ij v_j
ComplexVector ComplexMatrix::operator*(const DoubleVector & v) {
  ComplexMatrix m(*this);
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayCols() != v.displayEnd() || v.displayStart() != 1) {
    ostringstream ii;
    ii << "matrix times vector overload size problem\n" << m << "*" << v;
    throw ii.str();
  }
#endif
  
  ComplexVector temp(m.displayRows());
  int i, j;
  for(j=1; j<=m.displayRows(); j++) 
    for(i=1; i<=m.displayCols(); i++)
      temp(j) = temp(j) + m.display(j, i) * v.display(i);
  
  return temp;
}



// ---------------- diagonalisation algorithms --------------------
void diagonaliseSvd(DoubleMatrix & a, DoubleVector & w, DoubleMatrix & v)
{
  int n = a.displayCols(), m = a.displayRows();
  int flag, i, its, j, jj, k, l = 0, nm = 0; 
  double anorm, c, f, g, h, s, scale, x, y, z; 
  
  DoubleVector rv1(a.displayCols());
  
  g = scale = anorm = 0.0; 
  for (i = 1; i<= n; i++) {
    l = i + 1; 
    rv1(i) = scale * g; 
    g = s = scale = 0.0; 
    if (i <=  m) {
      for (k = i; k <= m; k++) scale += fabs(a(k, i)); 
      if (scale) {
	for (k = i; k <= m; k++) {
	  a(k, i) /=  scale; 
	  s +=  a(k, i) * a(k, i); 
	}
	f = a(i, i); 
	g = -sign(sqrt(s), f); 
	h = f * g - s; 
	a(i, i) = f-g; 
	for (j = l; j <= n; j++) {
	  for (s = 0.0, k = i; k <= m; k++) s +=  a(k, i) * a(k, j); 
	  f = s / h; 
	  for (k = i; k <= m; k++) a(k, j) +=  f * a(k, i); 
	}
	for (k = i; k <= m; k++) a(k, i)  *=  scale; 
      }
    }
    w(i) = scale * g; 
    g = s = scale = 0.0; 
    if (i  <=  m && i !=  n) {
      for (k = l; k <= n; k++) scale += fabs(a(i, k)); 
      if (scale) {
	for (k = l; k <= n; k++) {
	  a(i, k) /=  scale; 
	  s +=  a(i, k) * a(i, k); 
	}
	f = a(i, l); 
	g  =  -sign(sqrt(s), f); 
	h = f * g-s; 
	a(i, l) = f - g; 
	for (k = l; k <= n; k++) rv1(k) = a(i, k) / h; 
	for (j = l; j <= m; j++) {
	  for (s = 0.0, k = l; k <= n; k++) s +=  a(j, k) * a(i, k); 
	  for (k = l; k <= n; k++) a(j, k) +=  s * rv1(k); 
	}
	for (k = l; k <= n; k++) a(i, k)  *=  scale; 
      }
    }
    anorm = maximum(anorm, (fabs(w(i)) + fabs(rv1(i)))); 
  }
  for (i = n; i >= 1; i--) {
    if (i < n) {
      if (g) {
	for (j = l; j <= n; j++)
	  v(j, i) = (a(i, j) / a(i, l)) / g; 
	for (j = l; j <= n; j++) {
	  for (s = 0.0, k = l; k <= n; k++) s += a(i, k) * v(k, j); 
	  for (k = l; k <= n; k++) v(k, j) += s * v(k, i); 
	}
      }
      for (j = l; j <= n; j++) v(i, j) = v(j, i) = 0.0; 
    }
    v(i, i) = 1.0; 
    g = rv1(i); 
    l = i; 
  }
  for (i = minimum(m, n); i >= 1; i--) {
    l = i + 1; 
    g = w(i); 
    for (j = l; j <= n; j++) a(i, j) = 0.0; 
    if (g) {
      g = 1.0 / g; 
      for (j = l; j <= n; j++) {
	for (s = 0.0, k = l; k <= m; k++) s +=  a(k, i) * a(k, j); 
	f = (s / a(i, i)) * g; 
	for (k = i; k <= m; k++) a(k, j) +=  f * a(k, i); 
      }
      for (j = i; j <= m; j++) a(j, i) *=  g; 
    } else for (j = i; j <= m; j++) a(j, i) = 0.0; 
    ++a(i, i); 
  }
  for (k = n; k >= 1; k--) {
    for (its = 1; its <= 30; its++) {
      flag = 1; 
      for (l = k; l >= 1; l--) {
	nm = l - 1; 
	if ((double)(fabs(rv1(l))+anorm)  ==  anorm) {
	  flag = 0; 
	  break; 
	}
	if ((double)(fabs(w(nm))+anorm)  ==  anorm) break; 
      }
      if (flag) {
	c = 0.0; 
	s = 1.0; 
	for (i = l; i <= k; i++) {
	  f = s * rv1(i); 
	  rv1(i) = c * rv1(i); 
	  if ((double)(fabs(f)+anorm)  ==  anorm) break; 
	  g = w(i); 
	  h = pythagoras(f, g); 
	  w(i) = h; 
	  h = 1.0 / h; 
	  c = g * h; 
	  s  =  -f * h; 
	  for (j = 1; j <= m; j++) {
	    y = a(j, nm); 
	    z = a(j, i); 
	    a(j, nm) = y * c + z * s; 
	    a(j, i) = z * c - y * s; 
	  }
	}
      }
      z = w(k); 
      if (l  ==  k) {
	if (z < 0.0) {
	  w(k) = -z; 
	  for (j = 1; j <= n; j++) v(j, k) = -v(j, k); 
	}
	break; 
      }
      if (its  ==  30) 
	throw "no convergence in 30 diagonaliseSvd iterations"; 
      
      x = w(l); 
      nm = k - 1; 
      y = w(nm); 
      g = rv1(nm); 
      h = rv1(k); 
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y); 
      g = pythagoras(f, 1.0); 
      f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x; 
      c = s = 1.0; 
      for (j = l; j <= nm; j++) {
	i = j+1; 
	g = rv1(i); 
	y = w(i); 
	h = s * g; 
	g = c * g; 
	z = pythagoras(f, h); 
	rv1(j) = z; 
	c = f / z; 
	s = h / z; 
	f = x * c + g * s; 
	g  =  g * c - x * s; 
	h = y * s; 
	y  *=  c; 
	for (jj = 1; jj <= n; jj++) {
	  x = v(jj, j); 
	  z = v(jj, i); 
	  v(jj, j) = x * c + z * s; 
	  v(jj, i) = z * c - x * s; 
	}
	z = pythagoras(f, h); 
	w(j) = z; 
	if (z) {
	  z = 1.0 / z; 
	  c = f * z; 
	  s = h * z; 
	}
	f = c * g+s * y; 
	x = c * y-s * g; 
	for (jj = 1; jj <= m; jj++) {
	  y = a(jj, j); 
	  z = a(jj, i); 
	  a(jj, j) = y * c + z * s; 
	  a(jj, i) = z * c - y * s; 
	}
      }
      rv1(l) = 0.0; 
      rv1(k) = f; 
      w(k) = x; 
    }
  }
}

double pythagoras(double a,  double b) {
  double absa, absb; 
  absa = fabs(a); 
  absb = fabs(b); 
  if (absa > absb) return absa  *  sqrt(1.0 + sqr(absb / absa)); 
  else return (absb  == 0.0 ? 0.0 : absb * 
	       sqrt(1.0 + sqr(absa / absb))); 
}

inline void rotate(DoubleMatrix & a, int i, int j, int k, int l, double s, 
		   double tau) {
  double g = a(i, j);
  double h = a(k, l);
  a(i, j) = g - s * (h + g * tau);
  a(k, l) = h + s * (g - h * tau);
}

void diagonaliseJac(DoubleMatrix & a,  int n,  DoubleVector & d,  
		    DoubleMatrix & v,  int *nrot) {
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c; 
  
  DoubleVector b(n), z(n); 
  
  for (ip=1; ip<=n; ip++) {
    for (iq=1; iq<=n; iq++) v(ip,  iq) = 0.0; 
    v(ip,  ip) = 1.0; 
  }
  for (ip=1; ip<=n; ip++) {
    b(ip) = d(ip) = a(ip,  ip); 
    z(ip) = 0.0; 
  }
  *nrot = 0; 
  for (i=1; i<=50; i++) {
    sm = 0.0; 
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++)
	sm += fabs(a(ip,  iq)); 
    }
    if (sm == 0.0) return; 
    
    if (i < 4)
      tresh = 0.2  *sm/(n*n); 
    else
      tresh = 0.0; 
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++) {
	g=100.0 * fabs(a(ip,  iq)); 
	if (i > 4 && (double)(fabs(d(ip)) + g) == (double) fabs(d(ip))
	    && (double)(fabs(d(iq)) + g) == (double) fabs(d(iq)))
	  a(ip,  iq) = 0.0; 
	else if (fabs(a(ip,  iq)) > tresh) {
	  h = d(iq) - d(ip); 
	  if ((double)(fabs(h) + g) == (double) fabs(h))
	    t=(a(ip,  iq)) / h; 
	  else {
	    theta = 0.5 * h / (a(ip,  iq)); 
	    t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta)); 
	    if (theta < 0.0) t = -t; 
	  }
	  c = 1.0 / sqrt(1 + t * t); 
	  s = t * c; 
	  tau = s / (1.0 + c); 
	  h = t * a(ip,  iq); 
	  z(ip) -= h; 
	  z(iq) += h; 
	  d(ip) -= h; 
	  d(iq) += h; 
	  a(ip,  iq) = 0.0; 
	  for (j=1; j<=ip-1; j++) rotate(a, j, ip, j, iq, s, tau); 
	  for (j=ip+1; j<=iq-1; j++) rotate(a, ip, j, j, iq, s, tau);
	  for (j=iq+1; j<=n; j++) rotate(a, ip, j, iq, j, s, tau);
	  for (j=1; j<=n; j++) rotate(v, j, ip, j, iq, s, tau);
	  ++(*nrot); 
	}
      }
    }
    for (ip=1; ip<=n; ip++) {
      b(ip) += z(ip); 
      d(ip) = b(ip); 
      z(ip) = 0.0; 
    }
  }
  if (PRINTOUT) {
    cout << "Too many iterations in routine diagonaliseJac: diagonalising\n" 
	 << a << "to" << d << " with " << v << flush;
  }
}


Complex * ComplexVector::myComplexVector(const long nl, const long nh) {
  Complex *v;
  
  //  v=(Complex *) malloc((size_t) ((nh - nl + 2) * sizeof(Complex)));
  v = new Complex [(nh - nl + 2)];
  if (!v) throw "Allocation failure in vector()"; 
  
  return v - nl + 1;
}

void ComplexVector::setVec(const Complex *t, long st, long nd) { 
  x = myComplexVector(st, nd); int i; 
  for(i = st; i <= nd; i++) x[i] = t[i];
}

ComplexVector::ComplexVector(int e)
  : start(1), end(e) { 
  Complex * t = myComplexVector(1, e);
  int i; for(i=1; i<=e; i++) t[i] = Complex(0);
  setVec(t, 1, e);
  freeMyComplexVector(t, 1, e);
}

ComplexVector::ComplexVector(int s, int e)
  : start(s), end(e) {
  Complex * t = myComplexVector(s, e);
  int i; for(i=s; i<=e; i++) t[i] = Complex(0);
  setVec(t, s, e);
  freeMyComplexVector(t, s, e);
}

// Changes the length of a vector - copies as many elements of old one as
// possible, and fills any extra up with zeroes
void ComplexVector::setEnd(int e) {
  if (e < 1) {
    ostringstream ii;
    ii << "ComplexVector.setEnd called with incompatible length " << e <<
      endl; throw ii.str();
  }
  Complex * t = myComplexVector(start, e);
  int i; for(i=start; i<=e; i++) {
    if (i <= end) t[i] = x[i]; 
    else
      t[i] = Complex(0);
  }
  freeMyComplexVector(x, start, end);
  setVec(t, start, e);
  freeMyComplexVector(t, start, e);
  
  end = e;
}

ComplexVector & ComplexVector::operator=(const ComplexVector &v) {
  if (this == &v) return *this;
  // Are the vectors of equivalent size?
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.start || end != v.end) {
    ostringstream ii;
    ii << "ComplexVector = overload; incompatible lengths\n";
    ii << *this << "=\n" << v;
    throw ii.str();
  }
#endif
  freeMyComplexVector(x, start, end);
  setVec(v.x, v.start, v.end);
  
  return *this;
}

ComplexVector operator*(DoubleMatrix & m, ComplexVector & v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayCols() != v.displayEnd()) {
    ostringstream ii;
    ii << "Can't multiply " << m << "*" << v << ": incompatible sizes.";
    throw ii.str();
  }
#endif
  ComplexVector temp(m.displayCols()); 
  
  int i, j;
  for(j=1; j<=m.displayRows(); j++) 
    for(i=1; i<=m.displayCols(); i++)
      temp(j) = temp(j) + m.display(j, i) * v.display(i);
  
  return temp;
}

ComplexVector ComplexVector::operator+(const ComplexVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.displayStart() || end != v.displayEnd()) {
    ostringstream ii;
    ii << "ComplexVector + overload; incompatible lengths\n";
    ii << *this << "+\n" << v;
    throw ii.str();
  }
#endif
  ComplexVector temp(v);
  int i; for(i=start; i<=end; i++) temp(i) = v.display(i) + x[i];
  return temp;
}

ComplexVector ComplexVector::operator-(const ComplexVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.displayStart() || end != v.displayEnd()) {     
    ostringstream ii;
    ii << "ComplexVector - overload; incompatible lengths\n";
    ii << *this << "-\n" << v;
    throw ii.str();
  }
#endif
  ComplexVector temp(v);
  int i; for(i=start; i<=end; i++) temp(i) = x[i] - v.display(i);
  return temp;
}

// NOComplex dot product, but product of elements.
ComplexVector ComplexVector::operator*(const ComplexVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.displayStart() || end != v.displayEnd()) {
    ostringstream ii;
    ii << "ComplexVector * overload; incompatible lengths\n";
    ii << *this << "*\n" << v;
    throw ii.str();
  }
#endif
  ComplexVector temp(v);
  int i; for(i=start; i<=end; i++) temp(i) = v.display(i) * x[i];
  return temp;
}

ComplexVector ComplexVector::operator*(Complex f) {
  ComplexVector temp(*this);
  int i; for(i=start; i<=end; i++) temp(i) = temp(i) * f;
  return temp;
}

Complex ComplexVector::dot(const ComplexVector & v) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (v.displayStart() != start || v.displayEnd() != end) {
    ostringstream ii;
    ii << "Scalar product between incompatible vectors " << *this << 
      " and " << v;
    throw ii.str();
  }
#endif
  Complex f;
  int i; for (i=start; i<=end; i++) f = f + x[i] * v.display(i);
  return f;				   
}

ComplexVector ComplexVector::apply(Complex (*fn)(Complex)) {
  ComplexVector temp(displayStart(), displayEnd());
  int i;
  for (i=displayStart(); i<=displayEnd(); i++)
    temp(i) = Complex(fn(display(i)));
  return temp;
}

Complex ComplexVector::max() const {
  Complex m = Complex(0);
  int i; 
  for (i=displayStart(); i<=displayEnd();  i++)
    if (display(i).mod() > m.mod()) m = display(i);
  return m;
}

// Smallest absolute value of element
Complex ComplexVector::min(int & p) const {
  Complex m = display(displayStart());
  p = displayStart();
  int i; 
  for (i=displayStart() + 1; i<=displayEnd();  i++)
    if (display(i).mod() < m.mod()) 
      {
	m = display(i);
	p = i;
      }
  return m;
}

void ComplexVector::swap(int i, int j) {
  Complex m;
  ComplexVector temp(*this);
  m = temp.display(j);
  temp(j) = temp.display(i);
  temp(i) = m;
  *this = temp; 
}

ComplexVector operator*(Complex f, const ComplexVector & v) {
  ComplexVector temp(v);
  int i;
  for (i=v.displayStart(); i<=v.displayEnd(); i++)
    temp(i) = temp(i) * f;
  return temp;
} 

/*ComplexVector ComplexVector::operator*(double f) {
  DoubleVector temp(*this);
  int i; for(i=start; i<=end; i++) temp(i) = temp(i) * Complex(f);
  return temp;
  }*/


ostream & operator <<(ostream &left, const ComplexVector &v) {
  const int NUMDisplay = 5;
  left << "(" << v.displayStart() << "," << v.displayEnd() << "):\n";
  if (v.displayEnd() > NUMDisplay) 
    left << '\n';
  else
    left << " ";
  int i; for(i=v.displayStart(); i<=v.displayEnd(); i++) {
    if ((i % NUMDisplay) == 1) left << "(" << i << ")=";
    left << v.display(i);
    if ((v.displayEnd() % NUMDisplay) == 0 || (i == v.displayEnd())) 
      left << '\n';
    else
      left << " ";
  }
  return left;
}

istream & operator>>(istream & left, ComplexVector &v) {
  int i;
  ComplexVector x(v.displayStart(), v.displayEnd());
  for (i=v.displayStart(); i<=v.displayEnd(); i++)
    left >> x(i);
  v = x;
  
  return left;
}

/* -----------------------------------------------------*/
/* allocate a Complex matrix with subscript range m[nrl..nrh][ncl..nch] */
Complex ** ComplexMatrix::myComplexMatrix(long nrl, long nrh, long ncl, long nch) const {
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  Complex **m;
  
  /* allocate pointers to rows */
  //  m=(Complex **) malloc((size_t)((nrow + 1) * sizeof(Complex*)));
  m = new Complex* [nrow + 1]; 
  if (!m) { throw "allocation failure 1 in matrix()"; }
  m += 1;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  //  m[nrl]=(Complex *) 
  // malloc((size_t)((nrow * ncol + 1) * sizeof(Complex)));
  m[nrl] = new Complex [(nrow * ncol + 1)];
  if (!m[nrl]) { throw "allocation failure 2 in matrix()"; }
  m[nrl] += 1;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i] = m[i-1] + ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

Complex * ComplexMatrix::myComplexVector(const long nl, const long nh) const {
  Complex *v;
  
  //  v=(Complex *) malloc((size_t) ((nh - nl + 1 + 1) * sizeof(Complex)));
  v = new Complex [(nh - nl + 1 + 1)];
  if (!v) { throw "allocation failure in vector()"; }
  return v - nl + 1;
}

void ComplexMatrix::setMat(Complex **t, int r, int c) { 
  x = myComplexMatrix(1, r, 1, c); 
  int i,j; 
  for(i=1; i<=r; i++)
    for (j=1; j<=c; j++)
      x[i][j] = t[i][j];
}


ComplexMatrix::ComplexMatrix(int r, int c)
  : rows(r), cols(c) {
  Complex ** t = myComplexMatrix(1, r, 1, c);
  int i, j; 
  for(i=1; i<=r;i++)
    for(j=1; j<=c; j++)
      t[i][j] = Complex(0);
  setMat(t, r, c);
  freeMyComplexMatrix(t, 1, r, 1, c);
}

// Makes diagonal square matrix out of W
ComplexMatrix::ComplexMatrix(const ComplexVector &v)
  : rows(v.displayEnd()), cols(v.displayEnd()) { 
  int e = v.displayEnd();
  Complex ** t = myComplexMatrix(1, e, 1, e);
  int i, j;
  for (i=1; i<=e; i++)
    for (j=1; j<=e; j++)
      {	if (i==j) t[i][j] = v.display(i);
      if (i!=j) t[i][j] = Complex(0);
      }
  setMat(t, e, e); 
  freeMyComplexMatrix(t, 1, e, 1, e);
}

const ComplexMatrix & ComplexMatrix::operator=(const ComplexMatrix &v) {
  if (this == &v) return *this;
  // Are the vectors of equivalent size?
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != v.rows || cols != v.cols)  {
    ostringstream ii;
    ii << "ComplexMatrix = overload; incompatible sizes\n";
    ii << *this << "=\n" << v;
    throw ii.str();
  }
#endif
  freeMyComplexMatrix(x, 1, rows, 1, cols);
  setMat(v.x, v.rows, v.cols);
  return *this;
}

const ComplexMatrix & ComplexMatrix::operator=(const Complex &v) {
  // Are the vectors of equivalent size?
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "ComplexMatrix = Complex overload; must be square\n";
    ii << *this << "=\n" << v;
    throw ii.str();
  }
#endif
  freeMyComplexMatrix(x, 1, rows, 1, cols);
  ComplexMatrix temp(rows, cols);
  int i; for (i=1; i<=rows; i++) temp(i, i) = v;
  setMat(temp.x, temp.rows, temp.cols);
  return *this;
}

ComplexMatrix ComplexMatrix::operator+(const ComplexMatrix &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != v.displayRows() || cols != v.displayCols()) {
    ostringstream ii;
    ii << "ComplexMatrix + overload; incompatible lengths\n";
    ii << *this << "+\n" << v;
    throw ii.str();
  }
#endif
  ComplexMatrix temp(v);
  int i, j; 
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      temp(i, j) = v.display(i, j) + x[i][j];
  return temp;
}

ComplexMatrix ComplexMatrix::operator-(const ComplexMatrix &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != v.displayRows() || cols != v.displayCols()) {
    ostringstream ii;
    ii << "ComplexMatrix - overload; incompatible lengths\n";
    ii << *this << "-\n" << v;
    throw ii.str();
  }
#endif
  ComplexMatrix temp(v);
  int i, j; 
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      temp(i, j) = x[i][j] - v.display(i, j);
  return temp;
}

ComplexMatrix ComplexMatrix::operator+(Complex f) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols){
    ostringstream ii;
    ii << "ComplexMatrix + Complex overload; not square\n";
    ii << *this;
    throw ii.str();
  }
#endif
  ComplexMatrix temp(*this);
  int i; 
  for(i=1; i<=rows; i++) 
    temp(i, i) = x[i][i] + f;
  return temp;
}

ComplexMatrix ComplexMatrix::operator-(Complex f) {
  ComplexMatrix temp(*this);
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "Subtracting Complex to non-square matrix\n" << temp;
    throw ii.str();
  }
#endif
  int i; 
  for (i=1; i<=rows; i++)
    temp(i, i) = temp(i, i) - f;
  return temp;  
}

ComplexMatrix ComplexMatrix::operator*(const ComplexMatrix &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (cols != v.displayRows()){     
    ostringstream ii;
    ii << "ComplexMatrix * overload; incompatible sizes\n" << *this << "*\n" << v;
    throw ii.str();
  }
#endif
  ComplexMatrix temp(rows, v.displayCols());
  int i, j, k;
  for(i=1; i<=rows; i++) 
    for(j=1; j<=v.displayCols(); j++) {
      for (k=1; k<=cols; k++)
	temp(i, j) = temp(i, j) + x[i][k] * v.display(k, j);
    }
  return temp;
}

ComplexMatrix ComplexMatrix::operator*(const Complex f) {
  ComplexMatrix temp(*this);
  int i, j;
  for(i=1; i<=rows; i++) 
    for(j=1; j<=cols; j++)
      temp(i, j) = temp(i, j) * f;
  return temp;
}

Complex ComplexMatrix::trace() const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols)  {
    ostringstream ii;
    ii << "Complexrace of non-square matrix\n" << *this;
    throw ii.str();
  }
#endif
  Complex sum = Complex(0);
  // this-> after +=
  int i; for(i=1; i<=rows; i++) sum = sum + display(i, i);
  return sum;
}

ComplexMatrix ComplexMatrix::transpose() const {
  ComplexMatrix temp(*this);
  temp.rows = cols;
  temp.cols = rows;
  int i, j; 
  for (i=1; i<=rows; i++)
    for (j=1; j<=cols; j++)
      temp(j, i) = (*this).display(i, j);
  return temp;
}

ComplexMatrix ComplexMatrix::hermitianConjugate() const {
  ComplexMatrix temp(this->displayCols(), this->displayRows());
  temp.rows = cols;
  temp.cols = rows;
  int i, j; 
  for (i=1; i<=rows; i++)
    for (j=1; j<=cols; j++)
      temp(j, i) = (*this).display(i, j).conj();
  return temp;
}

ComplexMatrix ComplexMatrix::complexConjugate() const {
  ComplexMatrix temp(*this);
  temp.rows = cols;
  temp.cols = rows;
  int i, j; 
  for (i=1; i<=rows; i++)
    for (j=1; j<=cols; j++)
      temp(i, j) = (*this).display(i, j).conj();
  return temp;
}


void ComplexMatrix::swaprows(int i, int j) {
  ComplexMatrix temp(*this);
  for (int k=1; k<= cols; k++)
    {
      temp(i,k) = this->display(j,k);
      temp(j,k) = this->display(i,k);
    };
  *this = temp;
}

void ComplexMatrix::swapcols(int i,int j) {
  ComplexMatrix temp(*this);
  for (int k=1; k<= rows; k++)
    {
      temp(k, i) = this->display(k, j);
      temp(k, j) = this->display(k, i);
    };
  *this = temp;
}

void ComplexMatrix::symmetrise() {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "Error: symmetrising rectangular matrix " << *this;
    throw ii.str();
  }
#endif
  int i, j;
  for (i=2; i<=rows; i++)
    for (j=1; j<i; j++)
      x[i][j] = x[j][i];
}

// Gives sum of difference between two matrices
double ComplexMatrix::compare(const ComplexMatrix & a) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != a.displayRows() || displayCols() !=
      a.displayCols()) {
    ostringstream ii;
    ii << "Error: comparing matrices of different sizes" << *this << 
      " and " << a;
    throw ii.str();
  }
#endif
  double c = double(0);
  int i, j; 
  for (i=1; i<=displayRows(); i++)
    for (j=1; j<=displayCols(); j++)
      c = c + (a.display(i, j) - display(i, j)).mod();
  return c;
}

Complex ComplexMatrix::min(int & k, int & l) const {
  Complex m = display(1, 1);
  int i, j; 
  for (i=1; i<=displayRows();  i++)
    for (j=1; j<=displayCols(); j++)
      if (display(i, j).mod() < m.mod()) {
	m = display(i, j);
	k = i; l = j;
      }
  return m;
}

// You have to free memory for temp before this: dumps matrix in temp
void ComplexMatrix::displayMat(Complex ** temp) const {
  int i, j; 
  for(i=1; i<=displayRows(); i++)
    for(j=1; j<=displayCols(); j++)
      temp[i][j] = display(i, j);
}

// Assumes the Complex is multiplied by identity
ComplexMatrix operator*(Complex f, const ComplexMatrix &m) {
  ComplexMatrix temp(m);
  int i, j; 
  for (i=1; i<=temp.displayRows(); i++)
    for (j=1; j<=temp.displayCols(); j++)
      temp(i, j) = temp(i, j) * Complex(f);
  return temp;
}

ComplexMatrix operator*(double f, const ComplexMatrix &m) {
  ComplexMatrix temp(m);
  int i, j; 
  for (i=1; i<=temp.displayRows(); i++)
    for (j=1; j<=temp.displayCols(); j++)
      temp(i, j) = temp(i, j) * f;
  return temp;
}

// vI = MIj vJ
ComplexVector ComplexMatrix::operator*(const ComplexVector & v) {
  ComplexMatrix m(*this);
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayCols() != v.displayEnd() || v.displayStart() != 1) {
    ostringstream ii;
    ii << "matrix times vector overload size problem\n" << m << "*" << v;
    throw ii.str();
  }
#endif  
  
  ComplexVector temp(m.displayRows());
  int i, j;
  for(j=1; j<=m.displayRows(); j++) {
    temp(j) = Complex(0);
    for(i=1; i<=m.displayCols(); i++)
      temp(j) = temp(j) + m.display(j, i) * v.display(i);
  }
  return temp;
}

ComplexMatrix operator*(const DoubleMatrix & a, const ComplexMatrix &v)
{
#ifdef ARRAY_BOUNDS_CHECKING
  if (a.displayCols() != v.displayRows()) {
    ostringstream ii;
    ii << "ComplexMatrix * overload; incompatible sizes\n" << a <<
      "*\n" << v; 
    throw ii.str();
  }
#endif
  ComplexMatrix temp(a.displayRows(), v.displayCols());
  int i, j, k;
  for(i=1; i<=a.displayRows(); i++) 
    for(j=1; j<=v.displayCols(); j++) {
      for (k=1; k<=a.displayCols(); k++)
	temp(i, j) = temp(i, j) + a.display(i, k) * v.display(k, j);
    }
  return temp;
  
}

// v_i = v_j M_ji 
ComplexVector operator*(const ComplexVector & v, const ComplexMatrix & m) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayCols() != v.displayEnd() || v.displayStart() != 1) {
    ostringstream ii;
    ii << "ComplexMatrix times vector overload size problem\n" << m << "*"
       << v;
    throw ii.str();
  }
#endif
  ComplexVector temp(m.displayRows());
  int i, j;
  for(j=1; j<=m.displayRows(); j++) {
    temp(j) = Complex(0);
    for(i=1; i<=m.displayCols(); i++)
      temp(j) = temp(j) + v.display(i) * m.display(i, j);
  }
  return temp;
}

ComplexVector operator*(ComplexMatrix & m, DoubleVector & v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayCols() != v.displayEnd() || v.displayStart() != 1) {
    ostringstream ii;
    ii << "ComplexMatrix times vector overload size problem\n" << m << "*"
       << v;
    throw ii.str();
  }
#endif
  ComplexVector temp(m.displayRows());
  int i, j;
  for(j=1; j<=m.displayRows(); j++) {
    temp(j) = Complex(0.0);
    for(i=1; i<=m.displayCols(); i++)
      temp(j) = temp(j) + Complex(v.display(i)) * m.display(j, i);
  }
  return temp;
}


ComplexMatrix operator-(Complex f, const ComplexMatrix &m) {
  ComplexMatrix temp(m);
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayRows() != m.displayCols()) {
    ostringstream ii;
    ii << "Subtracting Complex from non-square matrix\n" << m;
    throw ii.str();
  }
#endif
  int i, j; 
  for (i=1; i<=temp.displayRows(); i++)
    for (j=1; j<=temp.displayCols(); j++) {
      if (i != j) temp(i, j) = Complex(0.0) - temp(i, j);
      if (i == j) temp(i, j) = f - temp(i, j);
    }
  return temp;
}

ComplexMatrix operator+(Complex f, const ComplexMatrix &m) {
  ComplexMatrix temp(m);
#ifdef ARRAY_BOUNDS_CHECKING
  if (m.displayRows() != m.displayCols()) {
    ostringstream ii;
    ii << "Subtracting Complex from non-square matrix\n" << m;
    throw ii.str();
  }
#endif
  int i; 
  for (i=1; i<=temp.displayRows(); i++)
    temp(i, i) = temp(i, i) + Complex(f);
  return temp;
}

// Mij = v1i v2j
ComplexMatrix outerProduct
(const ComplexVector & v1, const ComplexVector & v2){ 
  ComplexMatrix temp(v1.displayEnd(), v2.displayEnd());
  int i,j;
  for (i=v1.displayStart(); i<=v1.displayEnd(); i++)
    for (j=v2.displayStart(); j<=v2.displayEnd(); j++)
      temp(i, j) = v1.display(i) * v2.display(j);
  return temp;
}


istream & operator>>(istream & left, ComplexMatrix &m) {
  int i, j;
  ComplexMatrix x(m.displayRows(), m.displayCols());
  for (i=1; i<=m.displayRows(); i++)
    for (j=1; j<=m.displayCols(); j++)
      left >> x(i, j);
  m = x;
  
  return left;
}

ostream & operator <<(ostream &left, const ComplexMatrix &v) {
  left << "(" << v.displayRows() << "," <<
    v.displayCols() << "):\n";
  int i, j; 
  for(i=1; i<=v.displayRows(); i++) {
    for(j=1; j<=v.displayCols(); j++) {
      if (v.display(i, j).real() >= double(0.0))
	left << " ";  
      left << v.display(i, j) << " ";
    }
    left << '\n';
  }
  return left;
}

ComplexMatrix ComplexMatrix::operator*(const DoubleMatrix &v) {
#ifdef ARRAY_BOUNDS_CHECKING 
  if (cols != v.displayRows()) {
    ostringstream ii;
    ii << "ComplexMatrix * overload; incompatible sizes\n" << *this <<
      "*\n" << v; 
    throw ii.str();
  }
#endif
  ComplexMatrix temp(rows, v.displayCols());
  int i, j, k;
  for(i=1; i<=rows; i++) 
    for(j=1; j<=v.displayCols(); j++) {
      for (k=1; k<=cols; k++)
	temp(i, j) = temp(i, j) + x[i][k] * v.display(k, j);
    }
  return temp;
}

