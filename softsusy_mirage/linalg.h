
/** \file linalg.h
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach 
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
    - Description: DoubleVector and DoubleMatrix classes of doubles and 
                operations between them, complexified copies also

   $Log: linalg.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.3  2005/08/16 17:22:08  allanach
   Corrected electroweak sbottom corrections

   Revision 1.2  2005/07/26 10:51:22  allanach
   new method nmin added to DoubleVector and DoubleMatrix by Bernhardt

   Revision 1.20  2004/01/15 13:53:48  allanach
   New time-saving unary operators included.

   Revision 1.19  2003/08/04 16:44:22  allanach
   Added max function for matrix and vector of doubles

   Revision 1.18  2003/07/29 15:26:40  allanach
   Using new and delete rather than old malloc and free C-style routines

   Revision 1.17  2003/07/28 12:11:37  allanach
   More error trapping, and rearranging rpvsoftsusy to use correct Higgs VEV
   (which is sometimes called at MZ)

   Revision 1.16  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.15  2003/06/05 09:17:19  allanach
   Started coding Les Houches Discord

   Revision 1.14  2003/05/27 15:05:51  allanach
   Purely efficiency corrections: used variable rather than display() methods
   whenever possible

   Revision 1.13  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.12  2003/03/04 12:55:02  allanach
   Removed Slavich diagonalisation

   Revision 1.10  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.9  2002/09/23 18:13:50  allanach
   Eigenvalue order with non symmetric diagonalisation setx

   Revision 1.8  2002/07/19 15:54:06  allanach
   SOFTSUSY1.5 version

   Revision 1.7  2002/06/10 16:13:07  allanach
   Speeded things up by adding option for array bounds checking

   Revision 1.5  2002/04/30 16:23:59  allanach
   Added matrix inversion routine

   Revision 1.4  2001/10/03 13:34:16  allanach
   Changed name of complex.h to avoid conflict with STD libraries

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef LINALG_H
#define LINALG_H

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include "def.h"
#include "utils.h"
#include "mycomplex.h"

/// Sorts memory out for the vectors. Violates C-conventions, and starts at 1
/// rather than 0 by default
inline void freeMyDoubleVector(double *v, long nl, long nh) {
  //  free((char *)(v+nl-1)); 
  delete [] (v+nl-1);
}

/// DoubleVector is of variable length, and contains double precision
/// numbers. 
class DoubleVector {
private:
  int start, end; ///< vector dimensions
  double * x;     ///< array of values held in *x
  
  /// Memory handling routines: allocate a double matrix with subscript range
  /// m[nrl..nrh][ncl..nch] 
  double *myDoubleVector(const long nl, const long nh);
  void setVec(const double *t, long st, long nd);

public:
  DoubleVector(int e); ///< Default dimensions start at 1, end at e
  DoubleVector(int s, int e); ///< Dimensions start at s, end at e
  /// Sets new double vector equal to v
  inline DoubleVector(const DoubleVector &v); 
  ~DoubleVector(); ///< Destructor

  /// Changes the length of a vector - copies as many elements of old one as
  /// possible, and fills any extra up with zeroes
  void setEnd(int e); 
  /// Standard linear algebra follows
  const DoubleVector & operator=(const DoubleVector &v); 
  DoubleVector operator+(const DoubleVector &v); 
  DoubleVector operator-(const DoubleVector &v);
  /// NOT dot product, but product of elements.
  DoubleVector operator*(const DoubleVector &v);
  DoubleVector operator*(double f);
  double & operator() (int i); ///< Reference a single element
  double dot(const DoubleVector & v) const; ///< dot product

  /// Obvious display functions
  inline double display(int i) const { return x[i]; }; ///< display one element
  int displayStart() const { return start; }; /// start of dimension
  int displayEnd() const { return end; }; /// returns end of dimension
  DoubleVector display() const { return *this; }; /// returns whole thing
  /// Displays all of vector in *a. a is in old C convention ie index starts at
  /// ZERO! (For outputting to fortran)
  void display(double * a) const; 

  /// applies fn to every element in a vector
  DoubleVector apply(double (*fn)(double)); 
  void set(int i, double f); ///< sets ith element
  double max() const; ///< maximum element in vector
  double min(int & p) const; ///< minimum element in vector
  double nmin(int & p) const; ///< next-to-minimum element in vector
  
  void swap(int i, int j); ///< swaps ith and jth elements of a vector
};

/// float times vector
DoubleVector operator*(double f, const DoubleVector & V);
/// prints a vector out formatted, maximum of 5 elements per line
ostream & operator <<(ostream &left, const DoubleVector &V);
/// inputs a vector
istream & operator>>(istream & left, DoubleVector &V);

inline DoubleVector::~DoubleVector() { freeMyDoubleVector(x, start, end); } 

inline DoubleVector::DoubleVector(const DoubleVector &v)
  : start(v.start), end(v.end)
{ setVec(v.x, v.start, v.end); } 

inline double & DoubleVector::operator() (int i) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i>end || i<start) {
      cout << "Trying to access " << i << "th element of DoubleVector\n";
      cout << "start " << start << " end " << end;
      cout << *this;
    }
#endif
  return x[i];
}

inline void DoubleVector::set(int i, double f) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i < start || i > end) 
    {
      ostringstream ii;
      ii << "Cannot access " << i << "th element of vector " << *this;
      throw ii.str();
    }
#endif
  x[i] = f;
}


/** free a double matrix allocated by myDoubleMatrix() */
inline void freeMyDoubleMatrix(double **m, long nrl, long nrh, long ncl, long
		   nch) {
  delete [] (m[nrl]+ncl-1);
  delete [] (m+nrl-1);
  //  free((char *) (m[nrl]+ncl-1));
  //  free((char *) (m+nrl-1));
}

class DoubleMatrix; class ComplexMatrix;

/// Standard linear algebra
inline DoubleMatrix operator-(const DoubleMatrix & A);
/// Assumes the double is multiplied by identity
DoubleMatrix operator*(double f, const DoubleMatrix &M);
// \f$ v_i = v_j M_{ji} \f$
DoubleVector operator*(const DoubleVector & V, const DoubleMatrix & M);
DoubleMatrix operator-(double f, const DoubleMatrix &M);
DoubleMatrix operator+(double f, const DoubleMatrix &M);

/// Defines \f$ Mij = v1i v2j \f$
DoubleMatrix outerProduct(const DoubleVector & v1, const DoubleVector & v2);
/// Formatted input of a matrix: size already determined in M
istream & operator>>(istream & left, DoubleMatrix &M);
/// Formatted output of a matrix: maximum 5 elements per line
ostream & operator <<(ostream &left, const DoubleMatrix &V);
/// Returns a 2x2 orthogonal matrix of rotation by angle theta
// [ cos theta    sin theta ]
// [-sin theta    cos theta ]
DoubleMatrix rot2d(double theta);
/// Returns 2 by 2 orthogonal mixing matrix
// [ -sin(theta)  cos(theta) ]
// [  cos(theta)  sin(theta) ] --
DoubleMatrix rot2dTwist(double theta);

/// Matrix from 1..rows, 1..cols of double values
class DoubleMatrix {
private:
  int rows, cols;
  double **x; ///< Memory allocation goes in **x
  
  double **myDoubleMatrix(long nrl, long nrh, long ncl, long nch) const;
  double *myDoubleVector(const long nl, const long nh) const;
  void setMat(double **t, int r, int c); ///< Sets whole matrix equal to **t
  
public:
  DoubleMatrix(int r, int c); ///< Constructor for matrix of zeroes 1..r,1..c
  inline DoubleMatrix(const DoubleMatrix &m); ///< Constructor sets = m
  
  /// Makes diagonal square matrix: diagonal elements are from v, in order
  DoubleMatrix(const DoubleVector &v);  
  ~DoubleMatrix(); ///< Standard destructor
  
  /// Routines for outputting size of matrix
  int displayRows() const { return rows; };
  int displayCols() const { return cols; };

  /// Sets matrix equal to v
  const DoubleMatrix & operator=(const DoubleMatrix &v);  
  /// Sets diagonal entries equal to v, rest are 0
  const DoubleMatrix & operator=(const double &v);

  /// Standard Linear algebra calculations  
  DoubleMatrix operator+(const DoubleMatrix &v) const;  
  DoubleMatrix operator-(const DoubleMatrix &v) const;  
  DoubleMatrix operator+(double f) const;  
  DoubleMatrix operator-(double f) const;  
  DoubleMatrix operator*(const DoubleMatrix &v) const;  
  DoubleMatrix operator*(const double f) const;
  void operator*=(double f);
  void operator+=(DoubleMatrix & f);
  /// \f$ v_i = M_{ij} v_j \f$
  DoubleVector operator*(const DoubleVector &v) const;
  DoubleMatrix operator/(const double f) const;

  /// trace must only be performed on a square matrix
  double trace() const;  
  double & operator() (int i, int j); ///< to reference one element
  DoubleMatrix transpose() const; ///< can be any size

  /// Obvious elementary row/column operations
  void swaprows(int i, int j);
  void swapcols(int i,int j);
  /// Perform on a mixing matrix O such that
  /// diag = O^double M O. It will put the diag evals in abs ascending order 
  /// and change O accordingly
  void associateOrderAbs(DoubleVector &v);
  double display(int i, int j) const; ///< ijth element
  DoubleMatrix display() const { return *this; }; ///< whole matrix returned

  /// Returns sum of absolute values of all elements
  double sumElements() const;  
  /// fills in bottom-left hand corner of a matrix
  void symmetrise(); 
  /// Sums up the absolute value of the difference in each element of two
  /// matrices 
  double compare(const DoubleMatrix & a) const;
  double compare(const ComplexMatrix & a) const;
  /// Just returns a vector filled with the diagonal values 
  DoubleVector diagVals() const;
  /// Uses singular SVD algorithm:
  /// \f$ A = U.W.V^T \f$ where W is a matrix of the eigenvalues,
  /// therefore \f$ W = U^T A V \f$.
  double diagonalise(DoubleMatrix & u, DoubleMatrix & v, DoubleVector & w)
    const;
  /// For SYMMETRIC MATRICES ONLY!
  /// \f$ A = V.W.V^T \f$ where W is a matrix of the eigenvalue therefore 
  /// \f$ W = V^T A V \f$. 
  double diagonaliseSym(ComplexMatrix & v, DoubleVector & w) const;
  double diagonaliseSym(DoubleMatrix & v, DoubleVector & w) const;
  DoubleMatrix inverse() const; ///< returns inverse of a matrix

  /// Special case that's often used: eigenvalues of 2 by 2 symmetric matrix
  /// A(2,1) is assumed to be equal to A(1,2). theta is the mixing angle:
  /// m1 is NOT NECESSARILY > m2 
  // -[ cos theta    sin theta ]   A   [ cos theta -sin theta ] = diagonal
  // -[ -sin theta   cos theta ]       [ sin theta  cos theta ]
  DoubleVector sym2by2(double & theta) const; 
  /// Applies fn to every element of a matrix
  DoubleMatrix apply(double (*fn)(double));
  /// 2 by 2 asymmetric matrices
  // -[ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ] = diagonal
  // -[ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
  DoubleVector asy2by2(double & thetaL, double & thetaR) const;
  double min(int & k, int & l) const; ///< minimum element
  double nmin(int & k, int & l) const; ///< next-to-minimum element
  double max(int & k, int & l) const; ///< maximum element
  /// returns matrix contents in **temp. You have to free memory for temp
  /// before this 
  void displayMat(double ** temp) const; 
  // Perform on asymmetric matrices M such that
  /// diag = U^T double M V. It will put the diag evals in abs ascending order 
  /// and changes U,V accordingly
  void associateOrderAbs(DoubleMatrix & u, DoubleMatrix & v, DoubleVector &w) 
    const;
};

/// Redefines mixing matrices to be complex such that diagonal values are
/// positive for a 2 by 2: 
// [ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ]  = diag
// [ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
// as given by asy2by2!
/// \f$ u^* A v^+ \f$ = mdiagpositive
void positivise(double thetaL, double thetaR, const DoubleVector & diag,
		  ComplexMatrix & u, ComplexMatrix & v);
/// Diagonalisation routines
void diagonaliseSvd(DoubleMatrix & a, DoubleVector & w, DoubleMatrix & v);
double pythagoras(double a, double b);
void diagonaliseJac(DoubleMatrix & a,  int n,  DoubleVector & d, DoubleMatrix
		    & v,  int *nrot);

inline DoubleMatrix::DoubleMatrix(const DoubleMatrix &m) 
  : rows(m.rows), cols(m.cols) { setMat(m.x, rows, cols); }

inline DoubleMatrix::~DoubleMatrix() { freeMyDoubleMatrix(x, 1, rows, 1, cols); }

inline DoubleMatrix DoubleMatrix::operator/(const double f) const {
  return *this * (1.0 / f);
}

inline double & DoubleMatrix::operator() (int i, int j) { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i > rows || j > cols || i < 1 || j < 1) { 
    ostringstream ii;
    ii << "Trying to access " << i << "," << j << 
      "th element of DoubleMatrix\n" << *this;
    throw ii.str();
  }
#endif
  return x[i][j]; 
}  

inline double DoubleMatrix::display(int i, int j) const { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i < 1 || i > displayRows() || j < 1 || j > displayCols()) { 
    ostringstream ii;
    ii << "Error: Requested display (" << i << "," << j << 
	")th element of matrix " << *this; 
    throw ii.str();
    }
#endif
  return x[i][j]; 
}

inline DoubleMatrix operator-(const DoubleMatrix & A) { return -1.0 * A; }

/// Memory handling routine: dimension starts at 1 by default
inline void freeMyComplexVector(Complex *v, long nl, long nh) {
  //  free((char *)(v+nl-1)); 
  delete [] (v+nl-1);
}

/// Vector of double complex values
class ComplexVector {
private:
  int start, end; ///< Dimensions of vector
  Complex * x; ///< Memory storage of contents of vector
  
  /// Memory handling routines
  Complex *myComplexVector(const long nl, const long nh);
  void setVec(const Complex *t, long st, long nd);

public:
  ComplexVector(int e); ///< Default constructor sets e as dimension of vector
  ComplexVector(int s, int e); ///< Dimension starts at s and ends at e
  inline ComplexVector(const ComplexVector &v); ///< Sets contents = v
  ~ComplexVector(); ///< Destructor

  /// Changes the length of a vector - copies as many elements of old one as
  /// possible, and fills any extra up with zeroes
  void setEnd(int e); 
  
  ComplexVector & operator=(const ComplexVector &v);
  ComplexVector operator+(const ComplexVector &v);
  ComplexVector operator-(const ComplexVector &v);
  /// NOT Complex dot product, but product of elements.
  ComplexVector operator*(const ComplexVector &v);
  ComplexVector operator*(Complex f);
  Complex & operator() (int i); ///< reference one element
  Complex dot(const ComplexVector & v) const; ///< standard dot product
  inline Complex display(int i) const { return x[i]; }; ///<display ith element
  int displayStart() const { return start; }; ///< displays start of dimension
  int displayEnd() const { return end; }; ///< displays end of dimension
  ComplexVector display() const { return *this; }; ///< displays whole vector

  /// Apply fn to every element
  ComplexVector apply(Complex (*fn)(Complex)); 
  void set(int i, Complex f); ///< set ith element to f
  Complex max() const; ///< maximum absolute value
  Complex min(int & p) const; ///< smallest absolute element
  
  void swap(int i, int j); ///< swap ith and jth elements
};

ComplexVector operator*(double f, const ComplexVector & V);
ComplexVector operator*(Complex f, const ComplexVector & V);
/// Formatted output
ostream & operator<<(ostream &left, const ComplexVector &V); 
/// Formatted input
istream & operator>>(istream & left, ComplexVector &V);

inline ComplexVector::~ComplexVector() { freeMyComplexVector(x, start, end); }

inline ComplexVector::ComplexVector(const ComplexVector &v)
  : start(v.start), end(v.end)
{ setVec(v.x, v.start, v.end); } 

inline Complex & ComplexVector::operator() (int i) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i>end || i<start) {
      cout << "Trying to access " << i << "th element of DoubleVector\n";
      cout << "start " << start << " end " << end;
      cout << *this;
    }
#endif
  return x[i];
}

inline void ComplexVector::set(int i, Complex f) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i < start || i > end) {
      ostringstream ii;
      ii << "Cannot access " << i << "th element of vector " << *this;
      throw ii.str();
    }
#endif
  x[i] = f;
}


/// free a Complex matrix allocated by myComplexMatrix() 
inline void freeMyComplexMatrix(Complex **m, long nrl, long nrh, long ncl, long
		   nch) {
  //  free((char *) (m[nrl]+ncl-1));
  // free((char *) (m+nrl-1));
  delete [] (m[nrl]+ncl-1);
  delete [] (m+nrl-1);
}

class ComplexMatrix;
ComplexMatrix operator*(const ComplexMatrix & A, const ComplexMatrix & V);
/// Assumes the Complex is multiplied by identity
ComplexMatrix operator*(Complex f, const ComplexMatrix &M);
ComplexMatrix operator*(const DoubleMatrix & A, const ComplexMatrix &M);
ComplexMatrix operator*(double f, const ComplexMatrix &M);
ComplexVector operator*(const ComplexMatrix & M, const ComplexVector & V);
ComplexVector operator*(ComplexMatrix & M, DoubleVector & V);
/// \f$ v_i = v_j M_{ji} \f$
ComplexVector operator*(const ComplexVector & V, const ComplexMatrix & M);
ComplexMatrix operator-(Complex f, const ComplexMatrix &M);
ComplexMatrix operator+(Complex f, const ComplexMatrix &M);
/// \f$ M_{ij} = v1_i v2_j \f$
ComplexMatrix outerProduct(const ComplexVector & v1, const ComplexVector &
			   v2);

ComplexVector operator*(DoubleMatrix & m, ComplexVector & v);
/// Formatted input
istream & operator>>(istream & left, ComplexMatrix &M);
/// Formatted output
ostream & operator<<(ostream &left, const ComplexMatrix &V);

/// matrix of complex values, with standard linear algebra routines 
/// from 1..rows, 1..cols
class ComplexMatrix {
private:
  int rows, cols; ///< sets size of matrix
  Complex **x; ///< Memory stored in this variable
  /// allocate a Complex matrix with subscript range m[nrl..nrh][ncl..nch]  
  Complex **myComplexMatrix(long nrl, long nrh, long ncl, long nch) const;
  Complex *myComplexVector(const long nl, const long nh) const;
  /// set matrix according to elements of **t (size = r x c)
  void setMat(Complex **t, int r, int c); 
  
public:
  /// Default constructor: full of zeroes
  ComplexMatrix(int r, int c);
  /// Constructor sets matrix equal to m
  inline ComplexMatrix(const ComplexMatrix &m);
  /// Constructor sets square matrix's diagonal values equal to elements of v
  ComplexMatrix(const DoubleVector &v);
  /// Makes diagonal square matrix out of v
  ComplexMatrix(const ComplexVector &v);  
  ~ComplexMatrix();
  
  int displayRows() const { return rows; }; 
  int displayCols() const { return cols; };

  const ComplexMatrix & operator=(const ComplexMatrix &v);  
  const ComplexMatrix & operator=(const Complex &v);  
  ComplexMatrix operator+(const ComplexMatrix &v);  
  ComplexMatrix operator-(const ComplexMatrix &v);  
  /// Adds f times identity. Must be used on a square matrix
  ComplexMatrix operator+(Complex f);  
  /// Subtracts f times identity. Must be used on a square matrix
  ComplexMatrix operator-(Complex f); 
  ComplexMatrix operator*(const ComplexMatrix &v);  
  ComplexMatrix operator*(const Complex f);
  Complex trace() const;  
  /// \f$ v_i = M_{ij} v_j \f$
  ComplexVector operator*(const ComplexVector &);
  ComplexMatrix operator/(const Complex f);
  ComplexVector operator*(const DoubleVector & v);

  Complex & operator() (int i, int j); ///< Returns ijth element
  ComplexMatrix transpose() const;
  ComplexMatrix hermitianConjugate() const;
  ComplexMatrix complexConjugate() const;

  ComplexMatrix operator*(const DoubleMatrix &v);

  void swaprows(int i, int j); ///< Swaps row i with row j
  void swapcols(int i,int j); ///< Swaps column i with column j
  Complex display(int i, int j) const; ///< returns ijth element
  ComplexMatrix display() const { return *this; }; ///< returns whole matrix

  /// Fills in lower bottom half of a square matrix copying the top right
  void symmetrise();
  /// Returns the sum of the modulus of the difference of each element
  double compare(const ComplexMatrix & a) const;
  Complex min(int & k, int & l) const; ///< minimum element modulus value 
  /// whole matrix displayed in **temp. You have to free memory for temp
  /// before this 
  void displayMat(Complex ** temp) const;  
};

inline ComplexMatrix::ComplexMatrix(const ComplexMatrix &m) 
  : rows(m.rows), cols(m.cols) { setMat(m.x, rows, cols); }

inline ComplexMatrix::~ComplexMatrix() 
{ freeMyComplexMatrix(x, 1, rows, 1, cols); }

inline ComplexMatrix ComplexMatrix::operator/(const Complex f) {
  return *this * (1.0 / f);
}

inline Complex & ComplexMatrix::operator() (int i, int j) { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i > rows || j > cols || i < 0 || j < 0) { 
    ostringstream ii;
    ii << "Trying to access " << i << "," << j << 
      "th element of ComplexMatrix\n" << *this;
    throw ii.str();
  }
#endif
  return x[i][j]; 
}  

inline Complex ComplexMatrix::display(int i, int j) const { 
  if (i < 1 || i > displayRows() || j < 1 || j > displayCols()) { 
    ostringstream ii;
    ii << "Error: Requested display (" << i << "," << j << 
	")th element of matrix " << *this; 
    throw ii.str();
    }
  return x[i][j]; 
}

#endif
