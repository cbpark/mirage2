
/** \file mycomplex.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: complex numbers and operators between them

   $Log: mycomplex.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.5  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.4  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.3  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.3  2001/10/04 19:26:34  allanach
   New version deals with AMSB correctly

   Revision 1.2  2001/10/03 13:34:16  allanach
   Changed name of complex.h to avoid conflict with STD libraries
   Added proper header info
*/

#ifndef COMPLEX_H
#define COMPLEX_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "def.h"
#include "utils.h"

/// "Home-grown" complex class with some complex algebra included
class Complex {
private:
  double re, im; ///< Real and imaginary parts
public:
  Complex() : re(0.0), im(0.0) {}; ///< Empty (zero) constructor
  Complex(double r) : re(r), im(0.0) {}; ///< Real number constructor
  Complex(double r, double i) : re(r), im(i) {}; ///< General constructor
  /// Constructor, setting equal to another Complex number
  Complex(const Complex & cc) : re(cc.re), im(cc.im) {};
  /// Sets contents equal to another Complex number
  const Complex &operator=(const Complex & s);

  double real() const { return re; } ///< returns real part
  double imag() const { return im; } ///< returns imaginary part
  
  void setRe(double a) { re = a; } ///< sets real part
  void setIm(double a) { im = a; } ///< sets imaginary part

  double mod() const; ///< returns modulus of number
  double arg() const; ///< returns angle (in Argand diagram): theta=-pi->pi

  Complex conj() const; ///< Complex conjugate

  Complex operator/(double b) const;
  Complex operator*(const Complex & b) const;  
  Complex operator*(const double & a) const;
  /// Returns true if either real part or imaginary part is bigger than
  /// relevant bit of a
  bool operator>=(const Complex & a) const;
  Complex operator/(const Complex & a) const;
  Complex cc() const; ///< Complex conjugate
  Complex operator-(const Complex & a) const;
};

Complex log(const Complex &a); ///< Principal log
Complex exp(const Complex &a);
Complex sqrt(const Complex &a);

Complex operator*(double a, const Complex & b);
Complex operator+(const Complex & a, const Complex & b);
Complex operator/(double a, const Complex & b);
Complex operator-(double b, const Complex & a);
/// Formatted output
ostream & operator <<(ostream & , const Complex &);
/// Formatted input
inline istream & operator >> (istream & left, Complex & v){
  double r, i;
  left >> r >> i;
  v.setRe(r); v.setIm(i);
  return left;
}

#endif
