#ifndef POLY1D_H
#define POLY1D_H

#include <iostream>
#include <cmath>
#include <cstdarg>

using namespace std;

#include "Polynomial.h"

/** 
* @class Poly1D
*
* @brief This is a class for polynomial in one dimension.
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class Poly1D : public Polynomial {

public:
  /// Default constructor.	
  Poly1D();

  /// Constructor.
  Poly1D(const int);

  /// Copy constructor.
  Poly1D(const Poly1D&);

  /// Destructor.
  ~Poly1D();

  ///
  Poly1D& operator=(const Poly1D&);

  /// This is a function that set the degree of the polynomial.
  void set_degree(const int);

  /// This is a function that set the coefficients of the polynomial.
  void set_coeffs(...);

  /// This is a function that evaluates the polynomial at point x.
  double get_value(const double x) const;

  /// This is a function that evaluates the first order derivative of the polynomial with respect to x at point x.
  double get_dx_value(const double x) const;

  /// Second order derivative with respect to x.
  double get_dxx_value(const double x) const;

  /// This is a function that prints to the screen some information about the polynomial.
  void info() const;

};

#endif // POLY1D_H
