#ifndef POLY2D_H
#define POLY2D_H

#include <iostream>
#include <cmath>
#include <cstdarg>

using namespace std;

#include "Polynomial.h"

/** 
* @class Poly2D
*
* @brief This is a class for polynomial in two dimension.
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class Poly2D : public Polynomial {

public:
  /// Default constructor.
  Poly2D();

  /// Constructor.
  Poly2D(const int);

  /// Copy constructor.
  Poly2D(const Poly2D&);

  /// Destructor.
  ~Poly2D();

  ///
  Poly2D& operator=(const Poly2D&);

  /// This is a function that set the degree of the polynomial.
  void set_degree(const int);

  /// This is a function that set the coefficients of the polynomial.
  void set_coeffs(...);

  /// This is a function that evaluates the polynomial at point x, y.
  double get_value(const double x,
                   const double y) const;

  /// This is a function that evaluates the first order derivative of the polynomial with respect to x at point x, y.
  double get_dx_value(const double x,
                      const double y) const;

  /// First order derivative with respect to y.
  double get_dy_value(const double x,
                      const double y) const;

  /// Second order derivative with respect to x.
  double get_dxx_value(const double x,
                       const double y) const;

  /// Mixed second order derivative.
  double get_dxy_value(const double x,
                       const double y) const;

  /// Second order derivative with respect to y.
  double get_dyy_value(const double x,
                       const double y) const;

  /// This is a function that prints to the screen some information about the polynomial.
  void info() const;

private:
  /// Dimension (is equal to 2).
  int dim;

  /// Degree of the polynomial.
  int degree;

  /// Number of coefficients.
  int ncoeffs;

  /// The coefficients (lowest degree first, in the manner of Pascal's triangle).
  double* coeffs;

};

#endif // POLY2D_H
