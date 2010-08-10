#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <cmath>
#include <cstdarg>

using namespace std;

/** 
* @class Polynomial
*
* @brief This is an (abstract) class for polynomial.
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

int factorial(const int);

class Polynomial {

public:
  /// Default constructor.
  Polynomial();

  /// Constructor.
  Polynomial(const int);

  /// Constructor.
  Polynomial(const int,
             const int);

  /// Copy constructor.
  Polynomial(const Polynomial&);

  /// Destructor.
  ~Polynomial();

  Polynomial& operator=(const Polynomial&);

  ///
  virtual void set_degree(const int)
  {}

  ///
  virtual void set_coeffs(...)
  {}

  /// This is a function that evaluates the polynomial at point x.
  virtual double get_value(const double) const;

  /// This is a function that evaluates the first order derivative of the polynomial with respect to x at point x.
  virtual double get_dx_value(const double) const;

  /// Second order derivative with respect to x.
  virtual double get_dxx_value(const double) const;

  /// This is a function that evaluates the polynomial at point x, y.
  virtual double get_value(const double,
                           const double) const;

  /// This is a function that evaluates the first order derivative of the polynomial with respect to x at point x, y.
  virtual double get_dx_value(const double,
                              const double) const;

  /// First order derivative with respect to y.
  virtual double get_dy_value(const double,
                              const double) const;

  /// Second order derivative with respect to x.
  virtual double get_dxx_value(const double,
                               const double) const;

  /// Mixed second order derivative.
  virtual double get_dxy_value(const double,
                               const double) const;

  /// Second order derivative with respect to y.
  virtual double get_dyy_value(const double,
                               const double) const;

  /// This is a function that prints to the screen some information about the polynomial.
  virtual void info() const;

protected:
  /// Dimension of the polynomial.
  int dim;

  /// Degree of the polynomial.
  int degree; 

  /// Number of coefficients.
  int ncoeffs;

  /// The coefficients.
  double* coeffs; 

};

#endif // POLYNOMIAL_H
