#ifndef BASISFUNC_H
#define BASISFUNC_H

#include <iostream>
#include "Node.h"
#include "Quad.h"
#include "RefElement.h"

using namespace std;

/** 
* @class BasisFunc
*
* @brief This is a class for basis functions.
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class BasisFunc {

public:
  /// Constructor.
  BasisFunc(enum element_type);

  /// Copy constructor.
  BasisFunc(const BasisFunc&);

  /// Destructor.
  ~BasisFunc();

  ///
  BasisFunc& operator=(const BasisFunc&);

  ///
  int get_number_of_basis_functions() const;

  ///
  double get_value(const int i,
                   const double x) const;

  ///
  double get_dx_value(const int i,
                      const double x) const;

  ///
  double get_value(const int i,
                   const double x,
                   const double y) const;

  ///
  double get_dx_value(const int i,
                      const double x,
                      const double y) const;

  ///
  double get_dy_value(const int i,
                      const double x,
                      const double y) const;

  ///
  void info() const;

private:
  /// Type of element.
  enum element_type {
    /// Triangles with piecewise linear basis functions.
    TRIANGLE_LINEAR,
    /// Triangles with piecewise quadratic basis functions.
    TRIANGLE_QUADRATIC, 
    /// One-dimensional element with piecewise piecewise linear basis functions.
    ONED_LINEAR,
    /// One-dimensional element with piecewise piecewise quadratic basis functions.
    ONED_QUADRATIC
  } type;

  /// Number of basis function.
  int nfunc;

  /// The basis functions.
  Polynomial* basis_functions;
};

#endif // BASISFUNC_H

