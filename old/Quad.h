#ifndef QUAD_H
#define QUAD_H

#include <iostream>

/** 
* @class Quad
*
* @brief This is a class for the quadrature.
*
* ausdb
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class Quad {

public:
  /// Constructor.
  Quad(enum quad_type);
 
  /// Copy constructor.
  Quad(const Quad&);

  /// Destructor.  
  ~Quad();

  Quad& operator=(const Quad&);

  /// This function returns the number of quadrature points.
  int get_np() const;

  /// This function returns the ith weight.
  double get_w(int i) const;

  /// This functions returns the ith point.
  point get_lv(int i) const;

private:
  /// Quadrature type.  
  enum quad_type {
    /// Gaussian quadrature in 1D.
    TWOPT_GAUSS_ONED,
    /// Midpoint quadrature in 2D.
    MIDPOINT_TWOD
  } type;

  /// number of quadrature points
  int np;

  /// quadrature weights
  double *w;

  /// quadrature points in the local coordinates
  point *lx;

};

#endif // QUAD_H
